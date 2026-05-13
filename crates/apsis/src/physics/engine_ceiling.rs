//! Engine ceiling profiling harness.
//!
//! Lab notebook: `docs/experiments/2026-05-09-engine-ceiling.md`.
//!
//! Two opt-in tests:
//!
//! ```text
//! # Cell V (VelocityVerlet + BH octree, full 4-phase decomposition):
//! cargo test --release -p apsis engine_ceiling_v -- --ignored --nocapture
//!
//! # Cell I (IAS15, total step time only — escalate via --features
//! # ias15-profile if internal decomposition is needed):
//! cargo test --release -p apsis engine_ceiling_i -- --ignored --nocapture
//! ```
//!
//! Cell V uses a hand-rolled VV step that mirrors the production
//! [`VelocityVerlet`] integrator (force at t → kick(½dt) → drift(dt) →
//! force at t+dt → kick(½dt)) so the per-phase timing fits the protocol's
//! 4-phase split (`t_tree_build`, `t_bh_walk`, `t_integrator_overhead`,
//! `t_trail_record`) without contaminating the production integrator with
//! more cfg-feature timing flags. Force evaluation goes through
//! [`BarnesHutEngine::evaluate_profile`], which returns the aggregated
//! [`WalkCounters`] for the per-interaction analysis in §Tier 2.
//!
//! Cell I drives [`Ias15`] through its production [`Integrator`] trait;
//! per-phase decomposition is left to the existing `ias15-profile` feature
//! flag (`crates/apsis/Cargo.toml`) and is escalated only on demand. The
//! cell I path measures total step wall time and reports it together with
//! IAS15's adaptive `dt` so the SPS curve is interpretable.
//!
//! CSV output: `target/engine-ceiling/profile.csv`. Per-row schema in the
//! [`write_header`] body.
//!
//! Frozen variables (matching the lab notebook §Methodology):
//!
//! * Seed: `0x6E63696C` ("ncil")
//! * Body distribution: sphere log-normal mass (matches perf 2×2 family)
//! * VV dt: `1e-3`; IAS15 dt_hint: `1e-3`
//! * θ: `0.5`
//! * Warmup steps: 10 (discarded)
//! * Measured steps: 100
//! * Hardware identifier: recorded in §Results

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use crate::domain::body::Body;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::force_model::GravityForceModel;
use crate::physics::integrator::ias15::Ias15;
use crate::physics::integrator::traits::{Integrator, IntegratorContext};

const SEED: u64 = 0x6E63696C;
const WARMUP_STEPS: usize = 10;
const MEASURED_STEPS: usize = 100;
const DT: f64 = 1.0e-3;
const THETA: f64 = 0.5;

// ── Cell V: VelocityVerlet + BH octree ─────────────────────────────────────── //

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release engine_ceiling_v -- --ignored --nocapture"]
fn engine_ceiling_v() {
    let n_values = [100usize, 1_000, 5_000, 10_000, 50_000, 100_000];
    let trail_variants = [false, true];

    let out_dir = ceiling_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf output dir");
    let csv_path = out_dir.join("profile_v.csv");
    let mut writer = fs::File::create(&csv_path).expect("create profile_v.csv");
    write_header_v(&mut writer);

    eprintln!("[ceiling-V] cell=V (VV+BH), N={:?}, trail={:?}", n_values, trail_variants);
    let t_total = Instant::now();
    for &trail_on in &trail_variants {
        for &n in &n_values {
            let row = measure_v(n, trail_on);
            write_row_v(&mut writer, &row);
            print_row_v(&row);
        }
    }
    eprintln!("[ceiling-V] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[ceiling-V] wrote {}", csv_path.display());
}

// ── Cell I: IAS15 ──────────────────────────────────────────────────────────── //

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release engine_ceiling_i -- --ignored --nocapture"]
fn engine_ceiling_i() {
    let n_values = [100usize, 1_000, 10_000];

    let out_dir = ceiling_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf output dir");
    let csv_path = out_dir.join("profile_i.csv");
    let mut writer = fs::File::create(&csv_path).expect("create profile_i.csv");
    write_header_i(&mut writer);

    eprintln!("[ceiling-I] cell=I (IAS15 exact O(N^2)), N={:?}", n_values);
    let t_total = Instant::now();
    for &n in &n_values {
        let row = measure_i(n);
        write_row_i(&mut writer, &row);
        print_row_i(&row);
    }
    eprintln!("[ceiling-I] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[ceiling-I] wrote {}", csv_path.display());
}

// ── Cell V measurement ─────────────────────────────────────────────────────── //

#[derive(Debug, Clone, Copy)]
struct RowV {
    n: usize,
    trail_on: bool,
    t_total_ns: u64,
    t_tree_build_ns: u64,
    t_bh_walk_ns: u64,
    t_integrator_ns: u64,
    t_trail_ns: u64,
    n_node_visits: u64,
    n_bh_accepted: u64,
    n_leaf_interactions: u64,
}

fn measure_v(n: usize, trail_on: bool) -> RowV {
    let mut bodies = sphere_distribution_lognormal(n, SEED);
    let mut engine = BarnesHutEngine::new(16);
    let mut acc = vec![Vec3::ZERO; n];
    let mut arrays = BodyArrays::with_capacity(n);
    let mut trail: Vec<Vec3> =
        if trail_on { Vec::with_capacity(n * MEASURED_STEPS) } else { Vec::new() };

    // Warmup
    for _ in 0..WARMUP_STEPS {
        vv_step_untimed(&mut engine, &mut bodies, &mut arrays, &mut acc, &mut trail, trail_on);
    }

    // Measured
    let mut row = RowV {
        n,
        trail_on,
        t_total_ns: 0,
        t_tree_build_ns: 0,
        t_bh_walk_ns: 0,
        t_integrator_ns: 0,
        t_trail_ns: 0,
        n_node_visits: 0,
        n_bh_accepted: 0,
        n_leaf_interactions: 0,
    };
    let total_start = Instant::now();
    for _ in 0..MEASURED_STEPS {
        vv_step_timed(
            &mut engine,
            &mut bodies,
            &mut arrays,
            &mut acc,
            &mut trail,
            trail_on,
            &mut row,
        );
    }
    row.t_total_ns = total_start.elapsed().as_nanos() as u64;
    row
}

/// Hand-rolled VV step (mirrors `physics::integrator::velocity_verlet`):
/// force(t) → kick(½dt) → drift(dt) → force(t+dt) → kick(½dt).
/// Each phase timed; counters merged from `evaluate_profile`. Trail push
/// timed inside the same step boundary so the trail variant correctly
/// includes all per-step trail cost.
///
/// The SoA snapshot is repacked before each force eval (twice per step,
/// matching `GravityForceModel::compute`'s per-call pack pattern in
/// production). Pack cost lands inside `t_tree_build_ns` so the engine-
/// ceiling profile keeps a single tree-prep accounting bucket.
#[inline(always)]
fn vv_step_timed(
    engine: &mut BarnesHutEngine,
    bodies: &mut [Body],
    arrays: &mut BodyArrays,
    acc: &mut [Vec3],
    trail: &mut Vec<Vec3>,
    trail_on: bool,
    row: &mut RowV,
) {
    // Force(t): pack + build + walk
    let t = Instant::now();
    arrays.pack_from(bodies);
    engine.build(arrays);
    row.t_tree_build_ns += t.elapsed().as_nanos() as u64;

    let t = Instant::now();
    let (_, c) = engine.evaluate_profile(arrays, THETA, acc);
    row.t_bh_walk_ns += t.elapsed().as_nanos() as u64;
    row.n_node_visits += c.n_node_visits;
    row.n_bh_accepted += c.n_bh_accepted;
    row.n_leaf_interactions += c.n_leaf_interactions;

    // Integrator: kick(½dt) + drift(dt)
    let t = Instant::now();
    kick(bodies, acc, 0.5 * DT);
    drift(bodies, DT);
    row.t_integrator_ns += t.elapsed().as_nanos() as u64;

    // Force(t+dt): pack + build + walk
    let t = Instant::now();
    arrays.pack_from(bodies);
    engine.build(arrays);
    row.t_tree_build_ns += t.elapsed().as_nanos() as u64;

    let t = Instant::now();
    let (_, c) = engine.evaluate_profile(arrays, THETA, acc);
    row.t_bh_walk_ns += t.elapsed().as_nanos() as u64;
    row.n_node_visits += c.n_node_visits;
    row.n_bh_accepted += c.n_bh_accepted;
    row.n_leaf_interactions += c.n_leaf_interactions;

    // Final kick(½dt)
    let t = Instant::now();
    kick(bodies, acc, 0.5 * DT);
    row.t_integrator_ns += t.elapsed().as_nanos() as u64;

    if trail_on {
        let t = Instant::now();
        for b in bodies.iter() {
            trail.push(Vec3::new(b.pos_x, b.pos_y, b.pos_z));
        }
        row.t_trail_ns += t.elapsed().as_nanos() as u64;
    }
}

/// Same VV step without timing — used for warmup.
fn vv_step_untimed(
    engine: &mut BarnesHutEngine,
    bodies: &mut [Body],
    arrays: &mut BodyArrays,
    acc: &mut [Vec3],
    trail: &mut Vec<Vec3>,
    trail_on: bool,
) {
    arrays.pack_from(bodies);
    engine.build(arrays);
    let (_, _) = engine.evaluate_profile(arrays, THETA, acc);
    kick(bodies, acc, 0.5 * DT);
    drift(bodies, DT);
    arrays.pack_from(bodies);
    engine.build(arrays);
    let (_, _) = engine.evaluate_profile(arrays, THETA, acc);
    kick(bodies, acc, 0.5 * DT);
    if trail_on {
        for b in bodies.iter() {
            trail.push(Vec3::new(b.pos_x, b.pos_y, b.pos_z));
        }
    }
}

#[inline(always)]
fn kick(bodies: &mut [Body], acc: &[Vec3], half_dt: f64) {
    for (b, a) in bodies.iter_mut().zip(acc.iter()) {
        b.vel_x += a.x * half_dt;
        b.vel_y += a.y * half_dt;
        b.vel_z += a.z * half_dt;
    }
}

#[inline(always)]
fn drift(bodies: &mut [Body], dt: f64) {
    for b in bodies.iter_mut() {
        b.pos_x += b.vel_x * dt;
        b.pos_y += b.vel_y * dt;
        b.pos_z += b.vel_z * dt;
    }
}

// ── Cell I measurement ─────────────────────────────────────────────────────── //

#[derive(Debug, Clone, Copy)]
struct RowI {
    n: usize,
    t_total_ns: u64,
    n_steps: usize,
    final_dt: f64,
    sim_time: f64,
}

fn measure_i(n: usize) -> RowI {
    let bodies_init = sphere_distribution_lognormal(n, SEED);
    let mut bodies = bodies_init.clone();

    let mut force = GravityForceModel::new(THETA, 16);
    let mut perturbations: Vec<Box<dyn crate::physics::integrator::PerturbationForce>> = Vec::new();
    let mut acc: Vec<Vec3> = vec![Vec3::ZERO; n];
    let mut ias = Ias15::new();

    // Warmup: smaller count for IAS15 since each step is expensive.
    let warmup = WARMUP_STEPS.min(if n >= 10_000 { 2 } else { 5 });
    {
        let mut ctx = IntegratorContext {
            force: &mut force,
            g_factor: 1.0,
            perturbations: &perturbations,
            deadline: None,
        };
        for _ in 0..warmup {
            let _ = ias.step(&mut bodies, &mut ctx, DT, &mut acc);
        }
    }

    // Measured: shrink horizon for large N to keep runtime reasonable.
    let measured = if n >= 10_000 {
        5
    } else if n >= 1_000 {
        20
    } else {
        MEASURED_STEPS
    };

    let mut sim_time = 0.0_f64;
    let mut final_dt = DT;

    let total_start = Instant::now();
    {
        let mut ctx = IntegratorContext {
            force: &mut force,
            g_factor: 1.0,
            perturbations: &mut perturbations,
            deadline: None,
        };
        for _ in 0..measured {
            let r = ias.step(&mut bodies, &mut ctx, DT, &mut acc);
            sim_time += r.consumed_dt;
            final_dt = r.consumed_dt;
        }
    }
    let t_total_ns = total_start.elapsed().as_nanos() as u64;

    RowI { n, t_total_ns, n_steps: measured, final_dt, sim_time }
}

// ── CSV / printing ─────────────────────────────────────────────────────────── //

fn write_header_v(writer: &mut fs::File) {
    writeln!(
        writer,
        "n,trail_on,t_total_ns,t_tree_build_ns,t_bh_walk_ns,t_integrator_ns,t_trail_ns,\
         n_node_visits,n_bh_accepted,n_leaf_interactions,n_steps"
    )
    .unwrap();
}

fn write_row_v(writer: &mut fs::File, r: &RowV) {
    writeln!(
        writer,
        "{},{},{},{},{},{},{},{},{},{},{}",
        r.n,
        r.trail_on,
        r.t_total_ns,
        r.t_tree_build_ns,
        r.t_bh_walk_ns,
        r.t_integrator_ns,
        r.t_trail_ns,
        r.n_node_visits,
        r.n_bh_accepted,
        r.n_leaf_interactions,
        MEASURED_STEPS,
    )
    .unwrap();
}

fn write_header_i(writer: &mut fs::File) {
    writeln!(writer, "n,t_total_ns,n_steps,final_dt,sim_time").unwrap();
}

fn write_row_i(writer: &mut fs::File, r: &RowI) {
    writeln!(
        writer,
        "{},{},{},{:.6e},{:.6e}",
        r.n, r.t_total_ns, r.n_steps, r.final_dt, r.sim_time
    )
    .unwrap();
}

fn print_row_v(r: &RowV) {
    let total_ms = r.t_total_ns as f64 * 1e-6;
    let per_step_ms = total_ms / MEASURED_STEPS as f64;
    let sps = 1000.0 / per_step_ms;
    let phase_sum_ns = r.t_tree_build_ns + r.t_bh_walk_ns + r.t_integrator_ns + r.t_trail_ns;
    let sanity_pct = (phase_sum_ns as f64 / r.t_total_ns as f64) * 100.0;
    let total_interactions = r.n_bh_accepted + r.n_leaf_interactions;
    let t_per_interaction_ns = if total_interactions > 0 {
        r.t_bh_walk_ns as f64 / total_interactions as f64
    } else {
        0.0
    };
    let t_per_body_us = r.t_total_ns as f64 / 1_000.0 / (r.n as f64 * MEASURED_STEPS as f64);
    let bh_acceptance_ratio =
        if r.n_node_visits > 0 { r.n_bh_accepted as f64 / r.n_node_visits as f64 } else { 0.0 };
    let n_interactions_per_body = total_interactions as f64 / (r.n as f64 * MEASURED_STEPS as f64);

    eprintln!(
        "[ceiling-V]   N={:>6} trail={:>5} step={:>8.3}ms SPS={:>7.2} \
         build/walk/integ/trail={:>5.1}/{:>5.1}/{:>5.1}/{:>5.1}% sanity={:>5.1}% \
         t/int={:>6.1}ns t/body={:>6.2}\u{00b5}s int/body={:>6.0} accept={:>4.2}",
        r.n,
        r.trail_on,
        per_step_ms,
        sps,
        100.0 * r.t_tree_build_ns as f64 / r.t_total_ns as f64,
        100.0 * r.t_bh_walk_ns as f64 / r.t_total_ns as f64,
        100.0 * r.t_integrator_ns as f64 / r.t_total_ns as f64,
        100.0 * r.t_trail_ns as f64 / r.t_total_ns as f64,
        sanity_pct,
        t_per_interaction_ns,
        t_per_body_us,
        n_interactions_per_body,
        bh_acceptance_ratio,
    );
}

fn print_row_i(r: &RowI) {
    let per_step_ms = (r.t_total_ns as f64 * 1e-6) / r.n_steps as f64;
    let sps = 1000.0 / per_step_ms;
    let sim_per_wall = r.sim_time / (r.t_total_ns as f64 * 1e-9);
    eprintln!(
        "[ceiling-I]   N={:>5} step={:>10.3}ms SPS={:>7.3} \
         dt={:.2e} sim/wall={:.3} (n_steps={})",
        r.n, per_step_ms, sps, r.final_dt, sim_per_wall, r.n_steps
    );
}

fn ceiling_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/engine-ceiling")
}

// ── Body distribution (matches perf 2×2 sphere_distribution_lognormal) ─────── //

fn sphere_distribution_lognormal(n: usize, seed: u64) -> Vec<Body> {
    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut next_u64 = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        state
    };
    let mut next_unit = || (next_u64() >> 11) as f64 / (1u64 << 53) as f64;

    let mut bodies = Vec::with_capacity(n);
    while bodies.len() < n {
        let x = 2.0 * next_unit() - 1.0;
        let y = 2.0 * next_unit() - 1.0;
        let z = 2.0 * next_unit() - 1.0;
        if x * x + y * y + z * z > 1.0 {
            continue;
        }
        let u1 = next_unit().max(1e-12);
        let u2 = next_unit();
        let normal = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        let mass = normal.exp();
        let mut b = Body::rocky(mass).at(x, y).with_velocity(0.0, 0.0);
        b.pos_z = z;
        bodies.push(b);
    }
    bodies
}
