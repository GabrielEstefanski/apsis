//! Tree-maintenance perf harness — Tier 2 (build wall-time reduction),
//! Tier 3 (walk regression bound), and Tier 4 (chaotic-system regression).
//!
//! Lab notebook: `docs/experiments/2026-05-12-tree-incremental-updates.md`.
//!
//! Opt-in:
//!
//! ```text
//! cargo test --release -p apsis perf_tree_maintenance -- --ignored --nocapture
//! ```
//!
//! Sweeps the `(system × N × seed)` grid declared a priori in the lab
//! notebook §Methodology — stable system (sphere log-normal, gentle
//! drift) plus chaotic system (Plummer cluster, virial velocities) ×
//! `N ∈ {1_000, 5_000, 10_000, 50_000}` × three perf-canonical seeds —
//! and reports per-cell median build/maintain wall-time over 20
//! measured runs after 5 warmup runs. The Tier 1 correctness gate
//! lives next to the kernel itself (`tree.rs`).

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use crate::domain::body::Body;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::gravity::BarnesHutEngine;

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const N_VALUES: [usize; 4] = [1_000, 5_000, 10_000, 50_000];
const THETA: f64 = 0.5;
const DT: f64 = 1.0e-3;
const WARMUP_RUNS: usize = 5;
const MEASURED_RUNS: usize = 20;

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release perf_tree_maintenance -- --ignored --nocapture"]
fn perf_tree_maintenance() {
    let out_dir = perf_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf-tree output dir");
    let csv_path = out_dir.join("maintenance.csv");
    let mut writer = fs::File::create(&csv_path).expect("create maintenance.csv");
    writeln!(
        writer,
        "system,n,seed,t_build_ns_median,t_maintain_ns_median,t_walk_build_ns_median,\
         t_walk_maintain_ns_median,migrants_per_step_median"
    )
    .unwrap();

    eprintln!(
        "[perf-tree] Tier 2 + Tier 3 + Tier 4 — systems=[stable, chaotic] N={N_VALUES:?} seeds={SEEDS:#X?}"
    );
    let t_total = Instant::now();

    let mut rows: Vec<Row> = Vec::new();
    for &system in &[SystemKind::Stable, SystemKind::Chaotic] {
        for &n in &N_VALUES {
            for &seed in &SEEDS {
                let row = measure(system, n, seed);
                writeln!(
                    writer,
                    "{:?},{},0x{:X},{},{},{},{},{}",
                    row.system,
                    row.n,
                    row.seed,
                    row.t_build_ns_median,
                    row.t_maintain_ns_median,
                    row.t_walk_build_ns_median,
                    row.t_walk_maintain_ns_median,
                    row.migrants_per_step_median,
                )
                .unwrap();
                eprintln!(
                    "[perf-tree]  {:>8?}  N={:>5}  seed=0x{:X}  \
                     build={:>7.2}\u{00B5}s  maint={:>7.2}\u{00B5}s  ratio={:>5.2}\u{00D7}  \
                     walk_b={:>7.2}ms  walk_m={:>7.2}ms  migrants={:>5}",
                    row.system,
                    row.n,
                    row.seed,
                    row.t_build_ns_median as f64 * 1e-3,
                    row.t_maintain_ns_median as f64 * 1e-3,
                    row.t_maintain_ns_median as f64 / row.t_build_ns_median.max(1) as f64,
                    row.t_walk_build_ns_median as f64 * 1e-6,
                    row.t_walk_maintain_ns_median as f64 * 1e-6,
                    row.migrants_per_step_median,
                );
                rows.push(row);
            }
        }
    }

    print_summary_tables(&rows);

    eprintln!("[perf-tree] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[perf-tree] wrote {}", csv_path.display());
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SystemKind {
    Stable,
    Chaotic,
}

#[derive(Debug, Clone, Copy)]
struct Row {
    system: SystemKind,
    n: usize,
    seed: u64,
    t_build_ns_median: u64,
    t_maintain_ns_median: u64,
    t_walk_build_ns_median: u64,
    t_walk_maintain_ns_median: u64,
    migrants_per_step_median: u64,
}

fn measure(system: SystemKind, n: usize, seed: u64) -> Row {
    let bodies_init = match system {
        SystemKind::Stable => sphere_lognormal_drift(n, seed, 0.3),
        SystemKind::Chaotic => plummer_cluster(n, seed),
    };

    // Two parallel engines: one always rebuilds, one maintains.
    let mut bodies_b = bodies_init.clone();
    let mut bodies_m = bodies_init.clone();
    let mut arrays_b = BodyArrays::with_capacity(n);
    let mut arrays_m = BodyArrays::with_capacity(n);
    let mut engine_b = BarnesHutEngine::new(16);
    let mut engine_m = BarnesHutEngine::new(16);

    arrays_b.pack_from(&bodies_b);
    arrays_m.pack_from(&bodies_m);
    engine_b.build(&arrays_b);
    engine_m.build(&arrays_m);

    let mut acc_b = vec![Vec3::ZERO; n];
    let mut acc_m = vec![Vec3::ZERO; n];

    for _ in 0..WARMUP_RUNS {
        step_drift(&mut bodies_b);
        step_drift(&mut bodies_m);
        arrays_b.pack_from(&bodies_b);
        arrays_m.pack_from(&bodies_m);
        engine_b.build(&arrays_b);
        engine_m.maintain(&arrays_m);
        engine_b.evaluate(&arrays_b, THETA, &mut acc_b);
        engine_m.evaluate(&arrays_m, THETA, &mut acc_m);
    }

    let mut samples_build_ns = Vec::with_capacity(MEASURED_RUNS);
    let mut samples_maint_ns = Vec::with_capacity(MEASURED_RUNS);
    let mut samples_walk_build_ns = Vec::with_capacity(MEASURED_RUNS);
    let mut samples_walk_maint_ns = Vec::with_capacity(MEASURED_RUNS);
    let mut samples_migrants = Vec::with_capacity(MEASURED_RUNS);

    for _ in 0..MEASURED_RUNS {
        step_drift(&mut bodies_b);
        step_drift(&mut bodies_m);
        arrays_b.pack_from(&bodies_b);
        arrays_m.pack_from(&bodies_m);

        let t = Instant::now();
        engine_b.build(&arrays_b);
        samples_build_ns.push(t.elapsed().as_nanos() as u64);

        let cell_idx_before = engine_m.tree_cell_idx_snapshot();
        let t = Instant::now();
        engine_m.maintain(&arrays_m);
        samples_maint_ns.push(t.elapsed().as_nanos() as u64);
        let migrants = engine_m
            .tree_cell_idx_snapshot()
            .iter()
            .zip(cell_idx_before.iter())
            .filter(|(a, b)| a != b)
            .count() as u64;
        samples_migrants.push(migrants);

        let t = Instant::now();
        engine_b.evaluate(&arrays_b, THETA, &mut acc_b);
        samples_walk_build_ns.push(t.elapsed().as_nanos() as u64);

        let t = Instant::now();
        engine_m.evaluate(&arrays_m, THETA, &mut acc_m);
        samples_walk_maint_ns.push(t.elapsed().as_nanos() as u64);
    }

    samples_build_ns.sort_unstable();
    samples_maint_ns.sort_unstable();
    samples_walk_build_ns.sort_unstable();
    samples_walk_maint_ns.sort_unstable();
    samples_migrants.sort_unstable();

    Row {
        system,
        n,
        seed,
        t_build_ns_median: samples_build_ns[MEASURED_RUNS / 2],
        t_maintain_ns_median: samples_maint_ns[MEASURED_RUNS / 2],
        t_walk_build_ns_median: samples_walk_build_ns[MEASURED_RUNS / 2],
        t_walk_maintain_ns_median: samples_walk_maint_ns[MEASURED_RUNS / 2],
        migrants_per_step_median: samples_migrants[MEASURED_RUNS / 2],
    }
}

fn step_drift(bodies: &mut [Body]) {
    for b in bodies.iter_mut() {
        b.pos_x += DT * b.vel_x;
        b.pos_y += DT * b.vel_y;
        b.pos_z += DT * b.vel_z;
    }
}

fn print_summary_tables(rows: &[Row]) {
    eprintln!();
    eprintln!("[perf-tree] ── Tier 2: build/maintain wall-time ratio (per cell) ──");
    eprintln!(
        "[perf-tree] {:>8}  {:>5}  {:>16}  {:>16}  {:>16}",
        "system", "N", "seed=0x6F637472", "seed=0x71756164", "seed=0x6D6F7274"
    );
    for &system in &[SystemKind::Stable, SystemKind::Chaotic] {
        for &n in &N_VALUES {
            let ratios: Vec<String> = SEEDS
                .iter()
                .map(|&seed| {
                    let r = rows.iter().find(|r| r.system == system && r.n == n && r.seed == seed);
                    match r {
                        Some(r) if r.t_build_ns_median > 0 => format!(
                            "{:>14.3}\u{00D7}",
                            r.t_maintain_ns_median as f64 / r.t_build_ns_median as f64
                        ),
                        _ => "n/a".to_string(),
                    }
                })
                .collect();
            eprintln!(
                "[perf-tree] {:>8?}  {:>5}  {:>16}  {:>16}  {:>16}",
                system, n, ratios[0], ratios[1], ratios[2]
            );
        }
    }

    eprintln!();
    eprintln!("[perf-tree] ── Tier 3: walk regression (maintain vs build) ──");
    eprintln!(
        "[perf-tree] {:>8}  {:>5}  {:>16}  {:>16}  {:>16}  envelope",
        "system", "N", "seed=0x6F637472", "seed=0x71756164", "seed=0x6D6F7274"
    );
    for &system in &[SystemKind::Stable, SystemKind::Chaotic] {
        for &n in &N_VALUES {
            let ratios: Vec<String> = SEEDS
                .iter()
                .map(|&seed| {
                    let r = rows.iter().find(|r| r.system == system && r.n == n && r.seed == seed);
                    match r {
                        Some(r) if r.t_walk_build_ns_median > 0 => format!(
                            "{:>14.3}\u{00D7}",
                            r.t_walk_maintain_ns_median as f64 / r.t_walk_build_ns_median as f64
                        ),
                        _ => "n/a".to_string(),
                    }
                })
                .collect();
            eprintln!(
                "[perf-tree] {:>8?}  {:>5}  {:>16}  {:>16}  {:>16}  [0.95, 1.05]",
                system, n, ratios[0], ratios[1], ratios[2]
            );
        }
    }

    eprintln!();
    eprintln!("[perf-tree] ── Migrant rate (median bodies/step changing cell) ──");
    eprintln!(
        "[perf-tree] {:>8}  {:>5}  {:>16}  {:>16}  {:>16}",
        "system", "N", "seed=0x6F637472", "seed=0x71756164", "seed=0x6D6F7274"
    );
    for &system in &[SystemKind::Stable, SystemKind::Chaotic] {
        for &n in &N_VALUES {
            let mig: Vec<String> = SEEDS
                .iter()
                .map(|&seed| {
                    let r = rows.iter().find(|r| r.system == system && r.n == n && r.seed == seed);
                    match r {
                        Some(r) => {
                            let pct = 100.0 * r.migrants_per_step_median as f64 / r.n as f64;
                            format!("{:>5} ({:>5.2}%)", r.migrants_per_step_median, pct)
                        },
                        None => "n/a".to_string(),
                    }
                })
                .collect();
            eprintln!(
                "[perf-tree] {:>8?}  {:>5}  {:>16}  {:>16}  {:>16}",
                system, n, mig[0], mig[1], mig[2]
            );
        }
    }
}

// ── Distributions ─────────────────────────────────────────────────────────── //

/// Sphere log-normal distribution with caller-controlled velocity scale.
/// `vel_scale` controls per-body random velocity in each axis.
fn sphere_lognormal_drift(n: usize, seed: u64, vel_scale: f64) -> Vec<Body> {
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
        b.vel_x = vel_scale * (2.0 * next_unit() - 1.0);
        b.vel_y = vel_scale * (2.0 * next_unit() - 1.0);
        b.vel_z = vel_scale * (2.0 * next_unit() - 1.0);
        bodies.push(b);
    }
    bodies
}

/// Plummer cluster with isotropic Maxwellian velocities at virial scale.
/// Produces a chaotic system with frequent close encounters and high
/// migration rate.
fn plummer_cluster(n: usize, seed: u64) -> Vec<Body> {
    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut next_u64 = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        state
    };
    let mut next_unit = || (next_u64() >> 11) as f64 / (1u64 << 53) as f64;

    let scale_radius = 1.0_f64;
    let v_scale = 1.5_f64;

    let mut bodies = Vec::with_capacity(n);
    while bodies.len() < n {
        let u = next_unit().clamp(1e-9, 1.0 - 1e-9);
        let r = scale_radius * (u.powf(-2.0 / 3.0) - 1.0).powf(-0.5);
        if !r.is_finite() || r > 20.0 * scale_radius {
            continue;
        }
        let cos_theta = 2.0 * next_unit() - 1.0;
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        let phi = 2.0 * std::f64::consts::PI * next_unit();
        let x = r * sin_theta * phi.cos();
        let y = r * sin_theta * phi.sin();
        let z = r * cos_theta;

        let u1 = next_unit().max(1e-12);
        let u2 = next_unit();
        let nx = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        let u3 = next_unit().max(1e-12);
        let u4 = next_unit();
        let ny = (-2.0 * u3.ln()).sqrt() * (2.0 * std::f64::consts::PI * u4).cos();
        let u5 = next_unit().max(1e-12);
        let u6 = next_unit();
        let nz = (-2.0 * u5.ln()).sqrt() * (2.0 * std::f64::consts::PI * u6).cos();

        let mass = 1.0_f64;
        let mut b = Body::rocky(mass).at(x, y).with_velocity(0.0, 0.0);
        b.pos_z = z;
        b.vel_x = v_scale * nx;
        b.vel_y = v_scale * ny;
        b.vel_z = v_scale * nz;
        bodies.push(b);
    }
    bodies
}

fn perf_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-tree")
}
