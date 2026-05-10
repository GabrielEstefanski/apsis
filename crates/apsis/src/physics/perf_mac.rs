//! MAC comparison harness — Cells M0 (Classical) vs M1 (Barnes 1990).
//!
//! Lab notebook: `docs/experiments/2026-05-09-octree-mac.md`.
//!
//! Opt-in:
//!
//! ```text
//! cargo test --release -p apsis perf_mac_m0_vs_m1 -- --ignored --nocapture
//! ```
//!
//! ## Matched-accuracy framing
//!
//! The notebook's Tier 2 metric is wall-time at *matched accuracy*, not at
//! matched θ. Barnes 1990 is strictly more conservative than Classical at
//! fixed θ (it descends into more nodes), so a fixed-θ comparison would
//! penalise M1 for being more accurate. Instead:
//!
//! 1. Run M0 at θ₀ = 0.5; record its per-body p95 force error vs exact
//!    O(N²) reference (call this `p95_target`).
//! 2. Bisect θ for M1 in `[0.1, 1.5]` until M1's p95 force error matches
//!    `p95_target` within `MATCH_TOL` (5 %). Call the matched value `θ₁`.
//! 3. At (M0, θ₀) and (M1, θ₁), measure median wall-time of (build +
//!    evaluate) over `MEASURED_RUNS` independent calls.
//! 4. Report `t_M1 / t_M0` — values < 1 mean M1 wins at matched accuracy.
//!
//! Build cost is included in the wall-time because M1 adds a third
//! aggregation pass (`δ_max`) and we want the total cycle cost the
//! integrator would actually pay.
//!
//! ## Frozen variables
//!
//! * Seeds: `0x6F637472`, `0x71756164`, `0x6D6F7274` (the perf 2×2 / engine
//!   ceiling canonical set; cross-experiment comparability)
//! * Body distribution: sphere log-normal mass (matches perf 2×2 / engine
//!   ceiling)
//! * θ for M0: `0.5` (the production default tested elsewhere)
//! * Bisection range for M1: `[0.1, 1.5]`
//! * Match tolerance: 5 %
//! * Warmup runs (discarded): 3 per (cell, seed)
//! * Measured runs: 5 per (cell, seed), median wall-time within seed
//! * N grid: `{1_000, 5_000, 10_000}`
//!
//! CSV output: `target/perf-mac/profile.csv`. Per-row schema in
//! [`write_header`]; one row per (n, seed, mac).

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::gravity::tree::MacKind;

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const THETA0: f64 = 0.5;
const MATCH_TOL: f64 = 0.05;
const BISECT_THETA_LO: f64 = 0.1;
const BISECT_THETA_HI: f64 = 1.5;
const BISECT_MAX_ITERS: usize = 32;
const BISECT_THETA_TOL: f64 = 1e-3;
const WARMUP_RUNS: usize = 3;
const MEASURED_RUNS: usize = 5;

// ── Cell M0 vs M1 ──────────────────────────────────────────────────────────── //

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release perf_mac_m0_vs_m1 -- --ignored --nocapture"]
fn perf_mac_m0_vs_m1() {
    let n_values = [1_000usize, 5_000, 10_000];

    let out_dir = perf_mac_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf output dir");
    let csv_path = out_dir.join("profile.csv");
    let mut writer = fs::File::create(&csv_path).expect("create profile.csv");
    write_header(&mut writer);

    eprintln!("[perf-mac] M0 (Classical) vs M1 (Barnes 1990), matched-accuracy, N={:?}", n_values);
    eprintln!(
        "[perf-mac]   protocol: M0 at θ=0.5 → bisect M1 θ for p95 match (±{:.0}%)",
        MATCH_TOL * 100.0
    );

    let t_total = Instant::now();
    for &n in &n_values {
        for &seed in &SEEDS {
            let bodies = sphere_distribution_lognormal(n, seed);
            let exact_acc = compute_exact_acc(&bodies);

            let m0 = measure_cell(&bodies, MacKind::Classical, THETA0, &exact_acc);
            let theta_m1 = bisect_matched_theta(&bodies, m0.p95_force_err, &exact_acc);
            let m1 = measure_cell(&bodies, MacKind::Barnes1990, theta_m1, &exact_acc);

            write_row(&mut writer, n, seed, &m0);
            write_row(&mut writer, n, seed, &m1);
            print_pair(n, seed, &m0, &m1);
        }
    }
    eprintln!("[perf-mac] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[perf-mac] wrote {}", csv_path.display());
}

// ── Per-cell measurement ───────────────────────────────────────────────────── //

#[derive(Debug, Clone, Copy)]
struct Row {
    mac: MacKind,
    theta: f64,
    p95_force_err: f64,
    t_median_ms: f64,
    t_min_ms: f64,
    t_max_ms: f64,
    t_build_median_ms: f64,
    t_walk_median_ms: f64,
    n_node_visits: u64,
    n_bh_accepted: u64,
    n_leaf_interactions: u64,
}

fn measure_cell(bodies: &[Body], mac: MacKind, theta: f64, exact_acc: &[Vec3]) -> Row {
    let p95 = p95_force_err(bodies, mac, theta, exact_acc);

    let mut engine = BarnesHutEngine::new(16);
    engine.set_exact_threshold(1);
    engine.set_mac_kind(mac);
    let mut acc = vec![Vec3::ZERO; bodies.len()];

    // Warmup
    let mut bodies_buf = bodies.to_vec();
    for _ in 0..WARMUP_RUNS {
        engine.build(&bodies_buf);
        let _ = engine.evaluate_profile(&bodies_buf, theta, &mut acc);
        // Restore positions to the canonical input each iter — evaluate
        // does not mutate bodies, but cloning keeps the contract loose
        // in case future kernels do.
        bodies_buf.copy_from_slice(bodies);
    }

    let mut total_samples_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut build_samples_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut walk_samples_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut last_counters = (0u64, 0u64, 0u64);
    for _ in 0..MEASURED_RUNS {
        let t_total = Instant::now();
        let t_build = Instant::now();
        engine.build(&bodies_buf);
        let build_ms = t_build.elapsed().as_secs_f64() * 1000.0;

        let t_walk = Instant::now();
        let (_, c) = engine.evaluate_profile(&bodies_buf, theta, &mut acc);
        let walk_ms = t_walk.elapsed().as_secs_f64() * 1000.0;

        total_samples_ms.push(t_total.elapsed().as_secs_f64() * 1000.0);
        build_samples_ms.push(build_ms);
        walk_samples_ms.push(walk_ms);
        last_counters = (c.n_node_visits, c.n_bh_accepted, c.n_leaf_interactions);
        bodies_buf.copy_from_slice(bodies);
    }
    total_samples_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    build_samples_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    walk_samples_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let median = total_samples_ms[total_samples_ms.len() / 2];
    let min = total_samples_ms[0];
    let max = total_samples_ms[total_samples_ms.len() - 1];

    Row {
        mac,
        theta,
        p95_force_err: p95,
        t_median_ms: median,
        t_min_ms: min,
        t_max_ms: max,
        t_build_median_ms: build_samples_ms[build_samples_ms.len() / 2],
        t_walk_median_ms: walk_samples_ms[walk_samples_ms.len() / 2],
        n_node_visits: last_counters.0,
        n_bh_accepted: last_counters.1,
        n_leaf_interactions: last_counters.2,
    }
}

/// Per-body p95 relative force error against the exact O(N²) reference.
/// Force evaluation under a fixed (mac, θ) is deterministic, so a single
/// call suffices — no median needed.
fn p95_force_err(bodies: &[Body], mac: MacKind, theta: f64, exact_acc: &[Vec3]) -> f64 {
    let mut engine = BarnesHutEngine::new(16);
    engine.set_exact_threshold(1);
    engine.set_mac_kind(mac);
    engine.build(bodies);
    let mut acc = vec![Vec3::ZERO; bodies.len()];
    engine.evaluate(bodies, theta, &mut acc);

    let mut errs: Vec<f64> = acc
        .iter()
        .zip(exact_acc.iter())
        .map(|(a, e)| {
            let diff = (a.x - e.x, a.y - e.y, a.z - e.z);
            let num = (diff.0 * diff.0 + diff.1 * diff.1 + diff.2 * diff.2).sqrt();
            let den = (e.x * e.x + e.y * e.y + e.z * e.z).sqrt().max(1e-300);
            num / den
        })
        .collect();
    errs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    errs[(errs.len() as f64 * 0.95) as usize]
}

/// Bisect θ for `MacKind::Barnes1990` until `p95_force_err` lands within
/// `MATCH_TOL` of `target`.
///
/// Monotonicity assumption: p95 force error is non-decreasing in θ for a
/// fixed body distribution and MAC. This holds for both M0 and M1 because
/// a larger θ admits more node acceptances without recursion (each
/// acceptance contributes at most as much error as the recursive
/// expansion would).
fn bisect_matched_theta(bodies: &[Body], target: f64, exact_acc: &[Vec3]) -> f64 {
    let mut lo = BISECT_THETA_LO;
    let mut hi = BISECT_THETA_HI;
    let mut last_mid = 0.5 * (lo + hi);

    for it in 0..BISECT_MAX_ITERS {
        let mid = 0.5 * (lo + hi);
        last_mid = mid;
        let p95 = p95_force_err(bodies, MacKind::Barnes1990, mid, exact_acc);
        let rel = (p95 - target).abs() / target.max(1e-300);
        eprintln!(
            "[perf-mac]   bisect[N={}] iter={:>2} θ={:.4} p95={:.3e} target={:.3e} (rel={:.3})",
            bodies.len(),
            it,
            mid,
            p95,
            target,
            rel
        );
        if rel <= MATCH_TOL {
            return mid;
        }
        if p95 < target {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo) < BISECT_THETA_TOL {
            break;
        }
    }
    eprintln!(
        "[perf-mac]   bisect did not converge within {} iters; using θ={:.4}",
        BISECT_MAX_ITERS, last_mid
    );
    last_mid
}

fn compute_exact_acc(bodies: &[Body]) -> Vec<Vec3> {
    let mut engine = BarnesHutEngine::new(16);
    engine.set_exact_threshold(usize::MAX);
    engine.build(bodies);
    let mut acc = vec![Vec3::ZERO; bodies.len()];
    engine.evaluate(bodies, THETA0, &mut acc);
    acc
}

// ── CSV / printing ─────────────────────────────────────────────────────────── //

fn write_header(writer: &mut fs::File) {
    writeln!(
        writer,
        "n,seed,mac,theta,p95_force_err,t_median_ms,t_min_ms,t_max_ms,\
         t_build_median_ms,t_walk_median_ms,\
         n_node_visits,n_bh_accepted,n_leaf_interactions,\
         n_runs,warmup_runs"
    )
    .unwrap();
}

fn write_row(writer: &mut fs::File, n: usize, seed: u64, r: &Row) {
    let mac_str = match r.mac {
        MacKind::Classical => "M0_classical",
        MacKind::Barnes1990 => "M1_barnes1990",
    };
    writeln!(
        writer,
        "{},0x{:X},{},{:.6},{:.6e},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{},{},{}",
        n,
        seed,
        mac_str,
        r.theta,
        r.p95_force_err,
        r.t_median_ms,
        r.t_min_ms,
        r.t_max_ms,
        r.t_build_median_ms,
        r.t_walk_median_ms,
        r.n_node_visits,
        r.n_bh_accepted,
        r.n_leaf_interactions,
        MEASURED_RUNS,
        WARMUP_RUNS
    )
    .unwrap();
}

fn print_pair(n: usize, seed: u64, m0: &Row, m1: &Row) {
    let t_ratio = m1.t_median_ms / m0.t_median_ms;
    let p95_ratio = m1.p95_force_err / m0.p95_force_err.max(1e-300);
    let int0 = m0.n_bh_accepted + m0.n_leaf_interactions;
    let int1 = m1.n_bh_accepted + m1.n_leaf_interactions;
    let int_ratio = int1 as f64 / int0.max(1) as f64;

    eprintln!(
        "[perf-mac] N={:>5} seed=0x{:X}  M0(θ={:.3}): p95={:.3e} t={:>7.3}ms \
         (build={:.3} walk={:.3}) interactions={}",
        n,
        seed,
        m0.theta,
        m0.p95_force_err,
        m0.t_median_ms,
        m0.t_build_median_ms,
        m0.t_walk_median_ms,
        int0
    );
    eprintln!(
        "[perf-mac] N={:>5} seed=0x{:X}  M1(θ={:.3}): p95={:.3e} t={:>7.3}ms \
         (build={:.3} walk={:.3}) interactions={}",
        n,
        seed,
        m1.theta,
        m1.p95_force_err,
        m1.t_median_ms,
        m1.t_build_median_ms,
        m1.t_walk_median_ms,
        int1
    );
    eprintln!(
        "[perf-mac] N={:>5} seed=0x{:X}  ratios: p95={:.3}  t={:.3}  interactions={:.3}",
        n, seed, p95_ratio, t_ratio, int_ratio
    );
}

fn perf_mac_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-mac")
}

// ── Body distribution (matches engine_ceiling.rs / perf 2×2) ───────────────── //

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
        b.z = z;
        bodies.push(b);
    }
    bodies
}
