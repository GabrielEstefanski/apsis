//! Perf 2×2 Pareto-frontier harness for the octree multipole experiment.
//!
//! Lab notebook: `docs/experiments/2026-05-08-octree-perf-2x2.md`.
//!
//! Opt-in via:
//!
//! ```text
//! cargo test --release -p apsis perf_2x2_pareto_frontier -- --ignored --nocapture
//! ```
//!
//! Writes one CSV per seed to `target/perf-2x2/octree_pareto_<seed>.csv` with
//! columns `cell,N,theta,seed,max_rel_err,t_build_ms,t_walk_ms,t_eval_ms,
//! std_err_t_eval_ms`. The reference O(N²) force is computed once per
//! (seed, N) for `N ≤ N_REFERENCE_MAX`; above that, `max_rel_err` is empty
//! and only wall times are reported.
//!
//! Scope: PR-perf-1 covers cells A (mono) and C (quad) with Morton off in
//! both. PR-perf-2 adds cells B and D (Morton on) and writes §Decision from
//! the combined CSVs.
//!
//! Total runtime is dominated by the N = 100 000 wall-time loop and the
//! N = 10 000 reference; expect ≈ 5 minutes on the lab notebook hardware.

#![allow(dead_code)] // experiment-only harness; lifecycle bound to the perf 2x2 PRs

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::domain::body::Body;
use crate::math::Vec3;

use super::engine::BarnesHutEngine;
use super::tree::MultipoleOrder;

// ── Frozen-variable matrix (notebook §Methodology) ────────────────────────── //

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const NS: [usize; 3] = [1_000, 10_000, 100_000];
const N_REFERENCE_MAX: usize = 10_000;
const THETAS: [f64; 4] = [0.3, 0.5, 0.7, 0.9];
const WARMUP: usize = 1;
const MEASURED: usize = 10;

#[derive(Debug, Clone, Copy)]
struct Cell {
    name: &'static str,
    multipole: MultipoleOrder,
}

const CELLS_PR1: [Cell; 2] = [
    Cell { name: "A", multipole: MultipoleOrder::Monopole },
    Cell { name: "C", multipole: MultipoleOrder::Quadrupole },
];

// ── Harness ──────────────────────────────────────────────────────────────── //

#[test]
#[ignore = "perf harness; opt-in via cargo test --release perf_2x2_pareto_frontier -- --ignored --nocapture"]
fn perf_2x2_pareto_frontier() {
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-2x2");
    fs::create_dir_all(&out_dir).expect("create perf-2x2 output dir");

    eprintln!("[perf_2x2] CSV output dir: {}", out_dir.display());
    eprintln!(
        "[perf_2x2] cells = {:?}, N = {:?}, theta = {:?}, seeds = {}, warmup = {WARMUP}, measured = {MEASURED}",
        CELLS_PR1.iter().map(|c| c.name).collect::<Vec<_>>(),
        NS,
        THETAS,
        SEEDS.len(),
    );

    let t_total = Instant::now();
    for &seed in &SEEDS {
        run_seed(seed, &out_dir);
    }
    eprintln!("[perf_2x2] total runtime: {:.1}s", t_total.elapsed().as_secs_f64());
}

fn run_seed(seed: u64, out_dir: &Path) {
    let csv_path = out_dir.join(format!("octree_pareto_{seed:#x}.csv"));
    let mut writer = fs::File::create(&csv_path).expect("create csv");
    writeln!(
        writer,
        "cell,N,theta,seed,max_rel_err,t_build_ms,t_walk_ms,t_eval_ms,std_err_t_eval_ms"
    )
    .unwrap();

    eprintln!("[perf_2x2] === seed = {seed:#x} ===");

    for &n in &NS {
        let bodies = sphere_distribution_lognormal(n, seed);

        let reference = if n <= N_REFERENCE_MAX {
            let t = Instant::now();
            let mut bh_exact = BarnesHutEngine::new(16);
            bh_exact.set_exact_threshold(usize::MAX);
            bh_exact.build(&bodies);
            let mut acc = vec![Vec3::ZERO; n];
            let _ = bh_exact.evaluate(&bodies, 0.5, &mut acc);
            eprintln!("[perf_2x2]   N={n:>6} reference O(N²) in {:.2}s", t.elapsed().as_secs_f64());
            Some(acc)
        } else {
            eprintln!("[perf_2x2]   N={n:>6} reference skipped (above N_REFERENCE_MAX)");
            None
        };

        for &cell in &CELLS_PR1 {
            for &theta in &THETAS {
                let row = measure_cell(&bodies, cell, theta, reference.as_deref(), seed, n);
                write_row(&mut writer, &row);
                eprintln!(
                    "[perf_2x2]   cell={} N={n:>6} theta={theta} err={:>11} t_build={:>8.3}ms t_walk={:>8.3}ms t_eval={:>8.3}ms (sigma={:.3}ms)",
                    cell.name,
                    row.max_rel_err.map(|e| format!("{e:.4e}")).unwrap_or_else(|| "—".into()),
                    row.t_build_ms,
                    row.t_walk_ms,
                    row.t_eval_ms,
                    row.std_err_t_eval_ms,
                );
            }
        }
    }

    eprintln!("[perf_2x2]   wrote {}", csv_path.display());
}

#[derive(Debug)]
struct CsvRow {
    cell: &'static str,
    n: usize,
    theta: f64,
    seed: u64,
    max_rel_err: Option<f64>,
    t_build_ms: f64,
    t_walk_ms: f64,
    t_eval_ms: f64,
    std_err_t_eval_ms: f64,
}

fn measure_cell(
    bodies: &[Body],
    cell: Cell,
    theta: f64,
    reference: Option<&[Vec3]>,
    seed: u64,
    n: usize,
) -> CsvRow {
    let mut bh = BarnesHutEngine::new(16);
    bh.set_multipole_order(cell.multipole);

    let max_rel_err = reference.map(|ref_acc| {
        bh.build(bodies);
        let mut acc = vec![Vec3::ZERO; bodies.len()];
        bh.evaluate(bodies, theta, &mut acc);
        body_max_rel_error(&acc, ref_acc)
    });

    let mut t_build = Vec::with_capacity(MEASURED);
    for _ in 0..WARMUP {
        bh.build(bodies);
    }
    for _ in 0..MEASURED {
        let t = Instant::now();
        bh.build(bodies);
        t_build.push(t.elapsed().as_secs_f64() * 1000.0);
    }

    bh.build(bodies);
    let mut acc = vec![Vec3::ZERO; bodies.len()];
    let mut t_walk = Vec::with_capacity(MEASURED);
    for _ in 0..WARMUP {
        bh.evaluate(bodies, theta, &mut acc);
    }
    for _ in 0..MEASURED {
        let t = Instant::now();
        bh.evaluate(bodies, theta, &mut acc);
        t_walk.push(t.elapsed().as_secs_f64() * 1000.0);
    }

    let mut t_eval = Vec::with_capacity(MEASURED);
    for _ in 0..WARMUP {
        bh.build(bodies);
        bh.evaluate(bodies, theta, &mut acc);
    }
    for _ in 0..MEASURED {
        let t = Instant::now();
        bh.build(bodies);
        bh.evaluate(bodies, theta, &mut acc);
        t_eval.push(t.elapsed().as_secs_f64() * 1000.0);
    }

    CsvRow {
        cell: cell.name,
        n,
        theta,
        seed,
        max_rel_err,
        t_build_ms: median(&mut t_build),
        t_walk_ms: median(&mut t_walk),
        t_eval_ms: median(&mut t_eval),
        std_err_t_eval_ms: std_err(&t_eval),
    }
}

fn write_row(writer: &mut fs::File, row: &CsvRow) {
    let err = row.max_rel_err.map(|e| format!("{e:.6e}")).unwrap_or_default();
    writeln!(
        writer,
        "{},{},{},{:#x},{},{:.6},{:.6},{:.6},{:.6}",
        row.cell,
        row.n,
        row.theta,
        row.seed,
        err,
        row.t_build_ms,
        row.t_walk_ms,
        row.t_eval_ms,
        row.std_err_t_eval_ms
    )
    .unwrap();
}

// ── Stats ─────────────────────────────────────────────────────────────────── //

fn median(xs: &mut [f64]) -> f64 {
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xs[xs.len() / 2]
}

/// Standard error of the mean: σ / √n.
fn std_err(xs: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let mean = xs.iter().sum::<f64>() / n;
    let var = xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (var / n).sqrt()
}

fn body_max_rel_error(acc: &[Vec3], reference: &[Vec3]) -> f64 {
    acc.iter().zip(reference).fold(0.0_f64, |peak, (a, r)| {
        let r_mag = r.length().max(1e-30);
        let err = (*a - *r).length() / r_mag;
        peak.max(err)
    })
}

// ── Body distribution (matches engine.rs Tier 1 helper) ───────────────────── //

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
