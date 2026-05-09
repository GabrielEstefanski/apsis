//! Perf 2×2 Pareto-frontier harness for the octree multipole experiment.
//!
//! Lab notebook: `docs/experiments/2026-05-08-octree-perf-2x2.md`.
//!
//! Three opt-in tests live here:
//!
//! ```text
//! # Full Pareto sweep (3 seeds × 3 N × 4 θ × 4 cells), writes CSVs:
//! cargo test --release -p apsis perf_2x2_pareto_frontier -- --ignored --nocapture
//!
//! # Force-accuracy gates per notebook §Tier 1 (asserts on p50 / p95):
//! cargo test --release -p apsis tier1_perf_2x2_force_accuracy_gates -- --ignored --nocapture
//!
//! # Leaf-capacity sensitivity sweep (LEAF ∈ {4, 8, 16, 32} on cell A):
//! cargo test --release -p apsis leaf_capacity_sensitivity -- --ignored --nocapture
//! ```
//!
//! Per-body force error is measured against:
//! * an exact O(N²) reference for `N ≤ N_REFERENCE_FULL_MAX`, or
//! * a sampled reference (`K_SAMPLE = 512` randomly chosen bodies, exact
//!   pairwise force on each) for larger `N`.
//!
//! The error array is summarised as `{p50, p95, p99, max, mean, std}` and
//! exported per (cell, N, θ, seed). The notebook gates only `p50` and `p95`
//! (Salmon-Warren 1994 typical-body and tail bounds); `p99` and `max` are
//! recorded as informational signals — long-tail sensitivity to pathological
//! body configurations is structural to the BH approximation, not an
//! implementation defect.
//!
//! Scope: PR-perf-2 covers all four cells of the factorial:
//!
//! | | Morton off | Morton on |
//! | --- | --- | --- |
//! | **Monopole** | A | B |
//! | **Quadrupole** | C | D |
//!
//! Each cell is measured at the full grid (3 seeds × 3 N × 4 θ).

#![allow(dead_code)] // experiment-only harness; lifecycle bound to the perf 2x2 PRs

use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;

use crate::domain::body::Body;
use crate::math::Vec3;

use super::engine::{BarnesHutEngine, bh_eval_body};
use super::kernel::{G, Kernel, PlummerKernel, pair_eps2};
use super::tree::{MultipoleOrder, Octree};

// ── Frozen-variable matrix (notebook §Methodology) ────────────────────────── //

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const NS: [usize; 3] = [1_000, 10_000, 100_000];
const N_REFERENCE_FULL_MAX: usize = 10_000;
/// Sampled-reference size for `N > N_REFERENCE_FULL_MAX`. K = 512 puts ~25
/// samples in the p95 tail and ~5 in the p99 tail (vs. 13 / 2.5 at K = 256).
/// p99 is still informational under sampling — the notebook documents the
/// confidence bound — but K = 512 is the cheapest size where p95 estimation
/// is statistically solid (SE ≈ 1 %).
const K_SAMPLE: usize = 512;
const THETAS: [f64; 4] = [0.3, 0.5, 0.7, 0.9];
const WARMUP: usize = 1;
const MEASURED: usize = 10;

#[derive(Debug, Clone, Copy)]
struct Cell {
    name: &'static str,
    multipole: MultipoleOrder,
    morton: bool,
}

const CELLS: [Cell; 4] = [
    Cell { name: "A", multipole: MultipoleOrder::Monopole, morton: false },
    Cell { name: "B", multipole: MultipoleOrder::Monopole, morton: true },
    Cell { name: "C", multipole: MultipoleOrder::Quadrupole, morton: false },
    Cell { name: "D", multipole: MultipoleOrder::Quadrupole, morton: true },
];

// ── Pareto-frontier harness ──────────────────────────────────────────────── //

#[test]
#[ignore = "perf harness; opt-in via cargo test --release perf_2x2_pareto_frontier -- --ignored --nocapture"]
fn perf_2x2_pareto_frontier() {
    let out_dir = perf_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf-2x2 output dir");

    eprintln!("[perf_2x2] CSV output dir: {}", out_dir.display());
    eprintln!(
        "[perf_2x2] cells = {:?}, N = {:?}, theta = {:?}, seeds = {}, K_sample = {K_SAMPLE}",
        CELLS.iter().map(|c| c.name).collect::<Vec<_>>(),
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
        "cell,N,theta,seed,p50,p95,p99,max,mean,std,t_build_ms,t_walk_ms,t_eval_ms,std_err_t_eval_ms"
    )
    .unwrap();

    eprintln!("[perf_2x2] === seed = {seed:#x} ===");

    for &n in &NS {
        let bodies = sphere_distribution_lognormal(n, seed);
        let reference = build_reference(&bodies, seed);

        for &cell in &CELLS {
            for &theta in &THETAS {
                let row = measure_cell(&bodies, cell, theta, &reference, seed, n);
                write_row(&mut writer, &row);
                let e = &row.error;
                eprintln!(
                    "[perf_2x2]   cell={} N={n:>6} theta={theta} \
                     p50={:.2e} p95={:.2e} p99={:.2e} max={:.2e} \
                     t_build={:>7.2}ms t_walk={:>7.2}ms t_eval={:>7.2}ms",
                    cell.name,
                    e.p50,
                    e.p95,
                    e.p99,
                    e.max,
                    row.t_build_ms,
                    row.t_walk_ms,
                    row.t_eval_ms,
                );
            }
        }
    }

    eprintln!("[perf_2x2]   wrote {}", csv_path.display());
}

// ── Tier 1 force-accuracy gates ──────────────────────────────────────────── //

#[test]
#[ignore = "perf gate; opt-in via cargo test --release tier1_perf_2x2_force_accuracy_gates -- --ignored --nocapture"]
fn tier1_perf_2x2_force_accuracy_gates() {
    // Notebook §Tier 1 bounds at θ = 0.5. Reading: Salmon-Warren 1994's
    // "5 % per body at θ = 0.5" is the 95th-percentile bound of the per-body
    // error distribution, not the max. Hernquist & Katz 1989 reports the
    // same distribution shape with ≈ 10× scale reduction under quadrupole.
    // p99 and max are recorded but not gated — long-tail sensitivity to
    // pathological body configurations (boundary cells, anisotropy near the
    // opening criterion) is structural to BH and not an implementation
    // defect; gating max would conflate algorithmic behaviour with
    // implementation correctness.
    //
    // Bounds:
    //   Cell A (mono): p50 ≤ 1e-2 (typical), p95 ≤ 5e-2 (tail)
    //   Cell C (quad): p50 ≤ 1e-3 (10× mono), p95 ≤ 5e-3 (10× mono)
    let theta = 0.5;
    let bounds = [
        (Cell { name: "A", multipole: MultipoleOrder::Monopole, morton: false }, 1.0e-2, 5.0e-2),
        (Cell { name: "B", multipole: MultipoleOrder::Monopole, morton: true }, 1.0e-2, 5.0e-2),
        (Cell { name: "C", multipole: MultipoleOrder::Quadrupole, morton: false }, 1.0e-3, 5.0e-3),
        (Cell { name: "D", multipole: MultipoleOrder::Quadrupole, morton: true }, 1.0e-3, 5.0e-3),
    ];

    let mut violations: Vec<String> = Vec::new();

    for &n in &[1_000, 10_000] {
        for &seed in &SEEDS {
            let bodies = sphere_distribution_lognormal(n, seed);
            let reference = build_reference(&bodies, seed);

            for &(cell, p50_bound, p95_bound) in &bounds {
                let stats = measure_error_only(&bodies, cell, theta, &reference);

                eprintln!(
                    "[tier1] cell={} N={n:>5} theta={theta} seed={seed:#x} \
                     p50={:.3e} p95={:.3e} p99={:.3e} max={:.3e}",
                    cell.name, stats.p50, stats.p95, stats.p99, stats.max,
                );

                if stats.p50 > p50_bound {
                    violations.push(format!(
                        "  cell={} N={n} theta={theta} seed={seed:#x}: p50={:.3e} > bound {:.3e}",
                        cell.name, stats.p50, p50_bound,
                    ));
                }
                if stats.p95 > p95_bound {
                    violations.push(format!(
                        "  cell={} N={n} theta={theta} seed={seed:#x}: p95={:.3e} > bound {:.3e}",
                        cell.name, stats.p95, p95_bound,
                    ));
                }
            }
        }
    }

    assert!(violations.is_empty(), "[tier1] bound violations:\n{}", violations.join("\n"),);
}

// ── Leaf-capacity sensitivity sweep ──────────────────────────────────────── //

/// Sweep `LEAF` ∈ {4, 8, 16, 32} on cell A configuration (monopole, no
/// Morton) at a single (N, θ, seed) operating point, exporting per-LEAF
/// build / walk / accuracy stats to `target/perf-2x2/leaf_sensitivity.csv`.
///
/// Cell A is chosen as the cleanest baseline — both factorial axes
/// disabled, so any change observed is attributable to LEAF alone. Single
/// (seed, N, θ) keeps runtime under a minute and matches the way the
/// notebook §Threats item documents the LEAF dependency: a one-shot
/// characterisation, not a full grid.
///
/// Goes through the generic `bh_eval_body` directly rather than
/// `BarnesHutEngine` because the engine is pinned to `LEAF = DEFAULT_LEAF`
/// in production. The walk logic mirrors `BarnesHutEngine::evaluate`'s
/// natural-order branch (no Morton perm).
#[test]
#[ignore = "perf characterisation; opt-in via cargo test --release leaf_capacity_sensitivity -- --ignored --nocapture"]
fn leaf_capacity_sensitivity() {
    let n = 10_000;
    let theta = 0.5;
    let seed = SEEDS[0];
    let bodies = sphere_distribution_lognormal(n, seed);
    let reference = build_reference(&bodies, seed);

    let rows = [
        measure_leaf::<4>(&bodies, &reference, theta),
        measure_leaf::<8>(&bodies, &reference, theta),
        measure_leaf::<16>(&bodies, &reference, theta),
        measure_leaf::<32>(&bodies, &reference, theta),
    ];

    let out_dir = perf_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf-2x2 output dir");
    let csv_path = out_dir.join("leaf_sensitivity.csv");
    let mut writer = fs::File::create(&csv_path).expect("create csv");
    writeln!(writer, "leaf,N,theta,seed,p50,p95,p99,max,t_build_ms,t_walk_ms").unwrap();

    eprintln!("[leaf-sens] N={n} theta={theta} seed={seed:#x} cell=A (mono, no Morton)");
    for row in &rows {
        writeln!(
            writer,
            "{},{},{},{:#x},{:.6e},{:.6e},{:.6e},{:.6e},{:.6},{:.6}",
            row.leaf,
            n,
            theta,
            seed,
            row.error.p50,
            row.error.p95,
            row.error.p99,
            row.error.max,
            row.t_build_ms,
            row.t_walk_ms,
        )
        .unwrap();
        eprintln!(
            "[leaf-sens]   LEAF={:>2}  p50={:.2e}  p95={:.2e}  max={:.2e}  \
             t_build={:>7.3}ms  t_walk={:>7.3}ms",
            row.leaf, row.error.p50, row.error.p95, row.error.max, row.t_build_ms, row.t_walk_ms,
        );
    }
    eprintln!("[leaf-sens]   wrote {}", csv_path.display());
}

#[derive(Debug)]
struct LeafRow {
    leaf: usize,
    error: ErrorStats,
    t_build_ms: f64,
    t_walk_ms: f64,
}

fn measure_leaf<const LEAF: usize>(bodies: &[Body], reference: &Reference, theta: f64) -> LeafRow {
    let kernel = PlummerKernel::new();
    let mut tree: Octree<LEAF> = Octree::new(16);

    // Build phase timing.
    let mut t_build = Vec::with_capacity(MEASURED);
    for _ in 0..WARMUP {
        tree.build(bodies, MultipoleOrder::Monopole, false);
    }
    for _ in 0..MEASURED {
        let t = Instant::now();
        tree.build(bodies, MultipoleOrder::Monopole, false);
        t_build.push(t.elapsed().as_secs_f64() * 1000.0);
    }

    // Final build for the walk + error measurement.
    tree.build(bodies, MultipoleOrder::Monopole, false);
    let nodes = tree.nodes();

    let walk_once = |acc: &mut [Vec3]| {
        let results: Vec<(usize, Vec3, f64)> = (0..bodies.len())
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                let (a, phi) = bh_eval_body(
                    nodes,
                    i,
                    &bodies[i],
                    bodies,
                    theta,
                    &kernel,
                    MultipoleOrder::Monopole,
                    &mut stack,
                );
                (i, a, phi)
            })
            .collect();
        for (i, a, _) in results {
            acc[i] = a;
        }
    };

    let mut acc = vec![Vec3::ZERO; bodies.len()];
    let mut t_walk = Vec::with_capacity(MEASURED);
    for _ in 0..WARMUP {
        walk_once(&mut acc);
    }
    for _ in 0..MEASURED {
        let t = Instant::now();
        walk_once(&mut acc);
        t_walk.push(t.elapsed().as_secs_f64() * 1000.0);
    }

    let errors: Vec<f64> = match reference {
        Reference::Full { forces } => {
            acc.iter().zip(forces).map(|(a, r)| relative_force_error(*a, *r)).collect()
        },
        Reference::Sampled { indices, forces } => {
            indices.iter().zip(forces).map(|(&idx, r)| relative_force_error(acc[idx], *r)).collect()
        },
    };

    LeafRow {
        leaf: LEAF,
        error: distribution_stats(&errors),
        t_build_ms: median(&mut t_build),
        t_walk_ms: median(&mut t_walk),
    }
}

// ── Reference computation ────────────────────────────────────────────────── //

enum Reference {
    /// `forces[i]` is the exact acceleration on `bodies[i]`. Used for `N ≤
    /// N_REFERENCE_FULL_MAX`.
    Full { forces: Vec<Vec3> },
    /// `forces[k]` is the exact acceleration on `bodies[indices[k]]`.
    /// `indices` is a deterministic random sample of size `K_SAMPLE`.
    Sampled { indices: Vec<usize>, forces: Vec<Vec3> },
}

fn build_reference(bodies: &[Body], seed: u64) -> Reference {
    let n = bodies.len();
    if n <= N_REFERENCE_FULL_MAX {
        let all_indices: Vec<usize> = (0..n).collect();
        let t = Instant::now();
        let forces = exact_pairwise_forces(bodies, &all_indices);
        eprintln!(
            "[perf_2x2]   N={n:>6} reference (full O(N^2), independent kernel) in {:.2}s",
            t.elapsed().as_secs_f64()
        );
        Reference::Full { forces }
    } else {
        let indices = sample_body_indices(n, K_SAMPLE, seed);
        let t = Instant::now();
        let forces = exact_pairwise_forces(bodies, &indices);
        eprintln!(
            "[perf_2x2]   N={n:>6} reference (sampled K={K_SAMPLE}) in {:.2}s",
            t.elapsed().as_secs_f64()
        );
        Reference::Sampled { indices, forces }
    }
}

fn sample_body_indices(n: usize, k: usize, seed: u64) -> Vec<usize> {
    let target = k.min(n);
    let mut state = seed.wrapping_add(0xCAFE_BABE_F00D_BEEFu64);
    let mut next_u64 = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        state
    };

    let mut chosen: HashSet<usize> = HashSet::with_capacity(target);
    while chosen.len() < target {
        chosen.insert((next_u64() as usize) % n);
    }
    let mut v: Vec<usize> = chosen.into_iter().collect();
    v.sort_unstable();
    v
}

/// Exact O(N) per-target pairwise forces, parallelised over `indices`.
///
/// Hand-rolled — does not route through [`BarnesHutEngine`]. Sharing only
/// the [`PlummerKernel`] and `pair_eps2` primitives means a defect in the
/// engine's exact-mode accumulation cannot mask a BH defect by silently
/// agreeing with the reference. Used by both Full mode (`indices = 0..N`)
/// and Sampled mode (`indices` = a deterministic random subset).
fn exact_pairwise_forces(bodies: &[Body], indices: &[usize]) -> Vec<Vec3> {
    let kernel = PlummerKernel::new();
    indices
        .par_iter()
        .map(|&i| {
            let body_i = &bodies[i];
            let mut a = Vec3::ZERO;
            for (j, body_j) in bodies.iter().enumerate() {
                if j == i {
                    continue;
                }
                let dx = body_j.x - body_i.x;
                let dy = body_j.y - body_i.y;
                let dz = body_j.z - body_i.z;
                let eps2 = pair_eps2(body_i.softening, body_j.softening);
                let r_sq = dx * dx + dy * dy + dz * dz;
                let fac = G * body_j.mass * kernel.acceleration_factor(r_sq, eps2);
                a.x += dx * fac;
                a.y += dy * fac;
                a.z += dz * fac;
            }
            a
        })
        .collect()
}

// ── Per-cell measurement ─────────────────────────────────────────────────── //

#[derive(Debug, Clone, Copy)]
struct ErrorStats {
    p50: f64,
    p95: f64,
    p99: f64,
    max: f64,
    mean: f64,
    std: f64,
}

#[derive(Debug)]
struct CsvRow {
    cell: &'static str,
    n: usize,
    theta: f64,
    seed: u64,
    error: ErrorStats,
    t_build_ms: f64,
    t_walk_ms: f64,
    t_eval_ms: f64,
    std_err_t_eval_ms: f64,
}

fn measure_cell(
    bodies: &[Body],
    cell: Cell,
    theta: f64,
    reference: &Reference,
    seed: u64,
    n: usize,
) -> CsvRow {
    let mut bh = BarnesHutEngine::new(16);
    bh.set_multipole_order(cell.multipole);
    bh.set_morton_enabled(cell.morton);

    let error = error_stats(&mut bh, bodies, theta, reference);

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
        error,
        t_build_ms: median(&mut t_build),
        t_walk_ms: median(&mut t_walk),
        t_eval_ms: median(&mut t_eval),
        std_err_t_eval_ms: std_err(&t_eval),
    }
}

fn measure_error_only(
    bodies: &[Body],
    cell: Cell,
    theta: f64,
    reference: &Reference,
) -> ErrorStats {
    let mut bh = BarnesHutEngine::new(16);
    bh.set_multipole_order(cell.multipole);
    bh.set_morton_enabled(cell.morton);
    error_stats(&mut bh, bodies, theta, reference)
}

fn error_stats(
    bh: &mut BarnesHutEngine,
    bodies: &[Body],
    theta: f64,
    reference: &Reference,
) -> ErrorStats {
    bh.build(bodies);
    let mut acc = vec![Vec3::ZERO; bodies.len()];
    bh.evaluate(bodies, theta, &mut acc);

    let errors: Vec<f64> = match reference {
        Reference::Full { forces } => {
            acc.iter().zip(forces).map(|(a, r)| relative_force_error(*a, *r)).collect()
        },
        Reference::Sampled { indices, forces } => {
            indices.iter().zip(forces).map(|(&idx, r)| relative_force_error(acc[idx], *r)).collect()
        },
    };

    if let Some(signal) = small_force_signal(reference, &errors)
        && signal < 0.5
    {
        eprintln!(
            "[perf_2x2-diag]   worst-body |F|/median = {signal:.2} \
             (small-force outlier; max/p99 inflated by relative-error denominator)"
        );
    }

    distribution_stats(&errors)
}

/// Ratio `|F_worst| / median(|F|)` where `F_worst` is the reference force on
/// the body that hit the largest relative error. A value below ~0.5 means
/// the max-error body sits in a small-|F| pocket, so its relative error is
/// inflated by the denominator rather than by an actual BH defect — a known
/// structural artefact of the metric, not the algorithm. Returns `None` for
/// an empty error array.
fn small_force_signal(reference: &Reference, errors: &[f64]) -> Option<f64> {
    let (max_idx, _) = errors.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap())?;
    let ref_forces = match reference {
        Reference::Full { forces } => forces,
        Reference::Sampled { forces, .. } => forces,
    };
    let worst_f = ref_forces[max_idx].length();
    let mut mags: Vec<f64> = ref_forces.iter().map(|v| v.length()).collect();
    mags.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_f = mags[mags.len() / 2];
    Some(worst_f / median_f.max(1e-30))
}

fn relative_force_error(a: Vec3, r: Vec3) -> f64 {
    let r_mag = r.length().max(1e-30);
    (a - r).length() / r_mag
}

fn distribution_stats(values: &[f64]) -> ErrorStats {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();

    // Linear-interpolation-free percentile (nearest-rank). With n ≥ 100 the
    // discretisation error is below the per-body noise we care about; for
    // n = K_SAMPLE = 256, p95 falls on the 244th element and p99 on the
    // 254th — the latter has only 2 samples in the tail, hence the
    // notebook treats p99 as informational.
    let percentile = |q: f64| -> f64 {
        let idx = ((n as f64 - 1.0) * q).round() as usize;
        sorted[idx.min(n - 1)]
    };

    let mean = values.iter().sum::<f64>() / n as f64;
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;

    ErrorStats {
        p50: percentile(0.50),
        p95: percentile(0.95),
        p99: percentile(0.99),
        max: sorted[n - 1],
        mean,
        std: var.sqrt(),
    }
}

// ── CSV / paths ──────────────────────────────────────────────────────────── //

fn write_row(writer: &mut fs::File, row: &CsvRow) {
    let e = &row.error;
    writeln!(
        writer,
        "{},{},{},{:#x},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6},{:.6},{:.6},{:.6}",
        row.cell,
        row.n,
        row.theta,
        row.seed,
        e.p50,
        e.p95,
        e.p99,
        e.max,
        e.mean,
        e.std,
        row.t_build_ms,
        row.t_walk_ms,
        row.t_eval_ms,
        row.std_err_t_eval_ms,
    )
    .unwrap();
}

fn perf_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-2x2")
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
