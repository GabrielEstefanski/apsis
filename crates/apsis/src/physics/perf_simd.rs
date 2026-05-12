//! SIMD walk perf harness — Tier 3 (end-to-end walk wall-time speedup)
//! and Tier 4 (pack overhead) for the Plummer leaf-pair SIMD kernels.
//!
//! Lab notebook: `docs/experiments/2026-05-11-simd-kernel.md`.
//!
//! Opt-in:
//!
//! ```text
//! cargo test --release -p apsis perf_simd_walk -- --ignored --nocapture
//! ```
//!
//! Sweeps the `(N × seed × LeafPairKernel)` grid declared a priori in the
//! notebook §Methodology — `N ∈ {1_000, 5_000, 10_000}`, three seeds,
//! three dispatch paths (`Scalar`, `Avx2`, `Avx512` where available) —
//! and reports per-cell median walk wall-time over 5 measured runs after
//! 3 warmup runs. Speedup tables are computed against the `Scalar` cell
//! at the same `(N, seed)` and printed alongside the engine ceiling
//! envelope `[1.3, 2.0]× AVX2` / `[1.7, 2.7]× AVX-512` from the
//! notebook §Tier 3 bounds.
//!
//! Tier 0/1/2a gates live as unit tests next to the kernels themselves
//! (`crates/apsis/src/physics/gravity/engine.rs`); they are not re-run
//! by this harness.
//!
//! CSV output: `target/perf-simd/walk.csv`. The harness is deleted in
//! the bake commit per the perf-2 / perf-4 / perf-5 closure pattern.

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use crate::domain::body::Body;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::gravity::{BarnesHutEngine, LeafPairKernel};

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const N_VALUES: [usize; 3] = [1_000, 5_000, 10_000];
const THETA: f64 = 0.5;
const WARMUP_RUNS: usize = 3;
const MEASURED_RUNS: usize = 5;

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release perf_simd_walk -- --ignored --nocapture"]
fn perf_simd_walk() {
    let kernels = available_kernels();

    let out_dir = perf_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf-simd output dir");
    let csv_path = out_dir.join("walk.csv");
    let mut writer = fs::File::create(&csv_path).expect("create walk.csv");
    writeln!(
        writer,
        "kernel,n,seed,t_walk_ns_median,n_node_visits,n_bh_accepted,n_leaf_interactions"
    )
    .unwrap();

    eprintln!(
        "[perf-simd] Tier 3 walk + Tier 4 pack — kernels={kernels:?} N={N_VALUES:?} seeds={SEEDS:#X?}"
    );
    let t_total = Instant::now();

    let mut rows: Vec<Row> = Vec::new();
    for &kernel in &kernels {
        for &n in &N_VALUES {
            for &seed in &SEEDS {
                let row = measure_walk(kernel, n, seed);
                writeln!(
                    writer,
                    "{:?},{},0x{:X},{},{},{},{}",
                    row.kernel,
                    row.n,
                    row.seed,
                    row.t_walk_ns_median,
                    row.n_node_visits,
                    row.n_bh_accepted,
                    row.n_leaf_interactions,
                )
                .unwrap();
                eprintln!(
                    "[perf-simd]  {:>7?}  N={:>5}  seed=0x{:X}  t_walk={:>8.3}ms  \
                     accepted={:>6}  leafpair={:>7}",
                    row.kernel,
                    row.n,
                    row.seed,
                    row.t_walk_ns_median as f64 * 1e-6,
                    row.n_bh_accepted,
                    row.n_leaf_interactions,
                );
                rows.push(row);
            }
        }
    }

    print_speedup_table(&rows, &kernels);
    print_tier4_pack_overhead();

    eprintln!("[perf-simd] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[perf-simd] wrote {}", csv_path.display());
}

#[derive(Debug, Clone, Copy)]
struct Row {
    kernel: LeafPairKernel,
    n: usize,
    seed: u64,
    t_walk_ns_median: u64,
    n_node_visits: u64,
    n_bh_accepted: u64,
    n_leaf_interactions: u64,
}

fn measure_walk(kernel: LeafPairKernel, n: usize, seed: u64) -> Row {
    let bodies = sphere_distribution_lognormal(n, seed);
    let mut arrays = BodyArrays::with_capacity(n);
    arrays.pack_from(&bodies);

    let mut engine = BarnesHutEngine::new(16);
    engine.leaf_pair_kernel = kernel;
    engine.build(&arrays);

    let mut acc = vec![Vec3::ZERO; n];

    for _ in 0..WARMUP_RUNS {
        engine.evaluate(&arrays, THETA, &mut acc);
    }

    let mut samples_ns: Vec<u64> = Vec::with_capacity(MEASURED_RUNS);
    let mut last_counters = (0_u64, 0_u64, 0_u64);
    for _ in 0..MEASURED_RUNS {
        let t = Instant::now();
        let (_, c) = engine.evaluate_profile(&arrays, THETA, &mut acc);
        samples_ns.push(t.elapsed().as_nanos() as u64);
        last_counters = (c.n_node_visits, c.n_bh_accepted, c.n_leaf_interactions);
    }
    samples_ns.sort_unstable();
    let t_walk_ns_median = samples_ns[MEASURED_RUNS / 2];

    Row {
        kernel,
        n,
        seed,
        t_walk_ns_median,
        n_node_visits: last_counters.0,
        n_bh_accepted: last_counters.1,
        n_leaf_interactions: last_counters.2,
    }
}

fn print_speedup_table(rows: &[Row], kernels: &[LeafPairKernel]) {
    let scalar_baseline = |n: usize, seed: u64| -> Option<u64> {
        rows.iter()
            .find(|r| r.kernel == LeafPairKernel::Scalar && r.n == n && r.seed == seed)
            .map(|r| r.t_walk_ns_median)
    };

    eprintln!();
    eprintln!("[perf-simd] ── Tier 3 walk speedup (per cell, vs scalar) ──");
    eprintln!(
        "[perf-simd] {:>7}  {:>5}  {:>16}  {:>16}  {:>16}",
        "kernel", "N", "seed=0x6F637472", "seed=0x71756164", "seed=0x6D6F7274"
    );
    for &kernel in kernels {
        if kernel == LeafPairKernel::Scalar {
            continue;
        }
        for &n in &N_VALUES {
            let speedups: Vec<String> = SEEDS
                .iter()
                .map(|&seed| {
                    let t_scalar = scalar_baseline(n, seed);
                    let t_kernel = rows
                        .iter()
                        .find(|r| r.kernel == kernel && r.n == n && r.seed == seed)
                        .map(|r| r.t_walk_ns_median);
                    match (t_scalar, t_kernel) {
                        (Some(s), Some(k)) if k > 0 => format!("{:>14.3}×", s as f64 / k as f64),
                        _ => "n/a".to_string(),
                    }
                })
                .collect();
            eprintln!(
                "[perf-simd] {:>7?}  {:>5}  {:>16}  {:>16}  {:>16}",
                kernel, n, speedups[0], speedups[1], speedups[2]
            );
        }
    }
    eprintln!();
    eprintln!("[perf-simd] ── Tier 3 median speedup (across seeds, at each N) ──");
    eprintln!("[perf-simd] {:>7}  {:>5}  {:>10}  envelope", "kernel", "N", "median×");
    for &kernel in kernels {
        if kernel == LeafPairKernel::Scalar {
            continue;
        }
        let envelope = match kernel {
            #[cfg(target_arch = "x86_64")]
            LeafPairKernel::Avx2 => "[1.3, 2.0]× per notebook §Tier 3",
            #[cfg(target_arch = "x86_64")]
            LeafPairKernel::Avx512 => "[1.7, 2.7]× per notebook §Tier 3",
            _ => "—",
        };
        for &n in &N_VALUES {
            let mut ratios: Vec<f64> = SEEDS
                .iter()
                .filter_map(|&seed| {
                    let t_scalar = scalar_baseline(n, seed)?;
                    let t_kernel = rows
                        .iter()
                        .find(|r| r.kernel == kernel && r.n == n && r.seed == seed)?
                        .t_walk_ns_median;
                    if t_kernel == 0 {
                        return None;
                    }
                    Some(t_scalar as f64 / t_kernel as f64)
                })
                .collect();
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = ratios[ratios.len() / 2];
            eprintln!("[perf-simd] {:>7?}  {:>5}  {:>9.3}×  {}", kernel, n, median, envelope);
        }
    }
}

// ── Tier 4 — pack / compute overhead ──────────────────────────────────────── //

fn print_tier4_pack_overhead() {
    eprintln!();
    eprintln!("[perf-simd] ── Tier 4 pack / compute ratio ──");
    eprintln!(
        "[perf-simd] {:>5}  {:>16}  {:>16}  {:>10}  bound",
        "N", "t_pack_ns_median", "t_compute_ns_median", "ratio"
    );

    let kernel = pick_default_kernel();
    for &n in &N_VALUES {
        let bodies = sphere_distribution_lognormal(n, SEEDS[0]);
        let mut arrays = BodyArrays::with_capacity(n);
        let mut engine = BarnesHutEngine::new(16);
        engine.leaf_pair_kernel = kernel;
        let mut acc = vec![Vec3::ZERO; n];

        for _ in 0..WARMUP_RUNS {
            arrays.pack_from(&bodies);
            engine.build(&arrays);
            engine.evaluate(&arrays, THETA, &mut acc);
        }

        let mut pack_samples: Vec<u64> = Vec::with_capacity(MEASURED_RUNS);
        let mut compute_samples: Vec<u64> = Vec::with_capacity(MEASURED_RUNS);
        for _ in 0..MEASURED_RUNS {
            let t = Instant::now();
            arrays.pack_from(&bodies);
            pack_samples.push(t.elapsed().as_nanos() as u64);

            let t = Instant::now();
            engine.build(&arrays);
            engine.evaluate(&arrays, THETA, &mut acc);
            compute_samples.push(t.elapsed().as_nanos() as u64);
        }
        pack_samples.sort_unstable();
        compute_samples.sort_unstable();

        let t_pack = pack_samples[MEASURED_RUNS / 2];
        let t_compute = compute_samples[MEASURED_RUNS / 2];
        let ratio = t_pack as f64 / t_compute as f64;
        let bound_ok = ratio <= 0.01;
        eprintln!(
            "[perf-simd] {:>5}  {:>16}  {:>16}  {:>9.4}  {} (≤ 0.01)",
            n,
            t_pack,
            t_compute,
            ratio,
            if bound_ok { "PASS" } else { "FAIL" },
        );
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────── //

fn available_kernels() -> Vec<LeafPairKernel> {
    let mut out = vec![LeafPairKernel::Scalar];
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
            out.push(LeafPairKernel::Avx2);
        }
        if std::is_x86_feature_detected!("avx512f") {
            out.push(LeafPairKernel::Avx512);
        }
    }
    out
}

/// Pick the fastest available leaf-pair kernel for the host. Used by
/// Tier 4's pack-overhead measurement to mirror what production picks
/// at engine construction.
fn pick_default_kernel() -> LeafPairKernel {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx512f") {
            return LeafPairKernel::Avx512;
        }
        if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
            return LeafPairKernel::Avx2;
        }
    }
    LeafPairKernel::Scalar
}

fn perf_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-simd")
}

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
