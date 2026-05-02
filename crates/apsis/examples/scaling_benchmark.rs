//! Scaling benchmark — wall-time per step across N, integrators, and scenarios.
//!
//! Measures the empirical cost of a single integration step using the
//! scenario catalogue in [`common::scenarios`]. Each scenario declares
//! its own natural timestep, so the measured `ms/step` values are
//! comparable in units of "useful simulation advance" rather than
//! units of "arbitrary fixed dt" — see the module-level doc in
//! `common/scenarios.rs` for the rationale.
//!
//! ## Metrics reported per row
//!
//! - `median_ms`, `p95_ms` — per-step wall-time (directly tells you
//!   whether the configuration fits an interactive render budget).
//! - `sim/wall` — simulation-time advanced per wall-second. The
//!   scenario-independent cost metric: rolls `ms/step` and `dt` into
//!   a single number a researcher can plan a long-term run against.
//! - `tier` — `interactive` (≤ 33 ms/step = 30 FPS), `batch-realtime`
//!   (≤ 1 s/step), `batch-overnight` (≤ 60 s/step), `infeasible`.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example scaling_benchmark
//! ```
//!
//! Release mode is required: debug mode is roughly 10× slower and
//! would mislabel the tier boundaries.

use std::io::Write;
use std::time::{Duration, Instant};

use apsis::core::system::System;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

mod common;
use common::scenarios::{self, Scenario};

// ── Configuration ─────────────────────────────────────────────────────────── //

const N_VALUES: &[usize] = &[128, 256, 512, 1024, 2048, 4096, 8192, 16_384, 32_768, 65_536];
const WARMUP_STEPS: usize = 10;
const MEASURED_STEPS: usize = 50;
const MAX_PER_STEP: Duration = Duration::from_secs(10);
const SEED: u64 = 0x5EED;

// ── Measurement ───────────────────────────────────────────────────────────── //

#[derive(Clone, Copy)]
#[allow(dead_code)]
struct Row {
    scenario: &'static str,
    integrator: IntegratorKind,
    n: usize,
    dt_hint: f64,
    median_ms: f64,
    p95_ms: f64,
    steps_per_sec: f64,
    sim_rate: f64,
}

fn classify(median_ms: f64) -> &'static str {
    if median_ms <= 33.0 {
        "interactive"
    } else if median_ms <= 1_000.0 {
        "batch-realtime"
    } else if median_ms <= 60_000.0 {
        "batch-overnight"
    } else {
        "infeasible"
    }
}

fn measure(scenario: &Scenario, integrator: IntegratorKind, n: usize) -> Option<Row> {
    let bodies = (scenario.build)(n, SEED);
    let dt = (scenario.dt_hint)(n);
    let mut sys =
        System::new(bodies, UnitSystem::canonical()).with_integrator(integrator).with_dt(dt);

    // Warm-up.
    for _ in 0..WARMUP_STEPS {
        let t0 = Instant::now();
        sys.step();
        if t0.elapsed() > MAX_PER_STEP {
            return None;
        }
    }

    let mut samples: Vec<(Duration, f64)> = Vec::with_capacity(MEASURED_STEPS);
    for _ in 0..MEASURED_STEPS {
        let t_before = sys.t();
        let wall_t0 = Instant::now();
        sys.step();
        let wall_dt = wall_t0.elapsed();
        let sim_dt = sys.t() - t_before;
        if wall_dt > MAX_PER_STEP {
            return None;
        }
        samples.push((wall_dt, sim_dt));
    }

    let mut step_ms: Vec<f64> = samples.iter().map(|(w, _)| w.as_secs_f64() * 1_000.0).collect();
    step_ms.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_ms = step_ms[step_ms.len() / 2];
    let p95_idx = ((step_ms.len() as f64) * 0.95) as usize;
    let p95_ms = step_ms[p95_idx.min(step_ms.len() - 1)];

    let total_wall_s: f64 = samples.iter().map(|(w, _)| w.as_secs_f64()).sum();
    let total_sim: f64 = samples.iter().map(|(_, s)| *s).sum();
    let sim_rate = if total_wall_s > 0.0 { total_sim / total_wall_s } else { 0.0 };

    Some(Row {
        scenario: scenario.name,
        integrator,
        n,
        dt_hint: dt,
        median_ms,
        p95_ms,
        steps_per_sec: 1_000.0 / median_ms,
        sim_rate,
    })
}

/// Per-integrator N ceiling override. IAS15's per-step cost is a large
/// multiple of the fixed-step integrators' (adaptive Picard iterations
/// plus the deterministic-force requirement forcing direct O(N²)), so
/// its practical ceiling is around 4096 for this bench.
fn integrator_max_n(kind: IntegratorKind) -> usize {
    match kind {
        IntegratorKind::Ias15 => 4096,
        _ => usize::MAX,
    }
}

fn print_header(scenario: &Scenario) {
    println!();
    println!("## {} — {}", scenario.name, scenario.description);
    println!("   dt_hint at N=1024: {:.3e}", (scenario.dt_hint)(1024));
    println!();
    println!(
        "{:>12} {:>8} {:>11} {:>10} {:>10} {:>10} {:>12}  tier",
        "integrator", "N", "dt", "median_ms", "p95_ms", "steps/s", "sim/wall"
    );
    println!("{}", "-".repeat(100));
}

fn print_row(row: &Row) {
    println!(
        "{:>12} {:>8} {:>11.3e} {:>10.3} {:>10.3} {:>10.1} {:>12.3e}  {}",
        format!("{:?}", row.integrator),
        row.n,
        row.dt_hint,
        row.median_ms,
        row.p95_ms,
        row.steps_per_sec,
        row.sim_rate,
        classify(row.median_ms)
    );
}

fn print_infeasible(integrator: IntegratorKind, n: usize, dt: f64) {
    let dash = "—";
    println!(
        "{:>12} {:>8} {:>11.3e} {:>10} {:>10} {:>10} {:>12}  infeasible (> MAX_PER_STEP)",
        format!("{:?}", integrator),
        n,
        dt,
        dash,
        dash,
        dash,
        dash,
    );
}

fn run_integrator_sweep(scenario: &Scenario, kind: IntegratorKind) {
    let max_n = integrator_max_n(kind).min(scenario.max_n);
    for &n in N_VALUES {
        if n < scenario.min_n || n > max_n {
            continue;
        }
        print!("{:>12} {n:>8}  measuring...  \r", format!("{:?}", kind));
        std::io::stdout().flush().ok();

        match measure(scenario, kind, n) {
            Some(row) => print_row(&row),
            None => {
                print_infeasible(kind, n, (scenario.dt_hint)(n));
                break;
            },
        }
    }
}

fn main() {
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!("# apsis scaling benchmark");
    println!(
        "# build: {}, cores: {}",
        if cfg!(debug_assertions) { "debug (SLOW — run with --release)" } else { "release" },
        cores
    );
    println!(
        "# warm-up: {WARMUP_STEPS} steps, measured: {MEASURED_STEPS} steps, abort: {} s/step",
        MAX_PER_STEP.as_secs()
    );

    let integrators =
        [IntegratorKind::VelocityVerlet, IntegratorKind::Yoshida4, IntegratorKind::Ias15];

    for scenario in scenarios::all() {
        print_header(&scenario);
        for &kind in &integrators {
            run_integrator_sweep(&scenario, kind);
        }
    }

    println!();
    println!("# tier definitions:");
    println!("#   interactive      ≤ 33 ms/step (30 FPS render loop)");
    println!("#   batch-realtime   ≤ 1 s/step  (~10 orbits in minutes)");
    println!("#   batch-overnight  ≤ 60 s/step (~10³ steps in a working day)");
    println!("#   infeasible       otherwise");
    println!();
    println!("# metric: sim/wall  simulation-time advanced per wall-second");
    println!("#   fixed-step integrators: sim/wall = dt_hint / (step wall-time)");
    println!("#   IAS15 (adaptive):       sim/wall reflects the average dt actually taken");
}
