//! Scaling benchmark — wall-time vs N across scenarios and integrators.
//!
//! Measures the empirical cost of a single integration step across a
//! range of body counts, integrators, and physically-distinct scenarios.
//! Each scenario declares its own natural timestep so that the measured
//! wall-time-per-step is comparable in units of "useful simulation
//! advance" rather than in units of "arbitrary fixed dt." The latter is
//! a hidden bias: a hierarchical Solar-system-like configuration at
//! `dt = 10⁻³` is massively over-resolved relative to a dense cluster
//! at the same `dt`, and comparing them on `ms/step` alone is unfair to
//! the cluster scenario.
//!
//! ## Metrics reported per row
//!
//! - `median_ms`, `p95_ms` — per-step wall-time (compute-cost metric;
//!   directly tells you whether the configuration fits an interactive
//!   render budget).
//! - `sim_rate` — simulation-time advanced per wall-second. This is the
//!   honest cross-scenario metric: it rolls `ms/step` and `dt` into a
//!   single number a researcher can plan a run against
//!   ("how many orbits per wall-hour?").
//! - `tier` — `interactive` (≤ 33 ms/step = 30 FPS), `batch-realtime`
//!   (≤ 1 s/step), `batch-overnight` (≤ 60 s/step), `infeasible` otherwise.
//!
//! ## Scenarios
//!
//! - **friendly_cluster** — uniform 2D disk with virial-scale velocities.
//!   Smooth density, no sub-structure. The baseline regime.
//! - **hierarchical_kepler** — central 1 M☉ plus N−1 test bodies on
//!   log-uniform Kepler orbits from 0.3 to 30 AU. Exercises the
//!   Solar-system-like regime where Wisdom–Holman shines and the BH
//!   tree is sparsely populated in the outer region.
//! - **clustered_substructure** — a main Plummer-like cluster
//!   surrounded by several denser sub-clumps. Exercises tree-depth
//!   non-degeneracy and pushes the opening-criterion hot path.
//! - **multiple_binaries** — `N/2` isolated two-body systems placed
//!   far apart. Tests the tree's far-field behaviour when sources are
//!   sparse and well-separated.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example scaling_benchmark
//! ```
//!
//! Release mode is required: debug mode is roughly 10× slower and would
//! mislabel the tier boundaries.

use std::f64::consts::TAU;
use std::io::Write;
use std::time::{Duration, Instant};

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

// ── Configuration ─────────────────────────────────────────────────────────── //

/// Body counts to sweep. Log-spaced; each scenario filters to its valid range.
const N_VALUES: &[usize] = &[128, 256, 512, 1024, 2048, 4096, 8192, 16_384, 32_768, 65_536];

const WARMUP_STEPS: usize = 10;
const MEASURED_STEPS: usize = 50;
const MAX_PER_STEP: Duration = Duration::from_secs(10);
const SEED: u64 = 0x5EED;

// ── Scenario catalogue ────────────────────────────────────────────────────── //

/// A self-contained scenario specification. Pure data plus two function
/// pointers — [`Scenario::build`] materialises the bodies, [`Scenario::dt_hint`]
/// returns a natural timestep chosen from the scenario's shortest
/// dynamical timescale at the given `N`.
struct Scenario {
    name: &'static str,
    description: &'static str,
    build: fn(usize, u64) -> Vec<Body>,
    /// Natural timestep for the scenario at the given N. Picked per
    /// scenario to avoid the hidden bias of one-size-fits-all `dt`.
    dt_hint: fn(usize) -> f64,
    min_n: usize,
    max_n: usize,
}

fn all_scenarios() -> Vec<Scenario> {
    vec![
        Scenario {
            name: "friendly_cluster",
            description: "Uniform 2D disk, virial-scale velocities",
            build: build_friendly_cluster,
            dt_hint: dt_hint_cluster,
            min_n: 128,
            max_n: 65_536,
        },
        Scenario {
            name: "hierarchical_kepler",
            description: "1 M☉ + (N−1) test bodies on log-uniform Kepler orbits",
            build: build_hierarchical_kepler,
            dt_hint: dt_hint_hierarchical,
            min_n: 128,
            max_n: 65_536,
        },
        Scenario {
            name: "clustered_substructure",
            description: "Main Plummer-like cluster + sub-clumps",
            build: build_clustered_substructure,
            dt_hint: dt_hint_cluster,
            min_n: 128,
            max_n: 32_768,
        },
        Scenario {
            name: "multiple_binaries",
            description: "N/2 isolated two-body systems, far-separated",
            build: build_multiple_binaries,
            dt_hint: dt_hint_binaries,
            min_n: 128,
            max_n: 32_768,
        },
    ]
}

// ── dt_hint functions ────────────────────────────────────────────────────── //

/// For a uniform disk of radius R with total mass scaling as M ∝ N, the
/// characteristic dynamical time `t_dyn ~ √(R³/M)` decreases slowly with
/// N. `dt = t_dyn / 1000` gives ~1000 steps per crossing time: enough to
/// resolve a close encounter without over-resolving the bulk orbit.
fn dt_hint_cluster(n: usize) -> f64 {
    let r_disk = 10.0_f64;
    let m_total = (n as f64) * 1e-4;
    let t_dyn = (r_disk.powi(3) / m_total).sqrt();
    t_dyn / 1000.0
}

/// For a hierarchical Kepler system with innermost orbit at `a_min = 0.3`,
/// the shortest period is `T_inner = 2π·a_min^{3/2}` (with GM = 1). `dt`
/// is that period divided by 1000.
fn dt_hint_hierarchical(_n: usize) -> f64 {
    let a_min = 0.3_f64;
    let t_inner = TAU * a_min.powf(1.5);
    t_inner / 1000.0
}

/// For a collection of binaries with separation ~ 1 and mass ~ 0.5 each,
/// the orbital period is `T = 2π·√(a³/(M))` with `M = 2·0.5 = 1` → `T = 2π`.
fn dt_hint_binaries(_n: usize) -> f64 {
    TAU / 1000.0
}

// ── Scenario builders ─────────────────────────────────────────────────────── //

fn build_friendly_cluster(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let r_disk: f64 = 10.0;
    let m_total: f64 = (n as f64) * 1e-4;
    let m_each: f64 = m_total / n as f64;
    let v_scale = (m_total / r_disk).sqrt() * 0.2;

    (0..n)
        .map(|_| {
            let theta = rng.random::<f64>() * TAU;
            let r = r_disk * rng.random::<f64>().sqrt();
            let x = r * theta.cos();
            let y = r * theta.sin();
            let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_scale;
            let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_scale;
            Body::rocky(m_each).at(x, y).with_velocity(vx, vy)
        })
        .collect()
}

/// Sun at origin plus N−1 bodies on circular Kepler orbits with
/// log-uniform semi-major axes in [0.3, 30] AU.
fn build_hierarchical_kepler(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let mut bodies = Vec::with_capacity(n);
    bodies.push(Body::star(1.0).at(0.0, 0.0).unsoftened());

    let log_a_min = 0.3_f64.log10();
    let log_a_max = 30.0_f64.log10();
    let log_m_min = -6.0_f64;
    let log_m_max = -4.0_f64;

    for _ in 1..n {
        let a = 10f64.powf(log_a_min + (log_a_max - log_a_min) * rng.random::<f64>());
        let theta = rng.random::<f64>() * TAU;
        let x = a * theta.cos();
        let y = a * theta.sin();
        let v = (1.0 / a).sqrt();
        let vx = -v * theta.sin();
        let vy = v * theta.cos();
        let m = 10f64.powf(log_m_min + (log_m_max - log_m_min) * rng.random::<f64>());
        bodies.push(Body::rocky(m).at(x, y).with_velocity(vx, vy));
    }
    bodies
}

/// 80% of bodies form a diffuse main disk; 20% form `K_sub` dense
/// sub-clumps placed at random positions around the main centre.
fn build_clustered_substructure(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let r_main: f64 = 10.0;
    let r_sub: f64 = 1.0;
    let m_total: f64 = (n as f64) * 1e-4;
    let m_each: f64 = m_total / n as f64;
    let v_main = (m_total / r_main).sqrt() * 0.2;
    let v_sub = (m_total / r_sub).sqrt() * 0.2;

    let n_main = (n as f64 * 0.8) as usize;
    let n_sub_bodies = n - n_main;
    let k_sub = ((n_sub_bodies as f64).sqrt().ceil() as usize).max(2);
    let per_clump = n_sub_bodies / k_sub;

    let mut bodies = Vec::with_capacity(n);

    // Main disk.
    for _ in 0..n_main {
        let theta = rng.random::<f64>() * TAU;
        let r = r_main * rng.random::<f64>().sqrt();
        let (x, y) = (r * theta.cos(), r * theta.sin());
        let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_main;
        let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_main;
        bodies.push(Body::rocky(m_each).at(x, y).with_velocity(vx, vy));
    }

    // Sub-clumps.
    for clump in 0..k_sub {
        let theta_c = rng.random::<f64>() * TAU;
        let r_c = r_main * 0.5 * (1.0 + rng.random::<f64>());
        let (cx, cy) = (r_c * theta_c.cos(), r_c * theta_c.sin());
        let members = if clump == k_sub - 1 {
            n_sub_bodies - per_clump * clump
        } else {
            per_clump
        };
        for _ in 0..members {
            let theta = rng.random::<f64>() * TAU;
            let r = r_sub * rng.random::<f64>().sqrt();
            let (x, y) = (cx + r * theta.cos(), cy + r * theta.sin());
            let vx = (rng.random::<f64>() - 0.5) * 2.0 * v_sub;
            let vy = (rng.random::<f64>() - 0.5) * 2.0 * v_sub;
            bodies.push(Body::rocky(m_each).at(x, y).with_velocity(vx, vy));
        }
    }
    bodies
}

/// `N/2` equal-mass two-body systems, each with separation 1 AU and
/// circular velocity, placed at random positions in a large box so that
/// inter-binary distances dominate intra-binary distances.
fn build_multiple_binaries(n: usize, seed: u64) -> Vec<Body> {
    let mut rng = SmallRng::seed_from_u64(seed);
    let n_binaries = n / 2;
    let field_size = 100.0_f64;
    let sep = 1.0_f64;
    let m_each = 0.5_f64;
    // Circular orbit for the reduced mass: v_rel = √(GM/a) with M = 1, a = 1.
    let v_each = 0.5_f64;

    let mut bodies = Vec::with_capacity(n_binaries * 2);
    for _ in 0..n_binaries {
        let cx = (rng.random::<f64>() - 0.5) * field_size;
        let cy = (rng.random::<f64>() - 0.5) * field_size;
        let phi = rng.random::<f64>() * TAU;
        // Body 1 at (cx − sep/2·cosφ, cy − sep/2·sinφ), velocity ⊥ separation.
        bodies.push(
            Body::rocky(m_each)
                .at(cx - 0.5 * sep * phi.cos(), cy - 0.5 * sep * phi.sin())
                .with_velocity(v_each * phi.sin(), -v_each * phi.cos()),
        );
        bodies.push(
            Body::rocky(m_each)
                .at(cx + 0.5 * sep * phi.cos(), cy + 0.5 * sep * phi.sin())
                .with_velocity(-v_each * phi.sin(), v_each * phi.cos()),
        );
    }
    bodies
}

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
    /// Simulation time advanced per wall-clock second (units: sim-time / s).
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
    let mut sys = System::new(bodies).with_integrator(integrator).with_dt(dt);

    // Warm-up.
    for _ in 0..WARMUP_STEPS {
        let t0 = Instant::now();
        sys.step();
        if t0.elapsed() > MAX_PER_STEP {
            return None;
        }
    }

    // Measured window: record (wall_time, sim_advance) pairs.
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

    let mut step_ms: Vec<f64> =
        samples.iter().map(|(w, _)| w.as_secs_f64() * 1_000.0).collect();
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

// ── Main ─────────────────────────────────────────────────────────────────── //

fn print_header(scenario: &Scenario) {
    println!();
    println!("## {} — {}", scenario.name, scenario.description);
    println!("   dt_hint at N=1024: {:.3e}", (scenario.dt_hint)(1024));
    println!();
    println!(
        "{:>12} {:>8} {:>11} {:>10} {:>10} {:>10} {:>12}  {}",
        "integrator", "N", "dt", "median_ms", "p95_ms", "steps/s", "sim/wall", "tier"
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
    println!(
        "{:>12} {:>8} {:>11.3e} {:>10} {:>10} {:>10} {:>12}  {}",
        format!("{:?}", integrator),
        n,
        dt,
        "—",
        "—",
        "—",
        "—",
        "infeasible (> MAX_PER_STEP)"
    );
}

/// Per-integrator N ceiling override for this benchmark run.
///
/// IAS15's per-step cost is a large multiple of the fixed-step
/// integrators' (adaptive Picard iterations + the deterministic-force
/// requirement forcing direct O(N²)), so its practical ceiling is
/// around 4096. Running it beyond that inflates the total bench
/// wall-time without teaching us anything we do not already see at
/// the lower N.
fn integrator_max_n(kind: IntegratorKind) -> usize {
    match kind {
        IntegratorKind::Ias15 => 4096,
        _ => usize::MAX,
    }
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
            }
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

    let integrators = [
        IntegratorKind::VelocityVerlet,
        IntegratorKind::Yoshida4,
        IntegratorKind::Ias15,
    ];

    for scenario in all_scenarios() {
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
