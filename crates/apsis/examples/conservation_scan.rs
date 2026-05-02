//! Conservation scan — relative drift of energy and angular momentum
//! across N, integrators, and scenarios.
//!
//! Complements [`scaling_benchmark`]. Speed alone is not useful if the
//! integrator is degrading quality as N grows; this bench measures the
//! *other* half of the `N → (is-my-code-valid-here)` question.
//!
//! For each `(scenario, integrator, N)` cell the bench integrates for
//! [`PERIODS`] multiples of the scenario's characteristic timescale,
//! records the initial and final energy and z-angular-momentum, and
//! reports the relative drift `|dE/E|` along with the absolute
//! `|ΔL_z|` (reported as absolute because some scenarios have
//! `L_z ≈ 0` by construction and the relative quantity is undefined).
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example conservation_scan
//! ```
//!
//! Release is required — debug builds are ~10× slower and will push
//! total bench wall-time into tens of minutes. Typical release run on
//! a 12-core workstation: 5 to 15 minutes.

use std::io::Write;
use std::time::{Duration, Instant};

use apsis::core::system::System;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

mod common;
use common::scenarios::{self, Scenario};

// ── Configuration ─────────────────────────────────────────────────────────── //

/// Integer N values scanned. Conservatively capped — we integrate
/// `PERIODS × t_characteristic` per cell, so per-cell cost is much
/// larger than in `scaling_benchmark`. The largest N values here are
/// chosen so that VV/Y4 finish in under ~30 s per cell.
const N_VALUES: &[usize] = &[128, 512, 2048, 4096];

/// Number of characteristic timescales integrated per cell. 10 is
/// enough to see meaningful drift in non-symplectic or poorly-tuned
/// integrators without pushing the total bench past ~30 minutes.
const PERIODS: f64 = 10.0;

const SEED: u64 = 0x5EED;

/// Per-cell wall-time budget. If a cell exceeds this, it is recorded
/// as `aborted` with whatever sim-time has been advanced so far —
/// crucial because IAS15 in stiff scenarios can drive `dt` to its
/// floor (1e-12) and effectively stop progressing, which would
/// otherwise hang the entire bench. Set high enough that healthy
/// cells finish comfortably; low enough that pathological cells fail
/// fast.
const MAX_WALL_PER_CELL: Duration = Duration::from_secs(60);

// ── Per-integrator ceilings ──────────────────────────────────────────────── //

/// IAS15 hits the `dt` floor at N ≥ 1024 in dense random-cluster
/// scenarios. Beyond that, the integrator advances negligibly and
/// the per-cell wall guard would always abort with no useful
/// data. Cap IAS15 at the largest N where we have empirical
/// evidence of healthy convergence on every scenario in this set.
///
/// Yoshida-4 at the upper N values (4096) is roughly 30–60 s per
/// cell at PERIODS = 10; comfortable inside the wall guard.
fn integrator_max_n(kind: IntegratorKind) -> usize {
    match kind {
        IntegratorKind::Ias15 => 512,
        _ => usize::MAX,
    }
}

// ── Measurement ───────────────────────────────────────────────────────────── //

#[derive(Clone, Copy)]
#[allow(dead_code)]
struct Row {
    scenario: &'static str,
    integrator: IntegratorKind,
    n: usize,
    dt_hint: f64,
    t_total_target: f64,
    t_total_actual: f64,
    steps: u64,
    wall_s: f64,
    rel_de: f64,
    lz0: f64,
    abs_dlz: f64,
    aborted: bool,
}

fn measure(scenario: &Scenario, integrator: IntegratorKind, n: usize) -> Row {
    let bodies = (scenario.build)(n, SEED);
    let dt = (scenario.dt_hint)(n);
    let t_char = (scenario.t_characteristic)(n);
    let t_total = PERIODS * t_char;

    let mut sys =
        System::new(bodies, UnitSystem::canonical()).with_integrator(integrator).with_dt(dt);

    // One step populates the cached last_kinetic / last_potential so
    // that sys.energy() returns the real Hamiltonian rather than 0.
    sys.step();
    let e0 = sys.energy();
    let lz0 = sys.lz();
    let t_before = sys.t();
    let steps_before = sys.steps();

    // Manual integration loop with per-cell wall-time guard. We avoid
    // `sys.integrate_for(t_total)` precisely because it has no
    // mechanism to abort if the integrator's adaptive controller
    // hits a `dt` floor and stops advancing — that pathology would
    // otherwise hang the whole bench. `step()` returns at the
    // controller's own pace; we check wall-time after each one.
    let wall = Instant::now();
    let mut aborted = false;
    while sys.t() - t_before < t_total {
        sys.step();
        if wall.elapsed() >= MAX_WALL_PER_CELL {
            aborted = true;
            break;
        }
    }
    let wall_s = wall.elapsed().as_secs_f64();

    let t_actual = sys.t() - t_before;
    let steps = sys.steps() - steps_before;
    let e_f = sys.energy();
    let lz_f = sys.lz();

    let rel_de = (e_f - e0).abs() / e0.abs().max(1e-30);
    let abs_dlz = (lz_f - lz0).abs();

    Row {
        scenario: scenario.name,
        integrator,
        n,
        dt_hint: dt,
        t_total_target: t_total,
        t_total_actual: t_actual,
        steps,
        wall_s,
        rel_de,
        lz0,
        abs_dlz,
        aborted,
    }
}

// ── Presentation ──────────────────────────────────────────────────────────── //

fn print_header(scenario: &Scenario) {
    println!();
    println!("## {} — {}", scenario.name, scenario.description);
    println!(
        "   t_char at N=1024: {:.3e}   |   integrated: {:.1} × t_char",
        (scenario.t_characteristic)(1024),
        PERIODS
    );
    println!();
    println!(
        "{:>12} {:>8} {:>11} {:>9} {:>11} {:>11} {:>11} {:>9} {:>7}  status",
        "integrator", "N", "dt", "steps", "Lz_0", "|dE/E|", "|dLz|", "t_frac", "wall_s",
    );
    println!("{}", "-".repeat(115));
}

fn print_row(row: &Row) {
    let t_frac =
        if row.t_total_target > 0.0 { row.t_total_actual / row.t_total_target } else { 0.0 };
    let status = if row.aborted { "ABORTED" } else { "ok" };
    println!(
        "{:>12} {:>8} {:>11.3e} {:>9} {:>11.3e} {:>11.3e} {:>11.3e} {:>9.3} {:>7.2}  {}",
        format!("{:?}", row.integrator),
        row.n,
        row.dt_hint,
        row.steps,
        row.lz0,
        row.rel_de,
        row.abs_dlz,
        t_frac,
        row.wall_s,
        status,
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

        let row = measure(scenario, kind, n);
        print_row(&row);
    }
}

fn main() {
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!("# apsis conservation scan");
    println!(
        "# build: {}, cores: {}",
        if cfg!(debug_assertions) { "debug (SLOW — run with --release)" } else { "release" },
        cores
    );
    println!("# periods × t_characteristic per cell: {PERIODS}");
    println!("# metric notes:");
    println!("#   |dE/E|   relative energy drift from post-first-step baseline");
    println!("#   |dLz|    absolute L_z drift; use with |Lz_0| to read ratio when meaningful");
    println!("#   steps    steps taken during the integrated window (exclusive of warm-up)");

    let integrators =
        [IntegratorKind::VelocityVerlet, IntegratorKind::Yoshida4, IntegratorKind::Ias15];

    for scenario in scenarios::all() {
        print_header(&scenario);
        for &kind in &integrators {
            run_integrator_sweep(&scenario, kind);
        }
    }
}
