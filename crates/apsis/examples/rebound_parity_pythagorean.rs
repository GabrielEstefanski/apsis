//! REBOUND parity — Pythagorean three-body (Burrau 1913), apsis IAS15 side.
//!
//! Runs the canonical Burrau Pythagorean problem — masses 3, 4, 5 at the
//! vertices of a 3-4-5 right triangle, released from rest — under apsis
//! IAS15 to a fixed horizon (canonical t = 70), sampling state and total
//! energy at uniform cadence. Output is a CSV consumable by the matching
//! Python REBOUND harness for cross-implementation parity comparison.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_pythagorean -p apsis
//! cargo run --release --example rebound_parity_pythagorean -p apsis -- --output path/to/apsis.csv
//! cargo run --release --example rebound_parity_pythagorean -p apsis -- --horizon 200  # informational long-horizon stress
//! ```
//!
//! Default output path: `validation/rebound-parity/pythagorean/out/apsis.csv`
//! (relative to the workspace root).
//!
//! ## Protocol
//!
//! The full protocol — initial conditions, integrator settings, metrics,
//! tolerances declared *a priori*, and the metric tier hierarchy — is
//! specified in
//! [`docs/experiments/2026-04-30-rebound-parity-pythagorean.md`](../../../../docs/experiments/2026-04-30-rebound-parity-pythagorean.md).
//!
//! Constants in this file mirror the protocol's IC and run-parameter
//! values. Changes here are protocol changes — update the notebook in
//! lockstep.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Body masses in canonical units (G = 1). Burrau (1913) opposite-side
/// convention: the side opposite mass mᵢ has length mᵢ.
const MASSES: [f64; 3] = [3.0, 4.0, 5.0];

/// Burrau (1913) initial conditions. All position and velocity components
/// are integer-valued in canonical units, so the f64 bit pattern is
/// identical between Rust and Python implementations on any IEEE-754
/// platform — eliminating IC-construction divergence as a confound.
const R1: (f64, f64) = (1.0, 3.0);
const R2: (f64, f64) = (-2.0, -1.0);
const R3: (f64, f64) = (1.0, -1.0);
const V1: (f64, f64) = (0.0, 0.0);
const V2: (f64, f64) = (0.0, 0.0);
const V3: (f64, f64) = (0.0, 0.0);

/// Default integration horizon in canonical time units. Exceeds both
/// "completion" reference points in the literature (Aarseth 2003 §3 cites
/// t ≈ 46; Szebehely & Peters 1967 Fig. 5 extends through t ≈ 60) without
/// entering the regime where the post-ejection binary's rapid orbit
/// dominates substep selection.
const HORIZON: f64 = 70.0;

/// Dense analysis cadence: samples emitted per canonical time unit. Total
/// rows = HORIZON × SAMPLES_PER_TIME_UNIT + 1 = 2101 at the default horizon.
const SAMPLES_PER_TIME_UNIT: u64 = 30;

/// Initial timestep in canonical units. Matches the `three_body_pythagorean`
/// template preset's `suggested_dt`; the IAS15 controller grows and shrinks
/// from this seed.
const DT_INITIAL: f64 = 1.0e-3;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let cli = parse_cli();
    if let Some(parent) = cli.output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions ──────────────────────────────────────────────── //
    //
    // Burrau ICs satisfy Σ mᵢ 𝐫ᵢ = Σ mᵢ 𝐯ᵢ = Σ mᵢ 𝐫ᵢ × 𝐯ᵢ = 𝟎 exactly
    // (verified algebraically in the protocol notebook §Methodology). No
    // additional recentering: the ICs are integer-valued f64, and any
    // explicit shift would introduce an implementation-divergent
    // perturbation to a state already bit-identical between sides.
    let body1 = Body::rocky(MASSES[0]).at(R1.0, R1.1).with_velocity(V1.0, V1.1).unsoftened();
    let body2 = Body::rocky(MASSES[1]).at(R2.0, R2.1).with_velocity(V2.0, V2.1).unsoftened();
    let body3 = Body::rocky(MASSES[2]).at(R3.0, R3.1).with_velocity(V3.0, V3.1).unsoftened();

    // ── Integrator setup ────────────────────────────────────────────────── //
    let mut sys = System::new(vec![body1, body2, body3], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(DT_INITIAL);

    // ── CSV output ──────────────────────────────────────────────────────── //
    //
    // Wide format identical to `rebound_parity_figure8.rs`: one row per
    // dense sample, 15 columns (sample, t, per-body x/y/vx/vy × 3, total
    // energy). Initial state at sample = 0; subsequent samples at uniform
    // 1 / SAMPLES_PER_TIME_UNIT spacing. Total = HORIZON ×
    // SAMPLES_PER_TIME_UNIT + 1 rows.
    let horizon = cli.horizon;
    let total_samples = (horizon * SAMPLES_PER_TIME_UNIT as f64).round() as u64 + 1;
    let dt_sample = 1.0 / SAMPLES_PER_TIME_UNIT as f64;

    let file = File::create(&cli.output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity — Pythagorean three-body (Burrau 1913) — apsis IAS15 side")
        .unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-04-30-rebound-parity-pythagorean.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis)").unwrap();
    writeln!(w, "# units: canonical (G = 1)").unwrap();
    writeln!(
        w,
        "# masses=({m0},{m1},{m2}), horizon={horizon:.18e}",
        m0 = MASSES[0],
        m1 = MASSES[1],
        m2 = MASSES[2],
    )
    .unwrap();
    writeln!(w, "# samples_per_t_unit={SAMPLES_PER_TIME_UNIT}, dt0={dt:.18e}", dt = DT_INITIAL,)
        .unwrap();
    writeln!(w, "sample,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,x2,y2,vx2,vy2,e_total").unwrap();

    write_sample(&mut w, 0, &sys);
    for n in 1..total_samples {
        let t_target = (n as f64) * dt_sample;
        sys.integrate_until(t_target);
        write_sample(&mut w, n, &sys);
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {}", total_samples, cli.output_path.display());

    // Diagnostic counters surfaced for cross-implementation comparison.
    // Pythagorean is a stiff-mix scenario; the substep economy and
    // rejection-class breakdown characterise how the controller navigated
    // the close-encounter regime. Comparable values from REBOUND go to
    // its own stderr in `rebound_side.py`.
    eprintln!("apsis substeps total: {}", sys.steps());
    if let Some(stats) = sys.adaptive_stats() {
        eprintln!(
            "apsis adaptive stats: rejections={} (picard={}, truncation={}), \
             degraded={}, picard_iters={}, picard_stagnations={}, shrink_grow_cycles={}",
            stats.rejections,
            stats.rejections_picard,
            stats.rejections_truncation,
            stats.degraded,
            stats.picard_iters,
            stats.picard_stagnations,
            stats.shrink_grow_cycles,
        );
    }
}

// ── Output helper ───────────────────────────────────────────────────────── //

fn write_sample(w: &mut BufWriter<File>, sample: u64, sys: &System) {
    let bodies = sys.bodies();
    let b0 = &bodies[0];
    let b1 = &bodies[1];
    let b2 = &bodies[2];
    let e_total = total_energy(bodies);
    writeln!(
        w,
        "{sample},{t:.18e},{x0:.18e},{y0:.18e},{vx0:.18e},{vy0:.18e},{x1:.18e},{y1:.18e},{vx1:.18e},{vy1:.18e},{x2:.18e},{y2:.18e},{vx2:.18e},{vy2:.18e},{e:.18e}",
        t = sys.t(),
        x0 = b0.pos_x, y0 = b0.pos_y, vx0 = b0.vel_x, vy0 = b0.vel_y,
        x1 = b1.pos_x, y1 = b1.pos_y, vx1 = b1.vel_x, vy1 = b1.vel_y,
        x2 = b2.pos_x, y2 = b2.pos_y, vx2 = b2.vel_x, vy2 = b2.vel_y,
        e = e_total,
    )
    .unwrap();
}

/// Total mechanical energy, computed inline so the formula is visible at
/// the comparison site and matches REBOUND's `sim.energy()` convention
/// exactly: KE = ½ Σ mᵢ vᵢ², PE = −Σᵢ<ⱼ G mᵢ mⱼ / rᵢⱼ, with G = 1 and no
/// softening (verified by `Body::unsoftened()` on every body). The
/// expression is mass-distribution-independent — `b.mass` is per-body —
/// so the same code works for the Pythagorean (3, 4, 5) and the figure-8
/// (1, 1, 1) configurations without modification.
fn total_energy(bodies: &[Body]) -> f64 {
    let ke: f64 =
        bodies.iter().map(|b| 0.5 * b.mass * (b.vel_x * b.vel_x + b.vel_y * b.vel_y)).sum();
    let mut pe = 0.0;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let r = (dx * dx + dy * dy).sqrt();
            pe -= bodies[i].mass * bodies[j].mass / r;
        }
    }
    ke + pe
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

struct Cli {
    output_path: PathBuf,
    horizon: f64,
}

fn parse_cli() -> Cli {
    let mut output_path: Option<PathBuf> = None;
    let mut horizon: f64 = HORIZON;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" | "-o" => {
                output_path =
                    Some(PathBuf::from(args.next().expect("--output requires a path argument")));
            },
            "--horizon" => {
                horizon = args
                    .next()
                    .expect("--horizon requires a positive float")
                    .parse()
                    .expect("--horizon must be a positive float");
                assert!(horizon > 0.0, "--horizon must be > 0");
            },
            other => panic!("unknown argument: {other}"),
        }
    }

    Cli {
        output_path: output_path.unwrap_or_else(|| {
            PathBuf::from("validation/rebound-parity/pythagorean/out/apsis.csv")
        }),
        horizon,
    }
}
