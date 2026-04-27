//! REBOUND parity — Figure-8 choreography, apsis IAS15 side.
//!
//! Runs the canonical Chenciner–Montgomery figure-8 three-body orbit under
//! apsis IAS15 for `N_PERIODS` orbital periods, sampling state and total
//! energy at a dense cadence (`SAMPLES_PER_PERIOD` per period). Output is a
//! CSV consumable by the matching Python REBOUND harness for cross-
//! implementation parity comparison.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_figure8 -p apsis
//! cargo run --release --example rebound_parity_figure8 -p apsis -- --output path/to/apsis.csv
//! cargo run --release --example rebound_parity_figure8 -p apsis -- --periods 50  # 50T sanity run
//! ```
//!
//! Default output path: `validation/rebound-parity/figure8/out/apsis.csv`
//! (relative to the workspace root — `cargo run` should be invoked from the
//! repo root or from the orchestrator in that directory).
//!
//! ## Protocol
//!
//! The full protocol — initial conditions, integrator settings, metrics,
//! tolerances declared *a priori*, and the metric tier hierarchy — is
//! specified in
//! [`docs/experiments/2026-04-26-rebound-parity-figure8.md`](../../../../docs/experiments/2026-04-26-rebound-parity-figure8.md).
//!
//! Constants in this file (`R*`, `V*`, `MASS`, `PERIOD`, `N_PERIODS`,
//! `SAMPLES_PER_PERIOD`, `DT_FRACTION_OF_PERIOD`) mirror the protocol's IC
//! and run-parameter values. Changes here are protocol changes — update the
//! notebook in lockstep.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Equal mass for all three bodies (canonical units, G = 1).
const MASS: f64 = 1.0;

// Chenciner–Montgomery (2000) initial conditions, 8-digit literature form.
// Both implementations parse these exact string literals; the f64 bit
// pattern is therefore identical between sides on the same hardware.
const R1: (f64, f64) = (-0.97000436, 0.24308753);
const R2: (f64, f64) = (0.97000436, -0.24308753);
const R3: (f64, f64) = (0.0, 0.0);
const V1: (f64, f64) = (0.4662036850, 0.4323657300);
const V2: (f64, f64) = (0.4662036850, 0.4323657300);
const V3: (f64, f64) = (-0.93240737, -0.86473146);

/// Figure-8 orbital period in canonical units. Numeric value from
/// Chenciner & Montgomery (2000); refined values in the literature
/// (Simó 2002) differ in the 10th digit and have no effect on parity
/// metrics, only on the absolute time labelling.
const PERIOD: f64 = 6.3259_139_870;

/// Default number of orbital periods integrated for the gated baseline.
/// Overridable via `--periods` for the 50T informational extension.
const N_PERIODS: u64 = 10;

/// Dense analysis cadence: samples emitted per orbital period. The
/// comparator computes `max(·)` aggregates over all 200 × N_PERIODS + 1
/// samples; this is the resolution at which the gate is evaluated.
const SAMPLES_PER_PERIOD: u64 = 200;

/// Initial timestep, expressed as a fraction of the orbital period.
const DT_FRACTION_OF_PERIOD: f64 = 1.0e-3;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let cli = parse_cli();
    if let Some(parent) = cli.output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions ──────────────────────────────────────────────── //
    //
    // No explicit COM/momentum zeroing: the published ICs already satisfy
    // Σ mᵢ 𝐫ᵢ = Σ mᵢ 𝐯ᵢ = Σ mᵢ 𝐫ᵢ × 𝐯ᵢ = 𝟎 to the precision of the 8-digit
    // literals. Forcing additional correction here would introduce an
    // implementation-divergent f64 perturbation to ICs that should be
    // bit-identical between the apsis and REBOUND sides.
    let body1 = Body::rocky(MASS).at(R1.0, R1.1).with_velocity(V1.0, V1.1).unsoftened();
    let body2 = Body::rocky(MASS).at(R2.0, R2.1).with_velocity(V2.0, V2.1).unsoftened();
    let body3 = Body::rocky(MASS).at(R3.0, R3.1).with_velocity(V3.0, V3.1).unsoftened();

    // ── Integrator setup ────────────────────────────────────────────────── //
    let dt0 = PERIOD * DT_FRACTION_OF_PERIOD;
    let mut sys = System::new(vec![body1, body2, body3], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(dt0);

    // ── CSV output ──────────────────────────────────────────────────────── //
    //
    // Wide format: one row per dense sample, including all body state
    // and total energy. Initial state at sample = 0; subsequent samples
    // at uniform `T / SAMPLES_PER_PERIOD` spacing. Total = N_PERIODS *
    // SAMPLES_PER_PERIOD + 1 rows.
    let n_periods = cli.n_periods;
    let total_samples = n_periods * SAMPLES_PER_PERIOD + 1;
    let dt_sample = PERIOD / (SAMPLES_PER_PERIOD as f64);

    let file = File::create(&cli.output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity — Figure-8 choreography — apsis IAS15 side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-04-26-rebound-parity-figure8.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis)").unwrap();
    writeln!(w, "# units: canonical (G = 1)").unwrap();
    writeln!(w, "# mass={MASS}, period={PERIOD:.18e}").unwrap();
    writeln!(
        w,
        "# n_periods={n_periods}, samples_per_period={SAMPLES_PER_PERIOD}, dt0={dt0:.18e}"
    )
    .unwrap();
    writeln!(
        w,
        "sample,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,x2,y2,vx2,vy2,e_total"
    )
    .unwrap();

    write_sample(&mut w, 0, &sys);
    for n in 1..total_samples {
        let t_target = (n as f64) * dt_sample;
        sys.integrate_until(t_target);
        write_sample(&mut w, n, &sys);
    }

    w.flush().unwrap();
    eprintln!(
        "wrote {} samples to {}",
        total_samples,
        cli.output_path.display()
    );
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
        x0 = b0.x, y0 = b0.y, vx0 = b0.vx, vy0 = b0.vy,
        x1 = b1.x, y1 = b1.y, vx1 = b1.vx, vy1 = b1.vy,
        x2 = b2.x, y2 = b2.y, vx2 = b2.vx, vy2 = b2.vy,
        e = e_total,
    )
    .unwrap();
}

/// Total mechanical energy, computed inline so the formula is visible at
/// the comparison site and matches REBOUND's `sim.energy()` convention
/// exactly: KE = ½ Σ mᵢ vᵢ², PE = −Σᵢ<ⱼ G mᵢ mⱼ / rᵢⱼ, with G = 1 and
/// no softening (verified by `Body::unsoftened()` on every body).
fn total_energy(bodies: &[Body]) -> f64 {
    let ke: f64 = bodies.iter().map(|b| 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy)).sum();
    let mut pe = 0.0;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].x - bodies[j].x;
            let dy = bodies[i].y - bodies[j].y;
            let r = (dx * dx + dy * dy).sqrt();
            pe -= bodies[i].mass * bodies[j].mass / r;
        }
    }
    ke + pe
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

struct Cli {
    output_path: PathBuf,
    n_periods: u64,
}

fn parse_cli() -> Cli {
    let mut output_path: Option<PathBuf> = None;
    let mut n_periods: u64 = N_PERIODS;

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" | "-o" => {
                output_path = Some(PathBuf::from(
                    args.next().expect("--output requires a path argument"),
                ));
            }
            "--periods" => {
                n_periods = args
                    .next()
                    .expect("--periods requires a positive integer")
                    .parse()
                    .expect("--periods must be a positive integer");
                assert!(n_periods >= 1, "--periods must be ≥ 1");
            }
            other => panic!("unknown argument: {other}"),
        }
    }

    Cli {
        output_path: output_path
            .unwrap_or_else(|| PathBuf::from("validation/rebound-parity/figure8/out/apsis.csv")),
        n_periods,
    }
}
