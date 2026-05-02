//! REBOUND parity — Retrograde Kepler e=0.5, apsis IAS15 side.
//!
//! Mirror of [`rebound_parity_kepler.rs`] with exactly one IC sign-flip:
//! the secondary's tangential periapsis velocity is `-v_peri` instead of
//! `+v_peri`. Every other component, mass, length, and physical scale is
//! identical to the prograde test. This isolates sign-convention coverage
//! as the only experimental variable; the result tests whether the apsis
//! IAS15 inner loop, controller, and orbital-element bookkeeping handle
//! `L_z < 0` to the same precision as `L_z > 0`.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_retrograde -p apsis
//! cargo run --release --example rebound_parity_retrograde -p apsis -- --output path/to/apsis.csv
//! ```
//!
//! Default output path: `validation/rebound-parity/retrograde/out/apsis.csv`
//! (relative to the workspace root — `cargo run` should be invoked from
//! the repo root or from the orchestrator in that directory).
//!
//! ## Protocol
//!
//! The full protocol — initial conditions, integrator settings, two-horizon
//! design, magnitude/sign tier separation, and tolerances declared *a priori*
//! — is specified in
//! [`docs/experiments/2026-05-01-rebound-parity-retrograde.md`](../../../../docs/experiments/2026-05-01-rebound-parity-retrograde.md).
//!
//! Constants in this file (`A`, `E`, `M_PRIMARY`, `M_SECONDARY`, `N_ORBITS`)
//! are the protocol's IC values. Changes here are protocol changes — update
//! the notebook in lockstep.
//!
//! ## Two-horizon design
//!
//! `N_ORBITS = 10_000` is the long-horizon gate. The 100-orbit checkpoint
//! analysis is performed by the comparator on the first 101 samples of the
//! same CSV; no separate run is needed.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Semi-major axis (canonical units, G = 1).
const A: f64 = 1.0;
/// Eccentricity.
const E: f64 = 0.5;
/// Primary mass.
const M_PRIMARY: f64 = 1.0;
/// Secondary mass.
const M_SECONDARY: f64 = 1.0e-6;
/// Number of orbital periods integrated (long-horizon gate).
const N_ORBITS: u64 = 10_000;
/// Initial timestep, expressed as a fraction of the orbital period.
const DT_FRACTION_OF_PERIOD: f64 = 1.0e-3;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions ──────────────────────────────────────────────── //
    //
    // Place the secondary at periapsis with tangential velocity `-v_peri`
    // (retrograde direction); place the primary so the centre of mass sits
    // at the origin and the system's net momentum is zero. The single
    // change vs prograde is the sign of `secondary_vy` and consequently of
    // `primary_vy` (which carries the COM-zero counter-momentum).
    //
    // Sign-flipping a correctly-rounded f64 is an exact bit-level operation
    // (toggles the sign bit, IEEE-754), so `secondary_vy` is bit-identical
    // to the negation of its prograde counterpart on any IEEE-754 platform.
    let r_peri = A * (1.0 - E);
    let v_peri = ((1.0 + E) / (A * (1.0 - E))).sqrt();

    let m_total = M_PRIMARY + M_SECONDARY;
    let primary_x = -(M_SECONDARY / m_total) * r_peri;
    let primary_vy = (M_SECONDARY / m_total) * v_peri; // sign flipped vs prograde
    let secondary_x = (M_PRIMARY / m_total) * r_peri;
    let secondary_vy = -(M_PRIMARY / m_total) * v_peri; // sign flipped vs prograde

    let primary =
        Body::star(M_PRIMARY).at(primary_x, 0.0).with_velocity(0.0, primary_vy).unsoftened();
    let secondary =
        Body::rocky(M_SECONDARY).at(secondary_x, 0.0).with_velocity(0.0, secondary_vy).unsoftened();

    // ── Integrator setup ────────────────────────────────────────────────── //
    let period = 2.0 * std::f64::consts::PI;
    let dt0 = period * DT_FRACTION_OF_PERIOD;

    let mut sys = System::new(vec![primary, secondary], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(dt0);

    // ── CSV output ──────────────────────────────────────────────────────── //
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity — Retrograde Kepler e=0.5 — apsis IAS15 side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-05-01-rebound-parity-retrograde.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis)").unwrap();
    writeln!(w, "# units: canonical (G = 1)").unwrap();
    writeln!(w, "# a={A}, e={E}, m_primary={M_PRIMARY}, m_secondary={M_SECONDARY:e}").unwrap();
    writeln!(w, "# orbit_sense: retrograde (secondary_vy = -v_peri)").unwrap();
    writeln!(w, "# period={period:.18e}, dt0={dt0:.18e}, n_orbits={N_ORBITS}").unwrap();
    writeln!(w, "orbit,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,e_total").unwrap();

    write_sample(&mut w, 0, &sys);
    for orbit in 1..=N_ORBITS {
        let t_target = period * (orbit as f64);
        sys.integrate_until(t_target);
        write_sample(&mut w, orbit, &sys);

        // Periodic progress report — without this, a 10^4-orbit run is silent
        // for ~10 seconds. Same cadence as a noticeable progress signal but
        // fast enough to be ignorable in a log.
        if orbit % 1000 == 0 {
            eprintln!("  apsis: {orbit}/{N_ORBITS} orbits completed (t={:.3e})", sys.t());
        }
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {}", N_ORBITS + 1, output_path.display());
}

// ── Output helper ───────────────────────────────────────────────────────── //

fn write_sample(w: &mut BufWriter<File>, orbit: u64, sys: &System) {
    let bodies = sys.bodies();
    let b0 = &bodies[0];
    let b1 = &bodies[1];
    let e_total = total_energy(bodies);
    writeln!(
        w,
        "{orbit},{t:.18e},{x0:.18e},{y0:.18e},{vx0:.18e},{vy0:.18e},{x1:.18e},{y1:.18e},{vx1:.18e},{vy1:.18e},{e:.18e}",
        t = sys.t(),
        x0 = b0.x,
        y0 = b0.y,
        vx0 = b0.vx,
        vy0 = b0.vy,
        x1 = b1.x,
        y1 = b1.y,
        vx1 = b1.vx,
        vy1 = b1.vy,
        e = e_total,
    )
    .unwrap();
}

/// Total mechanical energy for the two-body unsoftened system. Mirrors the
/// Kepler-prograde harness's `total_energy`; identical formula since the
/// sign of `v` does not affect the kinetic-energy or potential-energy form.
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

fn parse_output_path() -> PathBuf {
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" || arg == "-o" {
            return PathBuf::from(args.next().expect("--output requires a path argument"));
        }
    }
    PathBuf::from("validation/rebound-parity/retrograde/out/apsis.csv")
}
