//! REBOUND parity — Kepler e=0.5, apsis IAS15 side.
//!
//! Runs a canonical Kepler two-body orbit at eccentricity 0.5 under apsis
//! IAS15 for 100 orbital periods, sampling state and total energy at the end
//! of each orbit. Output is a CSV consumable by the matching Python REBOUND
//! harness for cross-implementation parity comparison.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_kepler -p apsis
//! cargo run --release --example rebound_parity_kepler -p apsis -- --output path/to/apsis.csv
//! ```
//!
//! Default output path: `validation/rebound-parity/kepler/out/apsis.csv`
//! (relative to the workspace root — `cargo run` should be invoked from
//! the repo root or from the orchestrator in that directory).
//!
//! ## Protocol
//!
//! The full protocol — initial conditions, integrator settings, metrics, and
//! tolerances declared *a priori* — is specified in
//! [`docs/experiments/2026-04-25-rebound-parity-kepler.md`](../../../../docs/experiments/2026-04-25-rebound-parity-kepler.md).
//!
//! Constants in this file (`A`, `E`, `M_PRIMARY`, `M_SECONDARY`, `N_ORBITS`)
//! are the protocol's IC values. Changes here are protocol changes — update
//! the notebook in lockstep.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Semi-major axis (canonical units, G = 1).
const A: f64 = 1.0;
/// Eccentricity.
const E: f64 = 0.5;
/// Primary mass.
const M_PRIMARY: f64 = 1.0;
/// Secondary mass.
const M_SECONDARY: f64 = 1.0e-6;
/// Number of orbital periods integrated.
const N_ORBITS: u64 = 100;
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
    // Place the secondary at periapsis; place the primary so the centre of
    // mass sits at the origin and the system's net momentum is zero. This
    // eliminates COM drift as a source of `Δr` in the cross-implementation
    // comparison.
    //
    // The relative-motion vis-viva at periapsis is
    //
    //     v_rel = sqrt((1 + e) / (a (1 − e)))            (G M_total ≈ 1)
    //
    // matching the notebook's stated formula. Both implementations must use
    // the f64 evaluation of this exact expression to obtain bit-identical ICs.
    let r_peri = A * (1.0 - E);
    let v_peri = ((1.0 + E) / (A * (1.0 - E))).sqrt();

    let m_total = M_PRIMARY + M_SECONDARY;
    let primary_x = -(M_SECONDARY / m_total) * r_peri;
    let primary_vy = -(M_SECONDARY / m_total) * v_peri;
    let secondary_x = (M_PRIMARY / m_total) * r_peri;
    let secondary_vy = (M_PRIMARY / m_total) * v_peri;

    let primary =
        Body::star(M_PRIMARY).at(primary_x, 0.0).with_velocity(0.0, primary_vy).unsoftened();
    let secondary =
        Body::rocky(M_SECONDARY).at(secondary_x, 0.0).with_velocity(0.0, secondary_vy).unsoftened();

    // ── Integrator setup ────────────────────────────────────────────────── //
    //
    // Period in canonical units (G = 1, M ≈ 1, a = 1) is 2π. Initial dt is
    // a fixed fraction of the period; IAS15's adaptive controller takes over
    // from there. The integrator–force-model pairing rule auto-forces direct
    // O(N²) when IAS15 is selected (see ADR-003).
    let period = 2.0 * std::f64::consts::PI;
    let dt0 = period * DT_FRACTION_OF_PERIOD;

    let mut sys =
        System::new(vec![primary, secondary]).with_integrator(IntegratorKind::Ias15).with_dt(dt0);

    // ── CSV output ──────────────────────────────────────────────────────── //
    //
    // Wide format: one row per (orbit completion) sample, including all body
    // state and total energy. Initial state at orbit = 0; subsequent samples
    // at orbit completion. Total = N_ORBITS + 1 = 101 samples.
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity — Kepler e=0.5 — apsis IAS15 side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-04-25-rebound-parity-kepler.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis)").unwrap();
    writeln!(w, "# units: canonical (G = 1)").unwrap();
    writeln!(w, "# a={A}, e={E}, m_primary={M_PRIMARY}, m_secondary={M_SECONDARY:e}").unwrap();
    writeln!(w, "# period={period:.18e}, dt0={dt0:.18e}, n_orbits={N_ORBITS}").unwrap();
    writeln!(w, "orbit,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,e_total").unwrap();

    write_sample(&mut w, 0, &sys);
    for orbit in 1..=N_ORBITS {
        let t_target = period * (orbit as f64);
        sys.integrate_until(t_target);
        write_sample(&mut w, orbit, &sys);
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

/// Total mechanical energy for the two-body unsoftened system.
///
/// Computed inline rather than through `apsis::physics::energy` so the
/// formula is visible at the comparison site and matches REBOUND's
/// `sim.calculate_energy()` convention exactly: KE = ½ Σ mᵢ vᵢ², PE =
/// −Σᵢ<ⱼ G mᵢ mⱼ / rᵢⱼ, with G = 1 and no softening (verified by
/// `Body::unsoftened()` on every body in this configuration).
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
    PathBuf::from("validation/rebound-parity/kepler/out/apsis.csv")
}
