//! Long-horizon Mercury 1PN — apsis IAS15 + apsis-1pn side.
//!
//! Sun + Mercury (unsoftened) integrated under apsis IAS15 with the
//! `apsis-1pn` Schwarzschild test-particle GR correction registered, for
//! 1000 years (~4150 Mercury orbits) in the canonical-Hénon unit system
//! (`G = 1`, time = year / (2π), length = AU, mass = M☉). Output is one
//! row per completed Mercury orbit, recording state and osculating
//! orbital elements; consumed by `validation/mercury-1pn-long-horizon/`
//! comparator to extract cumulative `Δω` and gate it against the
//! closed-form Schwarzschild prediction `6πGM/(c²a(1−e²))` per orbit.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example mercury_1pn_long_horizon_ias15 -p apsis
//! cargo run --release --example mercury_1pn_long_horizon_ias15 -p apsis -- --output path/to/ias15.csv
//! ```
//!
//! Default output path: `validation/mercury-1pn-long-horizon/out/ias15.csv`
//! (relative to the workspace root).
//!
//! ## Protocol
//!
//! Initial conditions, integrator settings, metrics, and tolerances declared
//! *a priori* in
//! [`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`](../../../docs/experiments/2026-05-13-mercury-1pn-long-horizon.md).
//! Constants in this file (`A`, `E`, `M_MERCURY`, `N_ORBITS`, `DT0`) are
//! the protocol's IC values — changes here are protocol changes; update the
//! notebook in lockstep.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Mercury semi-major axis (AU; numerically the same in canonical Hénon).
const A: f64 = 0.387_098;
/// Mercury eccentricity.
const E: f64 = 0.205_63;
/// Mercury / Sun mass ratio.
const M_MERCURY: f64 = 1.660_114e-7;
/// Number of full Mercury orbits to integrate.
///
/// Mercury period in canonical Hénon (`mu = 1`): `T = 2π · √(a³) ≈ 1.513`.
/// Total integration time at `N_ORBITS = 4153` is `N · T ≈ 6283.19`
/// canonical units = 1000 physical years.
const N_ORBITS: u64 = 4153;
/// First-call seed for the IAS15 controller (canonical Hénon time units).
/// Matches `crates/apsis-1pn/tests/mercury_precession_gate.rs`.
const DT0: f64 = 1.0e-4;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions ──────────────────────────────────────────────── //
    //
    // Same construction as the existing 500-orbit gate
    // (`crates/apsis-1pn/tests/mercury_precession_gate.rs`): Mercury at
    // periapsis along +x with the vis-viva tangent velocity for
    // `mu = G · (M_sun + M_mercury) ≈ 1` in canonical units. The Sun
    // sits at the origin with zero velocity; the COM is then at
    // `m_M / (1 + m_M) · r_peri ≈ 5 × 10⁻⁸ AU`, which `IAS15` recenters
    // into the COM frame on its first sub-step via the COM-shift hook.
    let sun = Body::star(1.0).unsoftened();
    let r_peri = A * (1.0 - E);
    // vis-viva: v² = mu · (2/r − 1/a); mu = 1 in canonical Hénon.
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri).unsoftened();

    // ── System setup ────────────────────────────────────────────────────── //
    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(DT0);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(
        UnitSystem::solar_canonical(),
    )));

    // Snapshot the t=0 osculating state so the output schema can include
    // the initial sample (orbit = 0) before any integration takes place.
    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0)
        .expect("Mercury IC must produce a bound Keplerian orbit");
    let period = el0.period;

    // ── CSV output ──────────────────────────────────────────────────────── //
    //
    // Wide format: one row per completed Mercury orbit, including state
    // and osculating elements. Initial state at orbit = 0; subsequent
    // samples at each orbital-period boundary. Total = N_ORBITS + 1.
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# Long-horizon Mercury 1PN — apsis IAS15 side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-05-13-mercury-1pn-long-horizon.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis), perturbation: PostNewtonian1PN::for_units(UnitSystem::solar_canonical())")
        .unwrap();
    writeln!(w, "# units: canonical Hénon (G = 1)").unwrap();
    writeln!(w, "# a={A}, e={E}, m_mercury={M_MERCURY:e}").unwrap();
    writeln!(w, "# period={period:.18e}, dt0={DT0:.18e}, n_orbits={N_ORBITS}").unwrap();
    writeln!(w, "orbit,t,x,y,vx,vy,a_osc,e_osc,omega_osc").unwrap();

    write_sample(&mut w, 0, &sys, &el0);

    for orbit in 1..=N_ORBITS {
        let t_target = period * (orbit as f64);
        sys.integrate_until(t_target);
        let el = compute_elements(sys.bodies(), 1, 0, 1.0)
            .expect("Mercury orbit should remain bound under 1PN over 1000 yr");
        write_sample(&mut w, orbit, &sys, &el);
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {}", N_ORBITS + 1, output_path.display());
}

// ── Output helper ───────────────────────────────────────────────────────── //

fn write_sample(
    w: &mut BufWriter<File>,
    orbit: u64,
    sys: &System,
    el: &apsis::physics::orbital::OrbitalElements,
) {
    let mercury = &sys.bodies()[1];
    writeln!(
        w,
        "{orbit},{t:.18e},{x:.18e},{y:.18e},{vx:.18e},{vy:.18e},{a:.18e},{e:.18e},{omega:.18e}",
        t = sys.t(),
        x = mercury.pos_x,
        y = mercury.pos_y,
        vx = mercury.vel_x,
        vy = mercury.vel_y,
        a = el.a,
        e = el.e,
        omega = el.omega,
    )
    .unwrap();
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

fn parse_output_path() -> PathBuf {
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" || arg == "-o" {
            return PathBuf::from(args.next().expect("--output requires a path argument"));
        }
    }
    PathBuf::from("validation/mercury-1pn-long-horizon/out/ias15.csv")
}
