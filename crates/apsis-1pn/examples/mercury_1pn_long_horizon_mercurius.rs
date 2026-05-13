//! Long-horizon Mercury 1PN — apsis Mercurius + apsis-1pn side.
//!
//! Sun + Mercury (unsoftened) integrated under apsis Mercurius with the
//! `apsis-1pn` Schwarzschild test-particle GR correction registered, for
//! 1000 years (~4150 Mercury orbits) in the canonical-Hénon unit system.
//! Mirrors the IAS15 sibling example (`mercury_1pn_long_horizon_ias15.rs`)
//! and writes an identically-shaped CSV — the comparator pairs them
//! row-for-row by Mercury orbit.
//!
//! For `N = 2` (Sun + Mercury), Mercurius's encounter step never fires
//! (no planet-planet pairs at all), so the integrator reduces to its
//! WH-like outer step: K-weighted half-kick → Kepler drift → indirect
//! drift → K-weighted half-kick. The 1PN perturbation enters via
//! `Mercurius::interaction_step` (PR #86), which applies one full `dt`
//! of perturbation strength symmetrically split across the two τ/2
//! kicks. This is the same Strang-split position WH integrators have
//! used for smooth perturbations since 1991.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example mercury_1pn_long_horizon_mercurius -p apsis-1pn
//! cargo run --release --example mercury_1pn_long_horizon_mercurius -p apsis-1pn -- --output path/to/mercurius.csv
//! ```
//!
//! Default output path: `validation/mercury-1pn-long-horizon/out/mercurius.csv`.
//!
//! ## Protocol
//!
//! Initial conditions, integrator settings, metrics, and tolerances declared
//! *a priori* in
//! [`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`](../../../docs/experiments/2026-05-13-mercury-1pn-long-horizon.md).

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

// ── Protocol constants (mirrored in the notebook + IAS15 sibling) ───────── //

const A: f64 = 0.387_098;
const E: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;
const N_ORBITS: u64 = 4153;
/// Outer dt for the Mercurius fixed-step (canonical Hénon).
///
/// Mercury period is `T_M ≈ 1.513` canonical time units; `DT = 0.01`
/// gives ~151 outer steps per orbit, well-resolved for a 2nd-order
/// symplectic-class scheme. Total outer steps over 1000 yr: ~628 000.
const DT: f64 = 1.0e-2;

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions (identical to IAS15 sibling, bit-for-bit) ───── //
    let sun = Body::star(1.0).unsoftened();
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury =
        Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri).unsoftened();

    // ── System setup ────────────────────────────────────────────────────── //
    let mut sys = System::new(vec![sun, mercury], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Mercurius)
        .with_dt(DT);
    sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0)
        .expect("Mercury IC must produce a bound Keplerian orbit");
    let period = el0.period;

    // ── CSV output ──────────────────────────────────────────────────────── //
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# Long-horizon Mercury 1PN — apsis Mercurius side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-05-13-mercury-1pn-long-horizon.md").unwrap();
    writeln!(w, "# integrator: Mercurius (apsis), perturbation: PostNewtonian1PN::solar_units()").unwrap();
    writeln!(w, "# units: canonical Hénon (G = 1)").unwrap();
    writeln!(w, "# a={A}, e={E}, m_mercury={M_MERCURY:e}").unwrap();
    writeln!(w, "# period={period:.18e}, dt={DT:.18e}, n_orbits={N_ORBITS}").unwrap();
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

fn parse_output_path() -> PathBuf {
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" || arg == "-o" {
            return PathBuf::from(args.next().expect("--output requires a path argument"));
        }
    }
    PathBuf::from("validation/mercury-1pn-long-horizon/out/mercurius.csv")
}
