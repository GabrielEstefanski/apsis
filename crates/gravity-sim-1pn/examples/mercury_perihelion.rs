//! Mercury perihelion precession — the canonical 1PN test.
//!
//! Run with:
//!
//! ```text
//! cargo run --example mercury_perihelion --release -p gravity-sim-1pn
//! ```
//!
//! Integrates Sun + Mercury under Newtonian gravity + the 1PN correction
//! supplied by this crate, measures the drift of the argument of perihelion
//! over 500 Mercury orbits, and compares the rate against the GR prediction
//!
//! ```text
//!   Δω / orbit = 6π G M / (c² a (1 − e²))
//! ```
//!
//! which integrates to 43 arcseconds per century for the real Mercury.
//!
//! The softening length on both bodies is zeroed at construction: the
//! default material-scaled Plummer softening introduces a spurious
//! (non-GR) apsidal precession that swamps the 1PN signal by orders of
//! magnitude. A test of a general-relativistic correction needs the
//! exact 1/r potential as the baseline.

use gravity_sim_core::core::system::System;
use gravity_sim_core::domain::body::Body;
use gravity_sim_core::physics::integrator::IntegratorKind;
use gravity_sim_core::physics::orbital::compute_elements;
use gravity_sim_1pn::PostNewtonian1PN;

use std::f64::consts::PI;

/// Orbital parameters for Mercury in the simulator's canonical units.
const A_MERCURY: f64 = 0.387_098; // semi-major axis in AU
const E_MERCURY: f64 = 0.205_63; // eccentricity
const M_MERCURY: f64 = 1.660_114e-7; // Mercury / Sun mass ratio
const M_SUN: f64 = 1.0;

/// How many Mercury orbits to integrate. 500 gives ~52 arcsec of accumulated
/// precession — two decades above numerical noise while staying fast enough
/// for a one-minute release-mode run.
const N_ORBITS: u64 = 500;

fn main() {
    // ── Initial conditions ──────────────────────────────────────────────────
    let mut sun = Body::star(M_SUN);
    sun.softening = 0.0;

    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = (M_SUN * (2.0 / r_peri - 1.0 / A_MERCURY)).sqrt();
    let mut mercury = Body::rocky(M_MERCURY)
        .at(r_peri, 0.0)
        .with_velocity(0.0, v_peri);
    mercury.softening = 0.0;

    // ── Build the simulation ────────────────────────────────────────────────
    let mut sys = System::new(vec![sun, mercury])
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);

    // Attach the out-of-tree 1PN perturbation. Everything below this line
    // uses only the public API of `gravity-sim-core`; `gravity-sim-1pn` has
    // no other dependency on the workspace. This is the Phase 3 gate.
    sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));

    // ── Reference state at t = 0 ────────────────────────────────────────────
    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0)
        .expect("Mercury must produce bound elements at t = 0");
    let mercury_period = el0.period;

    // ── Integrate and sample ────────────────────────────────────────────────
    let t_end = mercury_period * (N_ORBITS as f64);
    println!("Mercury + Sun + 1PN @ IAS15");
    println!("  T_mercury      = {mercury_period:.6} sim units");
    println!("  integrating    = {N_ORBITS} orbits  →  t = {t_end:.2}");
    println!();

    println!("{:>6}  {:>14}  {:>14}  {:>14}", "orbit", "Δω (rad)", "Δω (arcsec)", "|δE/E|");

    let sample_every = N_ORBITS / 10;
    for k in 1..=N_ORBITS {
        sys.integrate_until(mercury_period * (k as f64));

        if k % sample_every == 0 || k == N_ORBITS {
            let el = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
            let d_omega = unwrap_angle(el.omega - el0.omega);
            let arcsec = d_omega.to_degrees() * 3600.0;
            println!(
                "{:>6}  {:>+14.6e}  {:>+14.4}  {:>14.3e}",
                k,
                d_omega,
                arcsec,
                sys.energy_delta().abs(),
            );
        }
    }

    // ── Compare against GR prediction ───────────────────────────────────────
    //
    // Schwarzschild perihelion advance per orbit:
    //     Δω = 6π G M / (c² a (1 − e²))
    let c = PostNewtonian1PN::solar_units().c();
    let predicted_per_orbit =
        6.0 * PI * M_SUN / (c * c * A_MERCURY * (1.0 - E_MERCURY * E_MERCURY));
    let predicted_total = predicted_per_orbit * (N_ORBITS as f64);

    let el_final = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
    let measured_total = unwrap_angle(el_final.omega - el0.omega);

    let predicted_arcsec = predicted_total.to_degrees() * 3600.0;
    let measured_arcsec = measured_total.to_degrees() * 3600.0;
    let rel_err = (measured_total - predicted_total) / predicted_total;

    // Observable rate: arcseconds per century.
    //   1 simulation year = 2π sim time units; 1 century = 200π sim time.
    let t_centuries = sys.t() / (200.0 * PI);
    let arcsec_per_century = measured_arcsec / t_centuries;

    println!();
    println!("── GR comparison over {N_ORBITS} orbits ──");
    println!("  predicted Δω      = {predicted_total:+.6e} rad  ({predicted_arcsec:+.4} arcsec)");
    println!("  measured  Δω      = {measured_total:+.6e} rad  ({measured_arcsec:+.4} arcsec)");
    println!("  relative error    = {rel_err:+.3e}");
    println!("  rate              = {arcsec_per_century:.3} arcsec/century  (GR expects 43)");
}

fn unwrap_angle(d: f64) -> f64 {
    let mut x = d;
    while x > PI {
        x -= 2.0 * PI;
    }
    while x <= -PI {
        x += 2.0 * PI;
    }
    x
}
