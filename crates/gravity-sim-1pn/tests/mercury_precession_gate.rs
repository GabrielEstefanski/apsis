//! Integration test — Phase 3 gate for the out-of-tree plugin claim.
//!
//! Runs the Mercury perihelion scenario at release-level fidelity and
//! asserts the measured precession matches GR within 1 % relative error.
//!
//! If this test fails, one of three things is true:
//!
//! 1. Someone changed the 1PN formula and broke the sign or coefficients.
//! 2. Someone broke the [`PerturbationForce`] contract in `gravity-sim-core`
//!    (e.g. the integrator stopped summing perturbation accelerations).
//! 3. The simulator's Newtonian 2-body integration regressed below
//!    machine-precision quality — in which case the GR signal gets
//!    swamped by numerical noise.
//!
//! All three failure modes are things a reviewer of the paper would want
//! caught automatically.

use std::f64::consts::PI;

use gravity_sim_1pn::PostNewtonian1PN;
use gravity_sim_core::core::system::System;
use gravity_sim_core::domain::body::Body;
use gravity_sim_core::physics::integrator::IntegratorKind;
use gravity_sim_core::physics::orbital::compute_elements;

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_matches_gr_within_one_percent() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 300;

    // Softening zeroed so the Newtonian baseline is bit-exact Keplerian.
    let mut sun = Body::star(1.0);
    sun.softening = 0.0;
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mut mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);
    mercury.softening = 0.0;

    let mut sys = System::new(vec![sun, mercury])
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
    sys.integrate_for(el0.period * (N_ORBITS as f64));
    let el_end = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();

    let measured = {
        let mut d = el_end.omega - el0.omega;
        while d > PI {
            d -= 2.0 * PI;
        }
        while d <= -PI {
            d += 2.0 * PI;
        }
        d
    };

    let c = PostNewtonian1PN::solar_units().c();
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (N_ORBITS as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel_err < 1e-2,
        "Mercury precession off by {rel_err:.3e} — measured {measured:.3e} rad vs predicted {predicted:.3e} rad"
    );
}

/// Sanity test — without the 1PN perturbation, the same integration
/// must produce zero precession up to machine precision. Locks in the
/// baseline that the PN measurement relies on.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn baseline_newtonian_kepler_is_closed() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 300;

    let mut sun = Body::star(1.0);
    sun.softening = 0.0;
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mut mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);
    mercury.softening = 0.0;

    let mut sys = System::new(vec![sun, mercury])
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    // No perturbation attached → pure Keplerian 2-body.

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
    sys.integrate_for(el0.period * (N_ORBITS as f64));
    let el_end = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();

    let drift = (el_end.omega - el0.omega).abs();
    assert!(
        drift < 1e-9,
        "pure Kepler 2-body drifted ω by {drift:.3e} rad over {N_ORBITS} orbits — \
         baseline integration degraded"
    );
}
