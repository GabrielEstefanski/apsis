//! Integration test — paper-grade gate on the headline 1PN validation
//! number.
//!
//! Runs the Mercury perihelion scenario at release-level fidelity and
//! asserts the measured precession matches GR within 10 ppm
//! (`rel_err < 10⁻⁵`). The threshold is the CI-protected number cited
//! in `README.md` and `paper.md` — a regression past 10 ppm degrades
//! the paper claim and must fail the build, not a soft warning.
//!
//! If this test fails, one of four things is true:
//!
//! 1. Someone changed the 1PN formula and broke the sign or coefficients.
//! 2. Someone broke the [`PerturbationForce`] contract in `apsis`
//!    (e.g. the integrator stopped summing perturbation accelerations).
//! 3. Someone regressed the IAS15 substep velocity prediction
//!    (`predict_v_ias15` in `crate::physics::integrator::dense`) —
//!    velocity-dependent perturbations integrate against stale `v` and
//!    accumulate `O(a · dt)` per-substep bias.
//!    See `docs/experiments/2026-04-28-ias15-velocity-prediction-bug.md`.
//! 4. The simulator's Newtonian 2-body integration regressed below
//!    machine-precision quality — the GR signal is swamped by
//!    numerical noise.
//!
//! All four failure modes are things a reviewer of the paper would want
//! caught automatically.

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_matches_gr_within_10ppm() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    // 500 orbits matches the README / paper.md regime. At shorter runs
    // the absolute error floor (f64 round-off in osculating-element
    // extraction) dominates the ratio; the prior 300-orbit form gave
    // `rel_err ≈ 10⁻⁴`, numerically above the tightened threshold.
    const N_ORBITS: u64 = 500;

    // Softening zeroed so the Newtonian baseline is bit-exact Keplerian.
    let sun = Body::star(1.0).unsoftened();
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri).unsoftened();

    let mut sys = System::new(vec![sun, mercury], UnitSystem::canonical())
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
    // 10 ppm. Margin: ~10× over the achieved ~1 ppm post-fix.
    assert!(
        rel_err < 1e-5,
        "Mercury precession off by {rel_err:.3e} — measured {measured:.3e} rad vs predicted {predicted:.3e} rad"
    );
}

/// Contract test — registering a 1PN perturbation into a softened
/// system must raise a warn-level diagnostic on the log bus. This is
/// the protection against the silent Plummer-swamps-GR failure mode
/// that tripped up the first end-to-end run.
#[test]
fn softened_system_triggers_diagnostic() {
    use apsis::core::log::{Event, Level, subscribe, unsubscribe};
    use std::sync::{Arc, Mutex};

    const MARKER: &str = "perturbation requires exact 1/r gravity";

    let captured: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = captured.clone();
    let id = subscribe(move |event: &Event| {
        if event.message.starts_with(MARKER) {
            sink.lock().unwrap().push(event.clone());
        }
    });

    // Default softening left in place — this is the trap.
    let mut sys = System::new(
        vec![Body::star(1.0), Body::rocky(1e-7).at(0.4, 0.0).with_velocity(0.0, 1.5)],
        UnitSystem::canonical(),
    )
    .with_integrator(IntegratorKind::Ias15);
    sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));

    let events = captured.lock().unwrap().clone();
    unsubscribe(id);

    assert_eq!(events.len(), 1, "exactly one warning expected");
    assert_eq!(events[0].level, Level::Warn);
    // Fields are present and describe the softening state.
    let field_names: Vec<&str> = events[0].fields.iter().map(|(k, _)| *k).collect();
    assert!(field_names.contains(&"softened_bodies"));
    assert!(field_names.contains(&"max_softening"));
}

/// Counterpart of the above: when the system is properly unsoftened,
/// registering the same perturbation must stay silent.
#[test]
fn exact_gravity_system_stays_silent() {
    use apsis::core::log::{Event, subscribe, unsubscribe};
    use std::sync::{Arc, Mutex};

    const MARKER: &str = "perturbation requires exact 1/r gravity";

    let captured: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
    let sink = captured.clone();
    let id = subscribe(move |event: &Event| {
        if event.message.starts_with(MARKER) {
            sink.lock().unwrap().push(event.clone());
        }
    });

    let mut sys = System::new(
        vec![Body::star(1.0), Body::rocky(1e-7).at(0.4, 0.0).with_velocity(0.0, 1.5)],
        UnitSystem::canonical(),
    )
    .with_exact_gravity()
    .with_integrator(IntegratorKind::Ias15);
    sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));

    let events = captured.lock().unwrap().clone();
    unsubscribe(id);

    assert!(
        events.is_empty(),
        "no warning expected for fully-unsoftened system, got {}",
        events.len()
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

    let mut sys = System::new(vec![sun, mercury], UnitSystem::canonical())
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
