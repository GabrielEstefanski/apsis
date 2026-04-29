//! Integration test — release-mode gate on the Mercury 1PN vs GR claim.
//!
//! Runs the Sun–Mercury 1PN scenario over 500 orbits and asserts the
//! measured perihelion precession matches the closed-form GR prediction
//! within 100 ppm (`rel_err < 10⁻⁴`). The threshold absorbs the
//! cross-platform f64 / LLVM / libm variance observed between
//! developer hardware (Windows MSVC, ~1 ppm) and the CI runner
//! (Linux glibc, ~30 ppm); both numbers sit at the f64 noise floor of
//! the test-particle 1PN approximation, but the floor itself is
//! platform-dependent at the ULP level. The headline figure cited in
//! `README.md` and `paper.md` (~1 ppm) is the developer-hardware
//! achievement; the gate is the portable lower bound — anything above
//! 100 ppm is a regression class, not a platform difference.
//!
//! If this test fails, one of four things is true:
//!
//! 1. Someone changed the 1PN formula and broke the sign or coefficients.
//! 2. Someone broke the [`PerturbationForce`] contract in `apsis`
//!    (e.g. the integrator stopped summing perturbation accelerations).
//! 3. Someone regressed the IAS15 substep velocity prediction
//!    (`predict_v_ias15` in `crate::physics::integrator::dense`) —
//!    velocity-dependent perturbations integrate against stale `v` and
//!    accumulate `O(a · dt)` per-substep bias. See
//!    `docs/experiments/2026-04-28-ias15-velocity-prediction-bug.md`.
//! 4. The Newtonian 2-body baseline regressed below machine-precision
//!    quality — the GR signal is swamped by numerical noise.

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_matches_gr_within_100ppm() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
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
    assert!(
        rel_err < 1e-4,
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
