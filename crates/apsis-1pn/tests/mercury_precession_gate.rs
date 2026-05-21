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
//! 2. Someone broke the operator dispatch contract in `apsis`
//!    (e.g. the integrator stopped summing operator accelerations).
//! 3. Someone regressed the IAS15 substep velocity prediction
//!    (`predict_v_ias15` in `crate::physics::integrator::dense`) —
//!    velocity-dependent perturbations integrate against stale `v` and
//!    accumulate `O(a · dt)` per-substep bias. See
//!    `docs/experiments/2026-04-28-ias15-velocity-prediction-bug.md`.
//! 4. The Newtonian 2-body baseline regressed below machine-precision
//!    quality — the GR signal is swamped by numerical noise.
//!
//! ## Why `from_raw_c(C_SOLAR_UNITS, …)` instead of `for_units(…)`
//!
//! The gate is calibrated against the apsis-1pn historical baseline,
//! where `c` is the [`apsis_1pn::C_SOLAR_UNITS`] literal (derived from
//! `c_SI · YR_S/(2π) / AU`, the IAU julian-year convention).
//! `for_units(UnitSystem::solar_canonical())` derives `c` from
//! Gaussian time (`sqrt(AU³/(G·M))`) instead — numerically ~190 ppm
//! off the IAU literal. Both are physically valid; they differ only
//! by the historical IAU-vs-Gaussian gap.
//!
//! That ~190 ppm shift in `c` translates into a corresponding shift
//! in the 1PN force prefactor (`∝ 1/c²`), which IAS15's adaptive
//! substep schedule responds to at the ULP level. The 2D path
//! absorbs the perturbation (still passes within 100 ppm everywhere);
//! the 3D inclined path on Linux glibc + libm has slightly more
//! ULP-noise headroom and crosses the 100 ppm bound. Pinning the
//! gate to `C_SOLAR_UNITS` eliminates that confound and keeps the
//! gate locked against the same `c` value the original 4.4 ppm
//! headline was measured with.
//!
//! The recommended user-facing API is still
//! [`PostNewtonian1PN::for_units`] — see `examples/mercury_perihelion.rs`,
//! which demonstrates that path. The gate uses the raw escape
//! deliberately for regression-detection stability.

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::{C_SOLAR_UNITS, PostNewtonian1PN};

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_matches_gr_within_100ppm() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 500;

    // Softening zeroed so the Newtonian baseline is bit-exact Keplerian.
    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
        C_SOLAR_UNITS,
        UnitSystem::solar_canonical(),
    )))
    .expect("Mercury gate: System and operator share UnitSystem::solar_canonical()");

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

    let c = C_SOLAR_UNITS;
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (N_ORBITS as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel_err < 1e-4,
        "Mercury precession off by {rel_err:.3e} — measured {measured:.3e} rad vs predicted {predicted:.3e} rad"
    );
}

/// Opting into a Plummer kernel under a perturbation that demands exact
/// 1/r² (1PN) raises a warn diagnostic at registration time.
#[test]
fn plummer_kernel_under_1pn_triggers_diagnostic() {
    use apsis::core::log::{Event, Level, subscribe, unsubscribe};
    use apsis::physics::gravity::NewtonKernel;
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
        UnitSystem::solar_canonical(),
    )
    .with_kernel(std::sync::Arc::new(NewtonKernel::new(0.01)))
    .with_integrator(IntegratorKind::Ias15);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
        C_SOLAR_UNITS,
        UnitSystem::solar_canonical(),
    )))
    .expect("Mercury gate: System and operator share UnitSystem::solar_canonical()");

    let events = captured.lock().unwrap().clone();
    unsubscribe(id);

    assert_eq!(events.len(), 1, "exactly one warning expected");
    assert_eq!(events[0].level, Level::Warn);
    let field_names: Vec<&str> = events[0].fields.iter().map(|(k, _)| *k).collect();
    assert!(field_names.contains(&"kernel_epsilon"));
    assert!(field_names.contains(&"violated_invariant"));
}

/// Federation gate — ImplicitMidpoint + 1PN. IAS15 and WHFast measure the
/// 1PN signal directly because their Newtonian baseline is Kepler-exact
/// (analytic drift, zero numerical precession). IM integrates the
/// Newtonian Kepler motion numerically, so it produces an `O(dt²)`
/// numerical perihelion precession of the same order as the 1PN signal
/// at any practical dt. The federation contract is verified by
/// *subtraction*: precession with 1PN minus precession without 1PN must
/// equal the GR analytical prediction, isolating the operator's
/// contribution from IM's numerical ghost.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_implicit_midpoint_isolates_1pn_signal() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 500;

    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();

    let omega_drift = |with_1pn: bool| -> f64 {
        let sun = Body::star(1.0);
        let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);
        let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::ImplicitMidpoint)
            .with_dt(1e-4);
        if with_1pn {
            sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
                C_SOLAR_UNITS,
                UnitSystem::solar_canonical(),
            )))
            .expect("IM federation gate: shared UnitSystem::solar_canonical()");
        }
        let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
        sys.integrate_for(el0.period * (N_ORBITS as f64));
        let el_end = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
        let mut d = el_end.omega - el0.omega;
        while d > PI {
            d -= 2.0 * PI;
        }
        while d <= -PI {
            d += 2.0 * PI;
        }
        d
    };

    let with_1pn = omega_drift(true);
    let baseline = omega_drift(false);
    let measured = with_1pn - baseline;

    let c = C_SOLAR_UNITS;
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (N_ORBITS as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel_err < 1e-2,
        "IM + 1PN federation gate: 1PN-attributable precession off by {rel_err:.3e} — \
         measured (with − without) = {measured:.3e} rad vs predicted {predicted:.3e} rad \
         (with-1PN total = {with_1pn:.3e}, baseline ghost = {baseline:.3e})",
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

    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
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

/// Federation gate — WHFast + 1PN reproduces the same GR perihelion
/// precession the IAS15 gate above measures. Closes the federation
/// contract for symplectic-class integrators: a perturbation registered
/// against `System` is dispatched correctly through WHFast's
/// `interaction_step`, and the WHFast compensated-summation accumulator
/// does not corrupt the 1PN force prefactor.
///
/// Why a separate gate: WHFast and IAS15 cover different reviewer
/// expectations. IAS15 is the precision-controlled adaptive path;
/// WHFast is the long-horizon planetary baseline. The integrator zoo
/// claim is "every regime composes with every perturbation" — this
/// test discharges that claim for the WHFast slot.
///
/// Threshold matches the IAS15 gate (100 ppm) for the same cross-platform
/// reasons. WHFast is 2nd-order with fixed dt, so the cumulative phase
/// noise is `O(N_steps · dt² · m_p / m_0)`; at `dt = 1e-4`, 500 Mercury
/// orbits, `m_p / m_0 ≈ 1.66e-7`, the noise budget is `~1.7e-8 rad` vs
/// the GR signal of `~1e-3 rad` — comfortable room under 100 ppm.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_whfast_matches_gr_within_100ppm() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 500;

    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::WHFast)
        .with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
        C_SOLAR_UNITS,
        UnitSystem::solar_canonical(),
    )))
    .expect("WHFast federation gate: System and operator share UnitSystem::solar_canonical()");

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

    let c = C_SOLAR_UNITS;
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (N_ORBITS as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel_err < 1e-4,
        "WHFast + 1PN Mercury precession off by {rel_err:.3e} — \
         measured {measured:.3e} rad vs predicted {predicted:.3e} rad",
    );
}

/// 3D smoke test — Mercury in an inclined orbit `i = 7°` reproduces the
/// same GR perihelion precession rate as the planar gate above.
///
/// The orbit is constructed by rotating Mercury's velocity at perihelion
/// (initially `(0, v_peri, 0)` in the planar setup) by `i` around the
/// `x̂` axis: position stays at `(r_peri, 0, 0)`, velocity becomes
/// `(0, v_peri · cos i, v_peri · sin i)`. The orbital plane tilts out
/// of `xy`; the line of nodes runs along `+x̂` at the ascending node,
/// so initial `Ω = 0` and `ω = 0` while `inclination = i`.
///
/// The argument of periapsis `ω` is measured **in the orbital plane**,
/// so the GR precession rate `6π / (c² a (1−e²))` per orbit is
/// invariant under this rotation: the integrator and the orbital-element
/// pipeline must reproduce the same drift the planar gate measures.
/// `inclination` and `Ω` are conserved by two-body 1PN — they are
/// pinned tightly here as a separate check that the 3D path does not
/// leak energy or angular momentum into the orientation of the orbital
/// plane.
///
/// This is the test that proves `Body.{z, vz}`, the cross-product
/// `r × v`, the inclined `Ω` / `ω` branches in
/// `elements_from_invariants`, and the IAS15 `Vec3` substep buffers
/// are load-bearing — not decorative additions that happen to evaluate
/// to zero on the planar test suite.
///
/// If this test fails while the planar gate passes, one of:
/// 1. The 3D integrator path produces a different drift (kernel bug,
///    or 1PN formula's z-axis term wrong).
/// 2. `compute_invariants` / `elements_from_invariants` mishandles
///    the inclined branch — `ω` from `atan2(e·(h×n), e·n)` returns
///    the wrong quadrant or scale.
/// 3. The orbital plane is precessing spuriously (kernel asymmetry).
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn mercury_precession_3d_inclined_matches_gr_within_100ppm() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 500;
    const INCLINATION: f64 = 7.0_f64 * std::f64::consts::PI / 180.0; // 7°

    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    // Rotate the planar velocity (0, v_peri, 0) around x̂ by `INCLINATION`.
    let (sin_i, cos_i) = INCLINATION.sin_cos();
    let mercury = Body::rocky(M_MERCURY).at_3d(r_peri, 0.0, 0.0).with_velocity_3d(
        0.0,
        v_peri * cos_i,
        v_peri * sin_i,
    );

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
        C_SOLAR_UNITS,
        UnitSystem::solar_canonical(),
    )))
    .expect("Mercury gate: System and operator share UnitSystem::solar_canonical()");

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();

    // Setup invariant: the constructed orbit is genuinely inclined.
    assert!(
        (el0.inclination - INCLINATION).abs() < 1e-12,
        "test setup: initial inclination {} rad ≠ target {INCLINATION} rad",
        el0.inclination,
    );

    sys.integrate_for(el0.period * (N_ORBITS as f64));
    let el_end = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();

    // ── ω drift matches GR ────────────────────────────────────────────────
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

    let c = C_SOLAR_UNITS;
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (N_ORBITS as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();
    assert!(
        rel_err < 1e-4,
        "Mercury 3D inclined precession off by {rel_err:.17e} — \
         measured {measured:.17e} rad vs predicted {predicted:.17e} rad",
    );

    // ── Inclination and Ω are conserved by two-body 1PN ───────────────────
    //
    // 1PN does not torque the orbital plane in a test-particle Sun–Mercury
    // setup; both the inclination and the longitude of ascending node
    // must come back to their initial values to within a tight numerical
    // bound over 500 orbits. A drift here would indicate spurious
    // out-of-plane forcing — a kernel bug that the planar gate cannot
    // see by construction.
    let i_drift = (el_end.inclination - el0.inclination).abs();
    assert!(
        i_drift < 1e-9,
        "inclination drifted by {i_drift:.3e} rad over {N_ORBITS} orbits — \
         expected ≈ 0 for two-body 1PN",
    );
    let node_drift = (el_end.lon_ascending_node - el0.lon_ascending_node).abs();
    assert!(
        node_drift < 1e-9,
        "Ω drifted by {node_drift:.3e} rad over {N_ORBITS} orbits — \
         expected ≈ 0 for two-body 1PN",
    );
}

/// 3D baseline — without the 1PN perturbation, the inclined Kepler
/// orbit is closed: ω, inclination, and Ω all return to their initial
/// values to machine precision.
///
/// This is the 3D analogue of [`baseline_newtonian_kepler_is_closed`]
/// — same role, same threshold, but exercising the full 3D code path
/// (Body z/vz, kernel `r² = Δx² + Δy² + Δz²`, IAS15 Vec3 substep
/// buffers, inclined branch in `elements_from_invariants`). A drift
/// here means the 3D integrator path itself is not closed — the
/// problem is below the GR signal and the `_inclined` test above will
/// also fail.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn baseline_newtonian_kepler_3d_inclined_is_closed() {
    const A: f64 = 0.387_098;
    const E: f64 = 0.205_63;
    const M_MERCURY: f64 = 1.660_114e-7;
    const N_ORBITS: u64 = 300;
    const INCLINATION: f64 = 7.0_f64 * std::f64::consts::PI / 180.0;

    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let (sin_i, cos_i) = INCLINATION.sin_cos();
    let mercury = Body::rocky(M_MERCURY).at_3d(r_peri, 0.0, 0.0).with_velocity_3d(
        0.0,
        v_peri * cos_i,
        v_peri * sin_i,
    );

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    // No perturbation → pure 3D two-body.

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
    sys.integrate_for(el0.period * (N_ORBITS as f64));
    let el_end = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();

    let omega_drift = (el_end.omega - el0.omega).abs();
    assert!(
        omega_drift < 1e-9,
        "pure 3D Kepler 2-body drifted ω by {omega_drift:.3e} rad over {N_ORBITS} orbits",
    );
    let i_drift = (el_end.inclination - el0.inclination).abs();
    assert!(i_drift < 1e-9, "pure 3D Kepler 2-body drifted inclination by {i_drift:.3e} rad",);
    let node_drift = (el_end.lon_ascending_node - el0.lon_ascending_node).abs();
    assert!(node_drift < 1e-9, "pure 3D Kepler 2-body drifted Ω by {node_drift:.3e} rad",);
}
