//! Validation gate: PR drag energy-loss rate vs Burns 1979 analytic.
//!
//! For a circular orbit at radius `a` around a source of mass `M`, with
//! receiver radiation-to-gravity ratio `β`, the Poynting–Robertson drag
//! force at the initial epoch is purely tangential (v_r = 0):
//!
//! ```text
//!   |a_PR| = β · G · M / (r² · c) · v
//! ```
//!
//! The instantaneous power dissipated per unit mass is `a_PR · v`, so
//! the specific orbital energy drifts at
//!
//! ```text
//!   dE/dt = -β · G · M · v² / (r² · c)        [circular initial state]
//! ```
//!
//! Integrating for `T` time units of approximately-circular motion:
//!
//! ```text
//!   ΔE_analytic ≈ -β · G · M · v² · T / (r² · c)
//! ```
//!
//! IAS15 + this operator must reproduce the analytic ΔE to better than
//! the loosest tolerance the orbit-decay correction allows over the
//! integration window. We pick a short window (10 orbits at high β)
//! where the orbit drifts only ~0.5 %; the analytic constant-r
//! approximation is then good to a few percent and IAS15's own error
//! is many orders below the signal.
//!
//! Reference: Burns, J. A., Lamy, P. L., & Soter, S. (1979). Radiation
//! forces on small particles in the solar system. Icarus 40, 1–48.

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_radiation::PoyntingRobertsonDrag;

/// Speed of light in solar canonical units (AU per Gaussian time unit).
/// Matches `apsis_1pn::C_SOLAR_UNITS` in derivation but recomputed
/// locally so this gate has no dependency on apsis-1pn. Built from
/// the operator's own `for_units` derivation so the test exercises
/// the same code path as production.
fn c_solar_gaussian() -> f64 {
    let units = UnitSystem::solar_canonical();
    const C_SI: f64 = 299_792_458.0;
    C_SI * units.time_scale_si() / units.length_scale_si()
}

/// Energy drift rate for a circular β=0.5 dust grain matches the Burns
/// 1979 analytic prediction within 5 % over 10 orbits.
///
/// 5 % is the regression bound. The analytic derivation assumes
/// constant r; the orbit drifts ~0.5 % inward over 10 orbits which
/// shifts the v² · GM / r² scaling by < 2 %. Empirical agreement
/// at the centre: 1.2 % (was 0.7 % before #133 — back-reaction
/// suppression at extreme mass ratios removed a partial cancellation
/// between spurious primary motion and the constant-r approximation).
#[test]
fn pr_drag_energy_drift_matches_burns_1979_to_5_percent() {
    let units = UnitSystem::solar_canonical();

    let beta = 0.5_f64;
    let r0 = 1.0_f64; // 1 AU
    let m_dust = 1e-15_f64; // mass irrelevant to per-unit-mass dynamics
    let v_k = 1.0_f64; // Keplerian velocity at 1 AU in G=M=1 units

    let sun = Body::star(1.0);
    let dust = Body::rocky(m_dust).at(r0, 0.0).with_velocity(0.0, v_k);

    let mut sys =
        System::new(vec![sun, dust], units).with_integrator(IntegratorKind::Ias15).with_dt(1e-3);

    sys.add_non_conservative_perturbation(Box::new(PoyntingRobertsonDrag::from_raw_betas(
        0,
        vec![0.0, beta],
        units,
    )))
    .expect("matching units; PR drag must register");

    // Populate the energy cache before measuring; sys.energy() reads
    // last_kinetic + last_potential which are 0 until the first step.
    sys.step();
    let e0 = sys.energy();
    let n_orbits = 10.0_f64;
    let t_orbit = 2.0 * std::f64::consts::PI; // Gaussian time per orbit at a = 1 AU
    let t_total = n_orbits * t_orbit;
    sys.integrate_for(t_total);
    let e1 = sys.energy();

    let de_observed = e1 - e0;
    // dE/dt = -β · G · M · v² · m_dust / (r² · c). Total ΔE multiplies
    // by t_total. G = M = 1 in canonical.
    let de_analytic = -beta * v_k * v_k * m_dust * t_total / (r0 * r0 * c_solar_gaussian());

    let rel = ((de_observed - de_analytic) / de_analytic).abs();
    eprintln!(
        "[burns-pr] observed = {de_observed:.6e}, analytic = {de_analytic:.6e}, rel = {rel:.4}"
    );
    assert!(
        rel < 0.05,
        "PR drag dE disagreement vs Burns 1979 analytic: \
         observed = {de_observed:.6e}, analytic = {de_analytic:.6e}, \
         relative error = {rel:.4} (gate: 0.05)",
    );
}

/// Sanity counter-test: with PR drag *not* registered, the same orbit
/// preserves energy at IAS15 floor. If this fails, the dE detected in
/// the gate above could be from integrator drift, not from PR drag.
#[test]
fn ias15_alone_preserves_energy_for_circular_baseline() {
    let units = UnitSystem::solar_canonical();
    let sun = Body::star(1.0);
    let dust = Body::rocky(1e-15).at(1.0, 0.0).with_velocity(0.0, 1.0);

    let mut sys =
        System::new(vec![sun, dust], units).with_integrator(IntegratorKind::Ias15).with_dt(1e-3);

    sys.step();
    let e0 = sys.energy();
    sys.integrate_for(10.0 * 2.0 * std::f64::consts::PI);
    let e1 = sys.energy();
    let rel = ((e1 - e0) / e0).abs();
    assert!(rel < 1e-12, "IAS15 must conserve energy on Keplerian baseline; got {rel:.2e}");
}
