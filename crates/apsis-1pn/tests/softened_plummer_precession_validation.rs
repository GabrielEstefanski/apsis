//! Validation harness for the softened-Plummer apsidal-precession
//! counter-test (paper §3.2). Pairs the satisfied 1PN case from
//! `mercury_precession_gate.rs` with the violated case here — the
//! exactness-requiring 1PN operator attached to a softened kernel.
//!
//! The measured drift is validated against the softened-Plummer apsidal
//! precession from a full-potential apsidal-angle quadrature (numerically
//! converged, independent of apsis; see the §3.2 derivation notebook).

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::NewtonKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::{C_SOLAR_UNITS, PostNewtonian1PN};

const A_MERCURY: f64 = 0.387_098;
const E_MERCURY: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;
const N_ORBITS: u64 = 500;
const EPSILON_VIOLATED: f64 = 0.02;

const T_MERCURY_DAYS: f64 = 87.969;
const DAYS_PER_CENTURY: f64 = 36_525.0;
const RAD_TO_ARCSEC: f64 = 180.0 * 3600.0 / PI;

fn orbits_per_century() -> f64 {
    DAYS_PER_CENTURY / T_MERCURY_DAYS
}

/// Leading-order closed form  Δϖ_orbit = -3π ε² / [a²(1-e²)²]  (paper
/// §3.2), scaled to arcseconds per Earth century. Reported alongside the
/// quadrature reference for context only: it is +2.66 % high here, a
/// quantified O(ε²) next-order term, so it is not the gate oracle.
fn predicted_softened_plummer_drift_arcsec_per_century(a_au: f64, e: f64, epsilon_au: f64) -> f64 {
    let per_orbit_rad = -3.0 * PI * epsilon_au.powi(2) / (a_au.powi(2) * (1.0 - e.powi(2)).powi(2));
    per_orbit_rad * RAD_TO_ARCSEC * orbits_per_century()
}

/// Softened-Plummer apsidal precession per Kepler period (rad) for this
/// orbit, from a full-potential apsidal-angle quadrature — numerically
/// converged (Gauss–Legendre), independent of apsis. The gate oracle: it
/// carries all orders in ε, where the closed form above is leading order
/// only. Derived in the §3.2 Plummer-apsidal notebook (its ε=0 self-check
/// returns 1.7e-11 rad, i.e. the Kepler orbit closes).
const QUADRATURE_DRIFT_PER_KEPLER_RAD: f64 = -2.671_786_35e-2;

/// Drift of the osculating periapsis longitude under a softened
/// `NewtonKernel(ε)` plus the 1PN operator, accumulated per Kepler period
/// over `N_ORBITS`, returned in arcsec/century.
///
/// 1PN is registered deliberately: this is the counter-test for attaching
/// the exactness-requiring 1PN operator to a softened kernel. The softening
/// artifact dominates the relativistic signal by ~5e4, so 1PN contributes
/// ~2e-5 of the measured drift (below the gate) — the comparison to the
/// pure-softened-Plummer quadrature oracle holds.
///
/// Per-orbit unwrap: the end-vs-initial mod-2π form would alias the true
/// ~ -13.7 rad cumulative drift to a fractional value; the per-orbit step
/// (~2.7e-2 rad) stays well inside ±π, so the running sum is lossless.
fn measure_drift_arcsec_per_century(epsilon: f64) -> f64 {
    let sun = Body::star(1.0);
    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = (2.0 / r_peri - 1.0 / A_MERCURY).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_kernel(std::sync::Arc::new(NewtonKernel::new(epsilon)))
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
        C_SOLAR_UNITS,
        UnitSystem::solar_canonical(),
    )))
    .expect("violated case: System and operator share UnitSystem::solar_canonical()");

    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0).unwrap();
    let period = el0.period;

    let mut omega_prev = el0.omega;
    let mut drift_rad = 0.0_f64;
    for k in 1..=N_ORBITS {
        sys.integrate_until(period * (k as f64));
        let el_k = compute_elements(sys.bodies(), 1, 0, 1.0)
            .expect("Mercury orbit must stay bound under softened kernel + 1PN");
        let mut step = el_k.omega - omega_prev;
        while step > PI {
            step -= 2.0 * PI;
        }
        while step <= -PI {
            step += 2.0 * PI;
        }
        drift_rad += step;
        omega_prev = el_k.omega;
    }

    drift_rad * RAD_TO_ARCSEC * orbits_per_century() / (N_ORBITS as f64)
}

#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn plummer_violated_case_drift_is_measured() {
    let measured = measure_drift_arcsec_per_century(EPSILON_VIOLATED);
    eprintln!(
        "[softened-precession] measured drift at eps={EPSILON_VIOLATED} AU: \
         {measured:.0} arcsec/century",
    );
    assert!(
        measured.abs() > 1_000.0,
        "violated case drift too small ({measured:.3e} arcsec/century) — \
         scenario or kernel setup regressed",
    );
}

/// Gate: apsis reproduces the quadrature apsidal precession to ~0.04 %
/// here (osculating-ω vs geometric-apsis residual at 500 orbits, IAS15
/// dt=1e-4). Bound at 0.5 % keeps ~10x headroom — far tighter than the
/// prior 5 % closed-form comparison, which only existed to absorb the
/// closed form's own +2.66 % leading-order truncation.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn plummer_drift_matches_quadrature() {
    let measured = measure_drift_arcsec_per_century(EPSILON_VIOLATED);
    let quadrature = QUADRATURE_DRIFT_PER_KEPLER_RAD * RAD_TO_ARCSEC * orbits_per_century();
    let closed_form =
        predicted_softened_plummer_drift_arcsec_per_century(A_MERCURY, E_MERCURY, EPSILON_VIOLATED);

    let rel_err = (measured - quadrature).abs() / quadrature.abs();
    let closed_pct = (closed_form / quadrature - 1.0) * 100.0;

    const ACCEPTANCE: f64 = 5e-3;

    eprintln!(
        "[softened-precession] eps={EPSILON_VIOLATED} AU:\n  \
         measured    = {measured:.6e} arcsec/century\n  \
         quadrature  = {quadrature:.6e} arcsec/century\n  \
         closed form = {closed_form:.6e} arcsec/century ({closed_pct:+.2} % vs quadrature, leading-order)\n  \
         rel_err     = {rel_err:.4e} (gate: {ACCEPTANCE})",
    );

    assert!(
        rel_err < ACCEPTANCE,
        "softened-Plummer apsidal drift mismatch vs quadrature: \
         measured {measured:.3e} arcsec/century vs quadrature {quadrature:.3e}; \
         rel_err = {rel_err:.3e} exceeds {ACCEPTANCE}",
    );
}
