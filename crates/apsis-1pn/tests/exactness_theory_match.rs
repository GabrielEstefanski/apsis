//! Theory-match harness for the Exactness counter-test.
//!
//! Pair to `mercury_precession_gate.rs`: that gate locks the *satisfied*
//! case ($\varepsilon = 0$, 1PN on) against the GR prediction. This file
//! locks the *violated* case ($\varepsilon \approx 0.02$ AU softened
//! Plummer, 1PN on) against the closed-form softened-Plummer apsidal-
//! precession prediction derived in
//! `paper/notebooks/2026-05-28-plummer-apsidal-derivation.md`.
//!
//! Both gates are release-mode integration tests; run with
//! `cargo test --release -p apsis-1pn --tests -- --ignored`. The
//! drift-measurement function samples the periapsis longitude per orbit
//! and accumulates the unwrapped step-by-step drift, avoiding the
//! end-vs-initial mod-$2\pi$ alias that hides the true ($\sim 13.7$ rad
//! over 500 orbits) Plummer-induced cumulative drift.

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

/// Mercury orbital period in Earth days (87.969 d). Used to scale per-orbit
/// drift to per-century. Hardcoded to the physical value rather than the
/// canonical-units literal so the comparison stays in observational units.
const T_MERCURY_DAYS: f64 = 87.969;
const DAYS_PER_CENTURY: f64 = 36_525.0;
const RAD_TO_ARCSEC: f64 = 180.0 * 3600.0 / PI;

fn orbits_per_century() -> f64 {
    DAYS_PER_CENTURY / T_MERCURY_DAYS
}

/// Closed-form prediction for the apsidal-precession rate induced by a
/// Plummer-softened pair potential `U(r) = -G M m / sqrt(r^2 + eps^2)`
/// at the given Sun–Mercury orbital parameters, expressed in arcseconds
/// per Earth century.
///
/// Derived in `paper/notebooks/2026-05-28-plummer-apsidal-derivation.md`
/// by two independent routes (near-circular frequency decomposition and
/// orbital-averaged disturbing-function Lagrange equation) and verified
/// against an independent scipy DOP853 integration to 3.2% relative
/// agreement at $\varepsilon^2/a^2 = 2.67\times 10^{-3}$:
///
/// $$\Delta\varpi_\text{orbit} = -\frac{3\pi\,\varepsilon^2}{a^2\,(1-e^2)^2}.$$
///
/// The sign is negative (retrograde) because Plummer softening shallows
/// the effective well at small $r$, delaying the next periapsis. The
/// leading-order $\varepsilon^2/a^2$ expansion is accurate to $\sim 3\%$
/// at Mercury parameters; higher-order corrections enter at $O(\varepsilon^4)$.
fn predicted_softened_plummer_drift_arcsec_per_century(a_au: f64, e: f64, epsilon_au: f64) -> f64 {
    let per_orbit_rad = -3.0 * PI * epsilon_au.powi(2) / (a_au.powi(2) * (1.0 - e.powi(2)).powi(2));
    per_orbit_rad * RAD_TO_ARCSEC * orbits_per_century()
}

/// Run the §3.2 violated case (Sun + Mercury, Plummer-softened kernel
/// at `epsilon` AU, 1PN registered) for 500 orbits and return the
/// cumulative apsidal drift in arcsec per Earth century.
///
/// **Per-orbit unwrap.** Samples the periapsis longitude at every orbit
/// and accumulates the unwrapped step-by-step drift. The end-vs-initial
/// `omega` difference modulo $2\pi$ aliases for large drifts: the
/// Plummer-violated cumulative is $\sim -13.7$ rad over 500 orbits,
/// well outside $(-\pi, \pi]$. Per-orbit unwrap is lossless because the
/// per-orbit step magnitude is bounded by the closed-form prediction
/// $3\pi\varepsilon^2 / [a^2(1-e^2)^2] \approx 2.7\times 10^{-2}$ rad,
/// safely inside $\pm\pi$.
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

/// Smoke gate — the violated case produces a drift large enough that
/// the upcoming theory-match comparison has a real signal to lock onto.
/// Records the measured value via `--nocapture` for cross-check against
/// the literal cited in `paper.md` §3.2.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn plummer_violated_case_drift_is_measured() {
    let measured = measure_drift_arcsec_per_century(EPSILON_VIOLATED);
    eprintln!(
        "[exactness-theory] measured drift at eps={EPSILON_VIOLATED} AU: \
         {measured:.0} arcsec/century",
    );
    assert!(
        measured.abs() > 1_000.0,
        "violated case drift too small ({measured:.3e} arcsec/century) — \
         scenario or kernel setup regressed",
    );
}

/// Theory-match gate — measured drift matches the closed-form prediction
/// from the lab notebook within the acceptance tolerance.
///
/// Tolerance is 5 %: the leading-order $\varepsilon^2/a^2$ expansion is
/// 3.2 % accurate at Mercury parameters per the scipy DOP853 cross-check
/// recorded in `paper/notebooks/2026-05-28-plummer-apsidal-derivation.md`,
/// and 5 % gives margin for the $O(\varepsilon^4)$ correction plus
/// any IAS15 substep-cadence noise.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn plummer_drift_matches_softened_theory() {
    let measured = measure_drift_arcsec_per_century(EPSILON_VIOLATED);
    let predicted =
        predicted_softened_plummer_drift_arcsec_per_century(A_MERCURY, E_MERCURY, EPSILON_VIOLATED);

    let rel_err = (measured - predicted).abs() / predicted.abs();

    const ACCEPTANCE: f64 = 5e-2;

    eprintln!(
        "[exactness-theory] eps={EPSILON_VIOLATED} AU:\n  \
         measured  = {measured:.6e} arcsec/century\n  \
         predicted = {predicted:.6e} arcsec/century\n  \
         rel_err   = {rel_err:.4e} (gate: {ACCEPTANCE})",
    );

    assert!(
        rel_err < ACCEPTANCE,
        "softened-Plummer apsidal drift mismatch: \
         measured {measured:.3e} arcsec/century vs predicted {predicted:.3e}; \
         rel_err = {rel_err:.3e} exceeds {ACCEPTANCE}",
    );
}
