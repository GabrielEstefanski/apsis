//! Synthetic round-trip gate for `CentralForce::from_apsidal_rate`.
//!
//! The observable-inversion constructor inverts a target apsidal precession rate
//! ω̇ into the coupling `A` that produces it on a near-circular orbit.
//! This gate closes the loop end-to-end: register at ω̇ = X, integrate
//! for N orbits, fit ω̇ from the trajectory, assert |ω̇_measured −
//! ω̇_input| / |ω̇_input| < tolerance.
//!
//! No reliance on any specific physics scenario (Mercury, Schwarzschild,
//! …) — the gate verifies that the *inversion arithmetic plus the
//! integrator* together reproduce whatever ω̇ the caller asks for.
//! That is what the observable-inversion constructor promises as a contract, independent of the
//! force's physical interpretation.

use std::f64::consts::TAU;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_central::CentralForce;

/// Tolerance for the synthetic round-trip. Sources of bias:
///
/// - **Tamayo formula approximation.** The inversion uses the
///   instantaneous separation `d` (perihelion = `a(1−e)` here) rather
///   than the secular semi-major axis `a`. For `e = 0.1` the
///   `d^(γ+2)` factor differs from `a^(γ+2)` by ≈ 11%, biasing the
///   recovered `A` by the same factor. Documented behaviour of the
///   formula; deviating would break the cross-implementation
///   round-trip.
/// - **IAS15 truncation** over 50 orbits: well below 10⁻⁵ on ω, way
///   under the gate.
/// - **`(1 − e²)` correction** in Murray-Dermott: ~1% at `e = 0.1`,
///   partially cancels the previous bias.
///
/// Empirical agreement (printed by the gate as `[central-roundtrip] rel`):
/// 2.38%. 5% is the regression bound, matching the tolerance pattern set by
/// `apsis-radiation`'s Burns 1979 gate.
const ROUND_TRIP_TOLERANCE: f64 = 0.05;

/// Observable-inversion contract: registering with `omega_dot = X` reproduces
/// `omega_dot = X` within 5 % over 50 orbits at γ = -3, e = 0.1.
/// Empirical agreement: 2.38 %.
#[test]
fn from_apsidal_rate_synthetic_round_trip() {
    let units = UnitSystem::solar_canonical();
    // Source + receiver at 1 AU on an e = 0.1 orbit at perihelion.
    // Tamayo's near-circular envelope holds well below e ~ 0.2; e = 0.1
    // gives an eccentricity vector with magnitude large enough that
    // its direction (= ω) is robust against per-step f64 noise. At
    // e = 1e-3 the angular noise on e_vec aliased catastrophically
    // with the per-orbit sampling, drowning the secular signal.
    let sun = Body::star(1.0);
    let bound_e = 0.1_f64;
    let a = 1.0_f64;
    let r0 = a * (1.0 - bound_e);
    let v0 = ((1.0 + bound_e) / (1.0 - bound_e)).sqrt(); // canonical: v_circ = 1 at a = 1
    let receiver = Body::rocky(1e-10).at(r0, 0.0).with_velocity(0.0, v0);
    let bodies = vec![sun, receiver];

    // Observable inversion: pick a target ω̇ and let the operator invert.
    let omega_dot_in = 1.5e-3_f64; // rad / Gaussian time unit — large enough to
    // accumulate ~ 0.4 rad over 50 orbits, well above the IAS15 noise floor.
    let gamma = -3.0_f64;
    let force = CentralForce::from_apsidal_rate(0, 1, omega_dot_in, gamma, &bodies, units)
        .expect("circular pair must invert");

    let mut sys = System::new(bodies, units).with_integrator(IntegratorKind::Ias15).with_dt(1e-3);
    sys.add_hamiltonian_perturbation(Box::new(force))
        .expect("CentralForce registration must succeed");

    // Sample ω at integer multiples of the unperturbed orbital
    // period T = 2π. Phase-locked sampling sees the same osculating
    // ω value plus the secular drift each time, so the per-orbit
    // oscillation cancels exactly out of the linear fit. Without
    // this, the osculating ω swings ~e radians per orbit and
    // aliases catastrophically against any non-period-locked
    // sampling cadence.
    let g_code = units.g();
    let n_orbits = 50_usize;
    let n_samples = n_orbits + 1;
    let t_orbit = TAU; // GM = a = 1 in canonical
    let mut ts = Vec::with_capacity(n_samples);
    let mut omegas = Vec::with_capacity(n_samples);

    sys.step(); // populate caches before measuring
    for k in 0..n_samples {
        let t_target = (k as f64) * t_orbit;
        sys.integrate_until(t_target);
        let elems = compute_elements(sys.bodies(), 1, 0, g_code).expect("orbit must stay bound");
        ts.push(sys.t());
        omegas.push(elems.omega);
    }

    // Unwrap ω from [-π, π] to a continuous function so the linear
    // fit sees the secular drift, not the principal-value jumps.
    let omegas_unwrapped = unwrap_angles(&omegas);

    // Linear regression: slope = Σ(t·ω) / Σ(t²) after centring on mean.
    let omega_dot_measured = linear_slope(&ts, &omegas_unwrapped);

    let rel = ((omega_dot_measured - omega_dot_in) / omega_dot_in).abs();
    eprintln!("[central-roundtrip] measured = {omega_dot_measured:.6e}, in = {omega_dot_in:.6e}, rel = {rel:.4}");
    assert!(
        rel < ROUND_TRIP_TOLERANCE,
        "Observable-inversion round-trip: input ω̇ = {omega_dot_in:.4e}, measured ω̇ = {omega_dot_measured:.4e}, \
         relative error = {rel:.4} (gate: {ROUND_TRIP_TOLERANCE})",
    );
}

/// Counter-test: with no operator registered, the same orbit should
/// *not* precess (Newtonian Keplerian closure). If this fails, the
/// drift detected in the headline test could be coming from IAS15 or
/// from the slight non-circularity, not from the operator.
#[test]
fn keplerian_baseline_does_not_precess() {
    let units = UnitSystem::solar_canonical();
    let sun = Body::star(1.0);
    let bound_e = 0.1_f64;
    let a = 1.0_f64;
    let r0 = a * (1.0 - bound_e);
    let v0 = ((1.0 + bound_e) / (1.0 - bound_e)).sqrt();
    let receiver = Body::rocky(1e-10).at(r0, 0.0).with_velocity(0.0, v0);
    let bodies = vec![sun, receiver];

    let mut sys = System::new(bodies, units).with_integrator(IntegratorKind::Ias15).with_dt(1e-3);

    let g_code = units.g();
    sys.step();
    let n_orbits = 50_usize;
    let n_samples = n_orbits + 1;
    let t_orbit = TAU;
    let mut ts = Vec::with_capacity(n_samples);
    let mut omegas = Vec::with_capacity(n_samples);
    for k in 0..n_samples {
        let t_target = (k as f64) * t_orbit;
        sys.integrate_until(t_target);
        let elems = compute_elements(sys.bodies(), 1, 0, g_code).expect("orbit bound");
        ts.push(sys.t());
        omegas.push(elems.omega);
    }
    let omegas_unwrapped = unwrap_angles(&omegas);
    let omega_dot_baseline = linear_slope(&ts, &omegas_unwrapped);

    eprintln!("[kepler-baseline] ω̇ = {omega_dot_baseline:.4e}");
    // Gate sits 6 orders below the operator-driven drift it must catch
    // (ω̇_in = 1.5e-3) and far above the single-platform slope here (~2e-17);
    // kept conservative for cross-platform robustness, not pinned to one
    // platform's f64 floor.
    assert!(
        omega_dot_baseline.abs() < 1e-9,
        "Keplerian baseline drift exceeds floor: ω̇ = {omega_dot_baseline:.4e} (bound 1e-9)",
    );
}

/// Unwrap a sequence of angles in [-π, π] to a continuous function by
/// folding 2π jumps. Standard signal-processing primitive; no
/// dependency on `numpy.unwrap`-equivalent in the apsis stack.
fn unwrap_angles(xs: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(xs.len());
    if xs.is_empty() {
        return out;
    }
    out.push(xs[0]);
    let mut accum = 0.0_f64;
    for i in 1..xs.len() {
        let delta = xs[i] - xs[i - 1];
        if delta > std::f64::consts::PI {
            accum -= TAU;
        } else if delta < -std::f64::consts::PI {
            accum += TAU;
        }
        out.push(xs[i] + accum);
    }
    out
}

/// Least-squares slope of `y` vs `x`. Centring on means avoids the
/// numerical conditioning issues of Σx²Σy − (Σx)(Σxy) at large `t`.
fn linear_slope(xs: &[f64], ys: &[f64]) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    let n = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (x, y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        num += dx * (y - mean_y);
        den += dx * dx;
    }
    num / den
}
