//! Conditioning-aware helpers for relative conservation metrics.

/// Minimum baseline magnitude for a well-defined relative drift
/// `(X − X₀) / |X₀|`. Below this the f64 round-off floor dominates
/// the signal: noise of order `ε` against denominator of order
/// `√ε` yields SNR ≈ 1.
pub const MIN_RELATIVE_DENOMINATOR: f64 = 1.4901161193847656e-8;
// f64::EPSILON.sqrt() — written as a literal because const fn sqrt
// is unstable.

/// Relative drift `delta / |baseline|`, or `None` when `|baseline|`
/// is below [`MIN_RELATIVE_DENOMINATOR`].
///
/// `delta` is the already-computed `current − baseline`. Taking the
/// pre-computed delta (rather than `current` and `baseline`
/// separately) prevents the helper from being miscalled as
/// `(value − baseline − baseline) / |baseline|`.
#[inline]
pub(crate) fn regime_aware_rel(delta: f64, baseline: f64) -> Option<f64> {
    if baseline.abs() < MIN_RELATIVE_DENOMINATOR { None } else { Some(delta / baseline.abs()) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn well_conditioned_baseline_returns_some() {
        let result = regime_aware_rel(1.0e-10, 0.5);
        assert_eq!(result, Some(2.0e-10));
    }

    #[test]
    fn baseline_at_threshold_returns_some() {
        let result = regime_aware_rel(1.0, MIN_RELATIVE_DENOMINATOR);
        assert!(result.is_some());
    }

    #[test]
    fn baseline_below_threshold_returns_none() {
        let result = regime_aware_rel(1.0, MIN_RELATIVE_DENOMINATOR / 2.0);
        assert_eq!(result, None);
    }

    #[test]
    fn dust_regime_baseline_returns_none() {
        // |E_initial| ≈ 5e-16 for radiation_dust.py
        let result = regime_aware_rel(1.0e-17, 5.0e-16);
        assert_eq!(result, None);
    }

    #[test]
    fn kepler_regime_baseline_returns_some() {
        // |E_initial| ≈ 0.5 for kepler_2body
        let result = regime_aware_rel(3.775e-15, -0.5);
        assert_eq!(result, Some(3.775e-15 / 0.5));
    }

    #[test]
    fn zero_baseline_returns_none() {
        let result = regime_aware_rel(1.0, 0.0);
        assert_eq!(result, None);
    }

    #[test]
    fn negative_baseline_uses_abs() {
        let result = regime_aware_rel(1.0, -2.0);
        assert_eq!(result, Some(0.5));
    }

    #[test]
    fn min_relative_denominator_matches_sqrt_epsilon() {
        let expected = f64::EPSILON.sqrt();
        assert_eq!(MIN_RELATIVE_DENOMINATOR, expected);
    }

    // ── Integration tests: end-to-end regime detection through System ──

    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::physics::integrator::IntegratorKind;
    use crate::units::UnitSystem;

    #[test]
    fn well_conditioned_kepler_returns_some_rel_drift() {
        let bodies = vec![
            Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1e-3).at(1.0, 0.0).with_velocity(0.0, 1.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(1e-3);
        for _ in 0..100 {
            sys.step();
        }
        assert!(sys.energy_delta().is_some(), "Kepler regime must report Some(rel)");
        assert!(sys.lz_delta().is_some(), "Kepler Lz regime must report Some(rel)");
    }

    #[test]
    fn precision_limited_dust_returns_none_rel_drift() {
        let bodies = vec![
            Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1e-15).at(1.0, 0.0).with_velocity(0.0, 1.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(1e-3);
        for _ in 0..10 {
            sys.step();
        }
        assert_eq!(sys.energy_delta(), None, "dust scenario |E_initial| ~ 1e-15 must report None");
        assert!(sys.abs_energy_drift().is_finite(), "abs drift remains finite");
    }
}
