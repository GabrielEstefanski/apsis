//! Operator regime-of-validity contract.
//!
//! Every operator's force law is derived under physical assumptions —
//! "test particle around a dominant body", "small eccentricity",
//! "v ≪ c", "two-body system", etc. Outside the assumed regime the
//! derivation no longer applies and the operator's output is no longer
//! physics, even though the integrator continues to integrate it as
//! if it were.
//!
//! `RegimeViolation` is the structured signal an operator returns
//! when the current body state crosses one of its assumed bounds.
//! [`crate::core::system::System`] runs
//! [`Operator::check_regime`](crate::physics::integrator::Operator::check_regime)
//! at registration (initial state) and periodically during integration
//! (at each operator's
//! [`regime_check_cadence`](crate::physics::integrator::Operator::regime_check_cadence)),
//! emitting one `warn_diag` per `(operator, bound)` per session.
//!
//! # Relationship to `KernelRequirements`
//!
//! `KernelRequirements` is the **numerical contract** — what the
//! gravitational kernel must guarantee for the operator's derivation
//! to be applicable (Exactness, Continuity). `RegimeBounds` is the
//! **physical contract** — what regime of body state the operator's
//! derivation was constructed for. Together they cover the full
//! precondition surface of a federated operator: kernel side checked
//! once at registration, regime side checked statically AND dynamically.

use crate::domain::body::Body;

/// How seriously to treat a regime violation. The integrator never
/// halts on its own — these are diagnostic levels for the structured
/// log bus.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Body state is within the operator's regime but approaching
    /// the bound. Heads-up for users running near the edge of
    /// validity (e.g. eccentricity creeping toward an averaging-
    /// theory limit).
    Approaching,
    /// Body state has crossed the bound. The operator's derivation
    /// no longer strictly applies; results are best-effort rather
    /// than the validated baseline. Integration continues.
    Exceeded,
    /// Body state is far outside the operator's envelope (e.g. a
    /// test-particle 1PN attached to an equal-mass binary). The
    /// operator's output is not physics; integration continues only
    /// because the system does not auto-halt, but the run is no
    /// longer scientifically defensible.
    Hard,
}

/// Diagnostic record produced by [`Operator::check_regime`] when the
/// current body state crosses a bound the operator's derivation
/// assumed.
///
/// The fields are structural — they identify the operator, the
/// specific bound, the value observed, the threshold the bound was
/// set at, and (optionally) which body triggered the check. The
/// structured log bus serialises them so consumers can pattern-match
/// or render for human display.
#[derive(Debug, Clone)]
pub struct RegimeViolation {
    /// Operator identifier from
    /// [`crate::physics::integrator::Operator::name`].
    pub operator: &'static str,
    /// Short tag identifying the bound that was crossed
    /// (e.g. `"max_secondary_to_primary_mass_ratio"`,
    /// `"max_v_over_c"`, `"max_eccentricity"`). Used as the
    /// dedup key alongside `operator` for warn-once.
    pub bound: &'static str,
    /// Current value of the quantity the bound constrains.
    pub value: f64,
    /// The bound the derivation was constructed against.
    pub threshold: f64,
    /// Severity of the violation — see [`Severity`].
    pub severity: Severity,
    /// Optional body index when the violation is per-body
    /// (mass-ratio against a specific body, periapse of a specific
    /// orbit). `None` when the violation is system-wide.
    pub body_index: Option<usize>,
    /// Free-text explanation of the underlying physics — what the
    /// derivation assumed and what breaks when the bound is crossed.
    /// One sentence, written for a researcher who knows the operator
    /// but may have forgotten the derivation's fine print.
    pub message: &'static str,
}

impl RegimeViolation {
    /// Stable dedup key for warn-once filtering: `(operator, bound)`.
    /// Two violations with the same key are considered the same
    /// alert; only the first one in a session is emitted.
    pub fn dedup_key(&self) -> (&'static str, &'static str) {
        (self.operator, self.bound)
    }
}

/// Helper for operators that compare a body's mass against a primary's:
/// returns the violation severity (and value/threshold pair) when the
/// ratio crosses the supplied warn/hard thresholds. Returns `None`
/// when within the regime.
///
/// `warn_threshold`: ratio at which to emit `Severity::Exceeded`.
/// `hard_threshold`: ratio at which to escalate to `Severity::Hard`.
/// Both are typically operator-specific physical bounds.
pub fn classify_mass_ratio(
    ratio: f64,
    warn_threshold: f64,
    hard_threshold: f64,
) -> Option<(Severity, f64)> {
    if ratio >= hard_threshold {
        Some((Severity::Hard, hard_threshold))
    } else if ratio >= warn_threshold {
        Some((Severity::Exceeded, warn_threshold))
    } else if ratio >= 0.9 * warn_threshold {
        Some((Severity::Approaching, warn_threshold))
    } else {
        None
    }
}

/// Pair-mass ratio between two bodies in the slice. Returns `None`
/// when either body has zero mass (no meaningful ratio).
pub fn mass_ratio(bodies: &[Body], primary: usize, secondary: usize) -> Option<f64> {
    let m_primary = bodies.get(primary)?.mass;
    let m_secondary = bodies.get(secondary)?.mass;
    if m_primary == 0.0 {
        return None;
    }
    Some(m_secondary / m_primary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_below_approaching_returns_none() {
        assert!(classify_mass_ratio(0.001, 0.01, 0.1).is_none());
    }

    #[test]
    fn classify_in_approaching_band_flags_approaching() {
        // 0.0095 is between 0.9 * 0.01 = 0.009 and 0.01 itself.
        let (sev, threshold) = classify_mass_ratio(0.0095, 0.01, 0.1).unwrap();
        assert_eq!(sev, Severity::Approaching);
        assert_eq!(threshold, 0.01);
    }

    #[test]
    fn classify_above_warn_returns_exceeded() {
        let (sev, threshold) = classify_mass_ratio(0.05, 0.01, 0.1).unwrap();
        assert_eq!(sev, Severity::Exceeded);
        assert_eq!(threshold, 0.01);
    }

    #[test]
    fn classify_above_hard_returns_hard() {
        let (sev, threshold) = classify_mass_ratio(0.5, 0.01, 0.1).unwrap();
        assert_eq!(sev, Severity::Hard);
        assert_eq!(threshold, 0.1);
    }

    #[test]
    fn dedup_key_pairs_operator_and_bound() {
        let v = RegimeViolation {
            operator: "TestOp",
            bound: "max_x",
            value: 1.0,
            threshold: 0.5,
            severity: Severity::Exceeded,
            body_index: None,
            message: "",
        };
        assert_eq!(v.dedup_key(), ("TestOp", "max_x"));
    }
}
