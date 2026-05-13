//! Parameter validation for operator constructors.
//!
//! Operators expose three kinds of constructors per the convention
//! documented in [`apsis::contract`](crate::contract):
//!
//! - **Named regime** (`OperatorName::for_xxx`, `OperatorName::with_xxx`) —
//!   encodes a physical regime in the name itself. No raw parameter.
//! - **Observable inversion** (`OperatorName::from_<observable>`) —
//!   inverts a desired observable to compute the operator's parameter.
//! - **Raw escape** (`OperatorName::from_raw_xxx`,
//!   `OperatorName::from_raw_xxx_validated`) — accepts the parameter
//!   directly; the `_validated` form cross-checks against a reference
//!   (typically a [`crate::units::UnitSystem`] or a representative
//!   [`crate::domain::body::Body`]) and returns
//!   [`ParameterValidationError`] when the cross-check fails.
//!
//! This module provides the error type used by every `_validated`
//! constructor across the perturbation surface.

/// Diagnostic produced when a raw constructor's value disagrees with
/// the reference it was validated against.
///
/// The error is structural: it names the operator, the parameter, the
/// value the caller passed, the value the validator derived from the
/// reference, and the relative tolerance the validator used. Consumers
/// can render it for human display via [`std::fmt::Display`] or pattern-
/// match on the fields for programmatic recovery.
#[derive(Debug, Clone)]
pub struct ParameterValidationError {
    /// Operator that produced the error (e.g. `"PostNewtonian1PN"`).
    pub operator: &'static str,
    /// Parameter name that failed validation (e.g. `"c"`).
    pub parameter: &'static str,
    /// Value the caller passed.
    pub got: f64,
    /// Value the validator derived from the reference.
    pub expected: f64,
    /// Relative-error tolerance band (e.g. `1e-3` = 0.1 %).
    pub tolerance: f64,
    /// Free-text context: which reference was used, what physics
    /// motivated the bound, suggested remediation. Used by `Display`.
    pub message: String,
}

impl std::fmt::Display for ParameterValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let rel_err = (self.got - self.expected).abs() / self.expected.abs().max(f64::EPSILON);
        write!(
            f,
            "{}::{}: got {:.6e}, expected {:.6e} (relative error {:.3e}, tolerance {:.3e}). {}",
            self.operator,
            self.parameter,
            self.got,
            self.expected,
            rel_err,
            self.tolerance,
            self.message,
        )
    }
}

impl std::error::Error for ParameterValidationError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_renders_all_fields() {
        let err = ParameterValidationError {
            operator: "TestOp",
            parameter: "x",
            got: 1.0,
            expected: 2.0,
            tolerance: 1e-3,
            message: "value is half of expected".to_string(),
        };
        let s = format!("{err}");
        assert!(s.contains("TestOp::x"));
        assert!(s.contains("1.0"));
        assert!(s.contains("2.0"));
        assert!(s.contains("half of expected"));
    }

    #[test]
    fn relative_error_includes_zero_expected_safely() {
        // The fmt path divides by expected.abs().max(EPSILON), so
        // expected=0 must not produce NaN/inf — it must still render.
        let err = ParameterValidationError {
            operator: "T",
            parameter: "p",
            got: 1.0,
            expected: 0.0,
            tolerance: 1e-3,
            message: String::new(),
        };
        let s = format!("{err}");
        assert!(!s.contains("NaN"));
        assert!(!s.contains("inf"));
    }
}
