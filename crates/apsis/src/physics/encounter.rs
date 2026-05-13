//! Close-encounter classification.
//!
//! The close-encounter regime is a first-class observable in apsis: it
//! gates the [`Mercurius`](crate::physics::integrator) hybrid integrator's
//! sub-integration, and a generic diagnostic surface lets *any* integrator
//! emit warnings when the system enters a regime where its truncation
//! analysis no longer applies.
//!
//! [`EncounterFlag`] grades the system-wide minimum pairwise separation
//! against a user-supplied threshold:
//!
//! | Variant | Condition | Meaning |
//! | --- | --- | --- |
//! | [`Far`](EncounterFlag::Far)         | `r_min ≥ threshold` | All pairs comfortably separated |
//! | [`Approaching`](EncounterFlag::Approaching) | `0.5 · threshold ≤ r_min < threshold` | A pair has entered the warning band |
//! | [`Close`](EncounterFlag::Close)     | `r_min < 0.5 · threshold` | A pair is in active close encounter |
//!
//! The half-threshold split for [`Approaching`](EncounterFlag::Approaching)
//! is a hysteresis hint, not load-bearing physics — it lets diagnostic
//! consumers (UI advisories, logging) distinguish "starting to get close"
//! from "actively close" without tracking transitions themselves.
//!
//! # When the threshold is `None`
//!
//! [`classify`](EncounterFlag::classify) returns [`Far`](EncounterFlag::Far)
//! for any input. This makes the unconfigured default a safe no-op: the
//! flag is observable but never escalates.

/// Close-encounter regime classification of the system-wide minimum
/// pairwise separation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EncounterFlag {
    /// `r_min ≥ threshold` (or threshold unset). No pair within the
    /// configured warning distance.
    Far,
    /// `0.5 · threshold ≤ r_min < threshold`. A pair has entered the
    /// warning band but is not yet in the close regime; integrators
    /// that switch behaviour on close encounters should not switch yet.
    Approaching,
    /// `r_min < 0.5 · threshold`. A pair is in active close encounter
    /// — the regime where Mercurius hands off to its IAS15
    /// sub-integration and where IAS15-only configurations expect
    /// elevated step-rejection rates.
    Close,
}

impl EncounterFlag {
    /// Classify a minimum pairwise separation against an optional
    /// threshold.
    ///
    /// Returns [`Far`](Self::Far) when `threshold` is `None` — the
    /// unconfigured default. This makes the diagnostic safe to read
    /// unconditionally on any system.
    pub fn classify(r_min: f64, threshold: Option<f64>) -> Self {
        match threshold {
            None => Self::Far,
            Some(t) if r_min >= t => Self::Far,
            Some(t) if r_min >= 0.5 * t => Self::Approaching,
            Some(_) => Self::Close,
        }
    }

    /// Short human-readable label, suitable for diagnostic output and UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::Far => "far",
            Self::Approaching => "approaching",
            Self::Close => "close",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::EncounterFlag;

    #[test]
    fn no_threshold_always_far() {
        assert_eq!(EncounterFlag::classify(0.0, None), EncounterFlag::Far);
        assert_eq!(EncounterFlag::classify(1e-30, None), EncounterFlag::Far);
        assert_eq!(EncounterFlag::classify(f64::MAX, None), EncounterFlag::Far);
    }

    #[test]
    fn at_threshold_is_far() {
        assert_eq!(EncounterFlag::classify(1.0, Some(1.0)), EncounterFlag::Far);
    }

    #[test]
    fn between_half_and_full_threshold_is_approaching() {
        assert_eq!(EncounterFlag::classify(0.75, Some(1.0)), EncounterFlag::Approaching);
        assert_eq!(EncounterFlag::classify(0.5, Some(1.0)), EncounterFlag::Approaching);
    }

    #[test]
    fn below_half_threshold_is_close() {
        assert_eq!(EncounterFlag::classify(0.4999, Some(1.0)), EncounterFlag::Close);
        assert_eq!(EncounterFlag::classify(0.0, Some(1.0)), EncounterFlag::Close);
    }

    #[test]
    fn label_round_trip() {
        for f in [EncounterFlag::Far, EncounterFlag::Approaching, EncounterFlag::Close] {
            assert!(!f.label().is_empty());
        }
    }
}
