//! `KernelProperties` and `KernelRequirements` — the type-level extension
//! contract for gravitational kernels.
//!
//! A kernel implementation reports the physical invariants it provides via
//! [`KernelProperties`]; a perturbation declares the invariants it requires
//! via [`KernelRequirements`]; the match between the two is computed by
//! [`KernelRequirements::check_against`] and surfaces as a vector of
//! [`RequirementViolation`] records.
//!
//! These invariants are not ad-hoc labels. Each field encodes a specific
//! statement about K(r) that appears in the derivation or analysis of one
//! or more perturbations:
//!
//! - [`Exactness`] is a statement about the Newtonian base potential:
//!   whether the library delivers K(r) = 1/r exactly, or a modified form.
//!   Perturbations derived from the Newtonian Hamiltonian (1PN, J2, tides)
//!   require the exact form; a softened or truncated base invalidates the
//!   derivation itself.
//! - [`Continuity`] is a statement about phase-space geometry: symplectic
//!   integration relies on the Hamiltonian flow preserving phase-space
//!   volume, which requires a smooth H. Force discontinuities produce
//!   impulsive accelerations that cannot be represented within any
//!   symplectic splitting scheme.

// ── KernelProperties ──────────────────────────────────────────────────────── //

/// Physical invariants a kernel satisfies.
///
/// Returned from [`Kernel::properties`](super::Kernel::properties). May
/// depend on runtime state: for example, a Plummer kernel dynamically
/// reports [`Exactness::Exact`] when every body in the system has softening
/// length zero, and [`Exactness::Softened`] otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KernelProperties {
    /// Structural class of K(r) — whether it is 1/r exact, softened, or
    /// some other deviation.
    pub exactness: Exactness,
    /// Smoothness class of K(r) — equivalently, the highest derivative
    /// order of the force −K'(r) that remains continuous on (0, ∞).
    pub continuity: Continuity,
}

/// Structural class of the pair potential K(r).
///
/// Ordered so that [`Exact`] is strictest and [`Modified`] weakest;
/// [`Exactness::satisfies`] uses the rank to test whether a provided
/// exactness satisfies a required one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Exactness {
    /// K(r) = 1/r exactly. Preserves Newtonian point-mass structure as
    /// assumed by any derivation from the Newtonian Hamiltonian.
    Exact,
    /// K(r) uses some form of softening. Specific softening parameters
    /// are not encoded at the type level; the variant communicates only
    /// that the kernel has a non-trivial ε-dependence.
    Softened,
    /// K(r) deviates from 1/r in a way that is neither exact nor simply
    /// Plummer-softened — e.g., truncated, tabulated, or user-defined.
    Modified,
}

impl Exactness {
    const fn rank(self) -> u8 {
        match self {
            Self::Modified => 0,
            Self::Softened => 1,
            Self::Exact => 2,
        }
    }

    /// Whether a kernel providing `self` satisfies an extension requiring
    /// `required`.
    ///
    /// Ordering: `Exact > Softened > Modified`. A kernel satisfies a
    /// requirement iff its rank is at least the required rank.
    #[inline]
    pub const fn satisfies(self, required: Self) -> bool {
        self.rank() >= required.rank()
    }
}

/// Smoothness class of the pair potential K(r).
///
/// Ordered so that [`Smooth`] is the strongest guarantee and [`C0`] the
/// weakest; [`Continuity::satisfies`] uses the rank to test whether a
/// provided continuity class satisfies a required minimum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Continuity {
    /// K ∈ C⁰ but not C¹. Force `F = −K'` may be discontinuous.
    C0,
    /// K ∈ C¹ but not C². Force continuous; force derivative may jump.
    C1,
    /// K ∈ C² but not C^∞.
    C2,
    /// K ∈ C^∞. Force and all derivatives continuous.
    Smooth,
}

impl Continuity {
    const fn rank(self) -> u8 {
        match self {
            Self::C0 => 0,
            Self::C1 => 1,
            Self::C2 => 2,
            Self::Smooth => 3,
        }
    }

    /// Whether a kernel with continuity `self` satisfies an extension
    /// requiring at least `required`.
    ///
    /// Ordering: `Smooth > C2 > C1 > C0`.
    #[inline]
    pub const fn satisfies(self, required: Self) -> bool {
        self.rank() >= required.rank()
    }
}

// ── KernelRequirements ────────────────────────────────────────────────────── //

/// Invariants an extension requires the active kernel to provide.
///
/// Each field is optional: `None` means the extension does not constrain
/// that invariant. The natural constructors are [`none`](Self::none) for
/// an extension with no kernel preconditions and
/// [`exact_and_smooth`](Self::exact_and_smooth) for the canonical
/// "Newtonian 1PN-compatible" set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct KernelRequirements {
    /// If `Some(x)`, the kernel's exactness must satisfy `x` (in the
    /// sense of [`Exactness::satisfies`]).
    pub required_exactness: Option<Exactness>,
    /// If `Some(c)`, the kernel's continuity must satisfy at least `c`
    /// (in the sense of [`Continuity::satisfies`]).
    pub min_continuity: Option<Continuity>,
}

impl KernelRequirements {
    /// No kernel requirements — the extension is compatible with any
    /// kernel the system might be configured with.
    #[inline]
    pub const fn none() -> Self {
        Self { required_exactness: None, min_continuity: None }
    }

    /// The canonical 1PN-compatible requirement set: exact 1/r base
    /// plus C^∞ smoothness.
    #[inline]
    pub const fn exact_and_smooth() -> Self {
        Self {
            required_exactness: Some(Exactness::Exact),
            min_continuity: Some(Continuity::Smooth),
        }
    }

    /// Match these requirements against a kernel's reported
    /// [`KernelProperties`]. Returns every violation identified — a
    /// single kernel change can violate multiple invariants at once.
    pub fn check_against(&self, props: &KernelProperties) -> Vec<RequirementViolation> {
        let mut out = Vec::new();

        if let Some(required) = self.required_exactness
            && !props.exactness.satisfies(required)
        {
            out.push(RequirementViolation::Exactness { required, provided: props.exactness });
        }

        if let Some(required) = self.min_continuity
            && !props.continuity.satisfies(required)
        {
            out.push(RequirementViolation::Continuity { required, provided: props.continuity });
        }

        out
    }
}

/// A single-invariant violation identified by
/// [`KernelRequirements::check_against`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RequirementViolation {
    Exactness { required: Exactness, provided: Exactness },
    Continuity { required: Continuity, provided: Continuity },
}

impl std::fmt::Display for RequirementViolation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Exactness { required, provided } => {
                write!(f, "kernel exactness {provided:?} does not satisfy required {required:?}",)
            },
            Self::Continuity { required, provided } => {
                write!(f, "kernel continuity {provided:?} does not satisfy minimum {required:?}",)
            },
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;

    // ── Exactness ordering ────────────────────────────────────────────── //

    #[test]
    fn exactness_exact_satisfies_exact() {
        assert!(Exactness::Exact.satisfies(Exactness::Exact));
    }

    #[test]
    fn exactness_exact_satisfies_softened_and_modified() {
        assert!(Exactness::Exact.satisfies(Exactness::Softened));
        assert!(Exactness::Exact.satisfies(Exactness::Modified));
    }

    #[test]
    fn exactness_softened_does_not_satisfy_exact() {
        assert!(!Exactness::Softened.satisfies(Exactness::Exact));
    }

    #[test]
    fn exactness_softened_satisfies_softened() {
        assert!(Exactness::Softened.satisfies(Exactness::Softened));
    }

    #[test]
    fn exactness_modified_only_satisfies_modified() {
        assert!(Exactness::Modified.satisfies(Exactness::Modified));
        assert!(!Exactness::Modified.satisfies(Exactness::Softened));
        assert!(!Exactness::Modified.satisfies(Exactness::Exact));
    }

    // ── Continuity ordering ───────────────────────────────────────────── //

    #[test]
    fn continuity_smooth_satisfies_all_weaker() {
        assert!(Continuity::Smooth.satisfies(Continuity::Smooth));
        assert!(Continuity::Smooth.satisfies(Continuity::C2));
        assert!(Continuity::Smooth.satisfies(Continuity::C1));
        assert!(Continuity::Smooth.satisfies(Continuity::C0));
    }

    #[test]
    fn continuity_c0_does_not_satisfy_higher_classes() {
        assert!(!Continuity::C0.satisfies(Continuity::C1));
        assert!(!Continuity::C0.satisfies(Continuity::C2));
        assert!(!Continuity::C0.satisfies(Continuity::Smooth));
    }

    #[test]
    fn continuity_c2_satisfies_c2_and_below() {
        assert!(Continuity::C2.satisfies(Continuity::C2));
        assert!(Continuity::C2.satisfies(Continuity::C1));
        assert!(Continuity::C2.satisfies(Continuity::C0));
        assert!(!Continuity::C2.satisfies(Continuity::Smooth));
    }

    // ── KernelRequirements matching ──────────────────────────────────── //

    #[test]
    fn empty_requirements_never_violate() {
        let req = KernelRequirements::none();
        let props = KernelProperties { exactness: Exactness::Modified, continuity: Continuity::C0 };
        assert!(req.check_against(&props).is_empty());
    }

    #[test]
    fn exact_smooth_requirements_satisfied_by_newton_smooth() {
        let req = KernelRequirements::exact_and_smooth();
        let props =
            KernelProperties { exactness: Exactness::Exact, continuity: Continuity::Smooth };
        assert!(req.check_against(&props).is_empty());
    }

    #[test]
    fn exact_smooth_requirements_report_exactness_violation_only() {
        let req = KernelRequirements::exact_and_smooth();
        let props =
            KernelProperties { exactness: Exactness::Softened, continuity: Continuity::Smooth };
        let violations = req.check_against(&props);
        assert_eq!(violations.len(), 1);
        assert!(matches!(
            violations[0],
            RequirementViolation::Exactness {
                required: Exactness::Exact,
                provided: Exactness::Softened
            }
        ));
    }

    #[test]
    fn truncated_kernel_triggers_both_exactness_and_continuity_violations() {
        // Representative of a `TruncatedPlummerKernel`: modified
        // exactness, C0 continuity.
        let req = KernelRequirements::exact_and_smooth();
        let props = KernelProperties { exactness: Exactness::Modified, continuity: Continuity::C0 };
        let violations = req.check_against(&props);
        assert_eq!(violations.len(), 2);
        assert!(violations.iter().any(|v| matches!(v, RequirementViolation::Exactness { .. })));
        assert!(violations.iter().any(|v| matches!(v, RequirementViolation::Continuity { .. })));
    }

    #[test]
    fn requirement_only_one_field_checked_independently() {
        let req =
            KernelRequirements { required_exactness: Some(Exactness::Exact), min_continuity: None };
        let props = KernelProperties { exactness: Exactness::Exact, continuity: Continuity::C0 };
        assert!(req.check_against(&props).is_empty());
    }

    // ── Display formatting ────────────────────────────────────────────── //

    #[test]
    fn violation_display_mentions_required_and_provided() {
        let v = RequirementViolation::Exactness {
            required: Exactness::Exact,
            provided: Exactness::Softened,
        };
        let msg = format!("{v}");
        assert!(msg.contains("Exact"));
        assert!(msg.contains("Softened"));
    }

    #[test]
    fn continuity_violation_display_mentions_both_classes() {
        let v = RequirementViolation::Continuity {
            required: Continuity::Smooth,
            provided: Continuity::C0,
        };
        let msg = format!("{v}");
        assert!(msg.contains("Smooth"));
        assert!(msg.contains("C0"));
    }
}
