//! System-wide conservation report derived from the registered
//! operator stack.
//!
//! `System::conservation_report()` returns one of these so the user
//! can answer "is my registered stack conservative, and if not, which
//! operator broke the invariant?" The report classifies the system as
//! a whole and lists each contributing operator with its role and
//! energy impact.
//!
//! # Attribution scope
//!
//! The report makes **structural** attribution claims only:
//!
//! - "operator X is non-conservative → it is a candidate source of
//!   energy drift"
//! - "operator Y is Hamiltonian but its V is `NotAvailable` → its
//!   contribution to `total_energy` is silently excluded"
//!
//! It does **not** make quantitative claims like "operator X
//! contributed ΔE = −1e-7 over this run." That kind of attribution
//! requires per-step integration of force · velocity per operator and
//! becomes meaningful only when a shadow-Hamiltonian tracker is
//! registered. Until then, the report is the first half of the audit
//! trail: which operators *can* break which invariant.

use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Potential,
};

// ── System-wide classification ───────────────────────────────────────────────

/// System-wide conservation property derived from the registered
/// operator stack. Composes one value per system, not per operator.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConservationClass {
    /// Every registered Hamiltonian operator exposes a closed-form V
    /// via [`Potential::Value`]. `total_energy` is exact at the
    /// integrator's truncation error floor — drift comes from the
    /// integrator, not from the registered stack.
    Hamiltonian,

    /// All registered operators are Hamiltonian-derived, but at least
    /// one returns [`Potential::NotAvailable`]. The integrator's
    /// conservation invariants still hold in derivation (no
    /// non-conservative term), but `total_energy` excludes the
    /// `NotAvailable` operators' contributions. Energy drift cannot
    /// be cleanly attributed to integrator vs. excluded V.
    HamiltonianForceOnly,

    /// At least one non-conservative operator is registered. Energy
    /// will drift at the dissipation rate of those operators by
    /// design. Symplectic integrators no longer guarantee long-term
    /// conservation invariants.
    Mixed,
}

// ── Per-operator role / status ───────────────────────────────────────────────

/// Which extension trait the operator implements. Implicit from the
/// `System` storage it sits in; surfaced in the report so a reader can
/// match each contributor to its conservation behaviour without
/// cross-referencing trait impls.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperatorRole {
    /// Force derives from a Hamiltonian.
    Hamiltonian,
    /// Force has no Hamiltonian; dissipative coupling.
    NonConservative,
}

/// Whether a Hamiltonian operator exposes a closed-form V.
/// `NotApplicable` for non-Hamiltonian operators, which do not have a
/// potential by construction.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PotentialStatus {
    /// Operator returns [`Potential::Value`]; contributes to
    /// `total_energy`.
    Available,
    /// Operator returns [`Potential::NotAvailable`]; excluded from
    /// `total_energy`.
    NotAvailable,
    /// Operator is non-conservative; no V exists by construction.
    NotApplicable,
}

/// Direction of the operator's effect on `total_energy`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyImpact {
    /// V is summed into `total_energy`; the operator's contribution
    /// is captured.
    Tracked,
    /// Operator is Hamiltonian-derived but V is not exposed; its
    /// contribution is silently excluded from `total_energy`.
    Excluded,
    /// Operator dissipates energy by design; `total_energy` drifts
    /// monotonically (or near-monotonically) in its presence.
    Dissipates,
}

/// One row of [`ConservationReport::contributors`]: one registered
/// operator and its conservation behaviour.
#[derive(Debug, Clone)]
pub struct ConservationContribution {
    /// Operator identifier from [`crate::physics::integrator::Operator::name`].
    pub operator: &'static str,
    /// Which extension trait the operator implements.
    pub role: OperatorRole,
    /// Whether closed-form V is exposed.
    pub potential_status: PotentialStatus,
    /// Direction of the operator's effect on `total_energy`.
    pub energy_impact: EnergyImpact,
}

/// Conservation property of the system as a whole, plus the per-operator
/// breakdown that produced it. Returned by
/// [`crate::core::system::System::conservation_report`].
#[derive(Debug, Clone)]
pub struct ConservationReport {
    /// System-wide classification.
    pub global: ConservationClass,
    /// Per-operator breakdown, in registration order. Hamiltonian
    /// operators precede non-conservative ones, matching dispatch
    /// order.
    pub contributors: Vec<ConservationContribution>,
}

impl ConservationReport {
    /// Build a report by classifying each registered operator and
    /// composing the system-wide property.
    ///
    /// The `bodies` slice is passed to [`HamiltonianOperator::potential`]
    /// at probe time to determine whether each operator returns a
    /// closed-form value at the current state. Operators whose
    /// `potential` is deterministically `NotAvailable` will return so
    /// regardless of body state; operators with state-dependent
    /// availability (rare; mostly hypothetical) are evaluated against
    /// the current configuration.
    pub fn build(
        bodies: &[crate::domain::body::Body],
        hamiltonian: &[Box<dyn HamiltonianOperator>],
        non_conservative: &[Box<dyn NonConservativeOperator>],
    ) -> Self {
        let mut contributors = Vec::with_capacity(hamiltonian.len() + non_conservative.len());
        let mut has_unavailable_potential = false;

        for op in hamiltonian {
            let (potential_status, energy_impact) = match op.potential(bodies) {
                Potential::Value(_) => (PotentialStatus::Available, EnergyImpact::Tracked),
                Potential::NotAvailable => {
                    has_unavailable_potential = true;
                    (PotentialStatus::NotAvailable, EnergyImpact::Excluded)
                },
            };
            contributors.push(ConservationContribution {
                operator: op.name(),
                role: OperatorRole::Hamiltonian,
                potential_status,
                energy_impact,
            });
        }

        let has_nc = !non_conservative.is_empty();
        for op in non_conservative {
            contributors.push(ConservationContribution {
                operator: op.name(),
                role: OperatorRole::NonConservative,
                potential_status: PotentialStatus::NotApplicable,
                energy_impact: EnergyImpact::Dissipates,
            });
        }

        let global = if has_nc {
            ConservationClass::Mixed
        } else if has_unavailable_potential {
            ConservationClass::HamiltonianForceOnly
        } else {
            ConservationClass::Hamiltonian
        };

        Self { global, contributors }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;
    use crate::math::Vec3;
    use crate::physics::integrator::operator::Operator;

    // ── Test fakes ────────────────────────────────────────────────────────────

    struct HamWithPotential;
    impl Operator for HamWithPotential {
        fn name(&self) -> &'static str {
            "HamWithPotential"
        }
    }
    impl HamiltonianOperator for HamWithPotential {
        fn accumulate_force(&self, _: &[Body], _: &mut [Vec3]) {}
        fn potential(&self, _: &[Body]) -> Potential {
            Potential::Value(-1.0)
        }
    }

    struct HamForceOnly;
    impl Operator for HamForceOnly {
        fn name(&self) -> &'static str {
            "HamForceOnly"
        }
    }
    impl HamiltonianOperator for HamForceOnly {
        fn accumulate_force(&self, _: &[Body], _: &mut [Vec3]) {}
        // potential() inherits default → NotAvailable
    }

    struct NcOperator;
    impl Operator for NcOperator {
        fn name(&self) -> &'static str {
            "NcOperator"
        }
    }
    impl NonConservativeOperator for NcOperator {
        fn accumulate_force(&self, _: &[Body], _: &mut [Vec3]) {}
    }

    // ── Classification ────────────────────────────────────────────────────────

    fn no_bodies() -> Vec<Body> {
        Vec::new()
    }

    #[test]
    fn empty_stack_is_hamiltonian() {
        let ham: Vec<Box<dyn HamiltonianOperator>> = Vec::new();
        let nc: Vec<Box<dyn NonConservativeOperator>> = Vec::new();
        let report = ConservationReport::build(&no_bodies(), &ham, &nc);
        assert_eq!(report.global, ConservationClass::Hamiltonian);
        assert!(report.contributors.is_empty());
    }

    #[test]
    fn all_ham_with_potential_is_hamiltonian() {
        let ham: Vec<Box<dyn HamiltonianOperator>> = vec![Box::new(HamWithPotential)];
        let nc: Vec<Box<dyn NonConservativeOperator>> = Vec::new();
        let report = ConservationReport::build(&no_bodies(), &ham, &nc);
        assert_eq!(report.global, ConservationClass::Hamiltonian);
        assert_eq!(report.contributors.len(), 1);
        assert_eq!(report.contributors[0].potential_status, PotentialStatus::Available);
        assert_eq!(report.contributors[0].energy_impact, EnergyImpact::Tracked);
    }

    #[test]
    fn unavailable_potential_degrades_to_hamiltonian_force_only() {
        let ham: Vec<Box<dyn HamiltonianOperator>> =
            vec![Box::new(HamWithPotential), Box::new(HamForceOnly)];
        let nc: Vec<Box<dyn NonConservativeOperator>> = Vec::new();
        let report = ConservationReport::build(&no_bodies(), &ham, &nc);
        assert_eq!(report.global, ConservationClass::HamiltonianForceOnly);

        let with_pot = &report.contributors[0];
        assert_eq!(with_pot.potential_status, PotentialStatus::Available);
        assert_eq!(with_pot.energy_impact, EnergyImpact::Tracked);

        let force_only = &report.contributors[1];
        assert_eq!(force_only.potential_status, PotentialStatus::NotAvailable);
        assert_eq!(force_only.energy_impact, EnergyImpact::Excluded);
    }

    #[test]
    fn any_nc_operator_makes_system_mixed() {
        let ham: Vec<Box<dyn HamiltonianOperator>> = vec![Box::new(HamWithPotential)];
        let nc: Vec<Box<dyn NonConservativeOperator>> = vec![Box::new(NcOperator)];
        let report = ConservationReport::build(&no_bodies(), &ham, &nc);
        assert_eq!(report.global, ConservationClass::Mixed);

        let nc_row = report
            .contributors
            .iter()
            .find(|c| c.role == OperatorRole::NonConservative)
            .expect("NC contributor must appear in report");
        assert_eq!(nc_row.potential_status, PotentialStatus::NotApplicable);
        assert_eq!(nc_row.energy_impact, EnergyImpact::Dissipates);
    }

    #[test]
    fn ham_operators_precede_nc_in_contributors() {
        let ham: Vec<Box<dyn HamiltonianOperator>> = vec![Box::new(HamWithPotential)];
        let nc: Vec<Box<dyn NonConservativeOperator>> = vec![Box::new(NcOperator)];
        let report = ConservationReport::build(&no_bodies(), &ham, &nc);
        assert_eq!(report.contributors[0].role, OperatorRole::Hamiltonian);
        assert_eq!(report.contributors[1].role, OperatorRole::NonConservative);
    }
}
