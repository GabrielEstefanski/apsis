//! Composable operators in the integration step.
//!
//! Three operator traits, distinguished by what they contribute:
//!
//! | Trait | Force | Hamiltonian potential |
//! |---|---|---|
//! | [`HamiltonianOperator`] | yes (= −∇V) | optional (V via [`potential`](HamiltonianOperator::potential)) |
//! | [`NonConservativeOperator`] | yes | no |
//! | [`Operator`] base | no | no |
//!
//! # Dispatch contract
//!
//! Symplectic-class integrators (`WisdomHolman`, `Mercurius`,
//! `Yoshida4`, `WHFast`) call operators at:
//!
//! 1. Pre-half-kick: sum of `accumulate_force` over Hamiltonian +
//!    non-conservative operators, applied as `v += dt/2 · acc`.
//! 2. Drift / Kepler: operators do not participate.
//! 3. Post-half-kick: same as 1.
//! 4. Step boundary, synchronized state: `observe` on every Operator.
//!
//! IAS15 calls `accumulate_force` once per Picard iteration inside
//! each adaptive sub-step; `observe` still fires once per outer
//! sub-step at synchronized state.
//!
//! Integrators route through [`crate::physics::integrator::operator_dispatch`].
//!
//! # Total energy
//!
//! ```text
//! E_total = T + V_kepler + Σᵢ HamiltonianOperator[i].potential
//! ```
//!
//! The sum is over Hamiltonian operators whose [`potential`] returns
//! [`Potential::Value`]. Operators that return [`Potential::NotAvailable`]
//! are silently excluded; [`crate::core::system::System::conservation_report`]
//! surfaces the exclusion so the omission is auditable.
//!
//! Non-conservative operators and observers contribute nothing to
//! `E_total`.
//!
//! [`potential`]: HamiltonianOperator::potential

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::gravity::kernel::KernelRequirements;
use crate::physics::integrator::citation::Citation;
use crate::physics::integrator::regime::RegimeViolation;
use crate::units::UnitSystem;

// ── Potential return type ────────────────────────────────────────────────────

/// Closed-form Hamiltonian potential V(bodies) — or an explicit signal
/// that the implementation does not expose one.
///
/// `#[non_exhaustive]` reserves space for future variants
/// (`Indeterminate(f64)` for stochastic terms with known bound,
/// `Pending` for explicit TODO markers) without breaking the API.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Potential {
    /// V evaluated in closed form for the given body configuration.
    Value(f64),

    /// V exists in principle for this operator's derivation but is not
    /// exposed by this implementation. Test-particle pairwise 1PN is
    /// the canonical case: the rigorous form is the full
    /// Einstein–Infeld–Hoffmann N-body Hamiltonian, which is out of
    /// scope for a test-particle approximation crate.
    /// [`crate::core::system::System::total_energy`] excludes these
    /// terms and emits one diagnostic when first observed.
    NotAvailable,
}

// ── Base trait ───────────────────────────────────────────────────────────────

/// Base trait for anything registrable in the integration loop.
///
/// Pure `Operator` implementations contribute no force and no energy —
/// they read body state at step boundaries via
/// [`observe`](Self::observe). Physics that contributes a force
/// implements [`HamiltonianOperator`] or [`NonConservativeOperator`].
pub trait Operator: Send + Sync {
    /// Human-readable identifier used by
    /// [`crate::core::system::System::conservation_report`] and structured
    /// diagnostics. Default: fully-qualified Rust type path (e.g.
    /// `apsis_1pn::PostNewtonian1PN`); override for a shorter or
    /// configuration-bearing string (e.g. `"PostNewtonian1PN(solar)"`).
    fn name(&self) -> &'static str {
        std::any::type_name::<Self>()
    }

    /// [`UnitSystem`] this operator was constructed for. Default:
    /// `None` — unit-system-agnostic (constant push, dimensionless
    /// coupling, …). Operators with dimensional parameters **must**
    /// override this: `Some(u)` is compared against the `System`'s
    /// units at registration, and a mismatch returns
    /// [`UnitSystemMismatch`] instead of silently integrating wrong
    /// physics.
    fn declared_units(&self) -> Option<UnitSystem> {
        None
    }

    /// Step-boundary observation. Called once per outer integration
    /// step, after body state is synchronised in the inertial frame.
    /// Default: no-op.
    ///
    /// `&mut self` lets the operator carry per-step state (rolling
    /// buffers, derivative estimates). Dispatch is serial within a
    /// step; cross-thread safety is the `Send + Sync` bound, not
    /// interior mutability.
    fn observe(&mut self, bodies: &[Body], t: f64, dt: f64) {
        let _ = (bodies, t, dt);
    }

    /// Kernel invariants the operator's derivation relies on. Default
    /// none. Violations against the active kernel emit one structured
    /// `warn_diag` per invariant at registration time.
    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::none()
    }

    /// Check the current body state against the operator's
    /// regime-of-validity bounds (mass-ratio, eccentricity, v/c,
    /// periapse, …). Returns one entry per crossed bound; an empty
    /// vector when the system is within the operator's envelope.
    /// Default: no check (operator is regime-agnostic, or has not yet
    /// declared its bounds).
    ///
    /// `System::add_*_perturbation` invokes this once at registration
    /// against the initial body state; `System::step` invokes it again
    /// at the cadence returned by [`regime_check_cadence`](Self::regime_check_cadence).
    /// One `warn_diag` is emitted per `(operator, bound)` pair per
    /// session — the same violation persisting across multiple checks
    /// does not respam the bus.
    ///
    /// See [`RegimeViolation`] for the structured fields and
    /// [`crate::physics::integrator::regime`] for the contract.
    fn check_regime(&self, bodies: &[Body], t: f64) -> Vec<RegimeViolation> {
        let _ = (bodies, t);
        Vec::new()
    }

    /// Outer steps between dynamic regime checks. Default 100; the
    /// dispatcher uses the minimum cadence across registered operators.
    /// Inert when [`check_regime`](Self::check_regime) is the no-op
    /// default.
    fn regime_check_cadence(&self) -> usize {
        100
    }

    /// Reference card: BibTeX of the source paper, DOI, and the
    /// implementing crate's name/version/commit. Default `None` (test
    /// fakes, internal tooling). Published perturbation crates should
    /// override it, using `env!("CARGO_PKG_NAME")` etc. so the captured
    /// values are the operator crate's, not apsis core's. See
    /// [`Citation`] for the field contract.
    fn citation(&self) -> Option<Citation> {
        None
    }
}

// ── Hamiltonian operator ─────────────────────────────────────────────────────

/// Operator derivable from a Hamiltonian: force = −∇V.
///
/// [`accumulate_force`](Self::accumulate_force) and
/// [`potential`](Self::potential) describe the same operator from two
/// angles. When both are provided, they must stay consistent — the
/// gradient of `potential` is `accumulate_force`, component-wise per
/// body. When `potential` returns [`Potential::NotAvailable`], only
/// the force half of the contract is honoured.
///
/// # Closed-form V — implementation notes
///
/// Operators with non-trivial closed-form `V` should handle singular
/// cases inside the method (e.g. logarithmic fallback at `γ = −1` for a
/// power-law central force `F ∝ r^γ`). Returning `Value(f64::NAN)` is a
/// contract violation; return [`Potential::NotAvailable`] instead when
/// the closed form is genuinely undefined or numerically unsafe in the
/// current regime.
pub trait HamiltonianOperator: Operator {
    /// Add this operator's force contribution to `acc[i]` for each
    /// body `i`. The integrator initialises `acc` before the dispatch
    /// loop; implementations must add, not overwrite.
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]);

    /// Closed-form Hamiltonian potential V(bodies).
    ///
    /// **Default: [`Potential::NotAvailable`].** Operators whose force
    /// derives from a Hamiltonian but whose closed-form V is not
    /// implemented in this crate (test-particle pairwise 1PN; custom
    /// researcher derivations with V deferred) inherit the default.
    /// Symplectic integrators do not depend on `potential` — they
    /// depend on the force derivation being conservation-friendly,
    /// which is the trait's promise. `System::total_energy` excludes
    /// `NotAvailable` contributions and surfaces the exclusion through
    /// `System::conservation_report`.
    fn potential(&self, bodies: &[Body]) -> Potential {
        let _ = bodies;
        Potential::NotAvailable
    }
}

// ── Non-conservative operator ────────────────────────────────────────────────

/// Operator with a force but no Hamiltonian — drag, radiation
/// reaction, dissipative coupling.
///
/// Symplectic integrators degrade silently when registered. The
/// integrator does not auto-degrade its splitting; the registration
/// site emits one `warn_diag` so the broken invariant is auditable.
pub trait NonConservativeOperator: Operator {
    /// Add this operator's force contribution to `acc`. Same dispatch
    /// position as [`HamiltonianOperator::accumulate_force`].
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]);
}

// ── Registration error ───────────────────────────────────────────────────────

/// Operator registration refused — the operator's
/// [`declared_units`](Operator::declared_units) disagree with the
/// `System`'s [`UnitSystem`].
///
/// Returned by [`crate::core::system::System::add_hamiltonian_perturbation`]
/// (and the non-conservative / observer counterparts) instead of panicking.
/// Caller decides the policy: propagate with `?`, log and skip, swap
/// the operator, fall back to a different unit system, etc. Two-tier
/// semantics on registration:
///
/// - **Hard** (`Result::Err(UnitSystemMismatch)`) — the operator's `c`
///   (or other dimensional parameter) would be interpreted in the
///   wrong frame; integration would silently produce wrong physics.
/// - **Soft** (`warn_diag`) — kernel-precondition violations,
///   regime-of-validity bounds, NC-on-symplectic. Integration can
///   continue with the user's choice; the bus carries the audit trail.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct UnitSystemMismatch {
    /// Operator identifier from
    /// [`Operator::name`].
    pub operator: &'static str,
    /// `UnitSystem` the operator was constructed for (its
    /// `declared_units`).
    pub operator_units: UnitSystem,
    /// `UnitSystem` the `System` was constructed for.
    pub system_units: UnitSystem,
}

impl std::fmt::Display for UnitSystemMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Unit-system mismatch on operator registration:\n  \
             operator '{}' was constructed for {}\n  \
             System was constructed for {}\n\
             Construct the operator with the same UnitSystem passed to \
             System::new(), or omit `declared_units` if the operator is \
             unit-system-agnostic.",
            self.operator, self.operator_units, self.system_units
        )
    }
}

impl std::error::Error for UnitSystemMismatch {}

// ── Descriptor for federation seam ───────────────────────────────────────────

/// Plugin metadata for downstream perturbation crates (`apsis-1pn`
/// and friends). UIs and headless runners collect descriptors into a
/// registry without learning concrete operator types.
///
/// `build` takes the target `UnitSystem` so the produced operator
/// carries the same units as the `System` that will register it —
/// the registration-time `declared_units` check then succeeds rather
/// than returning [`UnitSystemMismatch`] on the consumer. UI catalogs
/// query the active `System`'s units and pass them through.
pub trait HamiltonianOperatorDescriptor: Send + Sync {
    fn name(&self) -> &str;
    fn description(&self) -> &str;
    fn kernel_requirements(&self) -> KernelRequirements;
    fn build(&self, units: UnitSystem) -> Box<dyn HamiltonianOperator>;
}
