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
/// power-law central force `F ∝ r^γ`, following the REBOUNDx
/// `central_force` precedent). Returning `Value(f64::NAN)` is a
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

// ── Descriptor for federation seam ───────────────────────────────────────────

/// Plugin metadata for downstream perturbation crates (`apsis-1pn`
/// and friends). UIs and headless runners collect descriptors into a
/// registry without learning concrete operator types.
pub trait HamiltonianOperatorDescriptor: Send + Sync {
    fn name(&self) -> &str;
    fn description(&self) -> &str;
    fn kernel_requirements(&self) -> KernelRequirements;
    fn build(&self) -> Box<dyn HamiltonianOperator>;
}
