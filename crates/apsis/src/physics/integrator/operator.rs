//! Composable operators in the integration step.
//!
//! Three operator traits, distinguished by what they contribute:
//!
//! | Trait | Force | Hamiltonian |
//! |---|---|---|
//! | [`HamiltonianOperator`] | yes (= −∇V) | yes (V) |
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
//! E_total = T + V_kepler + Σᵢ HamiltonianOperator[i].energy_contribution
//! ```
//!
//! Non-conservative operators and observers contribute nothing to
//! `E_total`.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::gravity::kernel::KernelRequirements;

// ── Base trait ───────────────────────────────────────────────────────────────

/// Base trait for anything registrable in the integration loop.
///
/// Pure `Operator` implementations contribute no force and no energy —
/// they read body state at step boundaries via
/// [`observe`](Self::observe). Physics that contributes a force
/// implements [`HamiltonianOperator`] or [`NonConservativeOperator`].
pub trait Operator: Send + Sync {
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
/// `accumulate_force` and `energy_contribution` describe the same
/// operator from two angles and must stay consistent — the latter's
/// gradient is the former, component-wise per body.
pub trait HamiltonianOperator: Operator {
    /// Add this operator's force contribution to `acc[i]` for each
    /// body `i`. The integrator initialises `acc` before the dispatch
    /// loop; implementations must add, not overwrite.
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]);

    /// Hamiltonian term `V`. Summed by [`crate::core::system::System::total_energy`].
    fn energy_contribution(&self, bodies: &[Body]) -> f64;
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
