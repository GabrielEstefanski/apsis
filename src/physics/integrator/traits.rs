//! Core integrator abstraction: trait, context, result, and kind enum.
//!
//! # Architecture
//!
//! ```text
//! ┌────────────┐     ┌──────────────────┐     ┌───────────┐
//! │ Integrator │────▶│ IntegratorContext │────▶│ ForceModel│
//! │  (trait)   │     │  (force+params)  │     │  (trait)  │
//! └────────────┘     └──────────────────┘     └───────────┘
//!       │
//!       ▼
//! ┌────────────┐
//! │ StepResult │  ← returned after each integration step
//! └────────────┘
//! ```
//!
//! The [`Integrator`] trait replaces the old `Integrator` enum, enabling
//! new integration schemes to be added without touching the core.
//!
//! [`IntegratorKind`] is a plain enum retained for UI display, snapshot
//! serialisation, and `Metrics`.  It is **not** used for dispatch.

use crate::domain::body::Body;
use crate::physics::integrator::force_model::ForceModel;
use crate::physics::integrator::perturbation::PerturbationForce;

// ── IntegratorKind (serialisable enum) ────────────────────────────────────────

/// Identifies an integration algorithm without carrying behaviour.
///
/// Used for snapshot serialisation, UI combo-boxes, and `Metrics`.
/// The actual stepping logic lives in structs that implement [`Integrator`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegratorKind {
    VelocityVerlet,
    Yoshida4,
    WisdomHolman,
}

impl IntegratorKind {
    /// Short human-readable label used in the UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "Velocity Verlet (2nd)",
            Self::Yoshida4 => "Yoshida 4th-order",
            Self::WisdomHolman => "Wisdom–Holman (2nd, Keplerian)",
        }
    }

    /// Formal convergence order in the time step.
    pub fn order(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 2,
        }
    }

    /// Nominal number of force evaluations consumed per full time step.
    pub fn force_evals_per_step(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 1,
        }
    }

    /// One-line description shown in the UI tooltip.
    pub fn description(self) -> &'static str {
        match self {
            Self::VelocityVerlet => {
                "2nd-order symplectic leapfrog. Fast; energy oscillates around \
                 the initial value. Phase error ∝ dt². Good for real-time \
                 visualisation and short integrations."
            },
            Self::Yoshida4 => {
                "4th-order symplectic composition (Forest–Ruth). 4 force evals \
                 per step but phase error ∝ dt⁴ — allows 5–10× larger dt for \
                 the same energy conservation. Required for publication-quality \
                 long-term runs."
            },
            Self::WisdomHolman => {
                "Mixed-variable symplectic map. Keplerian two-body motion is \
                 solved analytically; perturbations are stepped numerically. \
                 Designed for hierarchical planetary systems."
            },
        }
    }

    /// All known variants, in the order shown in the UI combo-box.
    pub const ALL: [IntegratorKind; 3] = [
        IntegratorKind::VelocityVerlet,
        IntegratorKind::Yoshida4,
        IntegratorKind::WisdomHolman,
    ];
}

// ── IntegratorContext ─────────────────────────────────────────────────────────

/// Everything an integrator needs from the simulation besides bodies and dt.
///
/// Passed as `&mut` so the integrator can call `force.compute()` (which
/// requires `&mut self` for tree rebuilds, etc.).
///
/// # Design
///
/// This struct exists to **avoid coupling integrators to `System`**.
/// An integrator never sees the full `System`; it only sees this narrow
/// interface of force model + physical parameters.
pub struct IntegratorContext<'a> {
    /// The force model (e.g. Barnes-Hut gravity).
    pub force: &'a mut dyn ForceModel,

    /// Gravitational scaling factor: `G_eff = G₀ · g_factor`.
    pub g_factor: f64,

    /// Non-gravitational perturbation forces (radiation, drag, …).
    pub perturbations: &'a [Box<dyn PerturbationForce>],
}

// ── StepResult ────────────────────────────────────────────────────────────────

/// Output produced by a single integration step.
///
/// Centralises the physical diagnostics that `System` needs after each step,
/// so no integrator-specific logic leaks into the orchestrator.
pub struct StepResult {
    /// Gravitational potential energy at the end-of-step positions,
    /// **already scaled** by `g_factor`.
    pub potential_energy: f64,

    /// `true` if the integrator fell back to a different algorithm this step
    /// (e.g. Wisdom–Holman → Yoshida4 when the dominance criterion fails).
    pub used_fallback: bool,
}

// ── Integrator trait ──────────────────────────────────────────────────────────

/// A symplectic (or general) N-body integrator.
///
/// # Contract
///
/// - `step()` advances `bodies` by one time step `dt`.
/// - `step()` may call `ctx.force.compute()` one or more times.
/// - `step()` must leave `acc` consistent with the final body positions
///   (so that diagnostics can read it).
/// - `step()` must apply `ctx.g_factor` scaling and `ctx.perturbations`
///   at the appropriate points in the scheme.
///
/// # Mutability
///
/// `&mut self` is required because some integrators carry internal state
/// across steps (e.g. Wisdom–Holman's fallback integrator, IAS15's
/// predictor–corrector history).
pub trait Integrator: Send {
    /// Advance the system by one time step.
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<(f64, f64)>,
    ) -> StepResult;

    /// Which algorithm this integrator represents.
    fn kind(&self) -> IntegratorKind;
}
