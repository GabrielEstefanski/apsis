//! Immutable simulation view handed to hooks.
//!
//! Hooks receive a [`HookContext`] so they can observe state without mutating
//! it. All mutations must flow through [`Command`](super::commands::Command)s
//! returned by the hook, which are applied deterministically after dispatch.

use crate::domain::body::Body;

/// Read-only snapshot of simulation state at the moment a hook fires.
#[derive(Debug, Clone, Copy)]
pub struct HookPhase(pub HookPhaseKind);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HookPhaseKind {
    /// Before the integrator runs. Bodies still reflect the previous step.
    PreStep,
    /// After integration completes.
    PostStep,
    /// Lifecycle ending — hooks flush close-on-drop resources here.
    Finish,
}

#[derive(Clone)]
pub struct HookContext<'a> {
    pub bodies: &'a [Body],
    pub names: &'a [String],
    pub t: f64,
    pub dt: f64,
    pub steps: u64,
    pub rel_energy_error: Option<f64>,
    pub rel_angular_momentum_error: Option<f64>,
    pub phase: HookPhase,
    /// Serialised integrator scratch, populated by the orchestrator
    /// when at least one registered hook returns
    /// [`SimHook::wants_resume_state`](super::SimHook::wants_resume_state)`
    /// = true`. `None` otherwise — the field is per-step state, not a
    /// permanent capability.
    pub resume_state: Option<Vec<u8>>,
}
