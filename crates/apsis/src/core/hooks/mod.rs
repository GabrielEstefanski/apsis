//! Hook system — observer + command pattern for simulation extensions.
//!
//! # Design
//!
//! Hooks are user-supplied trait objects that observe the simulation and
//! request mutations through a restricted command channel. The three rules:
//!
//! 1. **Observation is immutable.** Hooks receive [`HookContext`], a read-only
//!    view. They cannot touch [`System`](crate::core::system::System) directly.
//! 2. **Mutation is deferred.** Hooks return [`Command`]s which the orchestrator
//!    applies after all hooks for the current phase have fired. This keeps
//!    dispatch order independent of side-effects.
//! 3. **Events are detected inside `step`, emitted afterwards.** Collision and
//!    escape detection run on the integrated state; hooks see the events only
//!    once the physics advance is complete. This preserves symplecticity —
//!    hooks cannot interfere with the integrator.
//!
//! # Ordering
//!
//! Hooks are stored in a [`HookRegistry`] as [`HookEntry`]s sorted by
//! `priority` (ascending — lower fires first). Ordering is stable across
//! insertions at equal priority.
//!
//! # Example
//!
//! ```ignore
//! struct EnergyLogger;
//! impl SimHook for EnergyLogger {
//!     fn post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
//!         eprintln!("t={:.3} ΔE/E={:.2e}", ctx.t, ctx.rel_energy_error);
//!         Vec::new()
//!     }
//! }
//!
//! system.hooks_mut().register(0, Box::new(EnergyLogger));
//! ```

pub mod commands;
pub mod context;
pub mod events;

pub use commands::Command;
pub use context::{HookContext, HookPhase, HookPhaseKind};
pub use events::{CollisionEvent, EscapeEvent};

/// A simulation hook. All methods default to no-ops so implementors only
/// override the phases they care about.
///
/// # Contract
///
/// - Hooks **must not** panic — a panicking hook aborts the simulation.
/// - Hooks **must not** retain references to `HookContext` beyond the call.
/// - Mutations are expressed solely via returned [`Command`]s.
pub trait SimHook: Send {
    /// Human-readable identifier used in logs and diagnostics.
    fn name(&self) -> &'static str {
        "unnamed"
    }

    /// Fired once before the integrator runs this step.
    fn pre_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired once after integration and event detection complete.
    fn post_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired for each [`CollisionEvent`] detected after the step.
    fn on_collision(&mut self, _event: &CollisionEvent, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired for each [`EscapeEvent`] detected after the step.
    fn on_escape(&mut self, _event: &EscapeEvent, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired once when the simulation lifecycle ends. Hooks holding
    /// open resources flush them here. Returned commands are ignored.
    fn on_finish(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Whether this hook reads [`HookContext::resume_state`]. When any
    /// registered hook returns `true`, the orchestrator calls
    /// [`Integrator::resume_state`](crate::physics::integrator::traits::Integrator::resume_state)
    /// before dispatch so the bytes are available; otherwise the field
    /// is `None` and the serialisation work is skipped.
    fn wants_resume_state(&self) -> bool {
        false
    }
}

/// A hook with its dispatch priority. Lower priorities fire first.
pub struct HookEntry {
    pub priority: i32,
    pub hook: Box<dyn SimHook>,
}

/// Container for all registered hooks.
///
/// Dispatch order within a phase is strictly by ascending `priority`. Insertion
/// order breaks ties (stable).
#[derive(Default)]
pub struct HookRegistry {
    entries: Vec<HookEntry>,
}

impl HookRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a hook at the given priority. Lower priorities fire first.
    pub fn register(&mut self, priority: i32, hook: Box<dyn SimHook>) {
        let pos = self.entries.partition_point(|e| e.priority <= priority);
        self.entries.insert(pos, HookEntry { priority, hook });
    }

    /// Remove all registered hooks.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// `true` when at least one registered hook reports
    /// [`SimHook::wants_resume_state`]. Read by the orchestrator before
    /// dispatch to decide whether to serialise the integrator's
    /// scratch state.
    pub fn any_wants_resume_state(&self) -> bool {
        self.entries.iter().any(|e| e.hook.wants_resume_state())
    }

    /// Dispatch `pre_step` to every hook and collect commands in order.
    pub fn dispatch_pre_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.pre_step(ctx));
        }
        out
    }

    /// Dispatch `post_step` to every hook and collect commands in order.
    pub fn dispatch_post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.post_step(ctx));
        }
        out
    }

    pub fn dispatch_collision(
        &mut self,
        event: &CollisionEvent,
        ctx: &HookContext<'_>,
    ) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.on_collision(event, ctx));
        }
        out
    }

    pub fn dispatch_escape(&mut self, event: &EscapeEvent, ctx: &HookContext<'_>) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.on_escape(event, ctx));
        }
        out
    }

    pub fn dispatch_finish(&mut self, ctx: &HookContext<'_>) {
        for entry in &mut self.entries {
            let _ = entry.hook.on_finish(ctx);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Recorder {
        label: &'static str,
        log: std::sync::Arc<std::sync::Mutex<Vec<&'static str>>>,
    }

    impl SimHook for Recorder {
        fn name(&self) -> &'static str {
            self.label
        }
        fn pre_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
            self.log.lock().unwrap().push(self.label);
            Vec::new()
        }
    }

    #[test]
    fn priority_order_is_ascending_with_stable_ties() {
        let log = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut reg = HookRegistry::new();

        reg.register(10, Box::new(Recorder { label: "late", log: log.clone() }));
        reg.register(0, Box::new(Recorder { label: "first", log: log.clone() }));
        reg.register(5, Box::new(Recorder { label: "mid_a", log: log.clone() }));
        reg.register(5, Box::new(Recorder { label: "mid_b", log: log.clone() }));

        let ctx = HookContext {
            bodies: &[],
            names: &[],
            t: 0.0,
            dt: 0.0,
            steps: 0,
            rel_energy_error: None,
            rel_angular_momentum_error: None,
            phase: HookPhase(HookPhaseKind::PreStep),
            resume_state: None,
        };
        reg.dispatch_pre_step(&ctx);

        assert_eq!(*log.lock().unwrap(), vec!["first", "mid_a", "mid_b", "late"]);
    }

    #[test]
    fn default_methods_are_no_ops() {
        struct Empty;
        impl SimHook for Empty {}

        let mut reg = HookRegistry::new();
        reg.register(0, Box::new(Empty));

        let ctx = HookContext {
            bodies: &[],
            names: &[],
            t: 0.0,
            dt: 0.0,
            steps: 0,
            rel_energy_error: None,
            rel_angular_momentum_error: None,
            phase: HookPhase(HookPhaseKind::PostStep),
            resume_state: None,
        };
        assert!(reg.dispatch_post_step(&ctx).is_empty());
    }
}
