//! Hook extension point — observers register against the `System` and
//! receive read-only [`HookContext`] views per phase. Mutations are
//! deferred through [`Command`] and applied after dispatch.
//!
//! ```ignore
//! struct EnergyLogger;
//! impl SimHook for EnergyLogger {
//!     fn post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
//!         eprintln!("t={:.3} ΔE/E={:?}", ctx.t, ctx.rel_energy_error);
//!         Vec::new()
//!     }
//! }
//!
//! system.hooks_mut().register(0, Box::new(EnergyLogger));
//! ```

use crate::domain::body::Body;

// ── Command ──────────────────────────────────────────────────────────────────

/// Action a hook may request of the orchestrator. Applied in insertion
/// order once all hooks for the current phase have fired.
#[derive(Debug, Clone)]
pub enum Command {
    /// Stop the main loop after the current step completes.
    Stop,
}

// ── Context ──────────────────────────────────────────────────────────────────

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

/// Read-only simulation view passed to a hook for the current phase.
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
    /// [`SimHook::wants_resume_state`]`= true`. `None` otherwise.
    pub resume_state: Option<Vec<u8>>,
}

// ── Trait ────────────────────────────────────────────────────────────────────

/// A simulation hook. Methods default to no-ops; implementors override
/// only the phases they care about.
pub trait SimHook: Send {
    /// Human-readable identifier used in logs and diagnostics.
    fn name(&self) -> &'static str {
        "unnamed"
    }

    /// Fired once before the integrator runs this step.
    fn pre_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired once after integration completes.
    fn post_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Fired once when the simulation lifecycle ends. Hooks holding open
    /// resources flush them here. Returned commands are ignored.
    fn on_finish(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
        Vec::new()
    }

    /// Whether this hook reads [`HookContext::resume_state`]. When any
    /// registered hook returns `true`, the orchestrator calls
    /// [`Integrator::resume_state`](crate::physics::integrator::traits::Integrator::resume_state)
    /// before dispatch so the bytes are available.
    fn wants_resume_state(&self) -> bool {
        false
    }
}

// ── Registry ─────────────────────────────────────────────────────────────────

/// A hook with its dispatch priority.
pub struct HookEntry {
    pub priority: i32,
    pub hook: Box<dyn SimHook>,
}

/// Container for registered hooks. Dispatch order is ascending priority
/// with stable ties (insertion order).
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
    /// [`SimHook::wants_resume_state`].
    pub fn any_wants_resume_state(&self) -> bool {
        self.entries.iter().any(|e| e.hook.wants_resume_state())
    }

    pub fn dispatch_pre_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.pre_step(ctx));
        }
        out
    }

    pub fn dispatch_post_step(&mut self, ctx: &HookContext<'_>) -> Vec<Command> {
        let mut out = Vec::new();
        for entry in &mut self.entries {
            out.extend(entry.hook.post_step(ctx));
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

    fn empty_ctx<'a>(phase: HookPhaseKind) -> HookContext<'a> {
        HookContext {
            bodies: &[],
            names: &[],
            t: 0.0,
            dt: 0.0,
            steps: 0,
            rel_energy_error: None,
            rel_angular_momentum_error: None,
            phase: HookPhase(phase),
            resume_state: None,
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

        reg.dispatch_pre_step(&empty_ctx(HookPhaseKind::PreStep));

        assert_eq!(*log.lock().unwrap(), vec!["first", "mid_a", "mid_b", "late"]);
    }

    #[test]
    fn default_methods_are_no_ops() {
        struct Empty;
        impl SimHook for Empty {}

        let mut reg = HookRegistry::new();
        reg.register(0, Box::new(Empty));

        assert!(reg.dispatch_post_step(&empty_ctx(HookPhaseKind::PostStep)).is_empty());
    }
}
