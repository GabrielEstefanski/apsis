//! Deferred mutation requested by hooks, applied after hook dispatch.

/// Action a hook may request of the orchestrator. Applied in insertion
/// order once all hooks for the current phase have fired.
#[derive(Debug, Clone)]
pub enum Command {
    /// Stop the main loop after the current step completes.
    Stop,
}
