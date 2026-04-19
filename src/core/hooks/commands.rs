//! Command pattern for deferred state mutation requested by hooks.
//!
//! Hooks never mutate [`System`](crate::core::system::System) directly. They
//! return a list of [`Command`]s that the orchestrator applies in deterministic
//! order once all hooks for the current phase have been dispatched. This keeps
//! symplectic invariants intact and makes hook-driven behaviour reproducible.

use crate::domain::body::{Body, NamedBody};

/// Mutation requested by a hook, applied after the current dispatch phase.
#[derive(Debug, Clone)]
pub enum Command {
    /// Remove the body at `index`.
    RemoveBody { index: usize },

    /// Append a new body (with optional explicit name).
    AddBody(NamedBody),

    /// Replace two bodies with a single merged body.
    ///
    /// `remove_a` and `remove_b` indices are dropped; the new body is added.
    /// Hooks should pre-compute the merge (mass-weighted position/velocity,
    /// momentum conservation) — the orchestrator trusts the payload.
    Merge { remove_a: usize, remove_b: usize, merged: Body, merged_name: Option<String> },

    /// Stop the main loop after the current step completes.
    Stop,
}

