//! Structured diagnostic events — the payload type flowing through
//! the [event bus](super::bus).
//!
//! The design keeps events cheap enough to construct on any code
//! path that needs them while carrying enough structure that UI
//! consumers can filter, group, and format without string-parsing.
//!
//! # Field ownership
//!
//! [`Event::message`] is `&'static str` because every call site
//! today is a literal. Field names are `&'static str` for the same
//! reason. Field values are owned `String` because values are
//! typically runtime data (numbers, formatted paths, debug dumps) —
//! borrowing would force every publisher to juggle lifetimes for
//! little gain at the emission rates we expect.
//!
//! # When to introduce a new [`Source`]
//!
//! Add a variant when a new *subsystem* starts emitting events that
//! UI consumers would want to filter on independently. Do not add a
//! variant per module; grouping by physical role of the emitter
//! keeps the filter list short and the mental model clean.

use std::time::SystemTime;

/// Severity of a diagnostic event. Maps directly to the three
/// canonical log levels the UI renders with distinct visual weight.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Level {
    /// Informational: a configuration change, a milestone in a long
    /// run, a state transition the user asked for.
    Info,
    /// Warning: something the user should notice (pairing auto-
    /// correction, scenario stiffness, scale advisory) but which
    /// does not stop execution.
    Warn,
    /// Error: execution cannot continue in its current form (failed
    /// checkpoint write, unreachable invariant, unrecoverable force
    /// evaluation). Emitting this does not by itself terminate —
    /// the publisher is responsible for follow-up action.
    Error,
}

impl Level {
    /// Short tag used in line-oriented stderr output. Kept stable
    /// for log-scraping continuity.
    pub fn tag(self) -> &'static str {
        match self {
            Level::Info => "INFO",
            Level::Warn => "WARN",
            Level::Error => "ERROR",
        }
    }
}

/// Subsystem that emitted the event. UI consumers filter on this.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Source {
    /// Integrator algorithm (IAS15 floor hits, Picard
    /// non-convergence, adaptive controller events).
    Integrator,
    /// Force-model configuration or pairing events.
    ForceModel,
    /// Physics thread orchestration (mode transitions, run
    /// start/complete, pause/abort lifecycle).
    PhysicsThread,
    /// System-level configuration and enforcement (integrator/force
    /// pairing auto-correction, parameter clamps).
    System,
    /// User-initiated action reaching the event bus for audit
    /// purposes (save, load, abort, commit).
    User,
}

impl Source {
    /// Short human-readable label. Used in UI filter chips and the
    /// line-oriented stderr bridge.
    pub fn label(self) -> &'static str {
        match self {
            Source::Integrator => "integrator",
            Source::ForceModel => "force-model",
            Source::PhysicsThread => "physics-thread",
            Source::System => "system",
            Source::User => "user",
        }
    }
}

/// One diagnostic event. Produced by the [`warn_diag!`](crate::warn_diag)
/// (and sibling) macros, delivered synchronously to every subscriber
/// registered on the [bus](super::bus).
#[derive(Debug, Clone)]
pub struct Event {
    /// Wall-clock time the event was created. `SystemTime` rather
    /// than `Instant` so UI consumers can render it against the
    /// user's locale clock.
    pub timestamp: SystemTime,
    /// Severity.
    pub level: Level,
    /// Emitting subsystem.
    pub source: Source,
    /// Primary message. Always a string literal so the event type
    /// stays cheap to clone and display.
    pub message: &'static str,
    /// Structured fields (`key`, `value`). Value is a pre-formatted
    /// string so consumers do not need to handle arbitrary types.
    pub fields: Vec<(&'static str, String)>,
    /// Optional coalescing key. Consumers that de-duplicate events
    /// (notification centers, toast queues) use this to merge
    /// runs of semantically-identical events into a single entry
    /// with a count. When `None`, each event stands alone.
    pub coalesce_key: Option<&'static str>,
}

impl Event {
    /// Construct an event with no fields and no coalesce key.
    /// Builder methods add the rest.
    pub fn new(level: Level, source: Source, message: &'static str) -> Self {
        Self {
            timestamp: SystemTime::now(),
            level,
            source,
            message,
            fields: Vec::new(),
            coalesce_key: None,
        }
    }

    /// Attach the coalescing key. See [`coalesce_key`](Self::coalesce_key).
    pub fn with_coalesce_key(mut self, key: &'static str) -> Self {
        self.coalesce_key = Some(key);
        self
    }

    /// Attach a single field. Prefer the [`warn_diag!`](crate::warn_diag)
    /// macro over calling this directly — the macro captures
    /// `stringify!(ident)` as the key automatically.
    pub fn with_field(mut self, key: &'static str, value: String) -> Self {
        self.fields.push((key, value));
        self
    }

    /// One-line string representation for the stderr bridge. Stable
    /// format: `[gravity-sim <TAG>] <message> { key1=value1 key2=value2 }`.
    /// Pre-existing log scrapers depend on the prefix token.
    pub fn format_single_line(&self) -> String {
        if self.fields.is_empty() {
            format!("[gravity-sim {}] {}", self.level.tag(), self.message)
        } else {
            use std::fmt::Write as _;
            let mut s = String::with_capacity(128);
            let _ = write!(&mut s, "[gravity-sim {}] {} {{", self.level.tag(), self.message);
            for (k, v) in &self.fields {
                // Always a leading space: the opening `{` has no
                // trailing space, and the `}` is added with a leading
                // space below. This keeps the separator symmetric
                // without per-iteration first/rest bookkeeping.
                let _ = write!(&mut s, " {}={}", k, v);
            }
            s.push_str(" }");
            s
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_single_line_no_fields() {
        let e = Event::new(Level::Warn, Source::System, "hello");
        let s = e.format_single_line();
        assert_eq!(s, "[gravity-sim WARN] hello");
    }

    #[test]
    fn format_single_line_with_fields() {
        let e = Event::new(Level::Warn, Source::Integrator, "floor reached")
            .with_field("dt", "1e-12".into())
            .with_field("count", "42".into());
        let s = e.format_single_line();
        assert_eq!(s, "[gravity-sim WARN] floor reached { dt=1e-12 count=42 }");
    }

    #[test]
    fn coalesce_key_is_preserved() {
        let e =
            Event::new(Level::Warn, Source::Integrator, "x").with_coalesce_key("ias15.floor_hit");
        assert_eq!(e.coalesce_key, Some("ias15.floor_hit"));
    }
}
