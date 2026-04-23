//! Trail recording — physics-side sampling, ring buffer, snapshot handoff.
//!
//! Trails are part of the **temporal evolution of the system**, not of the
//! viewer. The physics thread decides *when* to record a sample
//! ([`TrailSampler`]), the recorder captures the column ([`TrailRecorder`]),
//! and the ring buffer ([`TrailBuffer`]) holds the history that any
//! downstream consumer (renderer, CSV export, analysis script) can read.
//!
//! # Layering
//!
//! ```text
//!   core/trail (sampling + storage)        render (presentation)
//!   ───────────────────────────────        ──────────────────────
//!   TrailSampler  — when to record         TrailStyle    — how it looks
//!   TrailRecorder — side-effect of step    TrailRenderer — draw calls
//!   TrailBuffer   — what is stored
//! ```
//!
//! Read direction is `app/render → core`, never the reverse.

pub mod buffer;
pub mod recorder;
pub mod sampler;

pub use buffer::{TrailBuffer, adaptive_capacity};
pub use recorder::TrailRecorder;
pub use sampler::{ArcLengthSampler, StepSampler, TrailSampler, TrailSamplerKind};
