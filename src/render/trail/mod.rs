//! Trail visualization bounded context.
//!
//! Trails are a **rendering concern** — they record and display historical
//! body positions for the viewer, and know nothing about the physics domain
//! beyond consuming positions from published snapshots.
//!
//! # Layering
//!
//! ```text
//!   physics (domain)            render/trail (visualization)
//!   ─────────────────           ─────────────────────────────
//!   System                      TrailSampler   — when to record
//!   Body positions      ───→    TrailBuffer    — what is stored
//!                               TrailRenderer  — how it is drawn
//!                               TrailStyle     — how it looks
//! ```
//!
//! Phase 1 introduces [`TrailStyle`] (value object) and [`TrailSampler`]
//! (strategy trait) as the public abstractions. Phase 2 will migrate
//! [`crate::render::trail_buffer::TrailBuffer`] into this module and move
//! all sampling off the physics thread.

pub mod sampler;
pub mod style;

pub use sampler::{AdaptiveSampler, SampleDecision, TimeSampler, TrailSampler};
pub use style::{TrailStyle, TrailStylePreset};
