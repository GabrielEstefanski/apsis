//! Trail presentation — visual parameters consumed by the renderer.
//!
//! Sampling, recording, and storage live in
//! [`crate::core::trail`]; this module is strictly the "how it looks" half
//! of the trail pipeline.

pub mod style;

pub use style::{TrailStyle, TrailStylePreset};
