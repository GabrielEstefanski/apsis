//! Gravity simulator core — headless simulation engine.
//!
//! This crate is the physics + simulation stack with **zero UI
//! dependencies**. It is consumed by:
//!
//! * `gravity-sim-app` — interactive GUI (egui/wgpu).
//! * `benches/` — Criterion-driven IAS15 harness with versioned baselines.
//! * Out-of-tree client crates (e.g. perturbation-force plugins like
//!   `gravity-sim-1pn`) against the public API surfaced here.
//!
//! The `app/render → core` read direction is enforced by the workspace
//! split: nothing in this crate may import from the app crate.

pub mod core;
pub mod domain;
pub mod io;
pub mod physics;
pub mod templates;
