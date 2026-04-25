//! Shared infrastructure across benchmark examples.
//!
//! Any file in this directory is accessible to a sibling example via
//! `mod common;` — Cargo treats `examples/common/` as a plain source
//! directory, not as an example itself, because no `main.rs` lives at
//! the root.
//!
//! Each individual example consumes only the subset of this module it
//! needs, so Cargo's per-example compilation will flag the rest as
//! dead code when the example is built in isolation. The
//! `#[allow(dead_code)]` below suppresses those warnings without
//! weakening warnings in the rest of the crate; the shared
//! infrastructure as a whole is exercised across the set of examples.

#[allow(dead_code)]
pub mod scenarios;
