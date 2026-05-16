// ── Clippy allowances ─────────────────────────────────────────────────────────
//
// These are suppressed crate-wide rather than per-site because the pattern
// each one refers to is a deliberate style choice, not a bug to fix:
//
// * `needless_range_loop` — N-body physics code indexes `bodies` and parallel
//   accumulators by position; clippy's preferred `for b in bodies.iter()` loses
//   the index, which we need for the scratch-accumulator write. Switching to
//   `enumerate()` everywhere adds noise without changing performance.
// * `drop_non_drop` — borrow-checker dances use `drop(x)` to release a borrow
//   explicitly before the end of scope. Non-Drop types are still valid here;
//   the intent is scope-ending, not destructor invocation.
// * `too_many_arguments` — IAS15 and integrator internals take 7-8 arguments
//   deliberately: splitting into structs would either add indirection on the
//   hot loop or fragment closely-related state.
#![allow(clippy::needless_range_loop, clippy::drop_non_drop, clippy::too_many_arguments)]

//! Gravity simulator core — headless simulation engine.
//!
//! Zero UI dependencies, enforced by the `core-isolation` CI gate.
//! Federated operator crates and downstream consumers (Python
//! distribution, external visualisation shells, plugin authors)
//! depend on the public extension API surfaced here.

pub mod contract;
pub mod core;
pub mod domain;
pub mod io;
pub mod math;
pub mod physics;
pub mod templates;
pub mod units;
