//! Gravitational force evaluation for 3D N-body simulations.
//!
//! ## Module layout
//!
//! | Sub-module | Responsibility |
//! |---|---|
//! | [`kernel`] | [`Kernel`] trait and concrete implementations (Plummer default) |
//! | [`tree`] | Barnes-Hut octree data structure (algorithm) |
//! | [`engine`] | [`BarnesHutEngine`]: orchestrates tree + kernel (integration) |

mod engine;
pub mod kernel;
// `pub(crate)` so the in-flight MAC harness (`physics::perf_mac`) can name
// `tree::MacKind` directly. Reverts to private in the experiment's final
// commit when the toggle is removed.
pub(crate) mod tree;

pub use engine::BarnesHutEngine;
pub use kernel::{G, Kernel, PlummerKernel, pair_eps2};
