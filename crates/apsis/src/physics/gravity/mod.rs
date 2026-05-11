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
mod tree;

pub use engine::BarnesHutEngine;
pub use kernel::{G, Kernel, PlummerKernel, pair_eps2};

// Re-exported for the perf_soa harness's AoS shadow path; reverted with
// the harness when the SoA experiment closes.
#[cfg(test)]
pub(crate) use engine::WalkCounters;
#[cfg(test)]
pub(crate) use tree::{DEFAULT_LEAF, NO_CHILD, Node};
