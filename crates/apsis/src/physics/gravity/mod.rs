//! Gravitational force evaluation for 2-D N-body simulations.
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

#[cfg(test)]
mod perf_2x2;

pub use engine::BarnesHutEngine;
pub use kernel::{G, Kernel, PlummerKernel, pair_eps2};
