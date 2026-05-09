//! Gravitational force evaluation for 3D N-body simulations.
//!
//! ## Module layout
//!
//! | Sub-module | Responsibility |
//! |---|---|
//! | [`kernel`] | [`Kernel`] trait and concrete implementations (Plummer default) |
//! | [`tree`] | Barnes-Hut octree data structure (algorithm) |
//! | [`morton`] | Z-order spatial encoding for body insertion + walk locality |
//! | [`engine`] | [`BarnesHutEngine`]: orchestrates tree + kernel (integration) |

mod engine;
pub mod kernel;
mod morton;
mod tree;

#[cfg(test)]
mod perf_2x2;

pub use engine::BarnesHutEngine;
pub use kernel::{G, Kernel, PlummerKernel, pair_eps2};
