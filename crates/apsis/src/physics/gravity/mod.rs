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
#[cfg(target_arch = "x86_64")]
mod simd;
mod tree;

pub use engine::BarnesHutEngine;
#[cfg(test)]
pub(crate) use engine::LeafPairKernel;
pub use kernel::{G, Kernel, PlummerKernel, pair_eps2};
