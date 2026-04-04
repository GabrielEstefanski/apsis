//! Gravitational force evaluation for 2-D N-body simulations.
//!
//! ## Module layout
//!
//! | Sub-module | Responsibility |
//! |---|---|
//! | [`kernel`] | Pure Plummer-softened force/potential formulas (physics) |
//! | [`tree`] | Barnes-Hut quadtree data structure (algorithm) |
//! | [`engine`] | [`BarnesHutEngine`]: orchestrates tree + kernel (integration) |
//!
//! External code only needs to import [`G`] and [`BarnesHutEngine`]; the
//! sub-modules are implementation details and not part of the public API.

mod engine;
mod kernel;
mod tree;

pub use engine::BarnesHutEngine;
pub use kernel::G;
