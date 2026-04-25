//! Gravitational kernels — the pair potential K(r) defining pairwise
//! forces in the N-body simulation.
//!
//! ## Overview
//!
//! A kernel is characterised by a scalar function K: ℝ₊ → ℝ₊ such that
//! the pair potential energy between two bodies is
//!
//! ```text
//! U_ij = −G · m_i · m_j · K(r),     r = |x_i − x_j|
//! ```
//!
//! and the corresponding pair acceleration on body *i* is
//!
//! ```text
//! a_i = G · m_j · f(r) · (x_j − x_i),     f(r) = −K'(r) / r.
//! ```
//!
//! Different physical models of K correspond to different [`Kernel`]
//! implementations. The default is [`PlummerKernel`], which softens the
//! 1/r singularity with a spherically-symmetric Plummer sphere.
//!
//! ## Extension contract
//!
//! Kernels declare the physical invariants they provide via
//! [`KernelProperties`]; extensions (perturbations, integrators with
//! kernel preconditions) declare the invariants they require via
//! [`KernelRequirements`]. The match between the two is computed by
//! [`KernelRequirements::check_against`] and surfaces as
//! [`RequirementViolation`] records.
//!
//! ## Module layout
//!
//! | Sub-module | Responsibility |
//! |---|---|
//! | `trait_def` (private)   | The [`Kernel`] trait definition |
//! | `plummer` (private)     | [`PlummerKernel`] — the default impl |
//! | `properties` (private)  | Invariant types and matching logic |

mod plummer;
mod properties;
mod trait_def;
mod truncated;

pub use plummer::PlummerKernel;
pub use properties::{
    Continuity, Exactness, KernelProperties, KernelRequirements, RequirementViolation,
};
pub use trait_def::Kernel;
pub use truncated::{DEFAULT_TRUNCATED_OUTSIDE_SCALE, TruncatedPlummerKernel};

// ── Constants and helpers ─────────────────────────────────────────────────── //

/// Gravitational constant in simulation units.
///
/// All masses, lengths, and times in this simulation are expressed in a
/// unit system where G = 1. Physical results scale trivially: multiply
/// forces by `G_phys / 1` if real units are needed.
pub const G: f64 = 1.0;

/// Pairwise softening squared: ε²_ij = (ε²_i + ε²_j) / 2.
///
/// Specific to Plummer-style per-body softening. By averaging ε² rather
/// than ε, the kernel produces identical forces in both directions of the
/// pair, preserving Newton's 3rd law exactly.
#[inline]
pub fn pair_eps2(eps_i: f64, eps_j: f64) -> f64 {
    0.5 * (eps_i * eps_i + eps_j * eps_j)
}
