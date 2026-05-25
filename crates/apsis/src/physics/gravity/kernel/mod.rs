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
//! implementations. The default is [`NewtonKernel`], parameterised by a
//! single softening length `ε`: `ε = 0` is exact `1/r²` Newton; `ε > 0`
//! is the Plummer-softened regularisation. The `ε → 0` limit is
//! continuous, so Newton and Plummer share one impl.
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
//! | `newton` (private)      | [`NewtonKernel`] — `ε`-parameterised Newton/Plummer |
//! | `truncated` (private)   | [`TruncatedPlummerKernel`] — counter-test fixture |
//! | `properties` (private)  | Invariant types and matching logic |

mod newton;
mod properties;
mod trait_def;
mod truncated;

pub use newton::NewtonKernel;
pub use properties::{
    Continuity, Exactness, KernelProperties, KernelRequirements, RequirementViolation,
};
pub use trait_def::Kernel;
pub use truncated::{DEFAULT_TRUNCATED_OUTSIDE_SCALE, TruncatedPlummerKernel};

// ── Constants ─────────────────────────────────────────────────────────────── //

/// Gravitational constant in simulation units.
///
/// All masses, lengths, and times in this simulation are expressed in a
/// unit system where G = 1. Physical results scale trivially: multiply
/// forces by `G_phys / 1` if real units are needed.
pub const G: f64 = 1.0;
