//! The [`Kernel`] trait — interface for pair-potential evaluation.
//!
//! A gravitational kernel K(r) defines the pair potential
//!
//! ```text
//! U_ij = −G · m_i · m_j · K(r)
//! ```
//!
//! and the corresponding pair acceleration coefficient f(r) such that
//!
//! ```text
//! a_i = G · m_j · f · (x_j − x_i).
//! ```
//!
//! Implementations encode a specific physical model (exact 1/r,
//! Plummer-softened, truncated compact support, …) and report the
//! formal invariants they satisfy via [`Kernel::properties`].
//! Extensions match their [`KernelRequirements`] against those
//! properties at registration time.
//!
//! Softening, when used, is a property of the kernel itself, not of
//! individual bodies. Different physical models (Yukawa, MOND, …)
//! parameterise their state differently and ship as separate kernel
//! impls.

use super::properties::KernelProperties;

/// A scalar gravitational kernel.
///
/// The numeric methods take the squared separation `r_squared = |Δx|²`.
/// Any softening / truncation parameter the kernel uses is stored in
/// the kernel's own state and read directly by the implementation.
///
/// Implementations must be `Send + Sync` so engines sharing a kernel
/// through [`Arc`](std::sync::Arc) can be used under parallel traversal.
pub trait Kernel: Send + Sync {
    /// Pair potential factor K(r²) such that `U_ij = −G · m_i · m_j · K`.
    ///
    /// - Newton:  K = 1/√(r²)
    /// - Plummer: K = 1/√(r² + ε²)
    fn potential(&self, r_squared: f64) -> f64;

    /// Pair acceleration factor f(r²) such that
    /// `a_i = G · m_j · f · (x_j − x_i)`.
    ///
    /// Attractive convention: f ≥ 0 for gravitationally bound pairs.
    ///
    /// - Newton:  f = 1/(r²)^{3/2}
    /// - Plummer: f = 1/(r² + ε²)^{3/2}
    fn acceleration_factor(&self, r_squared: f64) -> f64;

    /// Physical invariants this kernel provides.
    ///
    /// Static — the kernel's state alone determines exactness and
    /// continuity. Bodies don't carry kernel parameters.
    fn properties(&self) -> KernelProperties;

    /// Returns `true` iff this kernel implements the standard
    /// Plummer-softened pair potential `K = 1/√(r² + ε²)`.
    ///
    /// Used by [`BarnesHutEngine`](super::super::BarnesHutEngine) to
    /// decide whether the leaf-pair phase of the BH walk can dispatch
    /// to the hand-vectorised SIMD path. Default `false` keeps custom
    /// kernels on the scalar dyn-dispatched path.
    fn is_plummer(&self) -> bool {
        false
    }

    /// Squared softening parameter, when the kernel has one.
    ///
    /// Returns `0.0` for kernels with no softening (Newton). Plummer
    /// returns `ε²` so the SIMD fast path can read it without a
    /// downcast.
    fn epsilon_squared(&self) -> f64 {
        0.0
    }
}
