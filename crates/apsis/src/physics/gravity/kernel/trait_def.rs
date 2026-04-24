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
//! Implementations encode a specific physical model (Plummer softening,
//! exact 1/r, truncated compact support, etc.). Future commits will attach
//! a `KernelProperties` record to each implementation that perturbations
//! can match against at registration time.

/// A scalar gravitational kernel.
///
/// All methods take the squared separation `r_squared = |Δx|²` and the
/// squared pairwise softening length `eps_squared`. Kernels that do not
/// use softening ignore the second parameter.
///
/// Implementations must be `Send + Sync` so that engines sharing a kernel
/// through [`Arc`](std::sync::Arc) can be used under parallel traversal.
pub trait Kernel: Send + Sync {
    /// Pair potential factor K(r², ε²) such that `U_ij = −G · m_i · m_j · K`.
    ///
    /// - Plummer:  K = 1/√(r² + ε²)
    /// - Newton:   K = 1/√(r²)            (ignores ε)
    fn potential(&self, r_squared: f64, eps_squared: f64) -> f64;

    /// Pair acceleration factor f(r², ε²) such that
    /// `a_i = G · m_j · f · (x_j − x_i)`.
    ///
    /// Attractive convention: f ≥ 0 for gravitationally bound pairs.
    ///
    /// - Plummer:  f = 1/(r² + ε²)^{3/2}
    /// - Newton:   f = 1/(r²)^{3/2}       (ignores ε)
    fn acceleration_factor(&self, r_squared: f64, eps_squared: f64) -> f64;
}
