//! Newton kernel — gravity parameterised by a single softening length.
//!
//! ```text
//! Φ(r) = −G m / √(r² + ε²)
//! a(r) = G m r / (r² + ε²)^{3/2}    (attractive; vanishes at r = 0)
//! ```
//!
//! `ε = 0` is exact `1/r²` Newton (default). `ε > 0` is the Plummer-
//! softened regularisation — same kernel family, just smoothed at small
//! separations. The `ε → 0` limit is continuous in both potential and
//! force, so there is no qualitative split between "Newton" and
//! "Plummer": they are the same kernel parameterised by ε.
//!
//! Used as the default for paper-grade planetary, two-body, and N-body
//! work (`ε = 0`). Cluster / cosmological work where particles
//! represent matter packets opts into `ε > 0` explicitly.
//!
//! ## References
//! - Plummer (1911). *Mon. Not. R. Astron. Soc.* 71, 460–470.
//! - Dehnen & Read (2011). *Eur. Phys. J. Plus* 126, 55. (softening review)

use super::Kernel;
use super::properties::{Continuity, Exactness, KernelProperties};

/// Newton kernel with a single softening length `ε`. Default `ε = 0`
/// gives exact `1/r²` gravity.
#[derive(Debug, Clone, Copy)]
pub struct NewtonKernel {
    pub epsilon: f64,
}

impl NewtonKernel {
    pub const fn new(epsilon: f64) -> Self {
        Self { epsilon }
    }

    /// Exact Newtonian (`ε = 0`).
    pub const fn exact() -> Self {
        Self { epsilon: 0.0 }
    }
}

impl Default for NewtonKernel {
    fn default() -> Self {
        Self::exact()
    }
}

impl Kernel for NewtonKernel {
    #[inline]
    fn potential(&self, r_squared: f64) -> f64 {
        (r_squared + self.epsilon * self.epsilon).sqrt().recip()
    }

    #[inline]
    fn acceleration_factor(&self, r_squared: f64) -> f64 {
        let inv_r = (r_squared + self.epsilon * self.epsilon).sqrt().recip();
        inv_r * inv_r * inv_r
    }

    fn properties(&self) -> KernelProperties {
        let exactness = if self.epsilon == 0.0 { Exactness::Exact } else { Exactness::Softened };
        KernelProperties { exactness, continuity: Continuity::Smooth }
    }

    #[inline]
    fn is_plummer(&self) -> bool {
        self.epsilon > 0.0
    }

    #[inline]
    fn epsilon_squared(&self) -> f64 {
        self.epsilon * self.epsilon
    }

    fn variant_name(&self) -> &'static str {
        "Newton"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_matches_inverse_cube() {
        let k = NewtonKernel::exact();
        assert!((k.acceleration_factor(1.0) - 1.0).abs() < 1e-15);
        assert!((k.acceleration_factor(4.0) - 0.125).abs() < 1e-15);
    }

    #[test]
    fn exact_potential_matches_inverse_r() {
        let k = NewtonKernel::exact();
        assert!((k.potential(4.0) - 0.5).abs() < 1e-15);
    }

    #[test]
    fn softened_acc_factor_is_finite_at_zero_separation() {
        let k = NewtonKernel::new(0.2);
        let f = k.acceleration_factor(0.0);
        assert!(f.is_finite() && f > 0.0);
        let inv_eps = 0.2_f64.recip();
        assert!((f - inv_eps.powi(3)).abs() < 1e-12);
    }

    #[test]
    fn properties_track_epsilon() {
        assert_eq!(NewtonKernel::exact().properties().exactness, Exactness::Exact);
        assert_eq!(NewtonKernel::new(0.1).properties().exactness, Exactness::Softened);
    }

    #[test]
    fn is_plummer_tracks_epsilon() {
        assert!(!NewtonKernel::exact().is_plummer());
        assert!(NewtonKernel::new(0.01).is_plummer());
    }
}
