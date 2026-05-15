//! Truncated Plummer kernel — counter-test demonstrator for the
//! `Continuity::C0` precondition class.
//!
//! Standard Plummer profile inside cutoff `R_c`, scaled Plummer outside.
//! `K(r)` is continuous at `R_c` but `dK/dr` is not, so the pair force
//! has a finite jump and a symplectic integrator emits an impulsive
//! energy-error event each crossing. Used in `apsis-1pn` tests to prove
//! the kernel-precondition system catches `Continuity::C0` violations
//! against operators that demand smooth gravity.
//!
//! ```text
//! K(r) = 1/√(r² + ε²)                       for r < R_c
//! K(r) = α · 1/√(r² + ε²) + (1 − α) · K_c   for r ≥ R_c
//! ```
//!
//! `K_c = 1/√(R_c² + ε²)` is the Plummer value at the cutoff and
//! `α ∈ [0, 1)` is the `outside_scale`. Default `α = 0.8` keeps a
//! canonical two-body orbit (equal masses, a = 1, e = 0.5) bound across
//! many crossings.

use super::Kernel;
use super::properties::{Continuity, Exactness, KernelProperties};

#[derive(Debug, Clone, Copy)]
pub struct TruncatedPlummerKernel {
    epsilon: f64,
    r_cut: f64,
    outside_scale: f64,
}

/// Default `outside_scale` — chosen so the canonical two-body counter-test
/// orbit (equal masses, a = 1, e = 0.5) stays bound across crossings.
pub const DEFAULT_TRUNCATED_OUTSIDE_SCALE: f64 = 0.8;

impl TruncatedPlummerKernel {
    pub const fn new(r_cut: f64) -> Self {
        Self { epsilon: 0.0, r_cut, outside_scale: DEFAULT_TRUNCATED_OUTSIDE_SCALE }
    }

    pub const fn with_outside_scale(r_cut: f64, outside_scale: f64) -> Self {
        Self { epsilon: 0.0, r_cut, outside_scale }
    }

    pub const fn with_epsilon(mut self, epsilon: f64) -> Self {
        self.epsilon = epsilon;
        self
    }

    pub const fn r_cut(&self) -> f64 {
        self.r_cut
    }

    pub const fn outside_scale(&self) -> f64 {
        self.outside_scale
    }

    pub const fn epsilon(&self) -> f64 {
        self.epsilon
    }
}

impl Kernel for TruncatedPlummerKernel {
    fn potential(&self, r_squared: f64) -> f64 {
        let eps_sq = self.epsilon * self.epsilon;
        let r_cut_sq = self.r_cut * self.r_cut;
        let inv_r = (r_squared + eps_sq).sqrt().recip();
        if r_squared < r_cut_sq {
            inv_r
        } else {
            let inv_r_cut = (r_cut_sq + eps_sq).sqrt().recip();
            self.outside_scale * inv_r + (1.0 - self.outside_scale) * inv_r_cut
        }
    }

    fn acceleration_factor(&self, r_squared: f64) -> f64 {
        let eps_sq = self.epsilon * self.epsilon;
        let r_cut_sq = self.r_cut * self.r_cut;
        let inv_r = (r_squared + eps_sq).sqrt().recip();
        let f_plummer = inv_r * inv_r * inv_r;
        if r_squared < r_cut_sq { f_plummer } else { self.outside_scale * f_plummer }
    }

    fn properties(&self) -> KernelProperties {
        KernelProperties { exactness: Exactness::Modified, continuity: Continuity::C0 }
    }

    #[inline]
    fn epsilon_squared(&self) -> f64 {
        self.epsilon * self.epsilon
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn potential_is_continuous_across_the_cutoff() {
        let k = TruncatedPlummerKernel::new(1.0);
        let inside = k.potential(0.999_999);
        let outside = k.potential(1.000_001);
        assert!(
            (inside - outside).abs() < 1e-5,
            "K has an unexpected jump at R_c: inside={inside}, outside={outside}"
        );
    }

    #[test]
    fn potential_matches_plummer_well_inside_cutoff() {
        let k = TruncatedPlummerKernel::new(10.0);
        let val = k.potential(1.0);
        assert!((val - 1.0).abs() < 1e-12);
    }

    #[test]
    fn acceleration_factor_has_finite_jump_at_cutoff() {
        let k = TruncatedPlummerKernel::new(1.0);
        let f_inside = k.acceleration_factor(0.999_999);
        let f_outside = k.acceleration_factor(1.000_001);
        let jump = f_inside - f_outside;
        let expected = 1.0 - DEFAULT_TRUNCATED_OUTSIDE_SCALE;
        assert!(
            (jump - expected).abs() < 1e-3,
            "expected jump of {expected} at R_c (α = {DEFAULT_TRUNCATED_OUTSIDE_SCALE}), got {jump}"
        );
    }

    #[test]
    fn acceleration_factor_is_scaled_plummer_outside_cutoff() {
        let k = TruncatedPlummerKernel::with_outside_scale(1.0, 0.3);
        let f = k.acceleration_factor(4.0);
        assert!((f - 0.0375).abs() < 1e-9);
    }

    #[test]
    fn properties_report_modified_and_c0() {
        let k = TruncatedPlummerKernel::new(1.0);
        let props = k.properties();
        assert_eq!(props.exactness, Exactness::Modified);
        assert_eq!(props.continuity, Continuity::C0);
    }

    #[test]
    fn default_outside_scale_matches_documented_constant() {
        let k = TruncatedPlummerKernel::new(1.0);
        assert_eq!(k.outside_scale(), DEFAULT_TRUNCATED_OUTSIDE_SCALE);
    }

    #[test]
    fn with_outside_scale_preserves_value() {
        let k = TruncatedPlummerKernel::with_outside_scale(2.0, 0.7);
        assert_eq!(k.r_cut(), 2.0);
        assert_eq!(k.outside_scale(), 0.7);
    }
}
