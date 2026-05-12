//! Plummer-softened gravitational kernel — the default [`Kernel`] impl.
//!
//! ## Physical model
//!
//! The Plummer softening regularises the 1/r singularity of point-mass
//! gravity by replacing the singular potential with a spherically-symmetric
//! Plummer sphere:
//!
//! ```text
//! Φ(r) = −G m / √(r² + ε²)
//! ```
//!
//! The corresponding acceleration (force per unit mass) is:
//!
//! ```text
//! a(r) = G m r / (r² + ε²)^(3/2)   [vector, attractive]
//! ```
//!
//! At `r ≫ ε` this recovers the exact inverse-square law.
//! At `r = 0` the force is zero (numerator vanishes) and the potential is
//! finite: `Φ(0) = −Gm/ε`.
//!
//! ## Pairwise softening
//!
//! When two bodies each carry a softening length (εᵢ, εⱼ), the combined
//! softening is computed (by the caller, via
//! [`pair_eps2`](super::pair_eps2)) as:
//!
//! ```text
//! ε²_ij = (ε²_i + ε²_j) / 2
//! ```
//!
//! Averaging ε² (not ε) ensures that the force on body i from body j uses
//! the same `ε²_ij` as the force on j from i — Newton's 3rd law is therefore
//! satisfied exactly at the kernel level, regardless of softening mismatch.
//!
//! ## References
//! - Plummer (1911). *Mon. Not. R. Astron. Soc.* 71, 460–470.
//! - Dehnen & Read (2011). *Eur. Phys. J. Plus* 126, 55. (softening review)

use crate::domain::body::Body;

use super::Kernel;
use super::properties::{Continuity, Exactness, KernelProperties};

/// Default Plummer-softened kernel.
///
/// Stateless: the pairwise softening `ε²_ij` is computed by the caller from
/// per-body softening lengths and passed into each trait method as
/// `eps_squared`. See [`pair_eps2`](super::pair_eps2).
#[derive(Debug, Default, Clone, Copy)]
pub struct PlummerKernel;

impl PlummerKernel {
    /// Construct a new Plummer kernel instance.
    pub const fn new() -> Self {
        Self
    }
}

impl Kernel for PlummerKernel {
    /// K(r², ε²) = 1/√(r² + ε²).
    #[inline]
    fn potential(&self, r_squared: f64, eps_squared: f64) -> f64 {
        (r_squared + eps_squared).sqrt().recip()
    }

    /// f(r², ε²) = 1/(r² + ε²)^{3/2}.
    ///
    /// Computed as `inv_r · inv_r · inv_r` rather than `inv_r.powi(3)` to
    /// match the exact floating-point sequence used in the pre-refactor
    /// implementation.
    #[inline]
    fn acceleration_factor(&self, r_squared: f64, eps_squared: f64) -> f64 {
        let inv_r = (r_squared + eps_squared).sqrt().recip();
        inv_r * inv_r * inv_r
    }

    /// Reports [`Exactness::Exact`] when every body has softening length
    /// zero, and [`Exactness::Softened`] otherwise.
    ///
    /// The rationale: a Plummer kernel multiplied by ε = 0 everywhere is
    /// mathematically indistinguishable from K(r) = 1/r, so an
    /// appropriately unsoftened configuration satisfies the Exactness
    /// invariant required by 1PN and similar extensions.
    ///
    /// Continuity is always [`Continuity::Smooth`]: Plummer's
    /// `1/√(r² + ε²)` is C^∞ on (0, ∞) regardless of ε.
    fn properties(&self, bodies: &[Body]) -> KernelProperties {
        let any_softened = bodies.iter().any(|b| b.softening > 0.0);
        KernelProperties {
            exactness: if any_softened { Exactness::Softened } else { Exactness::Exact },
            continuity: Continuity::Smooth,
        }
    }

    #[inline]
    fn is_plummer(&self) -> bool {
        true
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::super::{G, pair_eps2};
    use super::*;

    // ── pair_eps2 (helper; lives in parent module) ──────────────────────── //

    /// ε²_ij = (ε²_i + ε²_j) / 2 — explicit arithmetic check.
    #[test]
    fn pair_eps2_is_arithmetic_mean_of_squares() {
        let (eps_i, eps_j) = (0.1_f64, 0.3_f64);
        let expected = 0.5 * (0.01 + 0.09); // = 0.05
        assert!((pair_eps2(eps_i, eps_j) - expected).abs() < 1e-15);
    }

    /// ε²_ij = ε²_ji — symmetry required for Newton's 3rd law.
    #[test]
    fn pair_eps2_is_symmetric() {
        let (a, b) = (0.2_f64, 0.5_f64);
        assert_eq!(pair_eps2(a, b), pair_eps2(b, a));
    }

    // ── PlummerKernel::acceleration_factor ───────────────────────────────── //

    /// At ε = 0 and r = 1, acc factor = 1 (1/r³).
    #[test]
    fn acceleration_factor_is_one_at_unit_separation_no_softening() {
        let k = PlummerKernel::new();
        let f = k.acceleration_factor(1.0, 0.0);
        assert!((f - 1.0).abs() < 1e-12);
    }

    /// Factor is always positive (gravity is attractive).
    #[test]
    fn acceleration_factor_is_always_positive() {
        let k = PlummerKernel::new();
        assert!(k.acceleration_factor(1.0, 0.0) > 0.0);
        assert!(k.acceleration_factor(9.0, 0.01) > 0.0);
    }

    /// At r = 0 with ε > 0 the factor is finite = 1/ε³. The "zero
    /// acceleration at r = 0" property comes from multiplying the factor
    /// by Δx, which is zero in that case.
    #[test]
    fn acceleration_factor_is_finite_at_zero_separation() {
        let k = PlummerKernel::new();
        let eps2 = 0.04;
        let f = k.acceleration_factor(0.0, eps2);
        let inv_eps = eps2.sqrt().recip();
        assert!(f.is_finite());
        assert!((f - inv_eps * inv_eps * inv_eps).abs() < 1e-12);
    }

    /// Doubling r (r² → 4·r²) reduces the factor by 8× (inverse-cube
    /// scaling, ε = 0).
    #[test]
    fn acceleration_factor_scales_as_inverse_cube() {
        let k = PlummerKernel::new();
        let f_near = k.acceleration_factor(1.0, 0.0);
        let f_far = k.acceleration_factor(4.0, 0.0);
        assert!((f_near / f_far - 8.0).abs() < 1e-12);
    }

    /// a = G · m_j · f · Δx — verify full acceleration reconstruction.
    #[test]
    fn reconstructs_full_acceleration_formula() {
        let k = PlummerKernel::new();
        let (dx, dy, m, eps2) = (3.0, 4.0, 2.0, 0.25);
        let r_sq = dx * dx + dy * dy;
        let f = k.acceleration_factor(r_sq, eps2);
        let (ax, ay) = (G * m * f * dx, G * m * f * dy);
        // d² = 9 + 16 + 0.25 = 25.25
        let d2 = r_sq + eps2;
        let inv_r3 = {
            let inv_r = d2.sqrt().recip();
            inv_r * inv_r * inv_r
        };
        assert!((ax - G * m * dx * inv_r3).abs() < 1e-12);
        assert!((ay - G * m * dy * inv_r3).abs() < 1e-12);
    }

    // ── PlummerKernel::potential ────────────────────────────────────────── //

    /// K > 0 for every valid separation (sign of U is applied by caller).
    #[test]
    fn potential_factor_is_positive() {
        let k = PlummerKernel::new();
        assert!(k.potential(1.0, 0.0) > 0.0);
        assert!(k.potential(0.25, 0.01) > 0.0);
    }

    /// At r = 0 with ε > 0: K = 1/ε — finite, no singularity.
    #[test]
    fn potential_factor_is_finite_at_zero_separation() {
        let k = PlummerKernel::new();
        let eps = 0.5;
        let val = k.potential(0.0, eps * eps);
        assert!((val - 1.0 / eps).abs() < 1e-12);
    }

    /// U_ij = −G · m_i · m_j · K — verify full potential-energy
    /// reconstruction.
    #[test]
    fn reconstructs_full_potential_formula() {
        let k = PlummerKernel::new();
        let (dx, dy, m_i, m_j, eps2) = (3.0, 4.0, 1.0, 2.0, 0.25);
        let r_sq = dx * dx + dy * dy;
        let u = -G * m_i * m_j * k.potential(r_sq, eps2);
        let expected = -G * m_i * m_j / (r_sq + eps2).sqrt();
        assert!((u - expected).abs() < 1e-12);
    }

    // ── PlummerKernel::properties (dynamic, body-aware) ─────────────────── //

    #[test]
    fn properties_report_exact_for_empty_body_slice() {
        let k = PlummerKernel::new();
        let props = k.properties(&[]);
        assert_eq!(props.exactness, super::super::properties::Exactness::Exact);
        assert_eq!(props.continuity, super::super::properties::Continuity::Smooth);
    }

    #[test]
    fn properties_report_exact_when_all_bodies_unsoftened() {
        use crate::domain::body::Body;
        let k = PlummerKernel::new();
        let bodies = [
            Body::star(1.0).at(0.0, 0.0).unsoftened(),
            Body::rocky(1e-6).at(1.0, 0.0).unsoftened(),
        ];
        let props = k.properties(&bodies);
        assert_eq!(props.exactness, super::super::properties::Exactness::Exact);
    }

    #[test]
    fn properties_report_softened_when_any_body_has_positive_softening() {
        use crate::domain::body::Body;
        let k = PlummerKernel::new();
        let bodies = [
            Body::star(1.0).at(0.0, 0.0).unsoftened(), // ε = 0
            Body::rocky(1e-6).at(1.0, 0.0),            // default softening > 0
        ];
        let props = k.properties(&bodies);
        assert_eq!(props.exactness, super::super::properties::Exactness::Softened);
    }

    #[test]
    fn properties_continuity_is_always_smooth() {
        use crate::domain::body::Body;
        let k = PlummerKernel::new();
        for bodies in [
            vec![],
            vec![Body::rocky(1e-6).at(0.0, 0.0)],
            vec![Body::star(1.0).at(0.0, 0.0).unsoftened(), Body::rocky(1e-6).at(1.0, 0.0)],
        ] {
            let props = k.properties(&bodies);
            assert_eq!(
                props.continuity,
                super::super::properties::Continuity::Smooth,
                "Plummer continuity should be Smooth regardless of body state"
            );
        }
    }
}
