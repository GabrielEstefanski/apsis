//! Truncated Plummer kernel — a demonstrator for the `Continuity::C0` class
//! of precondition violation.
//!
//! The kernel follows the standard Plummer profile inside a cutoff radius
//! `R_c` and switches to a scaled Plummer profile outside, chosen so that
//! K itself is continuous at `R_c` but `dK/dr` is not. The pair force,
//! being proportional to `−dK/dr`, therefore has a finite jump at `R_c`,
//! and a symplectic integrator — whose derivation assumes the Hamiltonian
//! flow is smooth — produces an impulsive energy-error event each time a
//! trajectory crosses the discontinuity.
//!
//! ## Construction
//!
//! For separation `r` and pairwise softening `ε`:
//!
//! ```text
//! K(r) = 1/√(r² + ε²)                       for r < R_c
//! K(r) = α · 1/√(r² + ε²) + (1 − α) · K_c   for r ≥ R_c
//! ```
//!
//! where `K_c = 1/√(R_c² + ε²)` is the Plummer value at the cutoff and
//! `α ∈ [0, 1)` is the `outside_scale`. At `r = R_c` the two branches
//! give `K = K_c` (continuous), while `dK/dr` jumps from
//! `−R_c/(R_c² + ε²)^{3/2}` to `−α · R_c/(R_c² + ε²)^{3/2}`, a finite
//! discontinuity of `(1 − α) · R_c/(R_c² + ε²)^{3/2}`.
//!
//! With `α = 0.8` (the default chosen for the counter-test in
//! `apsis-1pn`) the outside strength is enough to keep a canonical
//! two-body orbit (equal masses, a = 1, e = 0.5) reliably bound across
//! many `R_c` crossings, yielding a reproducible series of impulsive
//! energy-error events rather than a single escape. The value is
//! deliberately above the marginal-binding threshold (`α ≈ 0.5` for the
//! canonical configuration) so that small perturbations do not tip the
//! orbit into spiralling escape.
//!
//! ## Physical invariants
//!
//! Returned from [`Kernel::properties`]:
//!
//! - [`Exactness::Modified`]: the kernel deviates from 1/r for
//!   `r ≥ R_c` regardless of softening.
//! - [`Continuity::C0`]: K is continuous, the force `−dK/dr` is not.
//!
//! Registering this kernel into a system with a perturbation declaring
//! [`KernelRequirements::exact_and_smooth`] — the 1PN Schwarzschild
//! correction does — therefore triggers both an Exactness and a
//! Continuity violation on the single `System::add_perturbation` call,
//! demonstrating that the match is compositional.
//!
//! [`Kernel::properties`]: super::Kernel::properties
//! [`KernelRequirements::exact_and_smooth`]: super::KernelRequirements::exact_and_smooth

use crate::domain::body::Body;

use super::Kernel;
use super::properties::{Continuity, Exactness, KernelProperties};

/// Plummer-softened kernel with a derivative discontinuity at a finite
/// cutoff radius `R_c`.
///
/// See the [module-level documentation](self) for the precise K(r)
/// definition and the rationale behind the outside-scale parameter.
#[derive(Debug, Clone, Copy)]
pub struct TruncatedPlummerKernel {
    r_cut: f64,
    outside_scale: f64,
}

/// Default `outside_scale` used by [`TruncatedPlummerKernel::new`].
///
/// Chosen so the canonical counter-test orbit (equal masses, a = 1,
/// e = 0.5) remains reliably bound across many `R_c` crossings: with
/// `α = 0.8` and `R_c = 1`, the bound apoapse sits near `r ≈ 2.06`,
/// comfortably inside the marginal-binding threshold at `α ≈ 0.5`.
pub const DEFAULT_TRUNCATED_OUTSIDE_SCALE: f64 = 0.8;

impl TruncatedPlummerKernel {
    /// Construct with the default [`DEFAULT_TRUNCATED_OUTSIDE_SCALE`].
    pub const fn new(r_cut: f64) -> Self {
        Self { r_cut, outside_scale: DEFAULT_TRUNCATED_OUTSIDE_SCALE }
    }

    /// Construct with a caller-specified outside scale in `[0, 1)`.
    ///
    /// `outside_scale = 0` is the hard-cutoff case (no force beyond
    /// `R_c`) and is not usable for multi-crossing scenarios.
    /// `outside_scale` approaching `1` reduces the discontinuity
    /// amplitude and with it the magnitude of the impulsive events.
    pub const fn with_outside_scale(r_cut: f64, outside_scale: f64) -> Self {
        Self { r_cut, outside_scale }
    }

    /// Cutoff radius `R_c`.
    pub const fn r_cut(&self) -> f64 {
        self.r_cut
    }

    /// Outside-scale multiplier `α`.
    pub const fn outside_scale(&self) -> f64 {
        self.outside_scale
    }
}

impl Kernel for TruncatedPlummerKernel {
    fn potential(&self, r_squared: f64, eps_squared: f64) -> f64 {
        let r_cut_sq = self.r_cut * self.r_cut;
        let inv_r = (r_squared + eps_squared).sqrt().recip();
        if r_squared < r_cut_sq {
            inv_r
        } else {
            let inv_r_cut = (r_cut_sq + eps_squared).sqrt().recip();
            self.outside_scale * inv_r + (1.0 - self.outside_scale) * inv_r_cut
        }
    }

    fn acceleration_factor(&self, r_squared: f64, eps_squared: f64) -> f64 {
        let r_cut_sq = self.r_cut * self.r_cut;
        let inv_r = (r_squared + eps_squared).sqrt().recip();
        let f_plummer = inv_r * inv_r * inv_r;
        if r_squared < r_cut_sq { f_plummer } else { self.outside_scale * f_plummer }
    }

    fn properties(&self, _bodies: &[Body]) -> KernelProperties {
        KernelProperties { exactness: Exactness::Modified, continuity: Continuity::C0 }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;

    // ── Continuity of K at r = R_c ─────────────────────────────────────── //

    #[test]
    fn potential_is_continuous_across_the_cutoff() {
        let k = TruncatedPlummerKernel::new(1.0);
        let eps_sq = 0.0;
        // r_cut = 1.0, so r_cut² = 1.0.
        let inside = k.potential(0.999_999, eps_sq);
        let outside = k.potential(1.000_001, eps_sq);
        assert!(
            (inside - outside).abs() < 1e-5,
            "K has an unexpected jump at R_c: inside={inside}, outside={outside}"
        );
    }

    #[test]
    fn potential_matches_plummer_well_inside_cutoff() {
        let k = TruncatedPlummerKernel::new(10.0); // far away cutoff
        // r² = 1, ε² = 0 → Plummer K = 1.
        let val = k.potential(1.0, 0.0);
        assert!((val - 1.0).abs() < 1e-12);
    }

    // ── Discontinuity of f (force) at r = R_c ──────────────────────────── //

    #[test]
    fn acceleration_factor_has_finite_jump_at_cutoff() {
        let k = TruncatedPlummerKernel::new(1.0);
        let eps_sq = 0.0;
        let f_inside = k.acceleration_factor(0.999_999, eps_sq);
        let f_outside = k.acceleration_factor(1.000_001, eps_sq);
        // Jump magnitude: (1 − α) · 1/(R_c)³ at r = R_c. With α = 0.8
        // (DEFAULT_TRUNCATED_OUTSIDE_SCALE) and R_c = 1, the jump is 0.2.
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
        let eps_sq = 0.0;
        // At r² = 4 (r = 2 > R_c): f_plummer = 1/8, scaled = 0.3/8 = 0.0375.
        let f = k.acceleration_factor(4.0, eps_sq);
        assert!((f - 0.0375).abs() < 1e-9);
    }

    // ── Properties ─────────────────────────────────────────────────────── //

    #[test]
    fn properties_report_modified_and_c0() {
        let k = TruncatedPlummerKernel::new(1.0);
        let props = k.properties(&[]);
        assert_eq!(props.exactness, Exactness::Modified);
        assert_eq!(props.continuity, Continuity::C0);
    }

    #[test]
    fn properties_do_not_depend_on_body_state() {
        use crate::domain::body::Body;
        let k = TruncatedPlummerKernel::new(1.0);
        let props_empty = k.properties(&[]);
        let props_with_bodies = k.properties(&[Body::star(1.0).unsoftened(), Body::rocky(1e-6)]);
        assert_eq!(props_empty, props_with_bodies);
    }

    // ── Default outside_scale sanity ───────────────────────────────────── //

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
