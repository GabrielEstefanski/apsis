//! Plummer-softened gravitational kernel — pure physics, no simulation state.
//!
//! ## Physical model
//!
//! The Plummer softening regularises the 1/r singularity of point-mass gravity
//! by replacing the singular potential with a spherically-symmetric Plummer sphere:
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
//! softening is computed as:
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

// ── Constants ─────────────────────────────────────────────────────────────── //

/// Gravitational constant in simulation units.
///
/// All masses, lengths, and times in this simulation are expressed in a
/// unit system where G = 1.  Physical results scale trivially: multiply
/// forces by G_phys / 1 if real units are needed.
pub const G: f64 = 1.0;

// ── Kernel functions ──────────────────────────────────────────────────────── //

/// Pairwise softening squared: `ε²_ij = (ε²_i + ε²_j) / 2`.
///
/// By averaging ε² rather than ε, the kernel `F_ij(ε²_ij)` is identical for
/// both directions of the pair, preserving Newton's 3rd law exactly.
#[inline]
pub fn pair_eps2(eps_i: f64, eps_j: f64) -> f64 {
    0.5 * (eps_i * eps_i + eps_j * eps_j)
}

/// Plummer-softened acceleration of body i due to body j.
///
/// `(dx, dy) = (xⱼ − xᵢ, yⱼ − yᵢ)` — the vector **from** i **toward** j.
///
/// Returns `(aₓ, aᵧ)` in the direction of `(dx, dy)`: the force is attractive,
/// so a body always accelerates toward its neighbor.
///
/// Formula: **a** = G mⱼ Δ**r** / (|Δ**r**|² + ε²)^(3/2)
#[inline]
pub fn plummer_acc(dx: f64, dy: f64, mass_j: f64, eps2: f64) -> (f64, f64) {
    let d2 = dx * dx + dy * dy + eps2;
    let inv_r = d2.sqrt().recip();
    let fac = G * mass_j * inv_r * inv_r * inv_r;
    (dx * fac, dy * fac)
}

/// Plummer-softened gravitational potential at body i due to body j.
///
/// `(dx, dy) = (xⱼ − xᵢ, yⱼ − yᵢ)`.
///
/// Formula: `Φᵢⱼ = −G mⱼ / √(|Δr|² + ε²)`
///
/// Always negative for positive mass; finite at r = 0: `Φᵢⱼ(0) = −G mⱼ / ε`.
/// To obtain the total pair energy, multiply by mᵢ: `Eᵢⱼ = mᵢ Φᵢⱼ`.
#[inline]
pub fn plummer_phi(dx: f64, dy: f64, mass_j: f64, eps2: f64) -> f64 {
    let d2 = dx * dx + dy * dy + eps2;
    -G * mass_j * d2.sqrt().recip()
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;

    // ── pair_eps2 ───────────────────────────────────────────────────────── //

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

    // ── plummer_acc ─────────────────────────────────────────────────────── //

    /// At ε = 0 and unit separation, acceleration = G·m — exact inverse-square law.
    #[test]
    fn plummer_acc_is_exact_inverse_square_at_unit_separation() {
        let (ax, ay) = plummer_acc(1.0, 0.0, 2.0, 0.0);
        // a = G·m/r² = 1·2/1² = 2 in +x
        assert!((ax - 2.0).abs() < 1e-12);
        assert!(ay.abs() < 1e-15);
    }

    /// Force is attractive: dx > 0 ⟹ aₓ > 0 (body i accelerates toward j).
    #[test]
    fn plummer_acc_direction_is_attractive() {
        let (ax, _) = plummer_acc(3.0, 0.0, 1.0, 0.01);
        assert!(ax > 0.0);
    }

    /// At coincident positions (r = 0) with ε > 0, force is zero.
    /// The numerator of the Plummer kernel contains Δr, so a(0) = 0.
    #[test]
    fn plummer_acc_is_zero_at_coincident_positions() {
        let (ax, ay) = plummer_acc(0.0, 0.0, 5.0, 0.04);
        assert_eq!(ax, 0.0);
        assert_eq!(ay, 0.0);
    }

    /// Doubling the separation reduces the force by 4× (inverse-square law).
    /// Uses ε = 0 for exact 1/r² scaling.
    #[test]
    fn plummer_acc_scales_as_inverse_square_of_distance() {
        let (ax_near, _) = plummer_acc(1.0, 0.0, 1.0, 0.0);
        let (ax_far, _) = plummer_acc(2.0, 0.0, 1.0, 0.0);
        assert!((ax_near / ax_far - 4.0).abs() < 1e-12);
    }

    /// aₓ = G·mⱼ·dx / (dx²+dy²+ε²)^(3/2) — full analytic formula check.
    #[test]
    fn plummer_acc_matches_analytic_formula() {
        let (dx, dy, m, eps2) = (3.0, 4.0, 2.0, 0.25);
        let (ax, ay) = plummer_acc(dx, dy, m, eps2);
        // d² = 9 + 16 + 0.25 = 25.25
        let d2 = dx * dx + dy * dy + eps2;
        let inv_r3 = d2.sqrt().recip().powi(3);
        assert!((ax - G * m * dx * inv_r3).abs() < 1e-12);
        assert!((ay - G * m * dy * inv_r3).abs() < 1e-12);
    }

    // ── plummer_phi ─────────────────────────────────────────────────────── //

    /// Gravitational potential is always negative (attractive interaction).
    #[test]
    fn plummer_phi_is_always_negative() {
        assert!(plummer_phi(1.0, 0.0, 1.0, 0.0) < 0.0);
        assert!(plummer_phi(0.5, 0.5, 2.0, 0.1) < 0.0);
    }

    /// At r = 0 with ε > 0: Φ = −G·m/ε — finite, no singularity.
    #[test]
    fn plummer_phi_is_finite_at_zero_separation() {
        let (m, eps) = (3.0_f64, 0.5_f64);
        let phi = plummer_phi(0.0, 0.0, m, eps * eps);
        assert!((phi - (-G * m / eps)).abs() < 1e-12);
    }

    /// Φ = −G·mⱼ / √(dx²+dy²+ε²) — full analytic formula check.
    #[test]
    fn plummer_phi_matches_analytic_formula() {
        let (dx, dy, m, eps2) = (3.0, 4.0, 2.0, 0.25);
        let phi = plummer_phi(dx, dy, m, eps2);
        let expected = -G * m / (dx * dx + dy * dy + eps2).sqrt();
        assert!((phi - expected).abs() < 1e-12);
    }
}
