//! Conservative-quantity observables for an N-body system.
//!
//! All functions are **pure**: they accept a slice of [`Body`] values (or plain
//! scalars) and return a scalar — no simulation state is modified.
//!
//! ## Physical identities encoded
//!
//! | Function | Formula |
//! |---|---|
//! | [`kinetic_energy`] | KE = ½ Σᵢ mᵢ (vxᵢ² + vyᵢ²) |
//! | [`angular_momentum_z`] | Lz = Σᵢ mᵢ (xᵢ vyᵢ − yᵢ vxᵢ) |
//! | [`total_energy`] | E = KE + PE |
//! | [`center_of_mass_state`] | r_com = Σᵢ mᵢ rᵢ / M, v_com = Σᵢ mᵢ vᵢ / M |

use crate::domain::body::Body;

/// Total kinetic energy of the system: KE = ½ Σᵢ mᵢ (vxᵢ² + vyᵢ²).
///
/// Always non-negative; equals zero if and only if every body is at rest.
pub fn kinetic_energy(bodies: &[Body]) -> f64 {
    bodies
        .iter()
        .map(|b| 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy))
        .sum()
}

/// Z-component of the orbital angular momentum: Lz = Σᵢ mᵢ (xᵢ vyᵢ − yᵢ vxᵢ).
///
/// This is the 2-D cross product r × (m·v) projected onto the z-axis, summed
/// over all bodies.  Positive for counter-clockwise bulk rotation (right-hand
/// rule with z out of the screen), negative for clockwise.
pub fn angular_momentum_z(bodies: &[Body]) -> f64 {
    bodies
        .iter()
        .map(|b| b.mass * (b.x * b.vy - b.y * b.vx))
        .sum()
}

/// Total mechanical energy: E = KE + PE.
///
/// The potential energy `pe` must be supplied by the caller (evaluated by the
/// force engine) because computing it requires O(N²) pair sums.
pub fn total_energy(kinetic: f64, potential: f64) -> f64 {
    kinetic + potential
}

/// Center-of-mass position and velocity: `(x_com, y_com, vx_com, vy_com)`.
///
/// Uses mass-weighted averages:
/// - r_com = Σᵢ mᵢ rᵢ / M
/// - v_com = Σᵢ mᵢ vᵢ / M
///
/// Returns `(0, 0, 0, 0)` for an empty or zero-mass slice.
pub fn center_of_mass_state(bodies: &[Body]) -> (f64, f64, f64, f64) {
    let mut m = 0.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;

    for b in bodies {
        m += b.mass;
        x += b.mass * b.x;
        y += b.mass * b.y;
        vx += b.mass * b.vx;
        vy += b.mass * b.vy;
    }

    if m == 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }

    (x / m, y / m, vx / m, vy / m)
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    // ── Kinetic energy ──────────────────────────────────────────────────── //

    /// KE = ½ m v²  for a single body.
    /// v = (3, 4) → v² = 25; KE = 0.5 × 2 × 25 = 25.
    #[test]
    fn kinetic_energy_single_body() {
        let b = Body::new(0.0, 0.0, 3.0, 4.0, 2.0, crate::domain::materials::Material::Rocky);
        assert!((kinetic_energy(&[b]) - 25.0).abs() < 1e-12);
    }

    /// KE = 0 when all bodies are at rest (vx = vy = 0).
    #[test]
    fn kinetic_energy_at_rest_is_zero() {
        let b = Body::new(1.0, 2.0, 0.0, 0.0, 5.0, crate::domain::materials::Material::Rocky);
        assert_eq!(kinetic_energy(&[b]), 0.0);
    }

    /// KE ≥ 0 always (velocity components appear only as squares).
    #[test]
    fn kinetic_energy_is_nonnegative() {
        let b = Body::new(0.0, 0.0, -3.0, 4.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(kinetic_energy(&[b]) >= 0.0);
    }

    /// KE is additive over disjoint subsets: KE(A ∪ B) = KE(A) + KE(B).
    #[test]
    fn kinetic_energy_is_additive() {
        // KE₁ = ½·1·1² = 0.5;  KE₂ = ½·2·2² = 4.0
        let b1 = Body::new(0.0, 0.0, 1.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 0.0, 0.0, 2.0, 2.0, crate::domain::materials::Material::Rocky);
        assert!((kinetic_energy(&[b1, b2]) - 4.5).abs() < 1e-12);
    }

    // ── Angular momentum ────────────────────────────────────────────────── //

    /// Lz = m r v  for a body in a circular orbit in the xy-plane.
    /// Body at (r, 0) with velocity (0, v): Lz = m (r·v − 0) = m r v.
    #[test]
    fn angular_momentum_z_circular_orbit() {
        let (r, v, m) = (3.0, 2.0, 4.0);
        let b = Body::new(r, 0.0, 0.0, v, m, crate::domain::materials::Material::Rocky);
        assert!((angular_momentum_z(&[b]) - m * r * v).abs() < 1e-12);
    }

    /// Lz > 0 for counter-clockwise motion (right-hand rule, +z out of screen).
    /// Body at (1, 0) with velocity in +y: CCW orbit.
    #[test]
    fn angular_momentum_z_positive_for_ccw() {
        let b = Body::new(1.0, 0.0, 0.0, 1.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(angular_momentum_z(&[b]) > 0.0);
    }

    /// Lz < 0 for clockwise motion.
    /// Body at (1, 0) with velocity in −y: CW orbit.
    #[test]
    fn angular_momentum_z_negative_for_cw() {
        let b = Body::new(1.0, 0.0, 0.0, -1.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(angular_momentum_z(&[b]) < 0.0);
    }

    /// Lz is additive over the slice: total Lz = Σ Lzᵢ.
    ///
    /// b1 at (1,0) v=(0,1):   Lz₁ = 1·(1·1 − 0·0) = 1
    /// b2 at (0,2) v=(−1,0):  Lz₂ = 1·(0·0 − 2·(−1)) = 2
    /// Total = 3.
    #[test]
    fn angular_momentum_z_is_additive() {
        let b1 = Body::new(1.0, 0.0, 0.0, 1.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 2.0, -1.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!((angular_momentum_z(&[b1, b2]) - 3.0).abs() < 1e-12);
    }

    // ── Total energy ────────────────────────────────────────────────────── //

    /// E = KE + PE — a trivial identity; verify for both signs of PE.
    #[test]
    fn total_energy_is_sum_of_ke_and_pe() {
        assert!((total_energy(3.0, -5.0) - (-2.0)).abs() < 1e-12);
        assert!((total_energy(0.0, -7.0) - (-7.0)).abs() < 1e-12);
    }

    // ── Center of mass ──────────────────────────────────────────────────── //

    /// r_com = midpoint for two equal-mass bodies.
    #[test]
    fn com_position_is_midpoint_for_equal_masses() {
        let b1 = Body::new(0.0, 0.0, 0.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(4.0, 2.0, 0.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let (cx, cy, _, _) = center_of_mass_state(&[b1, b2]);
        assert!((cx - 2.0).abs() < 1e-12);
        assert!((cy - 1.0).abs() < 1e-12);
    }

    /// v_com = Σ mᵢ vᵢ / M (mass-weighted mean).
    /// m₁=1 v₁=(4,0), m₂=3 v₂=(0,0) → v_com = (1, 0).
    #[test]
    fn com_velocity_is_mass_weighted_mean() {
        let b1 = Body::new(0.0, 0.0, 4.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 0.0, 0.0, 0.0, 3.0, crate::domain::materials::Material::Rocky);
        let (_, _, vx, vy) = center_of_mass_state(&[b1, b2]);
        assert!((vx - 1.0).abs() < 1e-12);
        assert!(vy.abs() < 1e-12);
    }

    /// An empty slice returns (0, 0, 0, 0) — no division by zero.
    #[test]
    fn com_of_empty_slice_returns_zero() {
        assert_eq!(center_of_mass_state(&[]), (0.0, 0.0, 0.0, 0.0));
    }
}
