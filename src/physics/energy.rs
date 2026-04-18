//! Conservative-quantity observables for an N-body system.
//!
//! All functions are pure: they accept a slice of [`Body`] values (or plain
//! scalars) and return a scalar. No simulation state is modified.

use crate::domain::body::Body;
use crate::physics::gravity::{G, pair_eps2};

/// Total kinetic energy of the system: `KE = 1/2 sum(m v^2)`.
pub fn kinetic_energy(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy)).sum()
}

/// Z-component of the orbital angular momentum: `Lz = sum(m (x vy - y vx))`.
pub fn angular_momentum_z(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| b.mass * (b.x * b.vy - b.y * b.vx)).sum()
}

/// Total mechanical energy: `E = KE + PE`.
pub fn total_energy(kinetic: f64, potential: f64) -> f64 {
    kinetic + potential
}

/// Center-of-mass position and velocity: `(x_com, y_com, vx_com, vy_com)`.
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

/// Computes the gravitational potential energy contributed to each body.
///
/// For each body `i`, `pe[i]` receives half of the pairwise potential with
/// every other body `j`:
///
/// ```text
/// pe_i = ½ Σⱼ≠ᵢ  −G mᵢ mⱼ / √(rᵢⱼ² + ε²ᵢⱼ)
/// ```
///
/// Summing `pe` over all `i` recovers the total potential energy exactly.
/// This is the standard symmetric partition used in N-body diagnostics.
///
/// Complexity: O(N²). Intended for recording/export, **not** the inner step loop.
pub fn per_body_potential_energy(bodies: &[Body], g_factor: f64) -> Vec<f64> {
    let n = bodies.len();
    let mut pe = vec![0.0_f64; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = bodies[j].x - bodies[i].x;
            let dy = bodies[j].y - bodies[i].y;
            let eps2 = pair_eps2(bodies[i].softening, bodies[j].softening);
            let d2 = dx * dx + dy * dy + eps2;
            let phi = -G * g_factor * bodies[i].mass * bodies[j].mass / d2.sqrt();
            pe[i] += phi * 0.5;
            pe[j] += phi * 0.5;
        }
    }

    pe
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    #[test]
    fn kinetic_energy_single_body() {
        let b = Body::new(0.0, 0.0, 3.0, 4.0, 2.0, crate::domain::materials::Material::Rocky);
        assert!((kinetic_energy(&[b]) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn kinetic_energy_at_rest_is_zero() {
        let b = Body::new(1.0, 2.0, 0.0, 0.0, 5.0, crate::domain::materials::Material::Rocky);
        assert_eq!(kinetic_energy(&[b]), 0.0);
    }

    #[test]
    fn kinetic_energy_is_nonnegative() {
        let b = Body::new(0.0, 0.0, -3.0, 4.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(kinetic_energy(&[b]) >= 0.0);
    }

    #[test]
    fn kinetic_energy_is_additive() {
        let b1 = Body::new(0.0, 0.0, 1.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 0.0, 0.0, 2.0, 2.0, crate::domain::materials::Material::Rocky);
        assert!((kinetic_energy(&[b1, b2]) - 4.5).abs() < 1e-12);
    }

    #[test]
    fn angular_momentum_z_circular_orbit() {
        let (r, v, m) = (3.0, 2.0, 4.0);
        let b = Body::new(r, 0.0, 0.0, v, m, crate::domain::materials::Material::Rocky);
        assert!((angular_momentum_z(&[b]) - m * r * v).abs() < 1e-12);
    }

    #[test]
    fn angular_momentum_z_positive_for_ccw() {
        let b = Body::new(1.0, 0.0, 0.0, 1.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(angular_momentum_z(&[b]) > 0.0);
    }

    #[test]
    fn angular_momentum_z_negative_for_cw() {
        let b = Body::new(1.0, 0.0, 0.0, -1.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!(angular_momentum_z(&[b]) < 0.0);
    }

    #[test]
    fn angular_momentum_z_is_additive() {
        let b1 = Body::new(1.0, 0.0, 0.0, 1.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 2.0, -1.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        assert!((angular_momentum_z(&[b1, b2]) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn total_energy_is_sum_of_ke_and_pe() {
        assert!((total_energy(3.0, -5.0) - (-2.0)).abs() < 1e-12);
        assert!((total_energy(0.0, -7.0) - (-7.0)).abs() < 1e-12);
    }

    #[test]
    fn com_position_is_midpoint_for_equal_masses() {
        let b1 = Body::new(0.0, 0.0, 0.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(4.0, 2.0, 0.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let (cx, cy, _, _) = center_of_mass_state(&[b1, b2]);
        assert!((cx - 2.0).abs() < 1e-12);
        assert!((cy - 1.0).abs() < 1e-12);
    }

    #[test]
    fn com_velocity_is_mass_weighted_mean() {
        let b1 = Body::new(0.0, 0.0, 4.0, 0.0, 1.0, crate::domain::materials::Material::Rocky);
        let b2 = Body::new(0.0, 0.0, 0.0, 0.0, 3.0, crate::domain::materials::Material::Rocky);
        let (_, _, vx, vy) = center_of_mass_state(&[b1, b2]);
        assert!((vx - 1.0).abs() < 1e-12);
        assert!(vy.abs() < 1e-12);
    }

    #[test]
    fn com_of_empty_slice_returns_zero() {
        assert_eq!(center_of_mass_state(&[]), (0.0, 0.0, 0.0, 0.0));
    }
}
