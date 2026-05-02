//! Pure radiation force kernels.
//!
//! All functions are stateless and operate exclusively on scalar physical
//! quantities — no [`Body`], no [`System`], no simulation state.  This makes
//! them straightforward to unit-test and to reuse outside the N-body context.
//!
//! # Physical model
//!
//! A radiating point source emits an isotropic photon flux.  The momentum
//! carried by that flux exerts two distinct effects on a small body:
//!
//! ## Direct radiation pressure
//!
//! The flux at distance `r` from a source of luminosity `L` is
//!
//! ```text
//! F = L / (4π r²)
//! ```
//!
//! The force on a body with cross-section `A`, efficiency `Q_pr`, and mass `m`
//! is directed radially away from the source:
//!
//! ```text
//! a_rad = (F · Q_pr · A) / (m · c)  ·  r̂
//! ```
//!
//! ## Poynting–Robertson drag
//!
//! In the rest frame of the source an orbiting body sees an aberrated flux.
//! To first order in `v/c` this produces an additional velocity-dependent
//! term:
//!
//! ```text
//! a_PR = a_mag · (r̂  −  v_rel / c)
//! ```
//!
//! where `a_mag` is the scalar magnitude from direct pressure and
//! `v_rel = v_body − v_source`.  The `r̂` term is the direct pressure;
//! the `−v_rel/c` term contains:
//!
//! - a **radial** component that slightly reduces the effective push, and
//! - a **tangential** component that removes orbital energy and angular
//!   momentum, driving the secular inspiral known as the Poynting–Robertson
//!   effect.
//!
//! Setting `v_rel = 0` recovers pure radiation pressure exactly.
//!
//! # References
//!
//! - Robertson, H. P. (1937). *Mon. Not. R. Astron. Soc.* 97, 423–438.
//! - Burns, J. A., Lamy, P. L., & Soter, S. (1979). *Icarus* 40, 1–48.
//!   Equation (5) is the canonical form implemented here.

use std::f64::consts::PI;

use crate::math::Vec3;
use crate::physics::radiation::{RadiationParams, RadiationSource};

// ── Kernels ───────────────────────────────────────────────────────────────────

/// Computes the direct radiation pressure acceleration on a body.
///
/// The acceleration is directed radially away from the source and falls
/// off as the inverse square of the distance:
///
/// ```text
/// a_rad = (L · Q_pr · A) / (4π r² · c · m)  ·  r̂
/// ```
///
/// Returns `Vec3::ZERO` when the body is coincident with the source
/// (`r² < 1e-30`).
pub fn radiation_acceleration(
    pos: Vec3,
    params: &RadiationParams,
    source: &RadiationSource,
) -> Vec3 {
    let dx = pos.x - source.x;
    let dy = pos.y - source.y;
    let dz = pos.z - source.z;
    let r2 = dx * dx + dy * dy + dz * dz;

    if r2 < 1e-30 {
        return Vec3::ZERO;
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;
    let flux = source.luminosity / (4.0 * PI * r2);
    let a_mag = flux * params.q_pr * params.area / (params.mass * source.c);

    Vec3::new(dx * inv_r * a_mag, dy * inv_r * a_mag, dz * inv_r * a_mag)
}

/// Computes the full Poynting–Robertson acceleration, comprising both direct
/// radiation pressure and the velocity-dependent drag term.
///
/// The first-order relativistic expression in the rest frame of the source is
///
/// ```text
/// a_PR = a_mag · (r̂  −  v_rel / c)
/// ```
///
/// where `v_rel = v_body − v_source` (Burns et al. 1979, eq. 5).
///
/// Returns `Vec3::ZERO` near the source.
pub fn pr_drag_acceleration(
    pos: Vec3,
    vel: Vec3,
    params: &RadiationParams,
    source: &RadiationSource,
) -> Vec3 {
    let dx = pos.x - source.x;
    let dy = pos.y - source.y;
    let dz = pos.z - source.z;
    let dvx = vel.x - source.vx;
    let dvy = vel.y - source.vy;
    let dvz = vel.z - source.vz;
    let r2 = dx * dx + dy * dy + dz * dz;

    if r2 < 1e-30 {
        return Vec3::ZERO;
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;
    let inv_c = 1.0 / source.c;
    let flux = source.luminosity / (4.0 * PI * r2);
    let a_mag = flux * params.q_pr * params.area / (params.mass * source.c);

    // a_PR = a_mag · (r̂ − v_rel / c)
    Vec3::new(
        a_mag * (dx * inv_r - dvx * inv_c),
        a_mag * (dy * inv_r - dvy * inv_c),
        a_mag * (dz * inv_r - dvz * inv_c),
    )
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// A minimal source at the origin with unit luminosity and c = 1.
    fn unit_source() -> RadiationSource {
        RadiationSource {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            luminosity: 1.0,
            c: 1.0,
        }
    }

    /// Unit radiation parameters: area = 1, mass = 1, Q_pr = 1.
    fn unit_params() -> RadiationParams {
        RadiationParams { area: 1.0, mass: 1.0, q_pr: 1.0 }
    }

    // ── radiation_acceleration ────────────────────────────────────────────────

    #[test]
    fn rad_points_away_from_source() {
        let a = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &unit_source());
        assert!(a.x > 0.0, "radial component should be positive (away from source)");
        assert!(a.y.abs() < 1e-15, "no y component for body on x-axis");
        assert!(a.z.abs() < 1e-15, "no z component for body on x-axis");
    }

    #[test]
    fn rad_inverse_square_falloff() {
        let a1 = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &unit_source());
        let a2 = radiation_acceleration(Vec3::new(2.0, 0.0, 0.0), &unit_params(), &unit_source());
        let ratio = a1.x / a2.x;
        assert!((ratio - 4.0).abs() < 1e-12, "expected r⁻² falloff (ratio = 4), got {ratio}");
    }

    #[test]
    fn rad_scales_linearly_with_luminosity() {
        let s2 = RadiationSource { luminosity: 2.0, ..unit_source() };
        let a1 = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &unit_source());
        let a2 = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &s2);
        assert!((a2.x / a1.x - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rad_scales_linearly_with_q_pr() {
        let p2 = RadiationParams { q_pr: 2.0, ..unit_params() };
        let a1 = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &unit_source());
        let a2 = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &p2, &unit_source());
        assert!((a2.x / a1.x - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rad_zero_luminosity_gives_zero() {
        let s = RadiationSource { luminosity: 0.0, ..unit_source() };
        let a = radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &s);
        assert_eq!(a, Vec3::ZERO);
    }

    #[test]
    fn rad_coincident_body_returns_zero() {
        let a = radiation_acceleration(Vec3::ZERO, &unit_params(), &unit_source());
        assert_eq!(a, Vec3::ZERO);
    }

    #[test]
    fn rad_diagonal_direction_is_normalised() {
        // Body at (1, 1): force should point along (1, 1) / √2.
        let a = radiation_acceleration(Vec3::new(1.0, 1.0, 0.0), &unit_params(), &unit_source());
        let mag = a.length();
        assert!((a.x / mag - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
        assert!((a.y / mag - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn rad_3d_direction_is_normalised() {
        // Body at (1, 1, 1): force should point along (1, 1, 1) / √3.
        let a = radiation_acceleration(Vec3::new(1.0, 1.0, 1.0), &unit_params(), &unit_source());
        let mag = a.length();
        let expected = 1.0 / 3.0_f64.sqrt();
        assert!((a.x / mag - expected).abs() < 1e-12);
        assert!((a.y / mag - expected).abs() < 1e-12);
        assert!((a.z / mag - expected).abs() < 1e-12);
    }

    // ── pr_drag_acceleration ──────────────────────────────────────────────────

    #[test]
    fn pr_at_rest_equals_radiation_pressure() {
        // A body at rest relative to the source experiences only direct pressure.
        let a_rad =
            radiation_acceleration(Vec3::new(1.0, 0.0, 0.0), &unit_params(), &unit_source());
        let a_pr = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            &unit_params(),
            &unit_source(),
        );
        assert!((a_rad.x - a_pr.x).abs() < 1e-15);
        assert!((a_rad.y - a_pr.y).abs() < 1e-15);
        assert!((a_rad.z - a_pr.z).abs() < 1e-15);
    }

    #[test]
    fn pr_tangential_drag_opposes_orbital_velocity() {
        // Body at (1, 0, 0) moving in +y (circular orbit direction).
        // The tangential drag component should be in −y.
        let a = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            &unit_params(),
            &unit_source(),
        );
        assert!(a.y < 0.0, "PR drag should decelerate the tangential velocity");
    }

    #[test]
    fn pr_radial_infall_increases_radial_push() {
        // Body falling radially inward (vx < 0) should see a reduced net
        // outward push compared to a body at rest (aberration reduces flux).
        let a_rest = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            &unit_params(),
            &unit_source(),
        );
        let a_infall = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-0.1, 0.0, 0.0),
            &unit_params(),
            &unit_source(),
        );
        assert!(
            a_infall.x > a_rest.x,
            "inward radial velocity should increase apparent flux (aberration)"
        );
    }

    #[test]
    fn pr_galilean_invariance() {
        // The force depends only on the relative velocity v_rel = v_body − v_source.
        // Therefore these two scenarios must give identical accelerations:
        //   (a) body at vel = +0.5, source at rest
        //   (b) body at rest,        source at vel = −0.5
        let s_moving_neg = RadiationSource { vx: -0.5, ..unit_source() };
        let a_body_moving = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            &unit_params(),
            &unit_source(),
        );
        let a_source_moving_neg = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            &unit_params(),
            &s_moving_neg,
        );
        assert_eq!(
            a_body_moving, a_source_moving_neg,
            "PR force must depend only on v_rel: body@+v ≡ source@−v"
        );
    }

    #[test]
    fn pr_zero_luminosity_gives_zero() {
        let s = RadiationSource { luminosity: 0.0, ..unit_source() };
        let a = pr_drag_acceleration(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 0.3, 0.0),
            &unit_params(),
            &s,
        );
        assert_eq!(a, Vec3::ZERO);
    }

    #[test]
    fn pr_coincident_body_returns_zero() {
        let a = pr_drag_acceleration(
            Vec3::ZERO,
            Vec3::new(1.0, 1.0, 1.0),
            &unit_params(),
            &unit_source(),
        );
        assert_eq!(a, Vec3::ZERO);
    }
}
