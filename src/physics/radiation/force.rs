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
/// # Arguments
///
/// - `pos` — inertial position of the body `[x, y]`
/// - `params` — radiation interaction parameters of the body
/// - `source` — radiating source
///
/// # Returns
///
/// Acceleration `[ax, ay]` in internal units.
/// Returns `[0, 0]` when the body is coincident with the source (`r² < 1e-30`).
pub fn radiation_acceleration(
    pos: [f64; 2],
    params: &RadiationParams,
    source: &RadiationSource,
) -> [f64; 2] {
    let dx = pos[0] - source.x;
    let dy = pos[1] - source.y;
    let r2 = dx * dx + dy * dy;

    if r2 < 1e-30 {
        return [0.0, 0.0];
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;
    let flux = source.luminosity / (4.0 * PI * r2);
    let a_mag = flux * params.q_pr * params.area / (params.mass * source.c);

    [dx * inv_r * a_mag, dy * inv_r * a_mag]
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
/// # Arguments
///
/// - `pos` — inertial position of the body `[x, y]`
/// - `vel` — inertial velocity of the body `[vx, vy]`
/// - `params` — radiation interaction parameters of the body
/// - `source` — radiating source; must carry velocity for aberration correction
///
/// # Returns
///
/// Acceleration `[ax, ay]`. Returns `[0, 0]` near the source.
pub fn pr_drag_acceleration(
    pos: [f64; 2],
    vel: [f64; 2],
    params: &RadiationParams,
    source: &RadiationSource,
) -> [f64; 2] {
    let dx = pos[0] - source.x;
    let dy = pos[1] - source.y;
    let dvx = vel[0] - source.vx;
    let dvy = vel[1] - source.vy;
    let r2 = dx * dx + dy * dy;

    if r2 < 1e-30 {
        return [0.0, 0.0];
    }

    let r = r2.sqrt();
    let inv_r = 1.0 / r;
    let inv_c = 1.0 / source.c;
    let flux = source.luminosity / (4.0 * PI * r2);
    let a_mag = flux * params.q_pr * params.area / (params.mass * source.c);

    // a_PR = a_mag · (r̂ − v_rel / c)
    [
        a_mag * (dx * inv_r - dvx * inv_c),
        a_mag * (dy * inv_r - dvy * inv_c),
    ]
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
            vx: 0.0,
            vy: 0.0,
            luminosity: 1.0,
            c: 1.0,
        }
    }

    /// Unit radiation parameters: area = 1, mass = 1, Q_pr = 1.
    fn unit_params() -> RadiationParams {
        RadiationParams {
            area: 1.0,
            mass: 1.0,
            q_pr: 1.0,
        }
    }

    // ── radiation_acceleration ────────────────────────────────────────────────

    #[test]
    fn rad_points_away_from_source() {
        let a = radiation_acceleration([1.0, 0.0], &unit_params(), &unit_source());
        assert!(
            a[0] > 0.0,
            "radial component should be positive (away from source)"
        );
        assert!(
            a[1].abs() < 1e-15,
            "no tangential component for body on x-axis"
        );
    }

    #[test]
    fn rad_inverse_square_falloff() {
        let a1 = radiation_acceleration([1.0, 0.0], &unit_params(), &unit_source());
        let a2 = radiation_acceleration([2.0, 0.0], &unit_params(), &unit_source());
        let ratio = a1[0] / a2[0];
        assert!(
            (ratio - 4.0).abs() < 1e-12,
            "expected r⁻² falloff (ratio = 4), got {ratio}"
        );
    }

    #[test]
    fn rad_scales_linearly_with_luminosity() {
        let s2 = RadiationSource {
            luminosity: 2.0,
            ..unit_source()
        };
        let a1 = radiation_acceleration([1.0, 0.0], &unit_params(), &unit_source());
        let a2 = radiation_acceleration([1.0, 0.0], &unit_params(), &s2);
        assert!((a2[0] / a1[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rad_scales_linearly_with_q_pr() {
        let p2 = RadiationParams {
            q_pr: 2.0,
            ..unit_params()
        };
        let a1 = radiation_acceleration([1.0, 0.0], &unit_params(), &unit_source());
        let a2 = radiation_acceleration([1.0, 0.0], &p2, &unit_source());
        assert!((a2[0] / a1[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn rad_zero_luminosity_gives_zero() {
        let s = RadiationSource {
            luminosity: 0.0,
            ..unit_source()
        };
        let a = radiation_acceleration([1.0, 0.0], &unit_params(), &s);
        assert_eq!(a, [0.0, 0.0]);
    }

    #[test]
    fn rad_coincident_body_returns_zero() {
        let a = radiation_acceleration([0.0, 0.0], &unit_params(), &unit_source());
        assert_eq!(a, [0.0, 0.0]);
    }

    #[test]
    fn rad_diagonal_direction_is_normalised() {
        // Body at (1, 1): force should point along (1, 1) / √2.
        let a = radiation_acceleration([1.0, 1.0], &unit_params(), &unit_source());
        let mag = (a[0] * a[0] + a[1] * a[1]).sqrt();
        assert!((a[0] / mag - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
        assert!((a[1] / mag - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12);
    }

    // ── pr_drag_acceleration ──────────────────────────────────────────────────

    #[test]
    fn pr_at_rest_equals_radiation_pressure() {
        // A body at rest relative to the source experiences only direct pressure.
        let a_rad = radiation_acceleration([1.0, 0.0], &unit_params(), &unit_source());
        let a_pr = pr_drag_acceleration([1.0, 0.0], [0.0, 0.0], &unit_params(), &unit_source());
        assert!((a_rad[0] - a_pr[0]).abs() < 1e-15);
        assert!((a_rad[1] - a_pr[1]).abs() < 1e-15);
    }

    #[test]
    fn pr_tangential_drag_opposes_orbital_velocity() {
        // Body at (1, 0) moving in +y (circular orbit direction).
        // The tangential drag component should be in −y.
        let a = pr_drag_acceleration([1.0, 0.0], [0.0, 1.0], &unit_params(), &unit_source());
        assert!(
            a[1] < 0.0,
            "PR drag should decelerate the tangential velocity"
        );
    }

    #[test]
    fn pr_radial_infall_increases_radial_push() {
        // Body falling radially inward (vx < 0) should see a reduced net
        // outward push compared to a body at rest (aberration reduces flux).
        let a_rest = pr_drag_acceleration([1.0, 0.0], [0.0, 0.0], &unit_params(), &unit_source());
        let a_infall =
            pr_drag_acceleration([1.0, 0.0], [-0.1, 0.0], &unit_params(), &unit_source());
        assert!(
            a_infall[0] > a_rest[0],
            "inward radial velocity should increase apparent flux (aberration)"
        );
    }

    #[test]
    fn pr_galilean_invariance() {
        // The force depends only on the relative velocity v_rel = v_body − v_source.
        // Therefore these two scenarios must give identical accelerations:
        //   (a) body at vel = +0.5, source at rest
        //   (b) body at rest,        source at vel = −0.5
        // Both yield dvx = vel[0] − source.vx = 0.5 − 0 = 0 − (−0.5) = 0.5.
        let s_moving_neg = RadiationSource {
            vx: -0.5,
            ..unit_source()
        };
        let a_body_moving =
            pr_drag_acceleration([1.0, 0.0], [0.5, 0.0], &unit_params(), &unit_source());
        let a_source_moving_neg =
            pr_drag_acceleration([1.0, 0.0], [0.0, 0.0], &unit_params(), &s_moving_neg);
        // Both paths compute the same dvx; the results must be bit-for-bit equal.
        assert_eq!(
            a_body_moving, a_source_moving_neg,
            "PR force must depend only on v_rel: body@+v ≡ source@−v"
        );
    }

    #[test]
    fn pr_zero_luminosity_gives_zero() {
        let s = RadiationSource {
            luminosity: 0.0,
            ..unit_source()
        };
        let a = pr_drag_acceleration([1.0, 0.0], [0.5, 0.3], &unit_params(), &s);
        assert_eq!(a, [0.0, 0.0]);
    }

    #[test]
    fn pr_coincident_body_returns_zero() {
        let a = pr_drag_acceleration([0.0, 0.0], [1.0, 1.0], &unit_params(), &unit_source());
        assert_eq!(a, [0.0, 0.0]);
    }
}
