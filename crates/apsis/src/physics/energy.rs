//! Conservative-quantity observables for an N-body system.
//!
//! Pure functions that compute scalar (or `Vec3`) summaries from a body
//! slice. No simulation state is mutated; no allocator is touched on the
//! hot path beyond `per_body_potential_energy`'s explicit `Vec<f64>`.
//!
//! ## Dependency direction (load-bearing)
//!
//! This module sits in the **domain physics** layer. It MUST NOT import
//! from:
//!
//! * [`crate::core::metrics`] — the `Metrics` DTO is a downstream
//!   consumer of these scalars,
//! * [`crate::io`] — CSV columns and binary snapshots read; they do
//!   not feed back,
//! * [`crate::domain::field`] — colour-bar samplers consume per-body
//!   versions of these quantities; the mapping happens at the field
//!   layer, not here.
//!
//! A `use crate::core::...` or `use crate::io::...` line in this file
//! is a review red flag: it inverts the layering and contaminates the
//! domain with output concerns. Reverse imports of this module by the
//! layers above are correct and expected.
//!
//! ## Read-only contract
//!
//! Every public function takes `bodies: &[Body]` (immutable borrow) and
//! returns either a scalar or an owned value. Combined with the
//! [`Body: Copy`] guarantee asserted below, this means:
//!
//! * a function in this module cannot mutate body state, and
//! * a function in this module cannot influence the integrator's next
//!   step, even indirectly.
//!
//! Callers may invoke these functions freely between integration steps;
//! they must not be called from inside the integrator's hot loop.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::gravity::{G, pair_eps2};

// ── Structural contract: Body is Copy ────────────────────────────────────────
//
// The read-only contract above relies on `Body` carrying no interior
// mutability (no `Cell`, `RefCell`, `Mutex`, etc.). `Copy` is
// incompatible with all of those — `&[Body]` cannot then be used to
// mutate body state through any back-channel. Removing `Copy` from
// `Body` would break this guarantee, so we wedge it as a compile-time
// assertion: a future change that makes `Body` non-Copy fails the
// build here, surfacing the contract violation before the layering it
// protects can be undermined.
const _: () = {
    const fn assert_copy<T: Copy>() {}
    assert_copy::<Body>();
};

// ── Energy ───────────────────────────────────────────────────────────────────

/// Total kinetic energy of the system: `KE = ½ Σ mᵢ |vᵢ|²`.
///
/// `|vᵢ|² = vxᵢ² + vyᵢ² + vzᵢ²` is summed component-by-component in
/// fixed `(x, y, z)` order. Re-associating the inner sum (e.g. via
/// `Vec3::dot`) is mathematically equivalent but shifts ULPs and is
/// observable on the energy-conservation gates pinned in
/// `docs/experiments/2026-04-29-3d-port-baseline.md`.
pub fn kinetic_energy(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz)).sum()
}

/// Total mechanical energy: `E = KE + PE`.
#[inline]
pub fn total_energy(kinetic: f64, potential: f64) -> f64 {
    kinetic + potential
}

// ── Angular momentum ─────────────────────────────────────────────────────────

/// Total orbital angular momentum vector: `L = Σ mᵢ (rᵢ × vᵢ)`.
///
/// Each per-body cross product is computed component-wise in the
/// canonical right-handed convention:
///
/// ```text
/// (r × v).x = ry · vz − rz · vy
/// (r × v).y = rz · vx − rx · vz
/// (r × v).z = rx · vy − ry · vx
/// ```
///
/// For systems confined to the `z = 0` plane with `vz = 0`, the `x` and
/// `y` components are exactly zero and the `z` component matches
/// [`angular_momentum_z`] bit-for-bit.
pub fn angular_momentum(bodies: &[Body]) -> Vec3 {
    bodies.iter().fold(Vec3::ZERO, |acc, b| {
        acc + b.mass
            * Vec3::new(b.y * b.vz - b.z * b.vy, b.z * b.vx - b.x * b.vz, b.x * b.vy - b.y * b.vx)
    })
}

/// Z-component of the orbital angular momentum: `Lz = Σ mᵢ (xᵢ vyᵢ − yᵢ vxᵢ)`.
///
/// Scalar projection of [`angular_momentum`] onto `ẑ`. Numerically
/// identical to `angular_momentum(bodies).z` but computed directly so
/// the per-body reduction is a single `mass · (x·vy − y·vx)` term,
/// matching the form used by the energy-conservation gates and the
/// CSV recorder column.
pub fn angular_momentum_z(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| b.mass * (b.x * b.vy - b.y * b.vx)).sum()
}

// ── Centre of mass ───────────────────────────────────────────────────────────

/// Centre-of-mass position and velocity in the inertial frame.
///
/// Returns `(position, velocity)` as a pair of [`Vec3`]. An empty body
/// slice or zero total mass returns `(Vec3::ZERO, Vec3::ZERO)`; callers
/// that need to distinguish degenerate input from a genuine zero-COM
/// configuration should test `bodies.is_empty()` themselves.
pub fn center_of_mass_state(bodies: &[Body]) -> (Vec3, Vec3) {
    let mut m = 0.0;
    let mut pos = Vec3::ZERO;
    let mut vel = Vec3::ZERO;

    for b in bodies {
        m += b.mass;
        pos.x += b.mass * b.x;
        pos.y += b.mass * b.y;
        pos.z += b.mass * b.z;
        vel.x += b.mass * b.vx;
        vel.y += b.mass * b.vy;
        vel.z += b.mass * b.vz;
    }

    if m == 0.0 {
        return (Vec3::ZERO, Vec3::ZERO);
    }

    (pos / m, vel / m)
}

// ── Per-body potential ───────────────────────────────────────────────────────

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
/// **Cost**: O(N²). Intended for export and offline analysis (CSV
/// recorder at the per-record interval, parameter scans). NOT intended
/// for per-frame UI sampling — wiring this into a render path will
/// silently quadratic-cost the frame budget.
pub fn per_body_potential_energy(bodies: &[Body], g_factor: f64) -> Vec<f64> {
    let n = bodies.len();
    let mut pe = vec![0.0_f64; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = bodies[j].x - bodies[i].x;
            let dy = bodies[j].y - bodies[i].y;
            let dz = bodies[j].z - bodies[i].z;
            let eps2 = pair_eps2(bodies[i].softening, bodies[j].softening);
            let d2 = dx * dx + dy * dy + dz * dz + eps2;
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

    // ── kinetic_energy ────────────────────────────────────────────────────────

    #[test]
    fn kinetic_energy_single_body() {
        let b = Body::rocky(2.0).at(0.0, 0.0).with_velocity(3.0, 4.0);
        assert!((kinetic_energy(&[b]) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn kinetic_energy_at_rest_is_zero() {
        let b = Body::rocky(5.0).at(1.0, 2.0).with_velocity(0.0, 0.0);
        assert_eq!(kinetic_energy(&[b]), 0.0);
    }

    #[test]
    fn kinetic_energy_is_nonnegative() {
        let b = Body::rocky(1.0).at(0.0, 0.0).with_velocity(-3.0, 4.0);
        assert!(kinetic_energy(&[b]) >= 0.0);
    }

    #[test]
    fn kinetic_energy_is_additive() {
        let b1 = Body::rocky(1.0).at(0.0, 0.0).with_velocity(1.0, 0.0);
        let b2 = Body::rocky(2.0).at(0.0, 0.0).with_velocity(0.0, 2.0);
        assert!((kinetic_energy(&[b1, b2]) - 4.5).abs() < 1e-12);
    }

    #[test]
    fn kinetic_energy_includes_vz_component() {
        // KE = ½ · 2 · (3² + 4² + 12²) = ½ · 2 · 169 = 169
        let b = Body::rocky(2.0).at_3d(0.0, 0.0, 0.0).with_velocity_3d(3.0, 4.0, 12.0);
        assert!((kinetic_energy(&[b]) - 169.0).abs() < 1e-12);
    }

    #[test]
    fn kinetic_energy_planar_input_is_unchanged_by_3d_path() {
        // Adding `+ 0.0² + 0.0²` (vz = 0) to a 2D KE must be IEEE-754
        // exact. Value-equal, not within tolerance.
        let b = Body::rocky(2.0).at(0.0, 0.0).with_velocity(3.0, 4.0);
        let pre_3d = 0.5 * 2.0 * (3.0_f64 * 3.0 + 4.0_f64 * 4.0);
        assert_eq!(kinetic_energy(&[b]), pre_3d);
    }

    // ── angular_momentum_z ────────────────────────────────────────────────────

    #[test]
    fn angular_momentum_z_circular_orbit() {
        let (r, v, m) = (3.0, 2.0, 4.0);
        let b = Body::rocky(m).at(r, 0.0).with_velocity(0.0, v);
        assert!((angular_momentum_z(&[b]) - m * r * v).abs() < 1e-12);
    }

    #[test]
    fn angular_momentum_z_positive_for_ccw() {
        let b = Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 1.0);
        assert!(angular_momentum_z(&[b]) > 0.0);
    }

    #[test]
    fn angular_momentum_z_negative_for_cw() {
        let b = Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, -1.0);
        assert!(angular_momentum_z(&[b]) < 0.0);
    }

    #[test]
    fn angular_momentum_z_is_additive() {
        let b1 = Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 1.0);
        let b2 = Body::rocky(1.0).at(0.0, 2.0).with_velocity(-1.0, 0.0);
        assert!((angular_momentum_z(&[b1, b2]) - 3.0).abs() < 1e-12);
    }

    // ── angular_momentum (Vec3) ───────────────────────────────────────────────

    #[test]
    fn angular_momentum_planar_orbit_has_only_z_component() {
        let b = Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 1.0);
        let l = angular_momentum(&[b]);
        assert_eq!(l.x, 0.0);
        assert_eq!(l.y, 0.0);
        assert_eq!(l.z, 1.0);
    }

    #[test]
    fn angular_momentum_z_matches_vector_z_component() {
        // Bit-exact: the scalar `_z` accessor and the vector function
        // must agree on the z-component for any input.
        let b1 = Body::rocky(1.0).at(1.0, 0.5).with_velocity(0.2, 1.3);
        let b2 = Body::rocky(2.5).at(-0.4, 2.0).with_velocity(-0.7, 0.9);
        let bodies = [b1, b2];
        assert_eq!(angular_momentum(&bodies).z, angular_momentum_z(&bodies));
    }

    #[test]
    fn angular_momentum_inclined_orbit_picks_up_x_y_components() {
        // Body at (1, 0, 0) moving in +z direction: L = m · (r × v)
        // = (1, 0, 0) × (0, 0, 1) = (0·1 − 0·0, 0·0 − 1·1, 1·0 − 0·0)
        // = (0, −1, 0).
        let b = Body::rocky(1.0).at_3d(1.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 1.0);
        let l = angular_momentum(&[b]);
        assert_eq!(l, Vec3::new(0.0, -1.0, 0.0));
    }

    // ── total_energy ──────────────────────────────────────────────────────────

    #[test]
    fn total_energy_is_sum_of_ke_and_pe() {
        assert!((total_energy(3.0, -5.0) - (-2.0)).abs() < 1e-12);
        assert!((total_energy(0.0, -7.0) - (-7.0)).abs() < 1e-12);
    }

    // ── center_of_mass_state ──────────────────────────────────────────────────

    #[test]
    fn com_position_is_midpoint_for_equal_masses() {
        let b1 = Body::rocky(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
        let b2 = Body::rocky(1.0).at(4.0, 2.0).with_velocity(0.0, 0.0);
        let (pos, _) = center_of_mass_state(&[b1, b2]);
        assert!((pos.x - 2.0).abs() < 1e-12);
        assert!((pos.y - 1.0).abs() < 1e-12);
        assert_eq!(pos.z, 0.0);
    }

    #[test]
    fn com_velocity_is_mass_weighted_mean() {
        let b1 = Body::rocky(1.0).at(0.0, 0.0).with_velocity(4.0, 0.0);
        let b2 = Body::rocky(3.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
        let (_, vel) = center_of_mass_state(&[b1, b2]);
        assert!((vel.x - 1.0).abs() < 1e-12);
        assert!(vel.y.abs() < 1e-12);
        assert_eq!(vel.z, 0.0);
    }

    #[test]
    fn com_z_components_track_3d_input() {
        let b1 = Body::rocky(1.0).at_3d(0.0, 0.0, 1.0).with_velocity_3d(0.0, 0.0, 4.0);
        let b2 = Body::rocky(3.0).at_3d(0.0, 0.0, 5.0).with_velocity_3d(0.0, 0.0, 0.0);
        let (pos, vel) = center_of_mass_state(&[b1, b2]);
        // z_com = (1·1 + 3·5) / 4 = 16/4 = 4
        assert!((pos.z - 4.0).abs() < 1e-12);
        // vz_com = (1·4 + 3·0) / 4 = 1
        assert!((vel.z - 1.0).abs() < 1e-12);
    }

    #[test]
    fn com_of_empty_slice_returns_zero() {
        assert_eq!(center_of_mass_state(&[]), (Vec3::ZERO, Vec3::ZERO));
    }
}
