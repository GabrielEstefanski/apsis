//! System-scale calibration and conservation-frame management.
//!
//! This module provides **pure, stateless functions** that operate on slices
//! of bodies and optionally trails.  No knowledge of `System` internals is
//! required; the caller decides when and how often to invoke each function.
//!
//! # Design rationale
//! Separating calibration from simulation orchestration keeps every function
//! single-purpose, trivially testable, and reusable across different integrators.

use crate::domain::body::{Body, default_moment_inertia, density_from_mass_radius};
use std::collections::VecDeque;

// ── Tuning constants ────────────────────────────────────────────────────────── //

/// Softening length as a fraction of the mean inter-particle spacing, for a
/// unit-mean-mass body.  1 % keeps the force accurate to ≈ 0.01 % at the
/// typical separation while still preventing singularities.
pub const SOFTENING_ETA: f64 = 0.01;

/// Collision radius as a fraction of the mean inter-particle spacing.  Set to
/// 0.2 % so that `radius ≤ softening × 0.2`, ensuring bodies are already deep
/// in the Plummer-softened region before surfaces touch.
pub const RADIUS_ETA: f64 = 0.002;

// ── Spatial scale ───────────────────────────────────────────────────────────── //

/// Approximate mean inter-particle spacing from the configuration bounding box.
///
/// Formula: `max(Δx, Δy) / sqrt(N − 1)`
///
/// This is translation-invariant (COM position irrelevant) and degrades
/// gracefully for collinear or near-coincident configurations.
///
/// Returns `0.0` when `bodies` has fewer than two elements.
pub fn system_length_scale(bodies: &[Body]) -> f64 {
    let n = bodies.len();
    if n < 2 {
        return 0.0;
    }

    let mut min_x = bodies[0].x;
    let mut max_x = bodies[0].x;
    let mut min_y = bodies[0].y;
    let mut max_y = bodies[0].y;

    for b in &bodies[1..] {
        min_x = min_x.min(b.x);
        max_x = max_x.max(b.x);
        min_y = min_y.min(b.y);
        max_y = max_y.max(b.y);
    }

    let extent = (max_x - min_x).max(max_y - min_y).max(1e-10);
    extent / ((n - 1) as f64).sqrt()
}

// ── Conservation-frame management ───────────────────────────────────────────── //

/// Remove the bulk centre-of-mass velocity so the system evolves in its own
/// rest frame.
///
/// # Returns
/// `true` if any correction was applied (i.e. the COM velocity was non-trivial).
/// The caller should reset the energy baseline when `true` is returned, because
/// removing bulk kinetic energy changes the conserved total energy.
pub fn zero_com_velocity(bodies: &mut [Body], total_mass: f64) -> bool {
    if total_mass <= 0.0 || bodies.is_empty() {
        return false;
    }

    let vx_cm: f64 = bodies.iter().map(|b| b.mass * b.vx).sum::<f64>() / total_mass;
    let vy_cm: f64 = bodies.iter().map(|b| b.mass * b.vy).sum::<f64>() / total_mass;

    if vx_cm.hypot(vy_cm) < 1e-15 * total_mass {
        return false;
    }

    for b in bodies.iter_mut() {
        b.vx -= vx_cm;
        b.vy -= vy_cm;
    }
    true
}

/// Shift all body positions (and stored trail points) so the centre of mass
/// lies exactly at the origin.
///
/// This is a **pure coordinate translation**: no velocities, forces, potential
/// energies, or relative separations change.  Performed periodically to
/// prevent floating-point drift from corrupting the angular-momentum
/// measurement `L_z = Σ m·(x·vy − y·vx)`, which is computed relative to the
/// origin.
pub fn recenter_com(bodies: &mut [Body], trails: &mut [VecDeque<(f64, f64)>], total_mass: f64) {
    if total_mass <= 0.0 || bodies.is_empty() {
        return;
    }

    let x_cm: f64 = bodies.iter().map(|b| b.mass * b.x).sum::<f64>() / total_mass;
    let y_cm: f64 = bodies.iter().map(|b| b.mass * b.y).sum::<f64>() / total_mass;

    if x_cm.hypot(y_cm) < 1e-14 {
        return;
    }

    for b in bodies.iter_mut() {
        b.x -= x_cm;
        b.y -= y_cm;
    }

    for trail in trails.iter_mut() {
        for (tx, ty) in trail.iter_mut() {
            *tx -= x_cm;
            *ty -= y_cm;
        }
    }
}

// ── Per-body parameter calibration ─────────────────────────────────────────── //

/// Set each body's gravitational softening length from the current system scale.
///
/// Formula: `ε_i = SOFTENING_ETA × (m_i / m_mean)^(1/3) × l_mean`
///
/// Scaling with `(m_i / m_mean)^(1/3)` keeps each body's Plummer-sphere volume
/// proportional to its mass (equal-mass softening criterion).
///
/// No-ops when fewer than two bodies are present (single-body default is kept).
pub fn calibrate_softening(bodies: &mut [Body], total_mass: f64) {
    let n = bodies.len();
    if n < 2 {
        return;
    }

    let l_mean = system_length_scale(bodies);
    if l_mean <= 0.0 {
        return;
    }

    let m_mean = total_mass / n as f64;
    for b in bodies.iter_mut() {
        b.softening = SOFTENING_ETA * (b.mass / m_mean).cbrt() * l_mean;
    }
}

/// Set each body's physical collision radius from the current system scale.
///
/// Formula: `r_i = RADIUS_ETA × (m_i / m_mean)^(1/3) × l_mean`
///
/// Enforces `r_i ≤ softening_i × 0.5` so the Plummer-force regime is always
/// entered before surfaces touch — preventing the "slingshot at contact"
/// numerical artefact.
///
/// Updates `moment_inertia` using the physical radius to preserve physical consistency.
pub fn calibrate_radii(bodies: &mut [Body], total_mass: f64) {
    let n = bodies.len();
    if n < 2 {
        return;
    }

    let l_mean = system_length_scale(bodies);
    if l_mean <= 0.0 {
        return;
    }

    let m_mean = total_mass / n as f64;
    for b in bodies.iter_mut() {
        let r = RADIUS_ETA * (b.mass / m_mean).cbrt() * l_mean;
        b.radius = r.min(b.softening * 0.5);
        b.moment_inertia = default_moment_inertia(b.mass, b.radius);
    }
}

// ── Unit tests ──────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;
    use std::collections::VecDeque;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::new(
            x,
            y,
            vx,
            vy,
            mass,
            crate::domain::materials::Material::Rocky,
        )
    }

    // ── zero_com_velocity ────────────────────────────────────────────────── //

    #[test]
    fn com_velocity_is_zero_after_correction() {
        let mut bodies = vec![
            body(-1.0, 0.0, 2.0, 1.0, 1.0),
            body(1.0, 0.0, 4.0, 3.0, 3.0),
        ];
        let m: f64 = bodies.iter().map(|b| b.mass).sum();
        zero_com_velocity(&mut bodies, m);

        let vcm_x: f64 = bodies.iter().map(|b| b.mass * b.vx).sum::<f64>() / m;
        let vcm_y: f64 = bodies.iter().map(|b| b.mass * b.vy).sum::<f64>() / m;
        assert!(vcm_x.abs() < 1e-12, "vx_cm must be 0 after correction");
        assert!(vcm_y.abs() < 1e-12, "vy_cm must be 0 after correction");
    }

    #[test]
    fn zero_com_velocity_preserves_relative_velocity() {
        let mut bodies = vec![body(0.0, 0.0, 1.0, 0.0, 1.0), body(1.0, 0.0, 3.0, 0.0, 1.0)];
        let m = 2.0;
        let dv_before = bodies[1].vx - bodies[0].vx;
        zero_com_velocity(&mut bodies, m);
        let dv_after = bodies[1].vx - bodies[0].vx;
        assert!(
            (dv_after - dv_before).abs() < 1e-12,
            "relative velocity must not change: only the bulk frame shifts"
        );
    }

    #[test]
    fn zero_com_velocity_returns_false_when_already_at_rest() {
        let mut bodies = vec![
            body(0.0, 0.0, 1.0, 0.0, 1.0),
            body(1.0, 0.0, -1.0, 0.0, 1.0),
        ];
        let m = 2.0;
        let changed = zero_com_velocity(&mut bodies, m);
        assert!(!changed, "no correction needed when v_cm ≈ 0");
    }

    // ── recenter_com ─────────────────────────────────────────────────────── //

    #[test]
    fn com_position_is_zero_after_recentering() {
        let mut bodies = vec![body(3.0, 1.0, 0.0, 0.0, 1.0), body(7.0, 5.0, 0.0, 0.0, 1.0)];
        let m = 2.0;
        let mut trails: Vec<VecDeque<(f64, f64)>> = vec![VecDeque::new(), VecDeque::new()];
        recenter_com(&mut bodies, &mut trails, m);

        let cx: f64 = bodies.iter().map(|b| b.mass * b.x).sum::<f64>() / m;
        let cy: f64 = bodies.iter().map(|b| b.mass * b.y).sum::<f64>() / m;
        assert!(cx.abs() < 1e-12, "x_cm must be 0 after recentering");
        assert!(cy.abs() < 1e-12, "y_cm must be 0 after recentering");
    }

    #[test]
    fn recenter_com_preserves_relative_positions() {
        let mut bodies = vec![
            body(100.0, 0.0, 0.0, 0.0, 2.0),
            body(104.0, 0.0, 0.0, 0.0, 1.0),
        ];
        let m = 3.0;
        let mut trails = vec![VecDeque::new(), VecDeque::new()];
        let dx_before = bodies[1].x - bodies[0].x;
        recenter_com(&mut bodies, &mut trails, m);
        let dx_after = bodies[1].x - bodies[0].x;
        assert!(
            (dx_after - dx_before).abs() < 1e-12,
            "separation must not change: recentering is a rigid translation"
        );
    }

    #[test]
    fn recenter_com_shifts_trail_points_consistently() {
        let mut bodies = vec![
            body(10.0, 0.0, 0.0, 0.0, 1.0),
            body(20.0, 0.0, 0.0, 0.0, 1.0),
        ];
        let m = 2.0;
        let mut trails: Vec<VecDeque<(f64, f64)>> = vec![
            VecDeque::from([(10.0_f64, 0.0_f64)]),
            VecDeque::from([(20.0_f64, 0.0_f64)]),
        ];
        recenter_com(&mut bodies, &mut trails, m);

        // Trail point for body 0 should have shifted by the same -x_cm
        let trail_dx = trails[1].front().unwrap().0 - trails[0].front().unwrap().0;
        assert!(
            (trail_dx - 10.0).abs() < 1e-12,
            "trail relative positions must match body relative positions after shift"
        );
    }

    // ── system_length_scale ──────────────────────────────────────────────── //

    #[test]
    fn length_scale_equals_separation_for_two_bodies() {
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, 1.0), body(6.0, 0.0, 0.0, 0.0, 1.0)];
        // extent = 6, N-1 = 1 → l = 6/sqrt(1) = 6
        let l = system_length_scale(&bodies);
        assert!((l - 6.0).abs() < 1e-12);
    }

    #[test]
    fn length_scale_returns_zero_for_single_body() {
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, 1.0)];
        assert_eq!(system_length_scale(&bodies), 0.0);
    }

    // ── calibrate_softening ──────────────────────────────────────────────── //

    #[test]
    fn heavier_body_gets_larger_softening() {
        let mut bodies = vec![
            body(0.0, 0.0, 0.0, 0.0, 1.0),
            body(5.0, 0.0, 0.0, 0.0, 8.0), // 8× heavier → ε scales as ∛8 = 2×
        ];
        let m: f64 = bodies.iter().map(|b| b.mass).sum();
        calibrate_softening(&mut bodies, m);
        assert!(
            bodies[1].softening > bodies[0].softening,
            "softening must scale with mass: heavier body needs larger ε"
        );
    }

    #[test]
    fn softening_scales_as_cube_root_of_mass_ratio() {
        let mut bodies = vec![body(0.0, 0.0, 0.0, 0.0, 1.0), body(5.0, 0.0, 0.0, 0.0, 8.0)];
        let m: f64 = bodies.iter().map(|b| b.mass).sum();
        calibrate_softening(&mut bodies, m);
        // ε_i = η × (m_i / m_mean)^(1/3) × l  →  ε_1/ε_0 = (m_1/m_0)^(1/3)
        let ratio = bodies[1].softening / bodies[0].softening;
        let expected = (8.0_f64 / 1.0_f64).cbrt();
        assert!(
            (ratio - expected).abs() < 1e-10,
            "softening ratio must equal (m_1/m_0)^(1/3)"
        );
    }

    // ── calibrate_radii ──────────────────────────────────────────────────── //

    #[test]
    fn radius_never_exceeds_half_softening() {
        let mut bodies = vec![
            body(0.0, 0.0, 0.0, 0.0, 0.001),
            body(5.0, 0.0, 0.0, 0.0, 1.0),
            body(10.0, 0.0, 0.0, 0.0, 100.0),
        ];
        let m: f64 = bodies.iter().map(|b| b.mass).sum();
        calibrate_softening(&mut bodies, m);
        calibrate_radii(&mut bodies, m);
        for b in &bodies {
            assert!(
                b.radius <= b.softening * 0.5 + 1e-15,
                "radius must satisfy r ≤ ε/2 to stay in Plummer flat core at contact"
            );
        }
    }

    #[test]
    fn moment_of_inertia_consistent_with_radius_after_calibration() {
        let mut bodies = vec![body(0.0, 0.0, 0.0, 0.0, 2.0), body(4.0, 0.0, 0.0, 0.0, 2.0)];
        let m = 4.0;
        calibrate_softening(&mut bodies, m);
        calibrate_radii(&mut bodies, m);
        for b in &bodies {
            let expected_i = 0.4 * b.mass * b.physical_radius * b.physical_radius;
            assert!(
                (b.moment_inertia - expected_i).abs() < 1e-15,
                "I_z must equal (2/5)·m·r² after radius calibration"
            );
        }
    }
}
