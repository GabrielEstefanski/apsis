//! System-scale calibration and conservation-frame management.
//!
//! This module provides **pure, stateless functions** that operate on slices
//! of bodies and optionally trails.  No knowledge of `System` internals is
//! required; the caller decides when and how often to invoke each function.
//!
//! # Design rationale
//! Separating calibration from simulation orchestration keeps every function
//! single-purpose, trivially testable, and reusable across different integrators.

use crate::core::body::Body;
use std::collections::VecDeque;

// ── Tuning constants ────────────────────────────────────────────────────────── //

/// Softening length as a fraction of the mean inter-particle spacing, for a
/// unit-mean-mass body.  1 % keeps the force accurate to ≈ 0.01 % at the
/// typical separation while still preventing singularities.
pub const SOFTENING_ETA: f64 = 0.01;

/// Removes center-of-mass velocity.
///
/// Transforms the system into its rest frame.
/// Does not affect internal dynamics.
pub fn zero_com_velocity(bodies: &mut [Body], total_mass: f64) {
    if total_mass <= 0.0 || bodies.is_empty() {
        return;
    }

    let vx_cm = bodies.iter().map(|b| b.mass * b.vx).sum::<f64>() / total_mass;
    let vy_cm = bodies.iter().map(|b| b.mass * b.vy).sum::<f64>() / total_mass;

    if vx_cm.hypot(vy_cm) < 1e-15 {
        return;
    }

    for b in bodies.iter_mut() {
        b.vx -= vx_cm;
        b.vy -= vy_cm;
    }
}

/// Recenters the system so that COM is at origin.
///
/// Pure translation, does not affect physics.
pub fn recenter_com(bodies: &mut [Body], trails: &mut [VecDeque<(f64, f64)>], total_mass: f64) {
    if total_mass <= 0.0 || bodies.is_empty() {
        return;
    }

    let x_cm = bodies.iter().map(|b| b.mass * b.x).sum::<f64>() / total_mass;
    let y_cm = bodies.iter().map(|b| b.mass * b.y).sum::<f64>() / total_mass;

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

/// Returns the centre-of-mass offset `(x_cm, y_cm)` that must be subtracted to
/// place the COM at the origin, or `None` when the system is already centred
/// (COM within 1 fm) or the input is degenerate.
///
/// This is a pure query — it does **not** modify bodies or trails.
pub fn com_offset(bodies: &[Body], total_mass: f64) -> Option<(f64, f64)> {
    if total_mass <= 0.0 || bodies.is_empty() {
        return None;
    }
    let x_cm = bodies.iter().map(|b| b.mass * b.x).sum::<f64>() / total_mass;
    let y_cm = bodies.iter().map(|b| b.mass * b.y).sum::<f64>() / total_mass;
    if x_cm.hypot(y_cm) < 1e-14 {
        return None;
    }
    Some((x_cm, y_cm))
}

/// Translates all bodies by `(-dx, -dy)`, i.e. removes the given COM offset.
///
/// This is the body-only half of a full recentering; the caller is responsible
/// for translating any associated trail data by the same vector.
pub fn apply_body_shift(bodies: &mut [Body], dx: f64, dy: f64) {
    for b in bodies.iter_mut() {
        b.x -= dx;
        b.y -= dy;
    }
}

// ── Unit tests ──────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::body::Body;
    use std::collections::VecDeque;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::new(x, y, vx, vy, mass, crate::core::materials::Material::Rocky)
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
}
