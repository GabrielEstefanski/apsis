//! System-scale calibration and conservation-frame management.
//!
//! This module provides **pure, stateless functions** that operate on slices
//! of bodies and optionally trails.  No knowledge of `System` internals is
//! required; the caller decides when and how often to invoke each function.
//!
//! # Design rationale
//! Separating calibration from simulation orchestration keeps every function
//! single-purpose, trivially testable, and reusable across different integrators.

use crate::domain::body::Body;

// ── Tuning constants ────────────────────────────────────────────────────────── //

/// Removes center-of-mass velocity.
///
/// Transforms the system into its rest frame.
/// Does not affect internal dynamics.
pub fn zero_com_velocity(bodies: &mut [Body], total_mass: f64) {
    if total_mass <= 0.0 || bodies.is_empty() {
        return;
    }

    let vx_cm = bodies.iter().map(|b| b.mass * b.vel_x).sum::<f64>() / total_mass;
    let vy_cm = bodies.iter().map(|b| b.mass * b.vel_y).sum::<f64>() / total_mass;

    if vx_cm.hypot(vy_cm) < 1e-15 {
        return;
    }

    for b in bodies.iter_mut() {
        b.vel_x -= vx_cm;
        b.vel_y -= vy_cm;
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
    let x_cm = bodies.iter().map(|b| b.mass * b.pos_x).sum::<f64>() / total_mass;
    let y_cm = bodies.iter().map(|b| b.mass * b.pos_y).sum::<f64>() / total_mass;
    if x_cm.hypot(y_cm) < 1e-14 {
        return None;
    }
    Some((x_cm, y_cm))
}

// ── Unit tests ──────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::rocky(mass).at(x, y).with_velocity(vx, vy)
    }

    // ── zero_com_velocity ────────────────────────────────────────────────── //

    #[test]
    fn com_velocity_is_zero_after_correction() {
        let mut bodies = vec![body(-1.0, 0.0, 2.0, 1.0, 1.0), body(1.0, 0.0, 4.0, 3.0, 3.0)];
        let m: f64 = bodies.iter().map(|b| b.mass).sum();
        zero_com_velocity(&mut bodies, m);

        let vcm_x: f64 = bodies.iter().map(|b| b.mass * b.vel_x).sum::<f64>() / m;
        let vcm_y: f64 = bodies.iter().map(|b| b.mass * b.vel_y).sum::<f64>() / m;
        assert!(vcm_x.abs() < 1e-12, "vx_cm must be 0 after correction");
        assert!(vcm_y.abs() < 1e-12, "vy_cm must be 0 after correction");
    }

    #[test]
    fn zero_com_velocity_preserves_relative_velocity() {
        let mut bodies = vec![body(0.0, 0.0, 1.0, 0.0, 1.0), body(1.0, 0.0, 3.0, 0.0, 1.0)];
        let m = 2.0;
        let dv_before = bodies[1].vel_x - bodies[0].vel_x;
        zero_com_velocity(&mut bodies, m);
        let dv_after = bodies[1].vel_x - bodies[0].vel_x;
        assert!(
            (dv_after - dv_before).abs() < 1e-12,
            "relative velocity must not change: only the bulk frame shifts"
        );
    }
}
