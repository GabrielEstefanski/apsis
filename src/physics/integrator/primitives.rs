//! Primitive kernels shared by all symplectic integrators.
//!
//! These are the three operators that every KDK-style scheme composes:
//!
//! | Operator | Effect |
//! |----------|--------|
//! | [`evaluate_accelerations`] | Rebuild gravity structure and fill `scratch_acc` |
//! | [`kick`]                   | `v += a · dt` — momentum update |
//! | [`drift`]                  | `x += v · dt` — position update |
//!
//! Kept here rather than on a specific integrator so that custom composition
//! schemes (Yoshida variants, SABA, user plugins) can reuse the same primitives
//! without pulling in the whole enum.

use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;

/// Rebuild the gravity structure and fill `scratch_acc` with accelerations.
///
/// Returns the raw (unscaled) gravitational potential energy.
pub fn evaluate_accelerations(
    bodies: &[Body],
    theta: f64,
    engine: &mut BarnesHutEngine,
    scratch_acc: &mut Vec<(f64, f64)>,
) -> f64 {
    if scratch_acc.len() != bodies.len() {
        scratch_acc.resize(bodies.len(), (0.0, 0.0));
    }
    engine.build(bodies);
    engine.evaluate(bodies, theta, scratch_acc)
}

/// Apply a velocity kick: `v += a · dt`.
///
/// Pass `0.5 * dt` for a half-kick (VV), or any scaled `w · dt` for
/// Yoshida sub-steps (including negative w for the middle sub-step).
pub fn kick(bodies: &mut [Body], acc: &[(f64, f64)], dt: f64) {
    for (body, &(ax, ay)) in bodies.iter_mut().zip(acc.iter()) {
        body.vx += ax * dt;
        body.vy += ay * dt;
    }
}

/// Advance all positions using the current velocities: `x += v · dt`.
pub fn drift(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}
