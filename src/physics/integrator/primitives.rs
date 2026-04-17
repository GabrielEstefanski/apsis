//! Primitive kernels shared by all symplectic integrators.
//!
//! These are the two stateless operators that every KDK-style scheme composes:
//!
//! | Operator | Effect |
//! |----------|--------|
//! | [`kick`]  | `v += a · dt` — momentum update |
//! | [`drift`] | `x += v · dt` — position update |
//!
//! Force evaluation is handled by the [`ForceModel`](super::force_model::ForceModel)
//! trait; integrators call it directly or via the thin
//! [`helpers::evaluate`](super::helpers::evaluate) wrapper.

use crate::domain::body::Body;

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
