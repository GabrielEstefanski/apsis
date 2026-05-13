//! Shared helper functions used by integrator implementations.
//!
//! These operate on raw acceleration buffers and are independent of any
//! specific integrator or force model. Perturbation dispatch lives in
//! [`crate::physics::integrator::operator_dispatch`] under the
//! operator-trait-aware helpers.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::force_model::ForceModel;

/// Ensure `acc` has the correct length, then evaluate forces via the model.
///
/// Returns the raw (unscaled) potential energy.
pub fn evaluate(bodies: &[Body], force: &mut dyn ForceModel, acc: &mut Vec<Vec3>) -> f64 {
    if acc.len() != bodies.len() {
        acc.resize(bodies.len(), Vec3::ZERO);
    }
    force.compute(bodies, acc)
}

/// Multiply every acceleration in `acc` by `g_factor` and return the scaled PE.
///
/// No-op (branchless fast path) when `g_factor == 1.0`.
pub fn scale_acc_and_pe(acc: &mut [Vec3], g_factor: f64, raw_pe: f64) -> f64 {
    if (g_factor - 1.0).abs() > 1e-15 {
        for a in acc.iter_mut() {
            *a *= g_factor;
        }
    }
    raw_pe * g_factor
}
