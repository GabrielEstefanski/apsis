//! Shared helper functions used by integrator implementations.
//!
//! These operate on raw acceleration buffers and are independent of any
//! specific integrator or force model.

use crate::domain::body::Body;
use crate::physics::integrator::force_model::ForceModel;
use crate::physics::integrator::perturbation::PerturbationForce;

/// Ensure `acc` has the correct length, then evaluate forces via the model.
///
/// Returns the raw (unscaled) potential energy.
pub fn evaluate(
    bodies: &[Body],
    force: &mut dyn ForceModel,
    acc: &mut Vec<(f64, f64)>,
) -> f64 {
    if acc.len() != bodies.len() {
        acc.resize(bodies.len(), (0.0, 0.0));
    }
    force.compute(bodies, acc)
}

/// Multiply every acceleration in `acc` by `g_factor` and return the scaled PE.
///
/// No-op (branchless fast path) when `g_factor == 1.0`.
pub fn scale_acc_and_pe(acc: &mut [(f64, f64)], g_factor: f64, raw_pe: f64) -> f64 {
    if (g_factor - 1.0).abs() > 1e-15 {
        for a in acc.iter_mut() {
            a.0 *= g_factor;
            a.1 *= g_factor;
        }
    }
    raw_pe * g_factor
}

/// Accumulate all registered perturbation forces into `acc`.
///
/// Perturbations are independent of `g_factor` — call this **after**
/// [`scale_acc_and_pe`].
pub fn apply_perturbations(
    bodies: &[Body],
    acc: &mut [(f64, f64)],
    perturbations: &[Box<dyn PerturbationForce>],
) {
    for p in perturbations {
        p.accumulate(bodies, acc);
    }
}

/// Variant of [`apply_perturbations`] for Wisdom–Holman sub-steps.
///
/// During WH the force tree is built from `bodies[1..]` only, so `acc`
/// has length `N − 1`.  This helper passes `bodies_planets` (= `bodies[1..]`)
/// and the correct `offset` to each perturbation.
pub fn apply_perturbations_planets(
    bodies_planets: &[Body],
    acc: &mut [(f64, f64)],
    perturbations: &[Box<dyn PerturbationForce>],
) {
    for p in perturbations {
        p.accumulate_offset(bodies_planets, acc, 1);
    }
}
