//! Public extension point for non-gravitational accelerations.
//!
//! Implementors of [`PerturbationForce`] can be attached to `System` to inject
//! additional forces (radiation pressure, drag, custom fields) without modifying
//! the core integrator.  This is the first trait in the `gravity-sim` public
//! extension API and is intended to remain stable across releases.
//!
//! # Contract
//!
//! Implementations must **add** their contribution to `scratch_acc`, not
//! overwrite it — multiple perturbations compose by accumulation.

use crate::domain::body::Body;

pub trait PerturbationForce: Send + Sync {
    /// Accumulates non-gravitational accelerations into `scratch_acc`.
    ///
    /// `scratch_acc[i]` corresponds to `bodies[i]`. Implementations must
    /// **add** to existing values, not overwrite, so multiple perturbations
    /// compose correctly.
    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)]);

    /// Accumulates accelerations for a sub-slice of bodies starting at
    /// global index `offset` within `System::bodies`.
    ///
    /// Used by [`crate::core::system::System::apply_perturbations_planets`]
    /// during the Wisdom–Holman sub-step, where `scratch_acc` covers only
    /// `bodies[1..]` and the global index of each entry is `local_index + offset`.
    ///
    /// The default implementation ignores `offset` and delegates to
    /// [`Self::accumulate`] — correct for perturbations that derive params
    /// dynamically from the body slice rather than from a pre-indexed vec.
    fn accumulate_offset(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)], offset: usize) {
        let _ = offset;
        self.accumulate(bodies, scratch_acc);
    }
}
