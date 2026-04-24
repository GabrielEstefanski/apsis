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

    /// Declares whether this perturbation needs the underlying gravity to be
    /// the exact `1/r` potential (no Plummer softening) to produce a
    /// physically meaningful result.
    ///
    /// Default: `false`. Perturbations whose magnitude is dominated by the
    /// Newtonian-force baseline (radiation pressure, drag, external fields
    /// on fuzzy N-body clusters) return `false` and are unaffected by
    /// softening.
    ///
    /// Corrections that measure *deviations* from `1/r` — general-relativistic
    /// post-Newtonian terms, frame-dragging, J2 oblateness — should return
    /// `true`. The simulator's default material-scaled softening (ε on the
    /// order of 10⁻² AU for a solar-mass body) introduces a purely numerical
    /// apsidal precession that, for Mercury, is ~2 × 10³ larger than the
    /// 43 arcsec/century GR effect. Registering a perturbation with
    /// `requires_exact_gravity() == true` into a system where any body still
    /// carries nonzero softening emits a diagnostic via
    /// [`crate::warn_diag!`] — the symptom would otherwise be silent.
    ///
    /// Call [`System::with_exact_gravity`](crate::core::system::System::with_exact_gravity)
    /// or [`Body::unsoftened`](crate::domain::body::Body::unsoftened) to
    /// silence the warning and make the perturbation measurement honest.
    fn requires_exact_gravity(&self) -> bool {
        false
    }
}
