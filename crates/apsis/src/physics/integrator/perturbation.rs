//! Public extension point for non-gravitational accelerations.
//!
//! Implementors of [`PerturbationForce`] can be attached to `System` to inject
//! additional forces (radiation pressure, drag, custom fields) without modifying
//! the core integrator.  This is the first trait in the `apsis` public
//! extension API and is intended to remain stable across releases.
//!
//! # Contract
//!
//! Implementations must **add** their contribution to `scratch_acc`, not
//! overwrite it — multiple perturbations compose by accumulation.
//!
//! # Kernel preconditions
//!
//! Perturbations whose derivation depends on structural invariants of the
//! gravitational kernel (exact 1/r base, smoothness, etc.) declare those
//! requirements through [`PerturbationForce::kernel_requirements`]. The
//! system checks the declared requirements against the active kernel's
//! [`KernelProperties`](crate::physics::gravity::kernel::KernelProperties)
//! at
//! [`System::add_perturbation`](crate::core::system::System::add_perturbation)
//! and emits a structured diagnostic for every invariant violation.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::gravity::kernel::{Exactness, KernelRequirements};

pub trait PerturbationForce: Send + Sync {
    /// Accumulates non-gravitational accelerations into `scratch_acc`.
    ///
    /// `scratch_acc[i]` corresponds to `bodies[i]`. Implementations must
    /// **add** to existing values, not overwrite, so multiple perturbations
    /// compose correctly.
    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [Vec3]);

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
    fn accumulate_offset(&self, bodies: &[Body], scratch_acc: &mut [Vec3], offset: usize) {
        let _ = offset;
        self.accumulate(bodies, scratch_acc);
    }

    /// Kernel invariants this perturbation's derivation relies on.
    ///
    /// Default: [`KernelRequirements::none`] — no constraints on the
    /// active kernel. Perturbations whose derivation assumes specific
    /// structural properties of K(r) should override this.
    ///
    /// A few representative cases:
    ///
    /// - A general-relativistic 1PN correction is derived from the
    ///   Newtonian Hamiltonian `H_N = p²/2m − GMm/r` and substituting a
    ///   softened potential invalidates the expansion itself; declare
    ///   [`KernelRequirements::exact_and_smooth`] (or at least
    ///   `required_exactness = Some(Exactness::Exact)`).
    /// - A perturbation whose derivation relies on smoothness of the
    ///   force (symplectic splitting schemes assume a smooth Hamiltonian)
    ///   but tolerates softening should declare only
    ///   `min_continuity = Some(Continuity::Smooth)`.
    /// - A perturbation that is physically meaningful against any
    ///   short-range force model (radiation pressure, Poynting-Robertson
    ///   drag, external fields on fuzzy N-body clusters) leaves the
    ///   default [`KernelRequirements::none`] in place.
    ///
    /// Registering a perturbation whose requirements the active kernel
    /// cannot satisfy emits one structured [`warn_diag!`](crate::warn_diag)
    /// per violated invariant, identifying the specific invariant, the
    /// value required, and the value the kernel provides. Dismiss the
    /// warning by adjusting the system configuration — for example,
    /// [`System::with_exact_gravity`](crate::core::system::System::with_exact_gravity)
    /// or per-body
    /// [`Body::unsoftened`](crate::domain::body::Body::unsoftened) for
    /// an exactness violation against the default Plummer kernel.
    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::none()
    }

    /// Legacy shorthand: whether the perturbation requires exact 1/r
    /// gravity.
    ///
    /// **Deprecated.** Declare kernel preconditions through
    /// [`kernel_requirements`](Self::kernel_requirements) instead — the
    /// richer record carries enough information to catch violations on
    /// invariants beyond Exactness (Continuity, etc.) and to emit
    /// diagnostics that name the specific invariant violated.
    ///
    /// Default implementation derives the answer from
    /// `kernel_requirements()` so that perturbations migrated to the new
    /// API continue to answer legacy callers correctly.
    #[deprecated(note = "override `kernel_requirements()` instead")]
    fn requires_exact_gravity(&self) -> bool {
        self.kernel_requirements().required_exactness == Some(Exactness::Exact)
    }
}
