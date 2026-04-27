//! Non-gravitational perturbation force registration.

use std::sync::Arc;

use crate::core::log::Source;
use crate::core::system::System;
use crate::physics::gravity::kernel::{Kernel, RequirementViolation};
use crate::physics::integrator::PerturbationForce;

impl System {
    /// Register a non-gravitational perturbation force.
    ///
    /// Applied at every subsequent integration step, additively with other
    /// registered perturbations. Remove with [`clear_perturbations`].
    ///
    /// # Kernel-precondition check
    ///
    /// The perturbation's declared
    /// [`kernel_requirements`](PerturbationForce::kernel_requirements) are
    /// matched against the active kernel's
    /// [`KernelProperties`](crate::physics::gravity::kernel::KernelProperties)
    /// computed from the current bodies. Every invariant violation emits a
    /// structured [`warn_diag!`](crate::warn_diag) naming the invariant,
    /// the value required, and the value the kernel provides.
    ///
    /// For the canonical case — a 1PN correction declaring
    /// `required_exactness = Exact` against a Plummer kernel with any
    /// softened body — the simulator emits an Exactness-violation warning
    /// that also carries the legacy `softened_bodies` / `max_softening`
    /// fields used by downstream log consumers. The numerical apsidal
    /// precession from Plummer softening would otherwise silently dominate
    /// the physical signal by orders of magnitude.
    ///
    /// Dismiss the warning by adjusting the configuration:
    ///
    /// ```ignore
    /// sys = sys.with_exact_gravity();                  // whole system
    /// // or per-body at construction:
    /// let sun = Body::star(1.0).unsoftened();
    /// ```
    pub fn add_perturbation(&mut self, p: Box<dyn PerturbationForce>) {
        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = p.kernel_requirements().check_against(&props);

        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        self.perturbations.push(p);
    }

    fn emit_kernel_requirement_violation(&self, v: &RequirementViolation) {
        match v {
            RequirementViolation::Exactness { required, provided } => {
                let softened = self.bodies.iter().filter(|b| b.softening != 0.0).count();
                let max_softening =
                    self.bodies.iter().map(|b| b.softening.abs()).fold(0.0_f64, f64::max);
                crate::warn_diag!(
                    Source::System,
                    "perturbation requires exact 1/r gravity but bodies are softened; \
                     call System::with_exact_gravity() or Body::unsoftened() — numerical \
                     apsidal precession from Plummer softening will otherwise swamp the signal",
                    softened_bodies = softened,
                    total_bodies = self.bodies.len(),
                    max_softening = max_softening,
                    violated_invariant = "Exactness",
                    kernel_exactness = format!("{provided:?}"),
                    required_exactness = format!("{required:?}"),
                );
            },
            RequirementViolation::Continuity { required, provided } => {
                crate::warn_diag!(
                    Source::System,
                    "perturbation requires a smooth kernel but the active kernel has \
                     force discontinuities — force discontinuities violate the smooth \
                     Hamiltonian assumption underlying symplectic integration and will \
                     produce impulsive energy-error events wherever the trajectory \
                     crosses the discontinuity",
                    violated_invariant = "Continuity",
                    kernel_continuity = format!("{provided:?}"),
                    required_continuity = format!("{required:?}"),
                );
            },
        }
    }

    /// Remove all registered perturbation forces.
    pub fn clear_perturbations(&mut self) {
        self.perturbations.clear();
    }

    pub fn perturbation_count(&self) -> usize {
        self.perturbations.len()
    }

    #[cfg(test)]
    pub(crate) fn softening_scale_value(&self) -> f64 {
        self.softening_scale
    }

    /// Zero the Plummer softening on every currently-registered body, and
    /// set the system's `softening_scale` to `0` so any body added later
    /// inherits exact `1/r` gravity too.
    ///
    /// Use at construction time for fine-physics experiments — post-Newtonian
    /// corrections, J2 oblateness, tidal dissipation — where the default
    /// material-scaled softening would contribute a numerical apsidal
    /// precession that dominates the measurement.
    ///
    /// ```ignore
    /// let mut sys = System::from_template(TemplateKind::SolarSystem)
    ///     .with_exact_gravity()
    ///     .with_integrator(IntegratorKind::Ias15);
    /// sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));
    /// ```
    #[must_use]
    pub fn with_exact_gravity(mut self) -> Self {
        for b in &mut self.bodies {
            b.softening = 0.0;
        }
        self.softening_scale = 0.0;
        self
    }

    /// Swap the gravitational kernel the system dispatches through.
    ///
    /// The default kernel is
    /// [`PlummerKernel`](crate::physics::gravity::PlummerKernel).
    /// Researchers use this builder to run experiments against
    /// non-default kernels — for example, the
    /// [`TruncatedPlummerKernel`](crate::physics::gravity::kernel::TruncatedPlummerKernel)
    /// that demonstrates the `Continuity::C0` precondition violation.
    ///
    /// The kernel affects both the force evaluation and the properties
    /// reported to the precondition-check inside
    /// [`System::add_perturbation`]: a perturbation whose
    /// [`KernelRequirements`](crate::physics::gravity::kernel::KernelRequirements)
    /// the new kernel cannot satisfy will emit one structured diagnostic
    /// per violated invariant.
    #[must_use]
    pub fn with_kernel(mut self, kernel: Arc<dyn Kernel>) -> Self {
        self.force_model.set_kernel(kernel);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;
    use crate::units::UnitSystem;

    #[test]
    fn with_exact_gravity_zeroes_existing_bodies() {
        let bodies = vec![Body::star(1.0), Body::rocky(3e-6)];
        // Pre-condition: material-scaled softening is nonzero.
        assert!(bodies.iter().all(|b| b.softening > 0.0));

        let sys = System::new(bodies, UnitSystem::canonical()).with_exact_gravity();
        assert!(sys.bodies().iter().all(|b| b.softening == 0.0));
        assert_eq!(sys.softening_scale_value(), 0.0);
    }

    #[test]
    fn with_exact_gravity_persists_for_later_added_bodies() {
        // Bodies added *after* `with_exact_gravity` must also end up
        // unsoftened — otherwise the guarantee is leaky.
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::canonical()).with_exact_gravity();
        sys.add_body(Body::rocky(3e-6));
        assert!(sys.bodies().iter().all(|b| b.softening == 0.0));
    }
}
