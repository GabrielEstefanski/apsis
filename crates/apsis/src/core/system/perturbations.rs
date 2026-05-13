//! Operator registration. See
//! [`crate::physics::integrator::operator`] for the trait split.

use std::sync::Arc;

use crate::core::log::Source;
use crate::core::system::System;
use crate::physics::gravity::kernel::{Kernel, RequirementViolation};
use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Operator,
};
use crate::physics::integrator::traits::IntegratorKind;

impl System {
    /// Register a Hamiltonian-class perturbation. Applied at every
    /// integration step; its `potential` (when [`Potential::Value`]) is
    /// summed into [`System::total_energy`]. Operators whose `potential`
    /// is `NotAvailable` contribute force but not energy, and surface
    /// as `HamiltonianForceOnly` in [`Self::conservation_report`].
    /// Kernel-precondition violations against the active kernel emit
    /// one structured `warn_diag` per invariant.
    ///
    /// [`Potential::Value`]: crate::physics::integrator::Potential::Value
    pub fn add_hamiltonian_perturbation(&mut self, p: Box<dyn HamiltonianOperator>) {
        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = p.kernel_requirements().check_against(&props);

        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        self.hamiltonian_perturbations.push(p);
    }

    /// Register a non-conservative perturbation (drag, radiation
    /// reaction). Symplectic integrators lose conservation invariants
    /// with one of these registered; a `warn_diag` fires at
    /// registration time when the active integrator is symplectic-class.
    pub fn add_non_conservative_perturbation(&mut self, p: Box<dyn NonConservativeOperator>) {
        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = p.kernel_requirements().check_against(&props);
        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        if integrator_is_symplectic(self.integrator.kind()) {
            crate::warn_diag!(
                Source::System,
                "non-conservative perturbation registered against a symplectic-class integrator; \
                 conservation invariants are no longer guaranteed — energy will drift at the \
                 dissipation rate of this operator",
                integrator = self.integrator.kind().slug(),
                hint = "switch to ias15 if exact conservation is required, or accept the drift",
            );
        }

        self.non_conservative_perturbations.push(p);
    }

    /// Register a pure observer. Called at synchronized step boundaries.
    pub fn register_observer(&mut self, o: Box<dyn Operator>) {
        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = o.kernel_requirements().check_against(&props);
        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        self.observers.push(o);
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

    /// Remove all registered perturbation operators (Hamiltonian and
    /// non-conservative) and observers. Use when the caller wants a
    /// fresh slate; for granular replacement see
    /// [`clear_hamiltonian_perturbations`](Self::clear_hamiltonian_perturbations),
    /// [`clear_non_conservative_perturbations`](Self::clear_non_conservative_perturbations),
    /// and [`clear_observers`](Self::clear_observers).
    pub fn clear_perturbations(&mut self) {
        self.hamiltonian_perturbations.clear();
        self.non_conservative_perturbations.clear();
        self.observers.clear();
    }

    /// Remove only the registered Hamiltonian-class perturbations,
    /// leaving non-conservative perturbations and observers untouched.
    /// Used by callers that want to atomically replace the Hamiltonian
    /// stack without disturbing dissipative coupling or diagnostic
    /// observers.
    pub fn clear_hamiltonian_perturbations(&mut self) {
        self.hamiltonian_perturbations.clear();
    }

    /// Remove only the registered non-conservative perturbations,
    /// leaving Hamiltonian operators and observers untouched.
    pub fn clear_non_conservative_perturbations(&mut self) {
        self.non_conservative_perturbations.clear();
    }

    /// Remove only the registered observers, leaving force-contributing
    /// operators untouched.
    pub fn clear_observers(&mut self) {
        self.observers.clear();
    }

    /// Total count of registered Hamiltonian + non-conservative
    /// perturbations (excludes observers).
    pub fn perturbation_count(&self) -> usize {
        self.hamiltonian_perturbations.len() + self.non_conservative_perturbations.len()
    }

    /// Count of registered observers.
    pub fn observer_count(&self) -> usize {
        self.observers.len()
    }

    /// Conservation property of the registered operator stack as a
    /// whole, plus the per-operator breakdown that produced it. See
    /// [`crate::physics::integrator::ConservationReport`] for the
    /// classification rules and attribution scope.
    pub fn conservation_report(&self) -> crate::physics::integrator::ConservationReport {
        crate::physics::integrator::ConservationReport::build(
            &self.bodies,
            &self.hamiltonian_perturbations,
            &self.non_conservative_perturbations,
        )
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
    /// sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::solar_units()));
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
    #[must_use]
    pub fn with_kernel(mut self, kernel: Arc<dyn Kernel>) -> Self {
        self.force_model.set_kernel(kernel);
        self
    }
}

fn integrator_is_symplectic(kind: IntegratorKind) -> bool {
    match kind {
        IntegratorKind::VelocityVerlet
        | IntegratorKind::Yoshida4
        | IntegratorKind::WisdomHolman
        | IntegratorKind::Mercurius => true,
        IntegratorKind::Ias15 => false,
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;
    use crate::domain::body::Body;
    use crate::physics::integrator::IntegratorKind;
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
        let mut sys =
            System::new(vec![Body::star(1.0)], UnitSystem::canonical()).with_exact_gravity();
        sys.add_body(Body::rocky(3e-6));
        assert!(sys.bodies().iter().all(|b| b.softening == 0.0));
    }

    /// A pure observer: contributes no force, no energy, just counts
    /// `observe` calls. Smoke-tests the dispatch contract that observers
    /// fire once per outer integration step at synchronized state.
    struct StepCounter(Arc<AtomicU64>);

    impl Operator for StepCounter {
        fn observe(&mut self, _bodies: &[Body], _t: f64, _dt: f64) {
            self.0.fetch_add(1, Ordering::Relaxed);
        }
    }

    #[test]
    fn registered_observer_fires_once_per_step() {
        const N_STEPS: u64 = 10;

        let counter = Arc::new(AtomicU64::new(0));
        let mut sys = System::new(
            vec![Body::star(1.0), Body::rocky(1e-6).at(1.0, 0.0)],
            UnitSystem::canonical(),
        )
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3);

        sys.register_observer(Box::new(StepCounter(Arc::clone(&counter))));
        assert_eq!(sys.observer_count(), 1);
        assert_eq!(sys.perturbation_count(), 0, "observers do not count as perturbations");

        for _ in 0..N_STEPS {
            sys.step();
        }

        assert_eq!(
            counter.load(Ordering::Relaxed),
            N_STEPS,
            "observer.observe() must fire exactly once per outer step",
        );
    }
}
