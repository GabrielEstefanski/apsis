//! Operator registration. See
//! [`crate::physics::integrator::operator`] for the trait split.

use std::sync::Arc;

use crate::core::log::Source;
use crate::core::system::System;
use crate::physics::gravity::kernel::{Kernel, RequirementViolation};
use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Operator, UnitSystemMismatch,
};
use crate::physics::integrator::traits::IntegratorKind;

impl System {
    /// Register a Hamiltonian-class perturbation. Applied at every
    /// integration step; its `potential` (when [`Potential::Value`]) is
    /// summed into [`System::total_energy`]. Operators whose `potential`
    /// is `NotAvailable` contribute force but not energy, and surface
    /// as `HamiltonianForceOnly` in [`Self::conservation_report`].
    ///
    /// # Errors
    ///
    /// Returns [`UnitSystemMismatch`] when the operator's
    /// [`declared_units`](crate::physics::integrator::Operator::declared_units)
    /// disagrees with the `System`'s own [`UnitSystem`]. The caller
    /// owns the policy: propagate with `?`, log and skip, swap the
    /// operator, fall back, or `.expect(...)` for end-of-line scripts
    /// that treat the mismatch as unrecoverable.
    ///
    /// On error the operator is **not** registered and no other
    /// side-effects (kernel-precondition warnings, regime checks,
    /// `hamiltonian_perturbations.push`) fire.
    ///
    /// Kernel-precondition violations against the active kernel still
    /// emit one structured `warn_diag` per invariant on success path
    /// (non-fatal). Same for regime-of-validity bounds. Two-tier:
    /// `UnitSystemMismatch` is `Err` because integration would be
    /// silently wrong; the others are warnings because integration
    /// proceeds with the user's choice.
    ///
    /// [`Potential::Value`]: crate::physics::integrator::Potential::Value
    /// [`UnitSystem`]: crate::units::UnitSystem
    pub fn add_hamiltonian_perturbation(
        &mut self,
        p: Box<dyn HamiltonianOperator>,
    ) -> Result<(), Box<UnitSystemMismatch>> {
        self.check_units_match(p.as_ref())?;
        self.run_regime_check_on_operator(p.as_ref());

        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = p.kernel_requirements().check_against(&props);

        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        self.hamiltonian_perturbations.push(p);
        Ok(())
    }

    /// Register a non-conservative perturbation (drag, radiation
    /// reaction). Symplectic integrators lose conservation invariants
    /// with one of these registered; a `warn_diag` fires at
    /// registration time when the active integrator is symplectic-class.
    ///
    /// # Errors
    ///
    /// Returns [`UnitSystemMismatch`] on `UnitSystem` mismatch — same
    /// semantics as [`add_hamiltonian_perturbation`](Self::add_hamiltonian_perturbation).
    pub fn add_non_conservative_perturbation(
        &mut self,
        p: Box<dyn NonConservativeOperator>,
    ) -> Result<(), Box<UnitSystemMismatch>> {
        self.check_units_match(p.as_ref())?;
        self.run_regime_check_on_operator(p.as_ref());

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
        Ok(())
    }

    /// Register a pure observer. Called at synchronized step boundaries.
    ///
    /// # Errors
    ///
    /// Returns [`UnitSystemMismatch`] on `UnitSystem` mismatch.
    pub fn register_observer(
        &mut self,
        o: Box<dyn Operator>,
    ) -> Result<(), Box<UnitSystemMismatch>> {
        self.check_units_match(o.as_ref())?;
        self.run_regime_check_on_operator(o.as_ref());

        let kernel = self.force_model.kernel();
        let props = kernel.properties(&self.bodies);
        let violations = o.kernel_requirements().check_against(&props);
        for v in &violations {
            self.emit_kernel_requirement_violation(v);
        }

        self.observers.push(o);
        Ok(())
    }

    /// Run an operator's [`check_regime`](Operator::check_regime)
    /// against the current body state and emit one `warn_diag` per
    /// new (operator, bound) violation. Idempotent — already-emitted
    /// violations are filtered via `regime_warnings_emitted`.
    pub(crate) fn run_regime_check_on_operator(&mut self, op: &dyn Operator) {
        let violations = op.check_regime(&self.bodies, self.t);
        for v in violations {
            self.emit_regime_violation_once(v);
        }
    }

    /// Run regime checks on every registered operator. Used by
    /// `System::step` at the cadence-min over all operators.
    pub(crate) fn run_regime_checks_all(&mut self) {
        let mut violations: Vec<crate::physics::integrator::RegimeViolation> = Vec::new();
        for op in &self.hamiltonian_perturbations {
            violations.extend(op.check_regime(&self.bodies, self.t));
        }
        for op in &self.non_conservative_perturbations {
            violations.extend(op.check_regime(&self.bodies, self.t));
        }
        for op in &self.observers {
            violations.extend(op.check_regime(&self.bodies, self.t));
        }
        for v in violations {
            self.emit_regime_violation_once(v);
        }
    }

    /// Smallest cadence across all registered operators. The dynamic
    /// regime check fires every `cadence` outer steps.
    pub(crate) fn regime_check_cadence_min(&self) -> usize {
        let mut min = usize::MAX;
        for op in &self.hamiltonian_perturbations {
            min = min.min(op.regime_check_cadence());
        }
        for op in &self.non_conservative_perturbations {
            min = min.min(op.regime_check_cadence());
        }
        for op in &self.observers {
            min = min.min(op.regime_check_cadence());
        }
        // No registered operators = no checks at all (caller short-
        // circuits on `min == usize::MAX`).
        min
    }

    /// Emit a `warn_diag` for the violation iff its `(operator, bound)`
    /// pair has not already been reported in this `System`'s lifetime.
    /// Subsequent violations of the same pair are silently dropped.
    fn emit_regime_violation_once(&mut self, v: crate::physics::integrator::RegimeViolation) {
        let key = v.dedup_key();
        if !self.regime_warnings_emitted.insert(key) {
            return;
        }
        let severity = format!("{:?}", v.severity);
        let body_field = v.body_index.map(|i| i as i64).unwrap_or(-1);
        crate::warn_diag!(
            crate::core::log::Source::System,
            "operator regime-of-validity bound crossed; \
             integration continues but the operator's derivation no \
             longer strictly applies",
            operator = v.operator,
            bound = v.bound,
            value = v.value,
            threshold = v.threshold,
            severity = severity,
            body_index = body_field,
            message = v.message,
        );
    }

    /// Clear the warn-once dedup state for regime-of-validity
    /// diagnostics. Future violations of any `(operator, bound)`
    /// pair will fire one fresh `warn_diag` again. Useful when the
    /// caller deliberately changes scenario (loaded a new snapshot,
    /// reset bodies) and wants the bus re-armed.
    pub fn reset_regime_warnings(&mut self) {
        self.regime_warnings_emitted.clear();
    }

    /// Check that the operator's
    /// [`declared_units`](Operator::declared_units) matches the
    /// `System`'s own `UnitSystem`. Returns `Ok(())` when the operator
    /// is unit-agnostic (`declared_units` returns `None`) or units
    /// match. Returns [`UnitSystemMismatch`] otherwise.
    fn check_units_match(&self, op: &dyn Operator) -> Result<(), Box<UnitSystemMismatch>> {
        let Some(op_units) = op.declared_units() else {
            return Ok(());
        };
        if op_units == self.units {
            return Ok(());
        }
        Err(Box::new(UnitSystemMismatch {
            operator: op.name(),
            operator_units: op_units,
            system_units: self.units,
        }))
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
    /// let units = UnitSystem::solar_canonical();
    /// let mut sys = System::from_template(TemplateKind::SolarSystem, units)
    ///     .with_exact_gravity()
    ///     .with_integrator(IntegratorKind::Ias15);
    /// sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(units)));
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

        sys.register_observer(Box::new(StepCounter(Arc::clone(&counter))))
            .expect("StepCounter is unit-agnostic; registration must succeed");
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

    /// A Hamiltonian operator pinned to a fixed `UnitSystem`. Used to
    /// exercise the registration-time unit-system check.
    struct UnitBoundOp(UnitSystem);

    impl Operator for UnitBoundOp {
        fn declared_units(&self) -> Option<UnitSystem> {
            Some(self.0)
        }
    }

    impl HamiltonianOperator for UnitBoundOp {
        fn accumulate_force(&self, _bodies: &[Body], _acc: &mut [crate::math::Vec3]) {}
    }

    /// Operator and System share the same `UnitSystem` → registration
    /// returns `Ok(())` and the operator is pushed onto the stack.
    #[test]
    fn registration_succeeds_when_units_match() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        let result =
            sys.add_hamiltonian_perturbation(Box::new(UnitBoundOp(UnitSystem::solar_canonical())));
        assert!(result.is_ok(), "matching units must register: {result:?}");
        assert_eq!(sys.perturbation_count(), 1);
    }

    /// Operator carries a `UnitSystem` distinct from the System's →
    /// registration returns `Err(UnitSystemMismatch)` carrying the
    /// operator name and both unit systems. Operator is **not**
    /// pushed onto the stack.
    #[test]
    fn registration_returns_err_on_unit_system_mismatch() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        // Operator built for IAU solar (year, G≈4π²); System uses
        // canonical solar (year/2π, G=1). Same length scale, different
        // time scale — silently produces wrong dynamics if not caught.
        let err = sys
            .add_hamiltonian_perturbation(Box::new(UnitBoundOp(UnitSystem::solar())))
            .expect_err("mismatched units must produce Err");
        // err is Box<UnitSystemMismatch> — deref or as_ref to access fields.
        assert_eq!(err.operator_units, UnitSystem::solar());
        assert_eq!(err.system_units, UnitSystem::solar_canonical());
        assert_eq!(
            sys.perturbation_count(),
            0,
            "operator must not be registered when units mismatch",
        );
    }

    /// Display impl is human-readable and names both unit systems.
    /// Locks the message contract (consumers may grep for the
    /// "Unit-system mismatch" prefix).
    #[test]
    fn unit_system_mismatch_display_includes_both_units() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        let err = sys
            .add_hamiltonian_perturbation(Box::new(UnitBoundOp(UnitSystem::solar())))
            .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.starts_with("Unit-system mismatch"));
        assert!(msg.contains("AU") || msg.contains("yr"));
    }

    /// Operators that return `None` from `declared_units` are
    /// unit-agnostic — registration succeeds regardless of the
    /// System's own unit system.
    #[test]
    fn registration_succeeds_for_unit_agnostic_operator() {
        struct AgnosticOp;
        impl Operator for AgnosticOp {}
        impl HamiltonianOperator for AgnosticOp {
            fn accumulate_force(&self, _bodies: &[Body], _acc: &mut [crate::math::Vec3]) {}
        }

        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(AgnosticOp))
            .expect("unit-agnostic operator must register");
        assert_eq!(sys.perturbation_count(), 1);
    }

    // ── Regime-of-validity tests ─────────────────────────────────────────────

    /// Test fake declaring a per-body mass-ratio bound. Mirrors the
    /// 1PN check shape so the System-side dedup / cadence behaviour
    /// can be exercised without depending on apsis-1pn.
    struct MassRatioBound {
        warn: f64,
        cadence: usize,
    }

    impl Operator for MassRatioBound {
        fn name(&self) -> &'static str {
            "MassRatioBound"
        }
        fn check_regime(
            &self,
            bodies: &[Body],
            _t: f64,
        ) -> Vec<crate::physics::integrator::RegimeViolation> {
            let mut violations = Vec::new();
            if bodies.len() < 2 {
                return violations;
            }
            let m_primary = bodies[0].mass;
            for (i, b) in bodies.iter().enumerate().skip(1) {
                let ratio = b.mass / m_primary;
                if ratio >= self.warn {
                    violations.push(crate::physics::integrator::RegimeViolation {
                        operator: "MassRatioBound",
                        bound: "max_secondary_to_primary_mass_ratio",
                        value: ratio,
                        threshold: self.warn,
                        severity: crate::physics::integrator::Severity::Exceeded,
                        body_index: Some(i),
                        message: "test fake",
                    });
                }
            }
            violations
        }
        fn regime_check_cadence(&self) -> usize {
            self.cadence
        }
    }

    impl HamiltonianOperator for MassRatioBound {
        fn accumulate_force(&self, _bodies: &[Body], _acc: &mut [crate::math::Vec3]) {}
    }

    /// Capture every `Warn` event whose `operator` field equals the
    /// supplied name during the closure. Serialised on a per-test
    /// basis via the bus mutex so concurrent tests cannot
    /// cross-contaminate.
    fn capture_regime_warnings(
        target_operator: &'static str,
        body: impl FnOnce(),
    ) -> Vec<crate::core::log::Event> {
        use crate::core::log::{Event, Level, subscribe, unsubscribe};
        use std::sync::Mutex;
        // Serialise across regime tests in this module — the bus is
        // process-global.
        static LOCK: Mutex<()> = Mutex::new(());
        let _guard = LOCK.lock().unwrap_or_else(|p| p.into_inner());

        let captured: Arc<Mutex<Vec<Event>>> = Arc::new(Mutex::new(Vec::new()));
        let sink = captured.clone();
        let id = subscribe(move |event: &Event| {
            if event.level != Level::Warn {
                return;
            }
            let matches_op = event
                .fields
                .iter()
                .any(|(k, v)| *k == "operator" && v.trim_matches('"') == target_operator);
            if matches_op {
                sink.lock().unwrap().push(event.clone());
            }
        });

        body();

        let events = captured.lock().unwrap().clone();
        unsubscribe(id);
        events
    }

    /// Static check at registration: registering against an
    /// out-of-regime body state fires one warning immediately.
    #[test]
    fn regime_check_fires_at_registration_when_initial_state_violates() {
        let warnings = capture_regime_warnings("MassRatioBound", || {
            let mut sys = System::new(
                // Equal-mass binary: ratio = 1.0, way over warn = 0.01
                vec![Body::star(1.0), Body::star(1.0)],
                UnitSystem::solar_canonical(),
            )
            .with_integrator(IntegratorKind::Ias15);
            sys.add_hamiltonian_perturbation(Box::new(MassRatioBound { warn: 0.01, cadence: 100 }))
                .expect("regime test fixture");
        });

        assert_eq!(
            warnings.len(),
            1,
            "expected exactly one regime warning at registration, got {}",
            warnings.len()
        );
    }

    /// Within-regime registration is silent — no false positives.
    #[test]
    fn regime_check_silent_when_initial_state_within_regime() {
        let warnings = capture_regime_warnings("MassRatioBound", || {
            let mut sys = System::new(
                // Sun + Mercury: ratio ≈ 1.7e-7, well inside the regime
                vec![Body::star(1.0), Body::rocky(1.66e-7)],
                UnitSystem::solar_canonical(),
            )
            .with_integrator(IntegratorKind::Ias15);
            sys.add_hamiltonian_perturbation(Box::new(MassRatioBound { warn: 0.01, cadence: 100 }))
                .expect("regime test fixture");
        });

        assert!(
            warnings.is_empty(),
            "expected no regime warnings for in-regime initial state, got {}",
            warnings.len()
        );
    }

    /// Warn-once dedup: a violation that persists across cadence
    /// boundaries fires exactly once, not once per check.
    #[test]
    fn regime_check_dedups_persistent_violation_across_steps() {
        let warnings = capture_regime_warnings("MassRatioBound", || {
            let mut sys =
                System::new(vec![Body::star(1.0), Body::star(1.0)], UnitSystem::solar_canonical())
                    .with_integrator(IntegratorKind::Ias15)
                    .with_dt(1e-3);
            sys.add_hamiltonian_perturbation(Box::new(MassRatioBound {
                warn: 0.01,
                cadence: 10, // check every 10 steps
            }))
            .expect("MassRatioBound is unit-agnostic; registration must succeed");
            // 50 steps → 5 cadence triggers; without dedup we'd see 6
            // warnings (1 at registration + 5 dynamic). With dedup: 1.
            for _ in 0..50 {
                sys.step();
            }
        });

        assert_eq!(
            warnings.len(),
            1,
            "warn-once dedup should fire exactly one warning across all checks, got {}",
            warnings.len()
        );
    }
}
