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
        let props = kernel.properties();
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
        let props = kernel.properties();
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
        let props = kernel.properties();
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

    /// Aggregate [`Citation`](crate::physics::integrator::Citation)
    /// entries from every registered operator (Hamiltonian +
    /// non-conservative + observers). Operators without a citation
    /// (default `None`) are silently skipped.
    ///
    /// Order: Hamiltonian operators first, in registration order;
    /// then non-conservative; then observers — same order as
    /// dispatch. Stable so a consumer can diff two `citations()`
    /// outputs to confirm the dependency graph stayed bit-equal.
    ///
    /// Integrator and kernel citations are not yet aggregated by
    /// this method; they live on different traits (`Integrator`,
    /// `Kernel`) that do not yet expose `citation()`. Future
    /// expansion will fold them in so the full reference list comes
    /// from one call.
    pub fn citations(&self) -> Vec<crate::physics::integrator::Citation> {
        let mut out = Vec::new();
        for op in &self.hamiltonian_perturbations {
            if let Some(c) = op.citation() {
                out.push(c);
            }
        }
        for op in &self.non_conservative_perturbations {
            if let Some(c) = op.citation() {
                out.push(c);
            }
        }
        for op in &self.observers {
            if let Some(c) = op.citation() {
                out.push(c);
            }
        }
        out
    }

    /// Render the registered operator stack's citations as a
    /// human-readable provenance block. Layout is the standard one
    /// from [`crate::physics::integrator::render_provenance`] —
    /// stable, diffable, suitable for paper supplementary material
    /// or for embedding in snapshot files.
    ///
    /// ```ignore
    /// println!("{}", sys.provenance());
    /// // Provenance (1 operator):
    /// //
    /// //   apsis-1pn 0.1.0 (commit f2d8e91)
    /// //     DOI: 10.1007/BF00769986
    /// //     @article{anderson1975, ...}
    /// ```
    pub fn provenance(&self) -> String {
        crate::physics::integrator::render_provenance(&self.citations())
    }

    /// Emit a BibTeX `@software` block for the registered operator
    /// stack — one entry per unique crate, deduped by `crate_name`,
    /// in registration order (Hamiltonian, then non-conservative,
    /// then observers). Each entry pins the crate against the
    /// workspace `Cargo.lock` blake3 hash and reports its declared
    /// `kernel_requirements`.
    ///
    /// Goes directly into a paper `.bib`; the upstream physics
    /// references (Anderson 1975, Burns 1979, Tamayo 2019, ...) stay
    /// reachable through [`Self::citations`] `[i].bibtex`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::records::provenance::ProvenanceError`] when
    /// the workspace `Cargo.lock` cannot be located or read — same
    /// lookup contract as `attach_record`.
    pub fn cite(&self) -> Result<String, crate::records::provenance::ProvenanceError> {
        let lock_hash = crate::records::provenance::lock_blake3(None)?;
        let mut entries = Vec::new();
        let mut seen: std::collections::HashSet<&'static str> = std::collections::HashSet::new();
        for op in &self.hamiltonian_perturbations {
            let Some(c) = op.citation() else { continue };
            if seen.insert(c.crate_name) {
                entries.push((c, op.kernel_requirements()));
            }
        }
        for op in &self.non_conservative_perturbations {
            let Some(c) = op.citation() else { continue };
            if seen.insert(c.crate_name) {
                entries.push((c, op.kernel_requirements()));
            }
        }
        for op in &self.observers {
            let Some(c) = op.citation() else { continue };
            if seen.insert(c.crate_name) {
                entries.push((c, op.kernel_requirements()));
            }
        }
        Ok(crate::physics::integrator::render_cite_block(&entries, &lock_hash))
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
                let kernel_eps = self.force_model.kernel().epsilon_squared().sqrt();
                crate::warn_diag!(
                    Source::System,
                    "perturbation requires exact 1/r gravity but the active kernel \
                     has ε > 0 — numerical apsidal precession from the softened \
                     kernel will otherwise swamp the signal; rebuild with \
                     NewtonKernel::exact() (or NewtonKernel::new(0.0)) to restore exact 1/r²",
                    kernel_epsilon = kernel_eps,
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
        | IntegratorKind::WHFast
        | IntegratorKind::Mercurius
        | IntegratorKind::ImplicitMidpoint => true,
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

    // ── Citation aggregation tests ───────────────────────────────────────────

    /// A Hamiltonian operator that publishes a fixed citation. Used to
    /// exercise `System::citations()` / `System::provenance()` without
    /// pulling in `apsis-1pn` (the workspace's only real citation
    /// publisher) — perturbations.rs lives in the core crate.
    struct CitedOp(crate::physics::integrator::Citation);

    impl Operator for CitedOp {
        fn name(&self) -> &'static str {
            "CitedOp"
        }
        fn citation(&self) -> Option<crate::physics::integrator::Citation> {
            Some(self.0)
        }
    }

    impl HamiltonianOperator for CitedOp {
        fn accumulate_force(&self, _bodies: &[Body], _acc: &mut [crate::math::Vec3]) {}
    }

    fn fake_citation(crate_name: &'static str) -> crate::physics::integrator::Citation {
        crate::physics::integrator::Citation {
            bibtex: "@article{fake, year={2026}}",
            doi: Some("10.0000/fake"),
            crate_name,
            crate_version: "0.1.0",
            commit_hash: None,
            description: Some("fake operator for tests"),
            url: Some("https://example.invalid/fake"),
        }
    }

    /// Empty stack — no citations, provenance string says so.
    #[test]
    fn citations_empty_when_no_operators_registered() {
        let sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical());
        assert!(sys.citations().is_empty());
        assert!(sys.provenance().contains("no operators"));
    }

    /// Operators that don't override `citation()` (default `None`) are
    /// silently skipped. Mixed stacks return only the publishers.
    #[test]
    fn citations_skip_operators_without_citation() {
        struct UncitedOp;
        impl Operator for UncitedOp {}
        impl HamiltonianOperator for UncitedOp {
            fn accumulate_force(&self, _b: &[Body], _a: &mut [crate::math::Vec3]) {}
        }

        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(UncitedOp)).expect("UncitedOp is unit-agnostic");
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-fake"))))
            .expect("CitedOp is unit-agnostic");

        let cites = sys.citations();
        assert_eq!(cites.len(), 1, "only the publishing operator should appear");
        assert_eq!(cites[0].crate_name, "apsis-fake");
    }

    /// Registration order is preserved — consumers diff `provenance()`
    /// across runs to confirm the operator stack stayed bit-equal.
    #[test]
    fn citations_preserve_registration_order() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-a"))))
            .expect("CitedOp is unit-agnostic");
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-b"))))
            .expect("CitedOp is unit-agnostic");

        let cites = sys.citations();
        assert_eq!(cites.len(), 2);
        assert_eq!(cites[0].crate_name, "apsis-a");
        assert_eq!(cites[1].crate_name, "apsis-b");
    }

    /// `cite()` on an empty stack returns an empty BibTeX string —
    /// the lockfile is still read (so `ProvenanceError` propagation
    /// stays exercised in CI), but no entries are emitted.
    #[test]
    fn cite_empty_when_no_operators_registered() {
        let sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical());
        let block = sys.cite().expect("Cargo.lock must be readable from the workspace");
        assert!(block.is_empty(), "no operators → no @software entries; got {block:?}");
    }

    /// Two operators from different crates produce two `@software`
    /// entries in registration order. Locks the dedupe-set ordering
    /// the consumer relies on for stable paper.bib output.
    #[test]
    fn cite_emits_one_entry_per_unique_crate_in_registration_order() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-a"))))
            .expect("CitedOp is unit-agnostic");
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-b"))))
            .expect("CitedOp is unit-agnostic");
        let block = sys.cite().expect("Cargo.lock must be readable from the workspace");
        let a_at = block.find("@software{apsis-a_0.1.0,").expect("apsis-a entry");
        let b_at = block.find("@software{apsis-b_0.1.0,").expect("apsis-b entry");
        assert!(a_at < b_at, "registration order: apsis-a must come before apsis-b");
    }

    /// Two operators from the same crate collapse to one `@software`
    /// entry — apsis-radiation publishes both `RadiationPressure`
    /// (Hamiltonian) and `PoyntingRobertsonDrag` (non-conservative)
    /// and the paper.bib should not list it twice.
    #[test]
    fn cite_dedups_when_one_crate_publishes_multiple_operators() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-radiation"))))
            .expect("CitedOp is unit-agnostic");
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-radiation"))))
            .expect("CitedOp is unit-agnostic");
        let block = sys.cite().expect("Cargo.lock must be readable from the workspace");
        assert_eq!(
            block.matches("@software{apsis-radiation_").count(),
            1,
            "duplicate crate registration must produce exactly one entry; got {block}",
        );
    }

    /// `provenance()` runs the standard renderer over the aggregated
    /// citations. Sanity: the rendered block names every operator's
    /// crate.
    #[test]
    fn provenance_renders_every_registered_citation() {
        let mut sys = System::new(vec![Body::star(1.0)], UnitSystem::solar_canonical())
            .with_integrator(IntegratorKind::Ias15);
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-a"))))
            .expect("CitedOp is unit-agnostic");
        sys.add_hamiltonian_perturbation(Box::new(CitedOp(fake_citation("apsis-b"))))
            .expect("CitedOp is unit-agnostic");

        let block = sys.provenance();
        assert!(block.contains("(2 operators):"));
        assert!(block.contains("apsis-a 0.1.0"));
        assert!(block.contains("apsis-b 0.1.0"));
        assert!(block.contains("DOI: 10.0000/fake"));
    }

    /// Warn-once dedup: a violation that persists across cadence
    /// boundaries fires exactly once, not once per check.
    #[test]
    fn regime_check_dedups_persistent_violation_across_steps() {
        let warnings = capture_regime_warnings("MassRatioBound", || {
            let bodies = vec![
                Body::star(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
                Body::star(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
            ];
            let mut sys = System::new(bodies, UnitSystem::solar_canonical())
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
