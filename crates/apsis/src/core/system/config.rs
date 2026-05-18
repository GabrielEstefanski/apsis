//! Getters and setters for all System configuration parameters.

use crate::core::adaptive::DtMode;
use crate::core::hooks::HookRegistry;
use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::traits::ExecutionProfile;
use crate::physics::integrator::{ForceModel, IntegratorKind, make_integrator};

/// Threshold above which selecting a [`ExecutionProfile::Precision`]
/// integrator on the current system emits a scale advisory through
/// `warn_diag!`.
///
/// The value is a soft hint — the call still proceeds. It exists
/// because IAS15's per-step wall time grows quickly with body count,
/// and by the time the user notices the interactive stutter, the
/// cascade may already be hundreds of substeps deep. The advisory
/// gives a cheaper signal.
///
/// 200 is at the low end of the regime where IAS15+direct becomes
/// uncomfortably slow for a 60 Hz render loop (direct O(N²) at N=200
/// is ~40 000 pair-evaluations per force call, ~50 µs typical, and
/// IAS15 does ~14 force calls per accepted substep).
pub(crate) const PRECISION_BODY_SOFT_WARN: usize = 200;

impl System {
    /// Immutable slice of all bodies in the simulation.
    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    /// Current effective timestep (may differ from `user_dt` in adaptive mode).
    pub fn dt(&self) -> f64 {
        self.current_dt
    }

    /// Total simulated time elapsed.
    pub fn t(&self) -> f64 {
        self.t
    }

    /// Number of integration steps completed.
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// Immutable access to the hook registry.
    pub fn hooks(&self) -> &HookRegistry {
        &self.hooks
    }

    /// Mutable access to the hook registry — register and configure hooks
    /// here.
    pub fn hooks_mut(&mut self) -> &mut HookRegistry {
        &mut self.hooks
    }

    /// Whether a hook has requested the main loop to stop.
    pub fn stop_requested(&self) -> bool {
        self.stop_requested
    }

    /// Clear the stop flag (headless runners call this after honouring it).
    pub fn clear_stop_request(&mut self) {
        self.stop_requested = false;
    }

    /// Total mass of the system.
    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    // ── Unit system ────────────────────────────────────────────────────────────

    /// The system's unit system, frozen at construction. Returns a borrow
    /// so the absence of `&mut UnitSystem` is visible at the type level.
    pub fn units(&self) -> &crate::units::UnitSystem {
        &self.units
    }

    // ── Gravitational scaling ──────────────────────────────────────────────────

    /// Override the runtime `G` multiplier (GUI slider). The unit system stays
    /// frozen; this scales `g_factor` on top of `units().g()`.
    pub fn set_g_factor(&mut self, g: f64) {
        self.g_factor = g.max(0.0);
    }

    pub fn g_factor(&self) -> f64 {
        self.g_factor
    }

    // ── Timestep ───────────────────────────────────────────────────────────────

    /// Set the user-requested timestep and reset the adaptive controller.
    pub fn set_dt(&mut self, dt: f64) {
        self.user_dt = dt;
        self.current_dt = dt;
        self.dt_ctrl.reset();
    }

    /// User-requested timestep (baseline; may differ from effective `dt()` in
    /// adaptive mode).
    pub fn user_dt(&self) -> f64 {
        self.user_dt
    }

    /// Set the timestep management policy.
    ///
    /// Switching to `DtMode::Fixed` restores `current_dt` to `user_dt` and
    /// resets the controller's slew history.  `DtMode::Adaptive` breaks
    /// symplecticity — avoid for publication-quality runs.
    pub fn set_dt_mode(&mut self, mode: DtMode) {
        self.dt_mode = mode;
        if mode == DtMode::Fixed {
            self.current_dt = self.user_dt;
            self.dt_ctrl.reset();
        }
    }

    pub fn dt_mode(&self) -> DtMode {
        self.dt_mode
    }

    // ── Integrator ─────────────────────────────────────────────────────────────

    pub fn integrator_kind(&self) -> IntegratorKind {
        self.integrator.kind()
    }

    /// Active gravitational kernel (used by `apsis::records` to capture
    /// the kernel variant in the record header).
    pub fn kernel(&self) -> std::sync::Arc<dyn crate::physics::gravity::kernel::Kernel> {
        self.force_model.kernel()
    }

    /// Read-only access to the registered Hamiltonian perturbations.
    /// Used by `apsis::records` provenance to enumerate operators.
    pub fn hamiltonian_perturbations(
        &self,
    ) -> &[Box<dyn crate::physics::integrator::HamiltonianOperator>] {
        &self.hamiltonian_perturbations
    }

    /// Read-only access to the registered non-conservative perturbations.
    pub fn non_conservative_perturbations(
        &self,
    ) -> &[Box<dyn crate::physics::integrator::NonConservativeOperator>] {
        &self.non_conservative_perturbations
    }

    /// Switch the integration algorithm. Takes effect on the next `step()`.
    ///
    /// ## Integrator–force compatibility
    ///
    /// Some integrators (IAS15) require the force model to be a
    /// deterministic function of state — bit-reproducible across
    /// Picard iterations. Barnes-Hut's position-dependent tree
    /// rebuild violates that invariant and, at large N, cascades into
    /// controller rejections and `dt` collapse.
    ///
    /// This method is the **single enforcement point** for that rule.
    /// Two diagnostic hooks fire here:
    ///
    /// * **Compatibility auto-correction** — if the new integrator
    ///   requires a deterministic force and the current force model
    ///   is not configured deterministically, the force model is
    ///   auto-reconfigured (exact threshold raised so BH is bypassed)
    ///   and a `warn_diag!` event is emitted. Downstream code
    ///   (physics thread, UI) never needs to re-check the pairing.
    ///
    /// * **Scale advisory for Precision-profile integrators** — if
    ///   the new integrator's execution profile is `Precision` and
    ///   the current body count is above
    ///   [`PRECISION_BODY_SOFT_WARN`], a second event is emitted
    ///   warning that interactive playback will not be real-time.
    ///   This is a hint, not a block — the caller may proceed.
    pub fn set_integrator(&mut self, kind: IntegratorKind) {
        let integrator = make_integrator(kind);

        if integrator.requires_deterministic_force() && !self.force_model.is_deterministic() {
            let prev_threshold = self.force_model.exact_threshold();
            // `usize::MAX` saturates to the engine's clamp ceiling,
            // which is the canonical "direct mode" threshold. See
            // `BarnesHutEngine::is_direct_mode`.
            self.force_model.set_exact_threshold(usize::MAX);
            let new_threshold = self.force_model.exact_threshold();
            let kind_label = kind.slug();
            crate::warn_diag!(
                crate::core::log::Source::System,
                "integrator requires deterministic force; switching force model to direct O(N²)",
                integrator = kind_label,
                exact_threshold_before = prev_threshold,
                exact_threshold_after = new_threshold,
                hint = "select velocity_verlet or yoshida4 for real-time playback",
            );
        }

        if integrator.execution_profile() == ExecutionProfile::Precision
            && self.bodies.len() > PRECISION_BODY_SOFT_WARN
        {
            let kind_label = kind.slug();
            let n = self.bodies.len();
            crate::warn_diag!(
                crate::core::log::Source::System,
                "Precision-profile integrator selected with many bodies; per-step wall time may spike",
                integrator = kind_label,
                n_bodies = n,
                soft_warn_threshold = PRECISION_BODY_SOFT_WARN,
                hint = "consider yoshida4 for interactive playback; IAS15 is designed for offline precision runs",
            );
        }

        // Wisdom-Holman emits a regime diagnostic at integrator selection time,
        // mirroring the `KernelRequirements` mismatch warnings that perturbation
        // forces issue at registration. The warning is observability only — the
        // integrator does not refuse to operate on a non-hierarchical system,
        // matching the discipline that `apsis` surfaces regime mismatches and
        // lets the caller decide.
        if matches!(kind, IntegratorKind::WisdomHolman) {
            let masses: Vec<f64> = self.bodies.iter().map(|b| b.mass).collect();
            let signal = crate::physics::integrator::traits::HierarchySignal::classify(&masses);
            match signal {
                crate::physics::integrator::traits::HierarchySignal::Hierarchical => {
                    // Quiet: the WH derivation operates inside its validated regime.
                },
                crate::physics::integrator::traits::HierarchySignal::Borderline => {
                    let m0 = self.bodies.first().map(|b| b.mass).unwrap_or(0.0);
                    let m_rest: f64 = self.bodies.iter().skip(1).map(|b| b.mass).sum();
                    let ratio = if m_rest > 0.0 { m0 / m_rest } else { f64::INFINITY };
                    crate::warn_diag!(
                        crate::core::log::Source::System,
                        "Wisdom-Holman selected with marginal hierarchy",
                        regime = signal.label(),
                        dominance_ratio = ratio,
                        threshold = 10.0,
                        hint = "energy drift may exceed the WH 1991 published floor; \
                                see docs/experiments/2026-05-03-wh-refactor.md",
                    );
                },
                crate::physics::integrator::traits::HierarchySignal::Violated => {
                    let m0 = self.bodies.first().map(|b| b.mass).unwrap_or(0.0);
                    let m_rest: f64 = self.bodies.iter().skip(1).map(|b| b.mass).sum();
                    let ratio = if m_rest > 0.0 { m0 / m_rest } else { f64::INFINITY };
                    crate::warn_diag!(
                        crate::core::log::Source::System,
                        "Wisdom-Holman selected on non-hierarchical configuration",
                        regime = signal.label(),
                        dominance_ratio = ratio,
                        threshold = 10.0,
                        hint = "WH derivation does not apply outside hierarchical regime; \
                                consider yoshida4 or ias15 for non-hierarchical systems",
                    );
                },
            }
        }

        self.integrator = integrator;
    }

    /// Set the IAS15 error tolerance. No-op for other integrators.
    pub fn set_ias15_epsilon(&mut self, eps: f64) {
        self.integrator.set_epsilon(eps);
    }

    /// Returns the active IAS15 epsilon, or `None` for other integrators.
    pub fn ias15_epsilon(&self) -> Option<f64> {
        self.integrator.epsilon()
    }

    /// Set the Mercurius Hill-radius multiplier (`α`). No-op for other
    /// integrators.
    pub fn set_mercurius_alpha(&mut self, alpha: f64) {
        self.integrator.set_hill_factor(alpha);
    }

    /// Returns the active Mercurius `α`, or `None` for other
    /// integrators.
    pub fn mercurius_alpha(&self) -> Option<f64> {
        self.integrator.hill_factor()
    }

    /// Set a cooperative wall-clock deadline for subsequent [`System::step`]
    /// calls. Adaptive integrators (IAS15) use this to short-circuit retry
    /// spins when the surrounding batch loop has already exhausted its
    /// budget. Passing `None` clears the deadline. Fixed-step integrators
    /// ignore it.
    pub fn set_step_deadline(&mut self, deadline: Option<std::time::Instant>) {
        self.step_deadline = deadline;
    }

    /// `true` if the system satisfies the Wisdom-Holman dominance criterion.
    pub fn is_wh_suitable(&self) -> bool {
        crate::physics::integrator::wisdom_holman::WisdomHolman::is_suitable_for(&self.bodies)
    }

    // ── Force model ────────────────────────────────────────────────────────────

    /// Replace the force model at runtime.
    ///
    /// Enables swapping in a different gravity engine (direct O(N²), GPU,
    /// post-Newtonian, …) without recreating the simulation.
    pub fn set_force_model(&mut self, model: Box<dyn ForceModel>) {
        self.force_model = model;
    }

    // ── Barnes-Hut θ ──────────────────────────────────────────────────────────

    /// Returns the current Barnes-Hut opening angle θ (or the force model's
    /// default if the active model does not use a tree).
    pub fn theta(&self) -> f64 {
        self.force_model.theta()
    }

    /// Set the opening angle θ (clamped to [0.05, 1.5] by the force model).
    ///
    /// No-op for non-tree force models.  Also syncs the adaptive controller.
    pub fn set_theta(&mut self, theta: f64) {
        self.force_model.set_theta(theta);
        self.theta_ctrl.set(self.force_model.theta());
    }

    /// Enable or disable the adaptive Barnes-Hut θ controller.
    ///
    /// Automatically disabled (no-op) when the active force model does not
    /// expose a Barnes-Hut engine.
    pub fn set_adaptive_theta(&mut self, enabled: bool) {
        self.adaptive_theta = enabled;
        if !enabled {
            self.theta_ctrl.set(self.force_model.theta());
        }
    }

    pub fn adaptive_theta_enabled(&self) -> bool {
        self.adaptive_theta
    }

    /// Returns the underlying Barnes-Hut engine if the active force model
    /// exposes one, or `None` for other backends.
    pub fn bh_engine(&self) -> Option<&BarnesHutEngine> {
        self.force_model.bh_engine()
    }

    // ── Reproducibility seed ──────────────────────────────────────────────────

    /// Current reproducibility seed.
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Set the reproducibility seed.
    ///
    /// The new seed is used by the next preset instantiation or cluster spawn.
    /// Does not retroactively affect already-spawned bodies.
    pub fn set_seed(&mut self, seed: u64) {
        self.seed = seed;
    }

    // ── Exact-evaluation threshold ────────────────────────────────────────────

    /// N threshold below which exact O(N²) pairwise gravity is used.
    ///
    /// For N ≤ threshold: direct O(N²) sum (exact, no tree overhead).
    /// For N > threshold: Barnes-Hut O(N log N) traversal.
    pub fn exact_threshold(&self) -> usize {
        self.force_model.exact_threshold()
    }

    /// Set the exact-evaluation threshold (clamped to [1, 10_000]).
    pub fn set_exact_threshold(&mut self, n: usize) {
        self.force_model.set_exact_threshold(n);
    }

    pub fn test_particle_threshold(&self) -> f64 {
        self.force_model.test_particle_threshold()
    }

    /// Clamped to `[0.0, 1.0]`. `0.0` disables suppression.
    pub fn set_test_particle_threshold(&mut self, threshold: f64) {
        self.force_model.set_test_particle_threshold(threshold);
    }

    #[must_use]
    pub fn with_test_particle_threshold(mut self, threshold: f64) -> Self {
        self.set_test_particle_threshold(threshold);
        self
    }

    // ── Close-encounter advisory ─────────────────────────────────────────────

    /// Set the close-encounter advisory threshold.
    ///
    /// `Some(t)` enables the [`EncounterFlag`](crate::physics::encounter::EncounterFlag)
    /// classification of the system-wide minimum pairwise separation;
    /// `None` (the default) disables the diagnostic. Setting a new
    /// threshold also resets the per-step transition tracker so the
    /// next descent into the `Close` band emits a warning event.
    pub fn set_close_encounter_threshold(&mut self, threshold: Option<f64>) {
        self.close_encounter_threshold = threshold;
        self.last_encounter_flag = crate::physics::encounter::EncounterFlag::Far;
    }

    /// Current close-encounter advisory threshold.
    pub fn close_encounter_threshold(&self) -> Option<f64> {
        self.close_encounter_threshold
    }

    /// Most recent [`EncounterFlag`](crate::physics::encounter::EncounterFlag)
    /// classification of the system-wide minimum separation. Always
    /// [`Far`](crate::physics::encounter::EncounterFlag::Far) when the
    /// threshold is unset.
    pub fn last_encounter_flag(&self) -> crate::physics::encounter::EncounterFlag {
        self.last_encounter_flag
    }

    // ── COM shift for TrailRecorder ───────────────────────────────────────────

    /// Takes (and clears) the accumulated COM translation since the last call.
    ///
    /// The physics-side COM recentering shifts all body positions; the
    /// [`TrailRecorder`](crate::core::trail::TrailRecorder) must apply
    /// the same shift to stored trail positions to keep them aligned.
    pub fn take_com_shift(&mut self) -> (f32, f32) {
        std::mem::replace(&mut self.pending_com_shift, (0.0, 0.0))
    }
}
