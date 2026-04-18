//! Getters and setters for all System configuration parameters.

use crate::core::adaptive::DtMode;
use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{ForceModel, Integrator, IntegratorKind, make_integrator};
use crate::render::trail_buffer::TrailBuffer;

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

    /// Shared reference to the GPU-ready trail ring buffer.
    pub fn trail_buf(&self) -> &TrailBuffer {
        &self.trail_buf
    }

    /// Mutable reference to the trail ring buffer.
    ///
    /// Required by the trail renderer to drain dirty flags each frame.
    pub fn trail_buf_mut(&mut self) -> &mut TrailBuffer {
        &mut self.trail_buf
    }

    /// Total mass of the system.
    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    // ── Gravitational scaling ──────────────────────────────────────────────────

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

    /// Switch the integration algorithm. Takes effect on the next `step()`.
    pub fn set_integrator(&mut self, kind: IntegratorKind) {
        self.integrator = make_integrator(kind);
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

    // ── Softening ─────────────────────────────────────────────────────────────

    pub fn softening_scale(&self) -> f64 {
        self.softening_scale
    }

    /// Set the global Plummer softening scale (`ε = ε_default · scale`) and
    /// rescale all existing body softenings immediately.
    pub fn set_softening_scale(&mut self, scale: f64) {
        use crate::domain::body::default_softening;
        self.softening_scale = scale.max(0.0);
        for b in &mut self.bodies {
            b.softening = default_softening(b.mass) * self.softening_scale;
        }
    }

    // ── Trails ────────────────────────────────────────────────────────────────

    pub fn trail_every(&self) -> usize {
        self.trail_every
    }

    pub fn set_trail_every(&mut self, n: usize) {
        self.trail_every = n.max(1);
    }

    /// Record current body positions into the trail ring buffer.
    ///
    /// Call once per rendered frame, not per physics step.
    pub fn push_trail(&mut self) {
        self.trail_buf.push(&self.bodies);
    }
}
