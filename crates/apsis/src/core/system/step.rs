//! Core simulation step and conservation-law tracking.

use crate::core::adaptive::{AccelerationStats, DtMode};
use crate::core::calibration;
use crate::core::hooks::{Command, HookContext, HookPhase, HookPhaseKind, HookRegistry};
use crate::core::system::System;
use crate::core::system::helpers::compute_closeness;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::encounter::EncounterFlag;
use crate::physics::energy::{angular_momentum_z, kinetic_energy, total_energy};
use crate::physics::integrator::IntegratorContext;
use crate::physics::integrator::operator_dispatch::dispatch_observers;
use crate::physics::integrator::{DenseSnapshot, IntegratorKind};

impl System {
    /// Advance the simulation by one time step using the configured integrator.
    ///
    /// Hook dispatch order (see [`crate::core::hooks`]):
    /// 1. `pre_step` — observe pre-integration state, queue commands.
    /// 2. Apply pre-step commands.
    /// 3. Integrator advances bodies.
    /// 4. `post_step` — observe integrated state, queue commands.
    /// 5. Apply post-step commands in insertion order.
    pub fn step(&mut self) {
        // Prime the conservation baseline so the first hook fires with
        // `rel_*_error = Some(0.0)`, not the uninitialised `None`.
        if self.initial_energy.is_none() {
            self.refresh_energy_diagnostics();
        }

        // ── 1. pre_step hooks (observe pre-integration state) ────────────────
        let pre_cmds = if !self.hooks.is_empty() {
            let mut hooks = take_hooks(self);
            let resume_state = if hooks.any_wants_resume_state() {
                Some(self.integrator.resume_state())
            } else {
                None
            };
            let mut ctx = build_hook_context(self, HookPhaseKind::PreStep);
            ctx.resume_state = resume_state;
            let cmds = hooks.dispatch_pre_step(&ctx);
            drop(ctx);
            restore_hooks(self, hooks);
            cmds
        } else {
            Vec::new()
        };
        self.apply_commands(pre_cmds);

        // ── 2. Integrate ─────────────────────────────────────────────────────
        let dt = self.current_dt;
        let g_factor = self.g_factor;

        // Capture pre-step state for the Order-2 dense-output fallback (VV,
        // Y4, WH). Skipped when `scratch_acc` is empty (first step) or when
        // its length differs from `bodies.len()` — WH leaves it sized N−1
        // and a misaligned snapshot panics at the renderer. See PR #14 and
        // issue #16.
        let pre_x0: Vec<Vec3>;
        let pre_v0: Vec<Vec3>;
        let pre_a0: Vec<Vec3>;
        let need_order2 = !self.scratch_acc.is_empty()
            && self.integrator.kind() != IntegratorKind::Ias15
            && self.scratch_acc.len() == self.bodies.len();
        if need_order2 {
            pre_x0 = self.bodies.iter().map(|b| Vec3::new(b.pos_x, b.pos_y, b.pos_z)).collect();
            pre_v0 = self.bodies.iter().map(|b| Vec3::new(b.vel_x, b.vel_y, b.vel_z)).collect();
            pre_a0 = self.scratch_acc.clone();
        } else {
            pre_x0 = Vec::new();
            pre_v0 = Vec::new();
            pre_a0 = Vec::new();
        }

        let mut ctx = IntegratorContext {
            force: &mut *self.force_model,
            g_factor,
            hamiltonian_perturbations: &self.hamiltonian_perturbations,
            non_conservative_perturbations: &self.non_conservative_perturbations,
            observers: &mut self.observers,
        };
        let result = self.integrator.step(&mut self.bodies, &mut ctx, dt, &mut self.scratch_acc);
        self.last_potential = result.potential_energy;
        self.last_step_degraded = result.degraded;

        // Advance by the time the integrator actually consumed. For fixed-step
        // integrators (VV, Y4, WH) this always equals `dt`; for IAS15 it is the
        // adaptive sub-step size chosen by the error controller. Using
        // `consumed_dt` here is what prevents the "teleport" class of bug where
        // `System::t` drifts away from the physical body state when an adaptive
        // integrator accepts a sub-step smaller than the caller's budget
        // (see ADR-004; Rein & Spiegel 2015, §2.3).
        let consumed_dt = result.consumed_dt;
        self.steps += 1;
        self.t += consumed_dt;

        // Operator boundary observation. Bodies are at synchronised
        // post-step state and `t` is post-advance; observers see the
        // canonical step boundary.
        dispatch_observers(
            &self.bodies,
            self.t,
            consumed_dt,
            &mut self.hamiltonian_perturbations,
            &mut self.non_conservative_perturbations,
            &mut self.observers,
        );

        // Dynamic regime-of-validity check, gated by the smallest
        // cadence across registered operators. Each `(operator, bound)`
        // pair fires at most once per session via the warn-once dedup
        // in `emit_regime_violation_once`.
        let cadence = self.regime_check_cadence_min();
        if cadence != usize::MAX && cadence > 0 && self.steps.is_multiple_of(cadence as u64) {
            self.run_regime_checks_all();
        }

        // Produce the dense-output snapshot.  t0 = system.t() - snapshot.dt
        // works for both cases: IAS15 sub-steps use their own dt, Order-2 uses
        // the full system dt.
        self.last_dense_snapshot = if let Some(mut snap) = result.step_snapshot {
            // IAS15 path: snapshot already has x0, v0, a0, b filled.
            snap.t0 = self.t - snap.dt;
            Some(snap)
        } else if need_order2 {
            Some(DenseSnapshot {
                t0: self.t - consumed_dt,
                dt: consumed_dt,
                x0: pre_x0,
                v0: pre_v0,
                a0: pre_a0,
                b: Vec::new(),
                kind: self.integrator.kind(),
                wh_data: None,
            })
        } else {
            None
        };

        self.last_diag = self.diagnostics.compute(&self.scratch_acc, &self.bodies, consumed_dt);

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        // `current_dt` is the value passed to the integrator on the
        // *next* call as its `dt_hint`, and is surfaced to downstream
        // observers as the simulation's current step size. Three regimes:
        //
        //   * Self-adaptive integrator (IAS15) — adopt the controller's
        //     proposed next step via [`Integrator::proposed_next_dt`].
        //   * `DtMode::Adaptive` — external `DtController` computes the
        //     next step from energy error and acceleration statistics.
        //   * `DtMode::Fixed` with a non-self-adaptive integrator — pin
        //     `current_dt = user_dt` (symplectic schemes need this for
        //     measure preservation).
        self.current_dt = if self.integrator.controls_own_step_size() {
            self.integrator.proposed_next_dt().unwrap_or(self.user_dt)
        } else {
            match self.dt_mode {
                DtMode::Fixed => self.user_dt,
                DtMode::Adaptive => {
                    let stats = AccelerationStats::new(self.last_diag.max_acc, self.last_diag.jerk);
                    self.dt_ctrl.update(self.user_dt, self.rel_energy_error, stats)
                },
            }
        };

        if self.adaptive_theta
            && !self.bodies.is_empty()
            && let Some(engine) = self.force_model.bh_engine()
        {
            // theta_error_proxy reads body 0's position from the SoA
            // snapshot. Pack a transient buffer for this one call
            // (~40 µs at N = 10⁴; called once per step) — the
            // ForceModel's own snapshot was packed at the previous
            // compute() and may be stale relative to current bodies.
            let theta = self.force_model.theta();
            let mut probe_arrays = BodyArrays::with_capacity(self.bodies.len());
            probe_arrays.pack_from(&self.bodies);
            let e_theta = engine.theta_error_proxy(0, &probe_arrays, theta);
            let new_theta = self.theta_ctrl.update(e_theta, self.current_dt);
            self.force_model.set_theta(new_theta);
        }

        if self.steps.is_multiple_of(97)
            && let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass)
        {
            // Route through the integrator's `recenter_bodies` rather
            // than a bare per-body subtraction: this preserves the
            // per-body compensation accumulators that IAS15 uses to
            // bound round-off error to `O(ε)` over long horizons.
            // Fixed-step integrators inherit the trait default (bare
            // subtraction), so behaviour is unchanged for them.
            self.integrator.recenter_bodies(&mut self.bodies, dx, dy);
            self.pending_com_shift.0 += -dx as f32;
            self.pending_com_shift.1 += -dy as f32;
        }

        self.r_min = compute_closeness(&self.bodies);

        // Close-encounter advisory. When `close_encounter_threshold` is
        // unset the flag is always `Far` and this branch is a single
        // comparison; when set it grades `r_min` against the threshold
        // and emits a structured event the first step a Close descent
        // becomes visible. Edge-triggered on the previous flag — once
        // the system is Close the warning does not repeat until it has
        // climbed back out of the band.
        let new_flag = EncounterFlag::classify(self.r_min, self.close_encounter_threshold);
        if new_flag == EncounterFlag::Close
            && self.last_encounter_flag != EncounterFlag::Close
            && let Some(threshold) = self.close_encounter_threshold
        {
            crate::warn_diag!(
                crate::core::log::Source::System,
                "close encounter detected",
                r_min = self.r_min,
                threshold = threshold,
                step = self.steps,
                t = self.t,
                hint = "consider Mercurius integrator for hybrid close-encounter handling",
            );
        }
        self.last_encounter_flag = new_flag;

        // ── 3. Dispatch post-step hooks ──────────────────────────────────────
        if !self.hooks.is_empty() {
            let mut hooks = take_hooks(self);
            let resume_state = if hooks.any_wants_resume_state() {
                Some(self.integrator.resume_state())
            } else {
                None
            };

            let mut ctx = build_hook_context(self, HookPhaseKind::PostStep);
            ctx.resume_state = resume_state;
            let cmds = hooks.dispatch_post_step(&ctx);
            drop(ctx);

            restore_hooks(self, hooks);
            self.apply_commands(cmds);
        }
    }

    /// Recompute the energy / angular-momentum cache from the current
    /// body state without advancing time. One force evaluation; idempotent.
    /// Use before reading [`energy`](Self::energy) on a freshly-constructed
    /// system that hasn't run a step yet.
    pub fn refresh_energy_diagnostics(&mut self) {
        if self.bodies.is_empty() {
            self.last_kinetic = 0.0;
            self.last_potential = 0.0;
            return;
        }
        if self.scratch_acc.len() < self.bodies.len() {
            self.scratch_acc.resize(self.bodies.len(), Vec3::ZERO);
        }
        let raw_potential = self.force_model.compute(&self.bodies, &mut self.scratch_acc);
        self.last_potential = self.g_factor * raw_potential;
        self.update_energy_tracking();
        self.update_angular_momentum_tracking();
    }

    pub(crate) fn update_energy_tracking(&mut self) {
        let kinetic = kinetic_energy(&self.bodies);
        self.last_kinetic = kinetic;

        let total = total_energy(kinetic, self.last_potential);

        let baseline = match self.initial_energy {
            Some(v) => v,
            None => {
                self.initial_energy = Some(total);
                total
            },
        };

        let delta = total - baseline;
        self.abs_energy_error = delta;
        self.rel_energy_error = crate::core::system::regime::regime_aware_rel(delta, baseline);
    }

    pub(crate) fn update_angular_momentum_tracking(&mut self) {
        let lz = angular_momentum_z(&self.bodies);

        let baseline = match self.initial_angular_momentum {
            Some(v) => v,
            None => {
                self.initial_angular_momentum = Some(lz);
                lz
            },
        };

        let delta = lz - baseline;
        self.abs_angular_momentum_error = delta.abs();
        self.rel_angular_momentum_error =
            crate::core::system::regime::regime_aware_rel(delta, baseline);
    }

    /// Apply hook-produced commands in insertion order.
    pub(crate) fn apply_commands(&mut self, cmds: Vec<Command>) {
        for cmd in cmds {
            match cmd {
                Command::Stop => self.stop_requested = true,
            }
        }
    }

    // ── High-level run methods ────────────────────────────────────────────────

    /// Advance the simulation by `duration` time units.
    ///
    /// Steps until `self.t` has advanced by `duration` relative to its
    /// value at the start of the call, landing exactly on the target
    /// time (the final step is clipped; see
    /// [`set_exact_finish_time`](crate::core::system::System::set_exact_finish_time)
    /// for the opt-out and its symplectic-rhythm caveat). Adaptive
    /// integrators (IAS15) decide their own sub-cadence; fixed-step
    /// integrators (Yoshida, Verlet) take `ceil(duration / dt)` steps,
    /// the last one shortened to fit.
    ///
    /// Respects [`stop_requested`](crate::core::system::System::stop_requested)
    /// and exits early if set — callers who want unconditional progress
    /// should call [`clear_stop_request`](crate::core::system::System::clear_stop_request)
    /// first.
    ///
    /// Returns the number of `step()` calls actually performed.
    pub fn integrate_for(&mut self, duration: f64) -> u64 {
        let t_end = self.t + duration;
        self.integrate_until(t_end)
    }

    /// Advance the simulation until `self.t == t_end` (exact finish
    /// time, the default) or until the first step at or past `t_end`
    /// (opt-out; see
    /// [`set_exact_finish_time`](crate::core::system::System::set_exact_finish_time)).
    ///
    /// No-op if `t_end <= self.t`. Respects `stop_requested` and exits
    /// early. Returns the number of `step()` calls actually performed.
    ///
    /// Without exact finish time, fixed-time measurements sample the
    /// state up to one step past `t_end`. See ADR-015 for why that is
    /// a measurement error and not a rounding detail.
    pub fn integrate_until(&mut self, t_end: f64) -> u64 {
        let start_steps = self.steps;
        while self.t < t_end && !self.stop_requested {
            if self.exact_finish_time {
                let remaining = t_end - self.t;
                if remaining < self.current_dt {
                    if self.integrator.controls_own_step_size() {
                        self.integrator.cap_next_step(remaining);
                    } else {
                        self.current_dt = remaining;
                    }
                }
            }
            self.step();
        }
        // Collapse the clipped step's `t += dt` round-off so callers
        // measure at exactly `t_end`. The state is at `t_end` to within
        // one ULP by construction (every crossing step was clipped).
        if self.exact_finish_time && !self.stop_requested && self.steps > start_steps {
            self.t = t_end;
        }
        self.steps - start_steps
    }

    /// Toggle exact-finish-time semantics for
    /// [`integrate_until`](Self::integrate_until) /
    /// [`integrate_for`](Self::integrate_for). On by default. Disabling
    /// lets the loop run whole steps past the target — useful when a
    /// fixed-step symplectic rhythm matters more than the endpoint
    /// time, at the cost of sampling the state up to one step late.
    pub fn set_exact_finish_time(&mut self, exact: bool) {
        self.exact_finish_time = exact;
    }

    /// Current exact-finish-time setting.
    pub fn exact_finish_time(&self) -> bool {
        self.exact_finish_time
    }

    /// Fire `on_finish` on every registered hook. Idempotent. Invoked
    /// automatically by `Drop` as a safety net.
    pub fn finish(&mut self) {
        if self.finished || self.hooks.is_empty() {
            self.finished = true;
            return;
        }
        let mut hooks = take_hooks(self);
        let ctx = build_hook_context(self, HookPhaseKind::Finish);
        hooks.dispatch_finish(&ctx);
        drop(ctx);
        restore_hooks(self, hooks);
        self.finished = true;
    }
}

impl Drop for System {
    fn drop(&mut self) {
        // Safety net for callers that forget the explicit close. Hooks
        // holding open resources see this as their last chance to flush.
        if !self.finished {
            self.finish();
        }
    }
}

// ── Hook borrow helpers ───────────────────────────────────────────────────────
//
// `dispatch_*` needs `&mut HookRegistry`, but building `HookContext` needs
// `&System` (including `&self.bodies`). To avoid aliasing, we temporarily move
// the registry out of `System`, dispatch, then move it back.

fn build_hook_context(system: &System, phase: HookPhaseKind) -> HookContext<'_> {
    HookContext {
        bodies: &system.bodies,
        t: system.t,
        dt: system.current_dt,
        steps: system.steps,
        rel_energy_error: system.rel_energy_error,
        rel_angular_momentum_error: system.rel_angular_momentum_error,
        phase: HookPhase(phase),
        resume_state: None,
    }
}

fn take_hooks(system: &mut System) -> HookRegistry {
    std::mem::take(&mut system.hooks)
}

fn restore_hooks(system: &mut System, hooks: HookRegistry) {
    system.hooks = hooks;
}
