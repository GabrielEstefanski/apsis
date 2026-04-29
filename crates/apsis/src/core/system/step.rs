//! Core simulation step and conservation-law tracking.

use crate::core::adaptive::{AccelerationStats, DtMode};
use crate::core::calibration;
use crate::core::hooks::{
    CollisionEvent, Command, EscapeEvent, HookContext, HookPhase, HookPhaseKind, HookRegistry,
};
use crate::core::system::System;
use crate::core::system::helpers::compute_closeness;
use crate::physics::energy::{angular_momentum_z, kinetic_energy, total_energy};
use crate::physics::integrator::IntegratorContext;
use crate::physics::integrator::{DenseSnapshot, IntegratorKind};

impl System {
    /// Advance the simulation by one time step using the configured integrator.
    ///
    /// Hook dispatch order (see [`crate::core::hooks`]):
    /// 1. `pre_step` — observe pre-integration state, queue commands.
    /// 2. Apply pre-step commands.
    /// 3. Integrator advances bodies.
    /// 4. Detect events (collisions, escapes) on the integrated state.
    /// 5. Dispatch event hooks and `post_step`, collect commands.
    /// 6. Optional `heartbeat` tick when `steps % heartbeat_interval == 0`.
    /// 7. Apply post-step / event commands in insertion order.
    pub fn step(&mut self) {
        // ── 1. pre_step hooks (observe pre-integration state) ────────────────
        let pre_cmds = if !self.hooks.is_empty() {
            let mut hooks = take_hooks(self);
            let ctx = build_hook_context(self, HookPhaseKind::PreStep);
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

        // Capture pre-step kinematics for Order-2 dense-output fallback.
        // `scratch_acc` holds a(t₀) from the previous step's end-of-step force
        // evaluation, which equals the start-of-this-step acceleration for all
        // four integrators (VV, Y4, WH each end with a force eval; IAS15 does too).
        // Skipped on the very first step when scratch_acc is empty.
        let pre_x0: Vec<crate::math::Vec3>;
        let pre_v0: Vec<crate::math::Vec3>;
        let pre_a0: Vec<crate::math::Vec3>;
        let need_order2 =
            !self.scratch_acc.is_empty() && self.integrator.kind() != IntegratorKind::Ias15;
        if need_order2 {
            pre_x0 = self.bodies.iter().map(|b| crate::math::Vec3::new(b.x, b.y, b.z)).collect();
            pre_v0 = self.bodies.iter().map(|b| crate::math::Vec3::new(b.vx, b.vy, b.vz)).collect();
            pre_a0 = self.scratch_acc.clone();
        } else {
            pre_x0 = Vec::new();
            pre_v0 = Vec::new();
            pre_a0 = Vec::new();
        }

        let mut ctx = IntegratorContext {
            force: &mut *self.force_model,
            g_factor,
            perturbations: &self.perturbations,
            deadline: self.step_deadline,
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
            })
        } else {
            None
        };

        self.last_diag = self.diagnostics.compute(&self.scratch_acc, &self.bodies, consumed_dt);

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        // `current_dt` is the value passed to the integrator on the *next*
        // call as its `dt_hint`, and is also surfaced to the UI / headless
        // CSV as the simulation's "current step size". Three regimes:
        //
        //   * Self-adaptive integrator (IAS15) — the integrator's own
        //     controller has already chosen the next step. Reading it via
        //     [`Integrator::proposed_next_dt`] keeps `current_dt` honest
        //     about the cadence the simulation is actually running at,
        //     rather than reporting `user_dt` perpetually. The hint we
        //     pass on the next call is the same value we reported, but
        //     by trait contract the integrator is free to refine it
        //     against the controller's internal state.
        //
        //   * `DtMode::Adaptive` — the *external* `DtController` (used by
        //     fixed-step schemes that want adaptive cadence) computes the
        //     next step from energy error and acceleration statistics.
        //
        //   * `DtMode::Fixed` with a non-self-adaptive integrator — pin
        //     `current_dt = user_dt` so the next step uses exactly the
        //     user's chosen cadence (the symplectic schemes need this for
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

        if self.adaptive_theta && !self.bodies.is_empty() {
            if let Some(engine) = self.force_model.bh_engine() {
                let theta = self.force_model.theta();
                let e_theta = engine.theta_error_proxy(0, &self.bodies, theta);
                let new_theta = self.theta_ctrl.update(e_theta, self.current_dt);
                self.force_model.set_theta(new_theta);
            }
        }

        if self.steps % 97 == 0 {
            if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
                // Route through the integrator's `recenter_bodies` rather
                // than the bare `apply_body_shift`: this preserves the
                // per-body compensation accumulators that IAS15 uses to
                // bound round-off error to `O(ε)` over long horizons.
                // Fixed-step integrators inherit the trait default (bare
                // subtraction), so behaviour is unchanged for them.
                self.integrator.recenter_bodies(&mut self.bodies, dx, dy);
                self.pending_com_shift.0 += -dx as f32;
                self.pending_com_shift.1 += -dy as f32;
            }
        }

        let (r_min, soft_max) = compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = soft_max;

        // ── 3. Detect events and dispatch post-step hooks ────────────────────
        if !self.hooks.is_empty() {
            let collisions = self.detect_collisions();
            let escapes = self.detect_escapes();

            let mut hooks = take_hooks(self);
            let heartbeat_interval = hooks.heartbeat_interval;
            let fire_heartbeat = heartbeat_interval > 0 && self.steps % heartbeat_interval == 0;

            let ctx = build_hook_context(self, HookPhaseKind::PostStep);
            let mut cmds = Vec::new();

            let event_ctx = HookContext { phase: HookPhase(HookPhaseKind::Event), ..ctx.clone() };
            for ev in &collisions {
                cmds.extend(hooks.dispatch_collision(ev, &event_ctx));
            }
            for ev in &escapes {
                cmds.extend(hooks.dispatch_escape(ev, &event_ctx));
            }
            cmds.extend(hooks.dispatch_post_step(&ctx));

            if fire_heartbeat {
                let hb_ctx = HookContext { phase: HookPhase(HookPhaseKind::Heartbeat), ..ctx };
                cmds.extend(hooks.dispatch_heartbeat(&hb_ctx));
            } else {
                drop(ctx);
            }
            drop(event_ctx);

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
            self.scratch_acc.resize(self.bodies.len(), crate::math::Vec3::ZERO);
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

        let denom = baseline.abs().max(1e-12);
        self.rel_energy_error = (total - baseline) / denom;
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

        self.abs_angular_momentum_error = (lz - baseline).abs();

        let denom = baseline.abs().max(1e-12);
        self.rel_angular_momentum_error = (lz - baseline) / denom;
    }

    /// Event detection stub — collision detection will arrive with the basic
    /// merge model. Returns empty until a [`CollisionHandler`]-style component
    /// is wired in.
    fn detect_collisions(&self) -> Vec<CollisionEvent> {
        Vec::new()
    }

    /// Event detection stub — escape detection will arrive with the boundary
    /// condition component.
    fn detect_escapes(&self) -> Vec<EscapeEvent> {
        Vec::new()
    }

    /// Apply a batch of hook-produced commands in order.
    ///
    /// Removals and merges are re-sorted by index (descending) so `swap_remove`
    /// on earlier indices cannot corrupt later ones. Other command kinds run
    /// in insertion order.
    pub(crate) fn apply_commands(&mut self, cmds: Vec<Command>) {
        if cmds.is_empty() {
            return;
        }

        // Split into removal-style and additive commands, preserving order
        // within each class. Removals are applied last, sorted descending, so
        // hook-side indices stay valid until we touch them.
        let mut removals: Vec<usize> = Vec::new();
        let mut additions: Vec<crate::domain::body::NamedBody> = Vec::new();
        let mut merges: Vec<(usize, usize, crate::domain::body::Body, Option<String>)> = Vec::new();

        for cmd in cmds {
            match cmd {
                Command::RemoveBody { index } => removals.push(index),
                Command::AddBody(nb) => additions.push(nb),
                Command::Merge { remove_a, remove_b, merged, merged_name } => {
                    merges.push((remove_a, remove_b, merged, merged_name));
                },
                Command::Stop => self.stop_requested = true,
            }
        }

        // Merges: queue both indices for removal and add the merged body.
        for (a, b, merged, name) in merges {
            removals.push(a);
            removals.push(b);
            additions.push(crate::domain::body::NamedBody { body: merged, name });
        }

        // Remove in descending, deduplicated order.
        removals.sort_unstable_by(|a, b| b.cmp(a));
        removals.dedup();
        for idx in removals {
            self.remove_body(idx);
        }

        if !additions.is_empty() {
            self.add_named_bodies(additions);
        }
    }

    // ── High-level run methods ────────────────────────────────────────────────

    /// Advance the simulation by `duration` time units.
    ///
    /// Steps until `self.t` has advanced by at least `duration` relative to
    /// its value at the start of the call. Uses the currently configured
    /// timestep (`self.current_dt`) — so adaptive integrators (IAS15) decide
    /// their own sub-cadence within each logical step, while fixed-step
    /// integrators (Yoshida, Verlet) take exactly `ceil(duration / dt)` steps.
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

    /// Advance the simulation until `self.t >= t_end`.
    ///
    /// No-op if `t_end <= self.t`. Uses the currently configured timestep
    /// (`self.current_dt`). Respects `stop_requested` and exits early.
    ///
    /// Returns the number of `step()` calls actually performed.
    pub fn integrate_until(&mut self, t_end: f64) -> u64 {
        let start_steps = self.steps;
        while self.t < t_end && !self.stop_requested {
            self.step();
        }
        self.steps - start_steps
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
        names: &system.names,
        t: system.t,
        dt: system.current_dt,
        steps: system.steps,
        rel_energy_error: system.rel_energy_error,
        rel_angular_momentum_error: system.rel_angular_momentum_error,
        phase: HookPhase(phase),
    }
}

fn take_hooks(system: &mut System) -> HookRegistry {
    std::mem::take(&mut system.hooks)
}

fn restore_hooks(system: &mut System, hooks: HookRegistry) {
    system.hooks = hooks;
}
