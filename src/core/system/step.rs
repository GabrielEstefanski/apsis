//! Core simulation step and conservation-law tracking.

use crate::core::adaptive::{AccelerationStats, DtMode};
use crate::core::calibration;
use crate::core::hooks::{
    CollisionEvent, Command, EscapeEvent, HookContext, HookPhase, HookPhaseKind, HookRegistry,
};
use crate::core::system::System;
use crate::core::system::helpers::compute_closeness;
use crate::physics::energy::{angular_momentum_z, kinetic_energy, total_energy};
use crate::physics::integrator::{DenseSnapshot, IntegratorKind};
use crate::physics::integrator::IntegratorContext;

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
        let pre_x0: Vec<(f64, f64)>;
        let pre_v0: Vec<(f64, f64)>;
        let pre_a0: Vec<(f64, f64)>;
        let need_order2 = !self.scratch_acc.is_empty()
            && self.integrator.kind() != IntegratorKind::Ias15;
        if need_order2 {
            pre_x0 = self.bodies.iter().map(|b| (b.x, b.y)).collect();
            pre_v0 = self.bodies.iter().map(|b| (b.vx, b.vy)).collect();
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
        };
        let result = self.integrator.step(&mut self.bodies, &mut ctx, dt, &mut self.scratch_acc);
        self.last_potential = result.potential_energy;

        self.steps += 1;
        self.t += dt;

        // Produce the dense-output snapshot.  t0 = system.t() - snapshot.dt
        // works for both cases: IAS15 sub-steps use their own dt, Order-2 uses
        // the full system dt.
        self.last_dense_snapshot = if let Some(mut snap) = result.step_snapshot {
            // IAS15 path: snapshot already has x0, v0, a0, b filled.
            snap.t0 = self.t - snap.dt;
            Some(snap)
        } else if need_order2 {
            Some(DenseSnapshot {
                t0: self.t - dt,
                dt,
                x0: pre_x0,
                v0: pre_v0,
                a0: pre_a0,
                b: Vec::new(),
                kind: self.integrator.kind(),
            })
        } else {
            None
        };

        self.last_diag = self.diagnostics.compute(&self.scratch_acc, &self.bodies, dt);

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        self.current_dt = match self.dt_mode {
            DtMode::Fixed => self.user_dt,
            DtMode::Adaptive => {
                let stats = AccelerationStats::new(self.last_diag.max_acc, self.last_diag.jerk);
                self.dt_ctrl.update(self.user_dt, self.rel_energy_error, stats)
            },
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
                calibration::apply_body_shift(&mut self.bodies, dx, dy);
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
