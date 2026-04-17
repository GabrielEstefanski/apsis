//! Core simulation step and conservation-law tracking.

use crate::core::adaptive::{AccelerationStats, DtMode};
use crate::core::calibration;
use crate::core::system::helpers::compute_closeness;
use crate::core::system::System;
use crate::physics::energy::{angular_momentum_z, kinetic_energy, total_energy};
use crate::physics::integrator::IntegratorContext;

impl System {
    /// Advance the simulation by one time step using the configured integrator.
    pub fn step(&mut self) {
        let dt = self.current_dt;
        let g_factor = self.g_factor;

        let mut ctx = IntegratorContext {
            force: &mut *self.force_model,
            g_factor,
            perturbations: &self.perturbations,
        };
        let result =
            self.integrator.step(&mut self.bodies, &mut ctx, dt, &mut self.scratch_acc);
        self.last_potential = result.potential_energy;

        let dt = self.current_dt;
        self.steps += 1;
        self.t += dt;

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
                self.trail_buf.translate(-dx as f32, -dy as f32);
            }
        }

        let (r_min, soft_max) = compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = soft_max;
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
}
