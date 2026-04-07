//! Simulation orchestrator.

use crate::core::adaptive::{AccelerationStats, DtController, ThetaController};
use crate::core::calibration;
use crate::core::collision;
use crate::core::collision::ImpactEvent;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::metrics::Metrics;
use crate::domain::body::Body;
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{drift, evaluate_accelerations, half_kick};
use std::collections::VecDeque;

/// The central simulation state.
pub struct System {
    bodies: Vec<Body>,
    trails: Vec<VecDeque<(f64, f64)>>,
    trail_cap: usize,
    trail_every: usize,

    total_mass: f64,
    last_kinetic: f64,
    last_potential: f64,

    initial_energy: Option<f64>,
    initial_energy_scale: Option<f64>,
    rel_energy_error: f64,
    max_rel_energy_error: f64,

    initial_angular_momentum: Option<f64>,
    initial_angular_momentum_scale: Option<f64>,
    rel_angular_momentum_error: f64,
    max_rel_angular_momentum_error: f64,

    engine: BarnesHutEngine,
    scratch_acc: Vec<(f64, f64)>,

    theta_controller: ThetaController,
    dt_controller: DtController,

    diagnostics: DiagnosticsComputer,
    last_diag: SimulationDiagnostics,
    probe_interval: u64,

    steps: u64,
    current_dt: f64,
    last_proposed_dt: f64,

    theta_fixed_rel_error: f64,
    dt_fixed_rel_error: f64,
    last_theta_error_norm: f64,

    merges_this_step: usize,
    bounces_this_step: usize,
    near_miss_count: usize,
    fragments_spawned_this_step: usize,
    hit_and_runs_this_step: usize,
    /// Cumulative ejecta dust mass (too small to individually track).
    total_dust_mass: f64,

    /// Coefficient of restitution: 0.0 = perfectly inelastic (merge if bound),
    /// 1.0 = perfectly elastic.
    cor: f64,

    /// Gravitational strength multiplier.  1.0 = default (G₀ = 1).
    /// Scales all force evaluations and the orbital-energy bound test in the
    /// collision resolver.
    g_factor: f64,

    /// Collision events accumulated since the last call to `take_impact_events`.
    last_impact_events: Vec<ImpactEvent>,
}

impl System {
    pub fn new(
        trail_cap: usize,
        trail_every: usize,
        max_depth: usize,
        theta_controller: ThetaController,
        dt_controller: DtController,
        probe_interval: u64,
        bodies: Vec<Body>,
    ) -> Self {
        let trails = (0..bodies.len())
            .map(|_| VecDeque::with_capacity(trail_cap))
            .collect();
        let total_mass = bodies.iter().map(|b| b.mass).sum();

        Self {
            bodies,
            trails,
            trail_cap,
            trail_every: trail_every.max(1),
            total_mass,
            last_kinetic: 0.0,
            last_potential: 0.0,
            initial_energy: None,
            initial_energy_scale: None,
            rel_energy_error: 0.0,
            max_rel_energy_error: 0.0,
            initial_angular_momentum: None,
            initial_angular_momentum_scale: None,
            rel_angular_momentum_error: 0.0,
            max_rel_angular_momentum_error: 0.0,
            engine: BarnesHutEngine::new(max_depth),
            scratch_acc: Vec::new(),
            theta_controller,
            dt_controller,
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            probe_interval,
            steps: 0,
            current_dt: 0.0,
            last_proposed_dt: 0.0,
            theta_fixed_rel_error: 0.0,
            dt_fixed_rel_error: 0.0,
            last_theta_error_norm: 0.0,
            merges_this_step: 0,
            bounces_this_step: 0,
            near_miss_count: 0,
            fragments_spawned_this_step: 0,
            hit_and_runs_this_step: 0,
            total_dust_mass: 0.0,
            cor: 0.0,
            g_factor: 1.0,
            last_impact_events: Vec::new(),
        }
    }
}

impl System {
    fn refresh_collision_geometry(&mut self) {
        self.total_mass = self.bodies.iter().map(|b| b.mass).sum();
        calibration::calibrate_softening_and_radii(&mut self.bodies, self.total_mass);
    }

    pub fn step_adaptive(&mut self, proposed_dt: f64) -> f64 {
        self.last_proposed_dt = proposed_dt;

        let stats = AccelerationStats::new(self.last_diag.max_acc, self.last_diag.jerk);

        let safe_proposed = if self.rel_energy_error.abs()
            > self.dt_controller.config.target_rel_energy_error * 2.0
        {
            proposed_dt * 0.5
        } else {
            proposed_dt
        };

        let dt = self
            .dt_controller
            .update(safe_proposed, self.rel_energy_error, stats);

        self.step(dt);

        if self.rel_energy_error.abs() > self.dt_controller.config.target_rel_energy_error * 10.0 {
            self.last_proposed_dt = dt * 0.5;
        }

        dt
    }
    pub fn step(&mut self, dt: f64) {
        self.last_proposed_dt = dt;
        self.current_dt = dt;

        if self.scratch_acc.len() != self.bodies.len() {
            self.scratch_acc.resize(self.bodies.len(), (0.0, 0.0));
        }

        let theta = self.theta_controller.current();

        let dt_changed = self.current_dt > 0.0 && (dt / self.current_dt - 1.0).abs() > 0.5;
        let acc_stale = self.steps == 0 || dt_changed;

        if acc_stale {
            let raw_pe = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe);
        }

        half_kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);

        let collision_outcome = self.advance_with_ccd(dt, theta);
        self.merges_this_step = collision_outcome.merges;
        self.bounces_this_step = collision_outcome.bounces;
        self.near_miss_count = collision_outcome.near_misses;
        self.fragments_spawned_this_step = collision_outcome.fragments_spawned;
        self.hit_and_runs_this_step = collision_outcome.hit_and_runs;
        self.total_dust_mass += collision_outcome.total_dust_mass;
        self.last_impact_events
            .extend(collision_outcome.impact_events);

        let needs_final_acc = collision_outcome.merges == 0
            && collision_outcome.fragments_spawned == 0
            && collision_outcome.hit_and_runs == 0;

        if needs_final_acc {
            let raw_pe = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe);
        }

        half_kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);

        self.steps += 1;
        self.last_diag = self
            .diagnostics
            .compute(&self.scratch_acc, &self.bodies, dt);

        if self.steps % self.trail_every as u64 == 0 {
            for (i, b) in self.bodies.iter().enumerate() {
                let t = &mut self.trails[i];
                t.push_back((b.x, b.y));
                if t.len() > self.trail_cap {
                    t.pop_front();
                }
            }
        }

        if collision_outcome.merges > 0
            || collision_outcome.fragments_spawned > 0
            || collision_outcome.hit_and_runs > 0
        {
            self.reset_energy_baseline();
        }
        self.update_energy_tracking();

        if self.steps % self.probe_interval == 0 {
            self.update_separated_errors();
        }
        self.theta_controller
            .update(self.theta_fixed_rel_error, self.current_dt);

        if self.steps % 97 == 0 {
            calibration::recenter_com(&mut self.bodies, &mut self.trails, self.total_mass);
        }
    }

    fn advance_with_ccd(&mut self, dt: f64, theta: f64) -> collision::CollisionOutcome {
        let mut outcome = collision::CollisionOutcome::default();
        let mut remaining = dt.max(0.0);
        let max_iterations = 32;
        let mut iterations = 0;

        while remaining > 1e-8 && iterations < max_iterations {
            iterations += 1;

            let Some(event) =
                collision::find_earliest_contact(&self.bodies, &self.scratch_acc, remaining)
            else {
                drift(&mut self.bodies, remaining);
                return outcome;
            };

            if event.time > 0.0 {
                drift(&mut self.bodies, event.time);
                remaining -= event.time;
            }

            let event_outcome = collision::resolve_contact(
                &mut self.bodies,
                &mut self.trails,
                event.i,
                event.j,
                self.cor,
                self.g_factor,
            );

            outcome.merges += event_outcome.merges;
            outcome.bounces += event_outcome.bounces;
            outcome.near_misses += event_outcome.near_misses;
            outcome.fragments_spawned += event_outcome.fragments_spawned;
            outcome.hit_and_runs += event_outcome.hit_and_runs;
            outcome.total_dust_mass += event_outcome.total_dust_mass;
            outcome.impact_events.extend(event_outcome.impact_events);

            let topology_changed = event_outcome.merges > 0
                || event_outcome.fragments_spawned > 0
                || event_outcome.hit_and_runs > 0;

            let physics_changed = topology_changed || event_outcome.bounces > 0;

            if !physics_changed {
                let epsilon = remaining.min(1e-8);
                if epsilon <= 0.0 {
                    break;
                }
                drift(&mut self.bodies, epsilon);
                remaining -= epsilon;
                continue;
            }

            if topology_changed {
                self.refresh_collision_geometry();
                self.scratch_acc.resize(self.bodies.len(), (0.0, 0.0));
                self.update_separated_errors();
            }

            let raw_pe_mid = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe_mid);
        }

        if remaining > 1e-8 {
            drift(&mut self.bodies, remaining);
        }

        outcome
    }
}

impl System {
    /// Multiply every acceleration in `scratch_acc` and the raw potential
    /// by `g_factor`, then return the scaled potential.
    ///
    /// The engine always uses the hard-coded `G₀ = 1.0`; multiplying the
    /// output is equivalent to running with `G_eff = G₀ · g_factor`.
    fn scale_acc_and_pe(&mut self, raw_pe: f64) -> f64 {
        if (self.g_factor - 1.0).abs() > 1e-15 {
            for a in &mut self.scratch_acc {
                a.0 *= self.g_factor;
                a.1 *= self.g_factor;
            }
        }
        raw_pe * self.g_factor
    }
}

impl System {
    pub fn reset_energy_baseline(&mut self) {
        self.initial_energy = None;
        self.initial_energy_scale = None;
        self.rel_energy_error = 0.0;
        self.max_rel_energy_error = 0.0;
        self.initial_angular_momentum = None;
        self.initial_angular_momentum_scale = None;
        self.rel_angular_momentum_error = 0.0;
        self.max_rel_angular_momentum_error = 0.0;
    }

    fn update_energy_tracking(&mut self) {
        let kinetic = kinetic_energy(&self.bodies);
        self.last_kinetic = kinetic;
        let total = total_energy(kinetic, self.last_potential);

        let baseline = match self.initial_energy {
            Some(v) => v,
            None => {
                let scale = (kinetic.abs() + self.last_potential.abs()).max(1e-12);
                self.initial_energy = Some(total);
                self.initial_energy_scale = Some(scale);
                total
            }
        };

        let denom = self
            .initial_energy_scale
            .unwrap_or_else(|| baseline.abs().max(1e-12));
        self.rel_energy_error = (total - baseline) / denom;

        if self.rel_energy_error.abs() > self.max_rel_energy_error {
            self.max_rel_energy_error = self.rel_energy_error.abs();
        }

        // Angular momentum drift
        let lz = angular_momentum_z(&self.bodies);
        let lz_baseline = match self.initial_angular_momentum {
            Some(v) => v,
            None => {
                // Scale: use |Lz| if non-trivial, otherwise fall back to a
                // kinematic estimate Σ mᵢ |rᵢ||vᵢ| so near-zero-Lz systems
                // still get a meaningful relative error.
                let kinematic_scale: f64 = self
                    .bodies
                    .iter()
                    .map(|b| {
                        let r = (b.x * b.x + b.y * b.y).sqrt();
                        let v = (b.vx * b.vx + b.vy * b.vy).sqrt();
                        b.mass * r * v
                    })
                    .sum::<f64>()
                    .max(1e-12);
                let scale = lz.abs().max(kinematic_scale * 1e-3).max(1e-12);
                self.initial_angular_momentum = Some(lz);
                self.initial_angular_momentum_scale = Some(scale);
                lz
            }
        };

        let lz_denom = self
            .initial_angular_momentum_scale
            .unwrap_or_else(|| lz_baseline.abs().max(1e-12));
        self.rel_angular_momentum_error = (lz - lz_baseline) / lz_denom;

        if self.rel_angular_momentum_error.abs() > self.max_rel_angular_momentum_error {
            self.max_rel_angular_momentum_error = self.rel_angular_momentum_error.abs();
        }
    }

    fn update_separated_errors(&mut self) {
        let n = self.bodies.len();
        if n == 0 {
            return;
        }

        let base_energy = total_energy(self.last_kinetic, self.last_potential);
        let denom = self
            .initial_energy_scale
            .unwrap_or_else(|| base_energy.abs().max(1e-9));
        let theta = self.theta_controller.current();

        let k = ((n as f64).sqrt().ceil() as usize).min(n);
        let step_size = (n / k).max(1);
        let mut idx = (self.steps as usize) % n;
        let mut sum = 0.0_f64;

        for _ in 0..k {
            let e = self.engine.theta_error_proxy(idx, &self.bodies, theta);
            sum += e * e;
            idx = (idx + step_size) % n;
        }

        let raw = (sum / k as f64).sqrt();
        let alpha = (self.current_dt / (0.1 + self.current_dt)).clamp(0.05, 0.3);
        self.theta_fixed_rel_error = alpha * raw + (1.0 - alpha) * self.theta_fixed_rel_error;
        self.last_theta_error_norm =
            self.theta_fixed_rel_error / self.theta_controller.target_error;

        let dt = self.last_proposed_dt.clamp(
            self.dt_controller.config.min_dt,
            self.dt_controller.config.max_dt,
        );

        let specific_energy_scale = denom / self.total_mass.max(1e-12);
        let vel = self.last_diag.max_vel;
        let acc = self.last_diag.max_acc;
        let jerk = self.last_diag.jerk;
        self.dt_fixed_rel_error =
            (vel * acc * dt * dt + jerk * dt * dt * dt) / specific_energy_scale.max(1e-12);
    }
}

impl System {
    pub fn zero_com_velocity(&mut self) {
        if calibration::zero_com_velocity(&mut self.bodies, self.total_mass) {
            self.reset_energy_baseline();
        }
    }

    pub fn recenter_com(&mut self) {
        calibration::recenter_com(&mut self.bodies, &mut self.trails, self.total_mass);
    }
}

impl System {
    pub fn add_body(&mut self, mut body: Body) {
        body.sync_physical_properties();
        body.softening = body.softening.max(body.physical_radius * 2.0);
        body.radius = body.radius.min(body.softening * 0.5);
        self.bodies.push(body);
        self.trails.push(VecDeque::with_capacity(self.trail_cap));
        self.refresh_collision_geometry();
        self.initial_energy = None;
        self.initial_energy_scale = None;
    }

    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.trails.clear();
        self.scratch_acc.clear();
        self.total_mass = 0.0;

        for b in bodies {
            self.total_mass += b.mass;
            self.bodies.push(b);
            self.trails.push(VecDeque::with_capacity(self.trail_cap));
        }

        self.initial_energy = None;
        self.initial_energy_scale = None;
        self.rel_energy_error = 0.0;
        self.max_rel_energy_error = 0.0;
        self.initial_angular_momentum = None;
        self.initial_angular_momentum_scale = None;
        self.rel_angular_momentum_error = 0.0;
        self.max_rel_angular_momentum_error = 0.0;
        self.current_dt = 0.0;
        self.steps = 0;
        self.last_potential = 0.0;
        self.theta_fixed_rel_error = 0.0;
        self.dt_fixed_rel_error = 0.0;
        self.last_theta_error_norm = 0.0;
        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();

        self.zero_com_velocity();
        self.recenter_com();
        self.refresh_collision_geometry();
    }

    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            self.bodies.swap_remove(index);
            self.trails.swap_remove(index);
            self.refresh_collision_geometry();
            self.initial_energy = None;
            self.initial_energy_scale = None;
            self.rel_energy_error = 0.0;
            self.max_rel_energy_error = 0.0;
        }
    }

    pub fn update_body(&mut self, index: usize, body: Body) {
        if let Some(slot) = self.bodies.get_mut(index) {
            let mass_changed = (slot.mass - body.mass).abs() > 1e-15;
            let mut updated = body;
            updated.sync_physical_properties();
            updated.softening = updated.softening.max(updated.physical_radius * 2.0);
            updated.radius = updated.radius.min(updated.softening * 0.5);
            *slot = updated;
            self.refresh_collision_geometry();
            if mass_changed {
                self.initial_energy = None;
                self.initial_energy_scale = None;
            }
        }
    }
}

impl System {
    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    pub fn trail(&self, index: usize) -> Option<&VecDeque<(f64, f64)>> {
        self.trails.get(index)
    }

    pub fn trails(&self) -> &[VecDeque<(f64, f64)>] {
        &self.trails
    }

    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    pub fn set_theta(&mut self, theta: f64) {
        self.theta_controller.set(theta);
    }

    pub fn set_cor(&mut self, cor: f64) {
        self.cor = cor.clamp(0.0, 1.0);
    }

    pub fn set_g_factor(&mut self, g: f64) {
        self.g_factor = g.max(0.0);
    }

    pub fn g_factor(&self) -> f64 {
        self.g_factor
    }

    pub fn engine(&self) -> &BarnesHutEngine {
        &self.engine
    }

    pub fn metrics(&self) -> Metrics {
        let kinetic = self.last_kinetic;
        let potential = self.last_potential;
        let total = total_energy(kinetic, potential);

        let lz = angular_momentum_z(&self.bodies);
        let (com_x, com_y, com_vx, com_vy) = center_of_mass_state(&self.bodies);

        Metrics {
            kinetic,
            potential,
            total_energy: total,
            angular_momentum_z: lz,
            rel_angular_momentum_error: self.rel_angular_momentum_error,
            max_rel_angular_momentum_error: self.max_rel_angular_momentum_error,
            com_x,
            com_y,
            com_vx,
            com_vy,
            g_factor: self.g_factor,
            theta: self.theta_controller.current(),
            dt: self.current_dt,
            rel_energy_error: self.rel_energy_error,
            max_rel_energy_error: self.max_rel_energy_error,
            theta_fixed_rel_error: self.theta_fixed_rel_error,
            dt_fixed_rel_error: self.dt_fixed_rel_error,
            last_theta_error_norm: self.last_theta_error_norm,
            theta_error_smoothed_norm: self.theta_controller.error(),
            dt_controller_state: self.dt_controller.last_dt(),
            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,
            merges_this_step: self.merges_this_step,
            bounces_this_step: self.bounces_this_step,
            near_miss_count: self.near_miss_count,
            fragments_spawned_this_step: self.fragments_spawned_this_step,
            hit_and_runs_this_step: self.hit_and_runs_this_step,
            total_dust_mass: self.total_dust_mass,
        }
    }

    /// Drain and return all collision events since the last call.
    /// The internal buffer is cleared; events accumulate across `steps_per_frame` physics steps.
    pub fn take_impact_events(&mut self) -> Vec<ImpactEvent> {
        std::mem::take(&mut self.last_impact_events)
    }

    /// Accelerations computed during the last integration step.
    /// Each entry corresponds to `bodies()[i]`.
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }
}
