use crate::core::adaptive::{AccelerationStats, DtController, ThetaController};
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::domain::body::Body;
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{step, StepResult};
use std::collections::VecDeque;

#[derive(Debug, Clone, Copy)]
pub struct Metrics {
    pub kinetic: f64,
    pub potential: f64,
    pub total_energy: f64,
    pub angular_momentum_z: f64,
    pub com_x: f64,
    pub com_y: f64,
    pub com_vx: f64,
    pub com_vy: f64,
    pub theta: f64,
    pub dt: f64,
    pub rel_energy_error: f64,
    pub max_rel_energy_error: f64,
    pub theta_fixed_rel_error: f64,
    pub dt_fixed_rel_error: f64,
    pub last_theta_error_norm: f64,
    pub theta_error_smoothed_norm: f64,
    pub dt_controller_state: f64,
    pub max_acc: f64,
    pub jerk: f64,
}

pub struct System {
    bodies: Vec<Body>,
    trails: Vec<VecDeque<(f64, f64)>>,
    trail_cap: usize,
    trail_every: usize,

    last_potential: f64,
    last_kinetic: f64,
    total_mass: f64,
    steps: u64,

    engine: BarnesHutEngine,
    scratch_acc: Vec<(f64, f64)>,

    initial_energy: Option<f64>,
    initial_energy_scale: Option<f64>,
    rel_energy_error: f64,
    max_rel_energy_error: f64,

    theta_controller: ThetaController,
    dt_controller: DtController,

    diagnostics: DiagnosticsComputer,
    last_diag: SimulationDiagnostics,

    current_dt: f64,
    probe_interval: u64,

    theta_fixed_rel_error: f64,
    dt_fixed_rel_error: f64,

    last_proposed_dt: f64,
    last_theta_error_norm: f64,
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
            last_potential: 0.0,
            last_kinetic: 0.0,
            total_mass,
            steps: 0,
            engine: BarnesHutEngine::new(max_depth),
            scratch_acc: Vec::new(),
            initial_energy: None,
            initial_energy_scale: None,
            rel_energy_error: 0.0,
            max_rel_energy_error: 0.0,
            theta_controller,
            dt_controller,
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            current_dt: 0.0,
            probe_interval,
            theta_fixed_rel_error: 0.0,
            dt_fixed_rel_error: 0.0,
            last_proposed_dt: 0.0,
            last_theta_error_norm: 0.0,
        }
    }

    pub fn step_adaptive(&mut self, proposed_dt: f64) -> f64 {
        self.last_proposed_dt = proposed_dt;

        let stats = AccelerationStats::new(
            self.last_diag.max_acc,
            self.last_diag.jerk,
        );

        let dt = self.dt_controller.update(
            proposed_dt,
            self.rel_energy_error,
            stats,
        );

        self.step(dt);

        dt
    }

    pub fn step(&mut self, dt: f64) {
        self.last_proposed_dt = dt;
        self.current_dt = dt;

        let n = self.bodies.len();
        if self.scratch_acc.len() != n {
            self.scratch_acc.resize(n, (0.0, 0.0));
        }

        let theta = self.theta_controller.current();

        let StepResult { acc1, potential } = step(
            &mut self.bodies,
            dt,
            theta,
            &mut self.engine,
            &mut self.scratch_acc,
        );

        self.last_potential = potential;

        self.steps += 1;

        self.last_diag = self.diagnostics.compute(
            acc1,
            &self.bodies,
            dt,
        );

        if self.steps % self.trail_every as u64 == 0 {
            for (i, b) in self.bodies.iter().enumerate() {
                let t = &mut self.trails[i];
                t.push_back((b.x, b.y));
                if t.len() > self.trail_cap {
                    t.pop_front();
                }
            }
        }

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

        let denom = self.initial_energy_scale.unwrap_or_else(|| baseline.abs().max(1e-12));
        self.rel_energy_error = (total - baseline) / denom;

        if self.rel_energy_error.abs() > self.max_rel_energy_error {
            self.max_rel_energy_error = self.rel_energy_error.abs();
        }

        if self.steps % self.probe_interval == 0 {
            self.update_separated_errors();
        }

        self.adapt_theta();
    }

    pub fn reset_energy_baseline(&mut self) {
        self.initial_energy = None;
        self.initial_energy_scale = None;
        self.rel_energy_error = 0.0;
        self.max_rel_energy_error = 0.0;
    }

    pub fn zero_com_velocity(&mut self) {
        if self.total_mass <= 0.0 { return; }
        let com_vx: f64 = self.bodies.iter().map(|b| b.mass * b.vx).sum::<f64>() / self.total_mass;
        let com_vy: f64 = self.bodies.iter().map(|b| b.mass * b.vy).sum::<f64>() / self.total_mass;
        for b in &mut self.bodies {
            b.vx -= com_vx;
            b.vy -= com_vy;
        }
    }

    fn adapt_theta(&mut self) {
        self.theta_controller
            .update(self.theta_fixed_rel_error, self.current_dt);
    }

    fn update_separated_errors(&mut self) {
        let n = self.bodies.len();
        if n == 0 { return; }

        let base_energy = total_energy(self.last_kinetic, self.last_potential);
        let denom = self.initial_energy_scale.unwrap_or_else(|| base_energy.abs().max(1e-9));
        let theta = self.theta_controller.current();

        let k         = ((n as f64).sqrt().ceil() as usize).min(n);
        let step_size = (n / k).max(1);
        let mut idx   = (self.steps as usize) % n;
        let mut proxy_sum = 0.0_f64;

        for _ in 0..k {
            let e = self.engine.theta_error_proxy(idx, &self.bodies, theta);
            proxy_sum += e * e;
            idx = (idx + step_size) % n;
        }

        let raw_theta_error = (proxy_sum / k as f64).sqrt();
        let alpha = (self.current_dt / (0.1 + self.current_dt)).clamp(0.05, 0.3);
        self.theta_fixed_rel_error =
            alpha * raw_theta_error + (1.0 - alpha) * self.theta_fixed_rel_error;
        self.last_theta_error_norm =
            self.theta_fixed_rel_error / self.theta_controller.target_error;

        let dt = self.last_proposed_dt
            .clamp(self.dt_controller.config.min_dt, self.dt_controller.config.max_dt);

        let specific_energy_scale = denom / self.total_mass.max(1e-12);
        let vel  = self.last_diag.max_vel;
        let acc  = self.last_diag.max_acc;
        let jerk = self.last_diag.jerk;
        self.dt_fixed_rel_error =
            (vel * acc * dt * dt + jerk * dt * dt * dt) / specific_energy_scale.max(1e-12);
    }

    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    pub fn trail(&self, index: usize) -> Option<&VecDeque<(f64, f64)>> {
        self.trails.get(index)
    }

    pub fn trails(&self) -> &[VecDeque<(f64, f64)>] {
        &self.trails
    }

    pub fn add_body(&mut self, body: Body) {
        self.total_mass += body.mass;
        self.bodies.push(body);
        self.trails.push(VecDeque::with_capacity(self.trail_cap));
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
        self.current_dt = 0.0;
        self.steps = 0;
        self.last_potential = 0.0;
        self.theta_fixed_rel_error = 0.0;
        self.dt_fixed_rel_error = 0.0;
        self.last_theta_error_norm = 0.0;
        self.diagnostics = crate::core::diagnostics::DiagnosticsComputer::new();
        self.last_diag = crate::core::diagnostics::SimulationDiagnostics::default();
    }

    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            self.total_mass -= self.bodies[index].mass;
            self.bodies.swap_remove(index);
            self.trails.swap_remove(index);
            self.initial_energy = None;
            self.initial_energy_scale = None;
            self.rel_energy_error = 0.0;
            self.max_rel_energy_error = 0.0;
        }
    }

    pub fn set_theta(&mut self, theta: f64) {
        self.theta_controller.set(theta);
    }

    pub fn total_mass(&self) -> f64 { self.total_mass }

    pub fn metrics(&self) -> Metrics {
        let kinetic   = self.last_kinetic;
        let potential = self.last_potential;
        let total     = total_energy(kinetic, potential);

        let lz = angular_momentum_z(&self.bodies);
        let (com_x, com_y, com_vx, com_vy) =
            center_of_mass_state(&self.bodies);

        Metrics {
            kinetic,
            potential,
            total_energy: total,
            angular_momentum_z: lz,
            com_x,
            com_y,
            com_vx,
            com_vy,
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
        }
    }
}