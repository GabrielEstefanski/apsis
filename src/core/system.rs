//! Simulation orchestrator for an N-body gravitational system.
//!
//! This module defines the [`System`] type, responsible for advancing
//! the state of a set of massive bodies interacting via gravity.
//!
//! ## Design goals
//!
//! - Preserve physical consistency with a symplectic integrator (Velocity Verlet)
//! - Maintain deterministic and reproducible evolution
//! - Provide diagnostics for energy and angular momentum conservation
//! - Support Barnes–Hut acceleration for scalable simulations
//!
//! ## Important assumptions
//!
//! - Time step (`dt`) is constant (required for symplectic behavior)
//! - The force field is evaluated consistently at well-defined points
//! - No discrete events (e.g., collisions) are applied within integration steps
//!
//! ## Notes
//!
//! This system is intended for scientific and numerical experiments in
//! gravitational dynamics, not for general-purpose physics engines.

use crate::core::calibration;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::metrics::Metrics;
use crate::core::trail_buffer::{TrailBuffer, adaptive_capacity};
use crate::domain::body::Body;
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{drift, evaluate_accelerations, half_kick};

/// Central simulation state for an N-body gravitational system.
pub struct System {
    /// Bodies participating in the simulation.
    bodies: Vec<Body>,

    /// GPU-ready ring buffer of trail positions and colours.
    trail_buf: TrailBuffer,
    trail_every: usize,

    /// Total mass of the system (used for COM recentering).
    total_mass: f64,

    /// Last computed energies.
    last_kinetic: f64,
    last_potential: f64,

    /// Initial total energy (used as reference).
    initial_energy: Option<f64>,

    /// Relative energy error (diagnostic only).
    rel_energy_error: f64,

    /// Barnes–Hut engine for approximate force computation.
    engine: BarnesHutEngine,

    /// Scratch buffer for accelerations.
    scratch_acc: Vec<(f64, f64)>,

    /// Barnes–Hut opening angle parameter (θ).
    theta: f64,

    /// Diagnostics subsystem.
    diagnostics: DiagnosticsComputer,
    last_diag: SimulationDiagnostics,

    /// Step counter.
    steps: u64,

    /// Fixed time step.
    current_dt: f64,

    /// Gravitational scaling factor (G multiplier).
    g_factor: f64,

    /// Initial angular momentum (z-component) used as reference.
    initial_angular_momentum: Option<f64>,

    /// Relative angular momentum error (diagnostic only).
    rel_angular_momentum_error: f64,

    /// Absolute angular momentum error (always meaningful).
    abs_angular_momentum_error: f64,
}

impl System {
    /// Creates a new simulation system.
    ///
    /// # Parameters
    ///
    /// - `bodies`: Initial set of bodies
    /// - `theta`: Barnes–Hut opening angle (controls accuracy vs performance)
    /// - `dt`: Fixed time step
    /// - `max_depth`: Maximum tree depth for Barnes–Hut
    /// - `trail_every`: Sampling interval for trails (ring-buffer depth is
    ///   chosen automatically via [`adaptive_capacity`])
    ///
    /// # Notes
    ///
    /// - Smaller `theta` increases accuracy (approaches O(N²))
    /// - Smaller `dt` improves stability and energy conservation
    pub fn new(
        bodies: Vec<Body>,
        theta: f64,
        dt: f64,
        max_depth: usize,
        trail_every: usize,
    ) -> Self {
        let n = bodies.len();
        let mut trail_buf = TrailBuffer::new(n);
        trail_buf.update_colors(&bodies);

        let total_mass = bodies.iter().map(|b| b.mass).sum();

        Self {
            bodies,
            trail_buf,
            trail_every: trail_every.max(1),
            total_mass,
            last_kinetic: 0.0,
            last_potential: 0.0,
            initial_energy: None,
            rel_energy_error: 0.0,
            engine: BarnesHutEngine::new(max_depth),
            scratch_acc: Vec::new(),
            theta,
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            steps: 0,
            current_dt: dt,
            g_factor: 1.0,
            initial_angular_momentum: None,
            rel_angular_momentum_error: 0.0,
            abs_angular_momentum_error: 0.0,
        }
    }
}

impl System {
    pub fn step(&mut self) {
        let dt = self.current_dt;
        let theta = self.theta;

        let raw_pe =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);

        self.last_potential = self.scale_acc_and_pe(raw_pe);

        half_kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);

        drift(&mut self.bodies, dt);

        let raw_pe_after =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);
        self.last_potential = self.scale_acc_and_pe(raw_pe_after);

        half_kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt);

        self.steps += 1;

        self.last_diag = self
            .diagnostics
            .compute(&self.scratch_acc, &self.bodies, dt);

        if self.steps % self.trail_every as u64 == 0 {
            self.trail_buf.push(&self.bodies);
        }

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        // Periodically remove COM drift.  The trail buffer is translated by
        // the same vector so stored positions remain consistent.
        if self.steps % 97 == 0 {
            if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
                calibration::apply_body_shift(&mut self.bodies, dx, dy);
                self.trail_buf.translate(-dx as f32, -dy as f32);
            }
        }
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
    /// Updates energy diagnostics for the current simulation state.
    fn update_energy_tracking(&mut self) {
        let kinetic = kinetic_energy(&self.bodies);
        self.last_kinetic = kinetic;

        let total = total_energy(kinetic, self.last_potential);

        let baseline = match self.initial_energy {
            Some(v) => v,
            None => {
                self.initial_energy = Some(total);
                total
            }
        };

        let denom = baseline.abs().max(1e-12);
        self.rel_energy_error = (total - baseline) / denom;
    }

    /// Updates angular momentum diagnostics.
    fn update_angular_momentum_tracking(&mut self) {
        let lz = angular_momentum_z(&self.bodies);

        let baseline = match self.initial_angular_momentum {
            Some(v) => v,
            None => {
                self.initial_angular_momentum = Some(lz);
                lz
            }
        };

        self.abs_angular_momentum_error = (lz - baseline).abs();

        let denom = baseline.abs().max(1e-12);
        self.rel_angular_momentum_error = (lz - baseline) / denom;
    }

    /// Removes the centre-of-mass velocity so the system is in its rest frame.
    pub fn zero_com_velocity(&mut self) {
        calibration::zero_com_velocity(&mut self.bodies, self.total_mass);
    }

    /// Recenters the system so that the centre of mass is at the origin.
    ///
    /// The trail buffer is translated by the same vector so stored positions
    /// remain visually consistent.
    pub fn recenter_com(&mut self) {
        if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
            calibration::apply_body_shift(&mut self.bodies, dx, dy);
            self.trail_buf.translate(-dx as f32, -dy as f32);
        }
    }
}

impl System {
    /// Adds a new body to the simulation.
    ///
    /// The trail buffer is reset to accommodate the new body count; trail
    /// history is lost.  Energy baseline is reset because the system
    /// topology has changed.
    pub fn add_body(&mut self, mut body: Body) {
        body.sync_physical_properties();
        self.total_mass += body.mass;
        self.bodies.push(body);

        let n = self.bodies.len();
        self.trail_buf.reset(n, adaptive_capacity(n));
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
    }

    /// Replaces the entire set of bodies in the simulation.
    ///
    /// All previous state is cleared, the trail buffer is reset, and the
    /// system is normalised to its COM rest frame.
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.scratch_acc.clear();
        self.total_mass = 0.0;

        for mut b in bodies {
            b.sync_physical_properties();
            self.total_mass += b.mass;
            self.bodies.push(b);
        }

        let n = self.bodies.len();
        self.trail_buf.reset(n, adaptive_capacity(n));
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
        self.rel_energy_error = 0.0;
        self.steps = 0;
        self.last_potential = 0.0;
        self.last_kinetic = 0.0;
        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();

        self.zero_com_velocity();
        self.recenter_com();
    }

    /// Removes a body from the simulation.
    ///
    /// Uses `swap_remove` for O(1) removal.  The trail buffer is reset
    /// because body indices change.
    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            let removed = self.bodies.swap_remove(index);
            self.total_mass -= removed.mass;

            let n = self.bodies.len();
            self.trail_buf.reset(n, adaptive_capacity(n));
            self.trail_buf.update_colors(&self.bodies);

            self.initial_energy = None;
            self.rel_energy_error = 0.0;
        }
    }

    /// Updates a body in-place, recomputing derived physical properties.
    ///
    /// If the body colour changes, the trail colour buffer is re-uploaded on
    /// the next render frame.
    pub fn update_body(&mut self, index: usize, mut body: Body) {
        if let Some(slot) = self.bodies.get_mut(index) {
            let mass_changed = (slot.mass - body.mass).abs() > 1e-15;

            body.sync_physical_properties();

            if mass_changed {
                self.total_mass += body.mass - slot.mass;
            }

            *slot = body;

            if mass_changed {
                self.initial_energy = None;
                self.rel_energy_error = 0.0;
            }

            self.trail_buf.update_colors(&self.bodies);
        }
    }
}

impl System {
    /// Returns an immutable slice of all bodies in the simulation.
    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }

    pub fn dt(&self) -> f64 {
        self.current_dt
    }

    /// Returns a shared reference to the GPU-ready trail ring buffer.
    pub fn trail_buf(&self) -> &TrailBuffer {
        &self.trail_buf
    }

    /// Returns a mutable reference to the GPU-ready trail ring buffer.
    ///
    /// Required by the trail renderer to drain dirty flags each frame.
    pub fn trail_buf_mut(&mut self) -> &mut TrailBuffer {
        &mut self.trail_buf
    }

    /// Returns the total mass of the system.
    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    /// Sets the gravitational scaling factor.
    pub fn set_g_factor(&mut self, g: f64) {
        self.g_factor = g.max(0.0);
    }

    /// Returns the current gravitational scaling factor.
    pub fn g_factor(&self) -> f64 {
        self.g_factor
    }

    pub fn set_dt(&mut self, dt: f64) {
        self.current_dt = dt;
    }

    /// Returns a reference to the Barnes–Hut engine.
    pub fn engine(&self) -> &BarnesHutEngine {
        &self.engine
    }

    /// Returns diagnostic metrics for the current simulation state.
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
            rel_energy_error: self.rel_energy_error,

            angular_momentum_z: lz,
            rel_angular_momentum_error: self.rel_angular_momentum_error,
            abs_angular_momentum_error: self.abs_angular_momentum_error,

            com_x,
            com_y,
            com_vx,
            com_vy,

            g_factor: self.g_factor,
            theta: self.theta,
            dt: self.current_dt,

            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,
        }
    }

    /// Returns accelerations computed during the last integration step.
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }
}
