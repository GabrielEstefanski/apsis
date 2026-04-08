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
use crate::domain::body::Body;
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{drift, evaluate_accelerations, half_kick};
use std::collections::VecDeque;

/// Central simulation state for an N-body gravitational system.
pub struct System {
    /// Bodies participating in the simulation.
    bodies: Vec<Body>,

    /// Trajectory history for visualization/debugging.
    trails: Vec<VecDeque<(f64, f64)>>,
    trail_cap: usize,
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
    /// - `trail_cap`: Maximum stored points per trajectory
    /// - `trail_every`: Sampling interval for trails
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
        trail_cap: usize,
        trail_every: usize,
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
            for (i, b) in self.bodies.iter().enumerate() {
                let t = &mut self.trails[i];
                t.push_back((b.x, b.y));
                if t.len() > self.trail_cap {
                    t.pop_front();
                }
            }
        }

        self.update_energy_tracking();
        self.update_angular_momentum_tracking();

        if self.steps % 97 == 0 {
            calibration::recenter_com(&mut self.bodies, &mut self.trails, self.total_mass);
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
    ///
    /// # Computes
    ///
    /// - Kinetic energy
    /// - Total mechanical energy
    /// - Relative energy drift from the initial state
    ///
    /// # Notes
    ///
    /// - The baseline energy is initialized on first call
    /// - In a correct symplectic simulation:
    ///     - Energy should oscillate around the baseline
    ///     - No long-term drift should be observed
    fn update_energy_tracking(&mut self) {
        let kinetic = kinetic_energy(&self.bodies);
        self.last_kinetic = kinetic;

        let total = total_energy(kinetic, self.last_potential);

        // Initialize baseline once
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
    ///
    /// # Computes
    ///
    /// - Absolute angular momentum error
    /// - Relative angular momentum error (when meaningful)
    ///
    /// # Notes
    ///
    /// - Relative error becomes unstable when angular momentum ≈ 0
    /// - Absolute error is always reliable
    fn update_angular_momentum_tracking(&mut self) {
        let lz = angular_momentum_z(&self.bodies);

        let baseline = match self.initial_angular_momentum {
            Some(v) => v,
            None => {
                self.initial_angular_momentum = Some(lz);
                lz
            }
        };

        // Absolute error (always valid)
        self.abs_angular_momentum_error = (lz - baseline).abs();

        // Relative error (only meaningful when baseline is not near zero)
        let denom = baseline.abs().max(1e-12);
        self.rel_angular_momentum_error = (lz - baseline) / denom;
    }

    /// Resets the velocity of the system so that the center of mass is at rest.
    ///
    /// # Notes
    ///
    /// - This operation preserves relative motion
    /// - It removes global drift from numerical accumulation
    /// - Does not affect internal dynamics
    pub fn zero_com_velocity(&mut self) {
        calibration::zero_com_velocity(&mut self.bodies, self.total_mass);
    }

    /// Recenters the system so that the center of mass is at the origin.
    ///
    /// # Notes
    ///
    /// - Affects only the reference frame
    /// - Does not alter relative trajectories
    /// - Useful to avoid numerical drift over long simulations
    pub fn recenter_com(&mut self) {
        calibration::recenter_com(&mut self.bodies, &mut self.trails, self.total_mass);
    }
}

impl System {
    /// Adds a new body to the simulation.
    ///
    /// # Notes
    ///
    /// - Resets energy baseline since the system state changes
    /// - Does NOT apply any collision-related adjustments
    /// - Assumes the caller provides physically consistent values
    pub fn add_body(&mut self, mut body: Body) {
        body.sync_physical_properties();

        self.total_mass += body.mass;

        self.bodies.push(body);
        self.trails.push(VecDeque::with_capacity(self.trail_cap));

        // Reset energy baseline due to topology change
        self.initial_energy = None;
    }

    /// Replaces the entire set of bodies in the simulation.
    ///
    /// # Behavior
    ///
    /// - Clears all previous state
    /// - Recomputes total mass
    /// - Resets diagnostics and energy tracking
    ///
    /// # Notes
    ///
    /// - Center-of-mass is re-centered
    /// - System starts from a clean deterministic state
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.trails.clear();
        self.scratch_acc.clear();

        self.total_mass = 0.0;

        for mut b in bodies {
            b.sync_physical_properties();
            self.total_mass += b.mass;

            self.bodies.push(b);
            self.trails.push(VecDeque::with_capacity(self.trail_cap));
        }

        // Reset simulation state
        self.initial_energy = None;
        self.rel_energy_error = 0.0;

        self.steps = 0;
        self.last_potential = 0.0;
        self.last_kinetic = 0.0;

        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();

        // Normalize reference frame
        self.zero_com_velocity();
        self.recenter_com();
    }

    /// Removes a body from the system.
    ///
    /// # Notes
    ///
    /// - Uses `swap_remove` for O(1) removal
    /// - Resets energy baseline due to system change
    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            let removed = self.bodies.swap_remove(index);
            self.trails.swap_remove(index);

            self.total_mass -= removed.mass;

            self.initial_energy = None;
            self.rel_energy_error = 0.0;
        }
    }

    /// Updates a body in-place.
    ///
    /// # Behavior
    ///
    /// - Replaces the body at the given index
    /// - Recomputes derived physical properties
    ///
    /// # Notes
    ///
    /// - If mass changes, energy baseline is reset
    pub fn update_body(&mut self, index: usize, mut body: Body) {
        if let Some(slot) = self.bodies.get_mut(index) {
            let mass_changed = (slot.mass - body.mass).abs() > 1e-15;

            body.sync_physical_properties();

            // Update total mass if needed
            if mass_changed {
                self.total_mass += body.mass - slot.mass;
            }

            *slot = body;

            if mass_changed {
                self.initial_energy = None;
                self.rel_energy_error = 0.0;
            }
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

    /// Returns the trajectory (trail) of a specific body, if it exists.
    ///
    /// # Parameters
    /// - `index`: Index of the body
    ///
    /// # Returns
    /// - `Some(&trail)` if the body exists
    /// - `None` otherwise
    pub fn trail(&self, index: usize) -> Option<&VecDeque<(f64, f64)>> {
        self.trails.get(index)
    }

    /// Returns all stored trajectory trails.
    ///
    /// Each entry corresponds to one body.
    pub fn trails(&self) -> &[VecDeque<(f64, f64)>] {
        &self.trails
    }

    /// Returns the total mass of the system.
    ///
    /// This is maintained incrementally for performance.
    pub fn total_mass(&self) -> f64 {
        self.total_mass
    }

    /// Sets the gravitational scaling factor.
    ///
    /// # Notes
    ///
    /// - Effective gravitational constant becomes `G = g_factor`
    /// - Must be non-negative
    /// - Should remain constant during a simulation run
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
    ///
    /// Useful for inspection or advanced diagnostics.
    pub fn engine(&self) -> &BarnesHutEngine {
        &self.engine
    }

    /// Returns diagnostic metrics for the current simulation state.
    ///
    /// # Includes
    ///
    /// - Kinetic, potential, and total energy
    /// - Angular momentum (z-component)
    /// - Center-of-mass position and velocity
    /// - Relative energy error
    /// - Integration diagnostics (max acceleration, velocity, jerk)
    ///
    /// # Notes
    ///
    /// - All values correspond to the last completed integration step
    /// - No adaptive or heuristic metrics are included
    pub fn metrics(&self) -> Metrics {
        let kinetic = self.last_kinetic;
        let potential = self.last_potential;
        let total = total_energy(kinetic, potential);

        let lz = angular_momentum_z(&self.bodies);
        let (com_x, com_y, com_vx, com_vy) = center_of_mass_state(&self.bodies);

        Metrics {
            // ── Energetics ─────────────────────────────── //
            kinetic,
            potential,
            total_energy: total,
            rel_energy_error: self.rel_energy_error,

            // ── Angular momentum ───────────────────────── //
            angular_momentum_z: lz,
            rel_angular_momentum_error: self.rel_angular_momentum_error,
            abs_angular_momentum_error: self.abs_angular_momentum_error,

            // ── Center of mass ─────────────────────────── //
            com_x,
            com_y,
            com_vx,
            com_vy,

            // ── Simulation parameters ──────────────────── //
            g_factor: self.g_factor,
            theta: self.theta,
            dt: self.current_dt,

            // ── Diagnostics ────────────────────────────── //
            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,
        }
    }

    /// Returns accelerations computed during the last integration step.
    ///
    /// Each entry corresponds to `bodies()[i]`.
    ///
    /// # Notes
    ///
    /// - Values are updated during the last force evaluation
    /// - Useful for debugging or numerical analysis
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }
}
