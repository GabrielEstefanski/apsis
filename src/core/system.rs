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

/// Number of bodies that actually need individual trail rendering.
///
/// Belt members and sub-threshold bodies are excluded because their trails are
/// suppressed by the renderer anyway. Using this count for ring-buffer capacity
/// allocation keeps GPU memory proportional to what's actually rendered, not
/// to the total body count (which can be dominated by asteroid belts).
/// Generate an auto-name for a new body given existing names.
/// Counts existing names that start with the material prefix and appends N+1.
fn auto_name(material: crate::domain::materials::Material, existing: &[String]) -> String {
    let prefix = material.display_name();
    let count = existing.iter().filter(|n| n.starts_with(prefix)).count() + 1;
    format!("{prefix} {count}")
}

fn trail_body_count(bodies: &[Body]) -> usize {
    if bodies.is_empty() {
        return 0;
    }
    let max_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);
    if max_mass <= 0.0 {
        return bodies.len();
    }
    bodies
        .iter()
        .filter(|b| b.mass / max_mass > 1e-6)
        .count()
        .max(1)
}
use crate::physics::energy::{
    angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy,
};
use crate::physics::gravity::BarnesHutEngine;
use crate::physics::integrator::{Integrator, Y4_C, Y4_D, drift, evaluate_accelerations, kick};
use crate::physics::orbital::{self, OrbitalElements};

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

    /// Active integration algorithm.
    integrator: Integrator,

    /// Cached osculating orbital elements — one slot per body.
    /// Updated on demand via [`System::update_orbital_elements`], not every step.
    orbital_cache: Vec<Option<OrbitalElements>>,

    /// Global Plummer softening scale applied on top of the per-body
    /// mass-proportional default: `ε = EPS_BASE · m^(1/3) · softening_scale`.
    softening_scale: f64,

    /// Diagnostics subsystem.
    diagnostics: DiagnosticsComputer,
    last_diag: SimulationDiagnostics,

    /// Step counter.
    steps: u64,

    /// Total simulated time elapsed (t = steps × dt, but tracked as f64
    /// so it remains correct even if dt changes mid-run).
    t: f64,

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

    /// Human-readable label for each body, parallel to `bodies`.
    /// Kept separate because `Body` is `Copy` and cannot own a `String`.
    names: Vec<String>,

    /// Minimum pairwise separation cached from the most recent step.
    r_min: f64,

    /// Maximum effective pairwise softening length cached from the most recent step.
    softening_max: f64,
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
        let trail_n = trail_body_count(&bodies);
        let cap = adaptive_capacity(trail_n.max(1));
        let mut trail_buf = TrailBuffer::new_with_capacity(n, cap);
        trail_buf.update_colors(&bodies);

        let total_mass = bodies.iter().map(|b| b.mass).sum();
        let names = bodies.iter().map(|b| auto_name(b.material, &[])).collect::<Vec<_>>();
        // Re-generate with correct counters (so Star 1, Star 2 … instead of all "Star 1")
        let names = {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for b in &bodies {
                acc.push(auto_name(b.material, &acc));
            }
            acc
        };

        let (r_min, softening_max) = Self::compute_closeness(&bodies);

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
            integrator: Integrator::VelocityVerlet,
            orbital_cache: Vec::new(),
            softening_scale: 1.0,
            diagnostics: DiagnosticsComputer::new(),
            last_diag: SimulationDiagnostics::default(),
            steps: 0,
            t: 0.0,
            current_dt: dt,
            g_factor: 1.0,
            initial_angular_momentum: None,
            rel_angular_momentum_error: 0.0,
            abs_angular_momentum_error: 0.0,
            names,
            r_min,
            softening_max,
        }
    }
}

impl System {
    /// Advance the simulation by one time step using the configured integrator.
    pub fn step(&mut self) {
        match self.integrator {
            Integrator::VelocityVerlet => self.step_vv(),
            Integrator::Yoshida4 => self.step_yoshida4(),
        }

        let dt = self.current_dt;
        self.steps += 1;
        self.t += dt;

        self.last_diag = self
            .diagnostics
            .compute(&self.scratch_acc, &self.bodies, dt);

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

        // Update softening diagnostics every step (O(N²) but bounded by threshold).
        let (r_min, soft_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = soft_max;
    }

    // ── Velocity Verlet (KDK leapfrog, 2nd-order) ────────────────────────────

    fn step_vv(&mut self) {
        let dt = self.current_dt;
        let theta = self.theta;

        // Force at x(t)
        let raw_pe =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);
        self.last_potential = self.scale_acc_and_pe(raw_pe);

        kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt); // half-kick
        drift(&mut self.bodies, dt);

        // Force at x(t + dt)
        let raw_pe =
            evaluate_accelerations(&self.bodies, theta, &mut self.engine, &mut self.scratch_acc);
        self.last_potential = self.scale_acc_and_pe(raw_pe);

        kick(&mut self.bodies, &self.scratch_acc, 0.5 * dt); // half-kick
    }

    // ── Yoshida 4th-order (Forest–Ruth DKD, 4th-order) ───────────────────────
    //
    // Scheme: drift(c₀) → F → kick(d₀) → drift(c₁) → F → kick(d₁) → drift(c₂) → F → kick(d₂) → drift(c₃)
    //
    // d₁ = w₀ ≈ −1.70, so the middle sub-step is a *backward* kick — correct
    // and essential for the 4th-order cancellation of error terms.

    fn step_yoshida4(&mut self) {
        let dt = self.current_dt;
        let theta = self.theta;

        for i in 0..3 {
            drift(&mut self.bodies, Y4_C[i] * dt);

            let raw_pe = evaluate_accelerations(
                &self.bodies,
                theta,
                &mut self.engine,
                &mut self.scratch_acc,
            );
            self.last_potential = self.scale_acc_and_pe(raw_pe);

            kick(&mut self.bodies, &self.scratch_acc, Y4_D[i] * dt);
        }

        // Final drift to complete the DKD stencil
        drift(&mut self.bodies, Y4_C[3] * dt);
    }
}

impl System {
    /// Compute the minimum pairwise separation and maximum effective softening
    /// length over all body pairs.
    ///
    /// Skipped (returns sentinels) when N < 2 or N > [`N_CLOSENESS_THRESHOLD`],
    /// to keep overhead bounded for large asteroid-belt simulations.
    fn compute_closeness(bodies: &[Body]) -> (f64, f64) {
        const N_CLOSENESS_THRESHOLD: usize = 512;

        if bodies.len() < 2 || bodies.len() > N_CLOSENESS_THRESHOLD {
            return (f64::MAX, 0.0);
        }

        let mut r_min = f64::MAX;
        let mut soft_max = 0.0_f64;

        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                let dx = bodies[i].x - bodies[j].x;
                let dy = bodies[i].y - bodies[j].y;
                let r = (dx * dx + dy * dy).sqrt();
                if r < r_min {
                    r_min = r;
                }
                let eps2_ij = (bodies[i].softening * bodies[i].softening
                    + bodies[j].softening * bodies[j].softening)
                    * 0.5;
                let eps_ij = eps2_ij.sqrt();
                if eps_ij > soft_max {
                    soft_max = eps_ij;
                }
            }
        }

        (r_min, soft_max)
    }

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
        use crate::domain::body::default_softening;
        body.sync_physical_properties();
        if (self.softening_scale - 1.0).abs() > 1e-15 {
            body.softening = default_softening(body.mass) * self.softening_scale;
        }
        self.total_mass += body.mass;
        self.names.push(auto_name(body.material, &self.names));
        self.bodies.push(body);

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
    }

    /// Add multiple bodies in a single batch.
    ///
    /// More efficient than calling [`add_body`] in a loop: the trail buffer is
    /// reset only once and the energy baseline is invalidated once.
    pub fn add_bodies(&mut self, new_bodies: Vec<Body>) {
        use crate::domain::body::default_softening;
        for mut body in new_bodies {
            body.sync_physical_properties();
            if (self.softening_scale - 1.0).abs() > 1e-15 {
                body.softening = default_softening(body.mass) * self.softening_scale;
            }
            self.total_mass += body.mass;
            self.names.push(auto_name(body.material, &self.names));
            self.bodies.push(body);
        }

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);
        self.initial_energy = None;
    }

    /// Read the display name for body `idx`.
    pub fn name(&self, idx: usize) -> &str {
        self.names.get(idx).map(|s| s.as_str()).unwrap_or("")
    }

    /// All body names (parallel to `bodies()`).
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Rename body `idx`. Silently ignores out-of-range indices.
    pub fn set_name(&mut self, idx: usize, name: String) {
        if let Some(slot) = self.names.get_mut(idx) {
            *slot = name;
        }
    }

    /// Replaces the entire set of bodies in the simulation.
    ///
    /// All previous state is cleared, the trail buffer is reset, and the
    /// system is normalised to its COM rest frame.
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.scratch_acc.clear();
        self.names.clear();
        self.total_mass = 0.0;

        for mut b in bodies {
            b.sync_physical_properties();
            self.total_mass += b.mass;
            self.names.push(auto_name(b.material, &self.names));
            self.bodies.push(b);
        }

        let n = self.bodies.len();
        let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);

        self.initial_energy = None;
        self.rel_energy_error = 0.0;
        self.steps = 0;
        self.t = 0.0;
        self.last_potential = 0.0;
        self.last_kinetic = 0.0;
        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();

        self.zero_com_velocity();
        self.recenter_com();

        let (r_min, softening_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }

    /// Removes a body from the simulation.
    ///
    /// Uses `swap_remove` for O(1) removal.  The trail buffer is reset
    /// because body indices change.
    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            let removed = self.bodies.swap_remove(index);
            self.total_mass -= removed.mass;
            if index < self.names.len() {
                self.names.swap_remove(index);
            }

            let n = self.bodies.len();
            let cap = adaptive_capacity(trail_body_count(&self.bodies).max(1));
            self.trail_buf.reset(n, cap);
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

    /// Total simulated time elapsed.
    pub fn t(&self) -> f64 {
        self.t
    }

    /// Number of integration steps completed.
    pub fn steps(&self) -> u64 {
        self.steps
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

    /// Returns the active integrator.
    pub fn integrator(&self) -> Integrator {
        self.integrator
    }

    /// Switches the integration algorithm.  Takes effect on the next [`step`].
    pub fn set_integrator(&mut self, i: Integrator) {
        self.integrator = i;
    }

    /// Returns the current Barnes–Hut opening angle θ.
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Sets the Barnes–Hut opening angle θ (clamped to [0.1, 1.5]).
    ///
    /// Smaller θ → more accurate (approaches O(N²) as θ → 0).
    /// Larger θ → faster but less accurate.
    pub fn set_theta(&mut self, theta: f64) {
        self.theta = theta.clamp(0.05, 1.5);
    }

    /// Returns the current global softening scale factor.
    pub fn softening_scale(&self) -> f64 {
        self.softening_scale
    }

    /// Sets a global Plummer softening scale applied on top of the
    /// per-body mass-proportional default (`ε = ε_default · scale`).
    ///
    /// Also rescales all existing body softenings immediately.
    pub fn set_softening_scale(&mut self, scale: f64) {
        use crate::domain::body::default_softening;
        self.softening_scale = scale.max(0.0);
        for b in &mut self.bodies {
            b.softening = default_softening(b.mass) * self.softening_scale;
        }
    }

    pub fn trail_every(&self) -> usize {
        self.trail_every
    }

    pub fn set_trail_every(&mut self, n: usize) {
        self.trail_every = n.max(1);
    }

    /// Records the current body positions into the trail ring buffer.
    ///
    /// Call this **once per rendered frame** (not per physics step) so the
    /// trail density is proportional to the amount of simulated time per
    /// frame rather than to a fixed physics-step count.
    pub fn push_trail(&mut self) {
        self.trail_buf.push(&self.bodies);
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

            t: self.t,
            steps: self.steps,

            integrator: self.integrator,
            g_factor: self.g_factor,
            theta: self.theta,
            dt: self.current_dt,

            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,

            r_min: self.r_min,
            softening_max: self.softening_max,
        }
    }

    /// Returns accelerations computed during the last integration step.
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }

    // ── Orbital elements ─────────────────────────────────────────────────────

    /// Recomputes osculating orbital elements for all bodies and caches the result.
    ///
    /// This is O(N²) and should be called **once per rendered frame**, not every
    /// physics step. The result is available via [`orbital_elements`].
    pub fn update_orbital_elements(&mut self) {
        self.orbital_cache = orbital::compute_all(&self.bodies, self.g_factor);
    }

    /// Returns the cached osculating orbital elements (one slot per body).
    ///
    /// Call [`update_orbital_elements`] first to get fresh values.
    pub fn orbital_elements(&self) -> &[Option<OrbitalElements>] {
        &self.orbital_cache
    }

    // ── Snapshot (save / load) ───────────────────────────────────────────────

    /// Capture the minimal state required for deterministic reproduction.
    pub fn to_snapshot(&self) -> crate::core::snapshot::SimSnapshot {
        use crate::core::snapshot::{BodyRecord, SimSnapshot};
        SimSnapshot {
            save_id: 0, // caller sets this via new_id() or save_to_dir()
            t: self.t,
            steps: self.steps,
            dt: self.current_dt,
            theta: self.theta,
            softening_scale: self.softening_scale,
            g_factor: self.g_factor,
            integrator: self.integrator,
            trail_every: self.trail_every,
            sim_name: String::new(), // set by the app layer before saving
            seed: 0,                 // set by the app layer before saving
            trail: None,             // set by the app layer before saving
            bodies: self.bodies.iter().map(BodyRecord::from_body).collect(),
            names: self.names.clone(),
        }
    }

    /// Replace the current simulation state with a saved snapshot.
    ///
    /// The trail buffer is cleared (it is cosmetic and cannot be restored).
    /// Energy / angular-momentum references are reset so the first post-load
    /// step establishes new baselines.
    pub fn restore_from_snapshot(&mut self, snap: &crate::core::snapshot::SimSnapshot) {
        let bodies: Vec<Body> = snap.bodies.iter().map(|r| r.into_body()).collect();
        // Restore names: use saved names if present, else auto-generate
        self.names = if snap.names.len() == bodies.len() {
            snap.names.clone()
        } else {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for b in &bodies { acc.push(auto_name(b.material, &acc)); }
            acc
        };
        let n = bodies.len();

        self.bodies = bodies;
        self.total_mass = self.bodies.iter().map(|b| b.mass).sum();
        self.scratch_acc.clear();

        // Rebuild trail buffer — restore saved trail if dimensions match,
        // otherwise start empty.
        let cap =
            crate::core::trail_buffer::adaptive_capacity(trail_body_count(&self.bodies).max(1));
        self.trail_buf.reset(n, cap);
        self.trail_buf.update_colors(&self.bodies);
        if let Some(trail_snap) = &snap.trail {
            if trail_snap.n_bodies == n as u32
                && trail_snap.positions.len() == (trail_snap.n_bodies * trail_snap.capacity) as usize
            {
                self.trail_buf.restore_from_snapshot(trail_snap);
            }
        }

        // Restore physics parameters
        self.t = snap.t;
        self.steps = snap.steps;
        self.current_dt = snap.dt;
        self.theta = snap.theta;
        self.softening_scale = snap.softening_scale;
        self.g_factor = snap.g_factor;
        self.integrator = snap.integrator;
        self.trail_every = snap.trail_every.max(1);

        // Reset energy / angular-momentum baselines so next step establishes
        // fresh references relative to the restored state.
        self.initial_energy = None;
        self.initial_angular_momentum = None;
        self.rel_energy_error = 0.0;
        self.rel_angular_momentum_error = 0.0;
        self.abs_angular_momentum_error = 0.0;
        self.last_kinetic = 0.0;
        self.last_potential = 0.0;
        self.diagnostics = crate::core::diagnostics::DiagnosticsComputer::new();
        self.last_diag = crate::core::diagnostics::SimulationDiagnostics::default();
        self.orbital_cache.clear();

        let (r_min, softening_max) = Self::compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }
}
