//! Body management: add, remove, update, load, naming, COM calibration.

use crate::core::calibration;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::system::System;
use crate::core::system::helpers::{
    DEFAULT_NAME_PREFIX, auto_name, compute_closeness, resolved_name,
};
use crate::domain::body::{Body, NamedBody};

impl System {
    /// Adds a new body to the simulation.
    ///
    /// The trail buffer is reset to accommodate the new body count; trail
    /// history is lost. Energy baseline is reset because the system
    /// topology has changed. Also invalidates any template-rebuild source
    /// remembered by [`from_template`](System::from_template): a later
    /// [`with_seed`](System::with_seed) will no longer overwrite this
    /// manual addition.
    ///
    /// Auto-naming uses the generic `"Body N"` prefix because [`Body`]
    /// no longer carries a preset reference. Spawn UIs and template
    /// loaders that know which preset produced the body should pass
    /// the preset's `display_name` via [`add_named_body`] for
    /// `"Rocky 1"`-style names.
    pub fn add_body(&mut self, mut body: Body) {
        body.sync_physical_properties();
        self.total_mass += body.mass;
        self.names.push(auto_name(DEFAULT_NAME_PREFIX, &self.names));
        self.bodies.push(body);
        self.initial_energy = None;
        self.template_source = None;
    }

    /// Adds a single body while preserving an explicit display name when given.
    pub fn add_named_body(&mut self, named_body: NamedBody) {
        self.add_named_bodies(vec![named_body]);
    }

    /// Add multiple bodies in a single batch.
    ///
    /// More efficient than calling [`add_body`] in a loop: the trail buffer is
    /// reset only once and the energy baseline is invalidated once.
    pub fn add_bodies(&mut self, new_bodies: Vec<Body>) {
        for mut body in new_bodies {
            body.sync_physical_properties();
            self.total_mass += body.mass;
            self.names.push(auto_name(DEFAULT_NAME_PREFIX, &self.names));
            self.bodies.push(body);
        }

        self.initial_energy = None;
        self.template_source = None;
    }

    /// Add multiple bodies in a single batch while preserving explicit names.
    pub fn add_named_bodies(&mut self, new_bodies: Vec<NamedBody>) {
        for mut named_body in new_bodies {
            let body = named_body.body;
            let mut body = body;
            body.sync_physical_properties();
            self.total_mass += body.mass;
            let name = resolved_name(named_body.name.take(), DEFAULT_NAME_PREFIX, &self.names);
            self.names.push(name);
            self.bodies.push(body);
        }

        self.initial_energy = None;
        self.template_source = None;
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

            self.initial_energy = None;
            self.rel_energy_error = None;
            self.abs_energy_error = 0.0;
            self.template_source = None;
        }
    }

    /// Updates a body in-place, recomputing derived physical properties.
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
                self.rel_energy_error = None;
                self.abs_energy_error = 0.0;
            }
        }
    }

    /// Replaces the entire set of bodies in the simulation.
    ///
    /// All previous state is cleared, the trail buffer is reset, and the
    /// system is normalised to its COM rest frame. Bodies receive
    /// `"Body N"` auto-names; use [`load_named_bodies`](Self::load_named_bodies)
    /// to supply explicit names (e.g. preset display names from a template).
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.load_named_bodies(
            bodies.into_iter().map(|body| NamedBody { body, name: None }).collect(),
        );
    }

    /// Replaces the entire set of bodies in the simulation, preserving any
    /// explicit names attached to each body.
    ///
    /// Same state-reset semantics as [`load_bodies`](Self::load_bodies):
    /// previous bodies, scratch buffers, energy baselines, and integrator
    /// controllers are cleared, and the new system is normalised to its COM
    /// rest frame.
    pub fn load_named_bodies(&mut self, named_bodies: Vec<NamedBody>) {
        self.bodies.clear();
        self.scratch_acc.clear();
        self.names.clear();
        self.total_mass = 0.0;

        for mut named in named_bodies {
            let mut body = named.body;
            body.sync_physical_properties();
            self.total_mass += body.mass;
            let name = resolved_name(named.name.take(), DEFAULT_NAME_PREFIX, &self.names);
            self.names.push(name);
            self.bodies.push(body);
        }

        self.initial_energy = None;
        self.rel_energy_error = None;
        self.abs_energy_error = 0.0;
        self.steps = 0;
        self.t = 0.0;
        self.last_potential = 0.0;
        self.last_kinetic = 0.0;
        self.diagnostics = DiagnosticsComputer::new();
        self.last_diag = SimulationDiagnostics::default();
        self.last_step_degraded = false;
        self.dt_ctrl.reset();
        self.theta_ctrl.set(self.force_model.theta());

        self.zero_com_velocity();
        self.recenter_com();

        self.r_min = compute_closeness(&self.bodies);
        self.template_source = None;
    }

    /// Removes the centre-of-mass velocity so the system is in its rest frame.
    pub fn zero_com_velocity(&mut self) {
        calibration::zero_com_velocity(&mut self.bodies, self.total_mass);
    }

    /// Recenters the system so that the centre of mass is at the origin.
    ///
    /// The translation is routed through the active integrator's
    /// [`recenter_bodies`](crate::physics::integrator::Integrator::recenter_bodies)
    /// hook so any per-body compensation buffers (notably IAS15's `csx`)
    /// stay consistent with the post-shift positions. The fixed-step
    /// integrators have no such buffers and inherit the trait default
    /// (bare subtraction).
    pub fn recenter_com(&mut self) {
        if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
            self.integrator.recenter_bodies(&mut self.bodies, dx, dy);
            // Notify the render-side TrailRecorder of the shift so it can
            // keep stored trail positions aligned with the new body coordinates.
            self.pending_com_shift.0 += -dx as f32;
            self.pending_com_shift.1 += -dy as f32;
        }
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
}
