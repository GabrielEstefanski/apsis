//! Body management: add, remove, update, load, naming, COM calibration.

use crate::core::calibration;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::system::System;
use crate::core::system::helpers::{DEFAULT_NAME_PREFIX, compute_closeness, resolved_name};
use crate::domain::body::Body;

impl System {
    /// Adds a new body to the simulation.
    ///
    /// If `body.name` is `None` the system fills it with a stable
    /// `"Body N"` placeholder before pushing; an explicit name on the
    /// incoming `Body` (via [`Body::with_name`]) rides through unchanged.
    /// Trail buffer is reset to accommodate the new body count; trail
    /// history is lost. Energy baseline is reset because the system
    /// topology has changed. Also invalidates any template-rebuild
    /// source remembered by [`from_template`](System::from_template).
    pub fn add_body(&mut self, mut body: Body) {
        body.sync_physical_properties();
        self.total_mass += body.mass;
        let existing = self.existing_names();
        body.name = Some(resolved_name(body.name.take(), DEFAULT_NAME_PREFIX, &existing));
        self.bodies.push(body);
        self.initial_energy = None;
        self.template_source = None;
    }

    /// Add multiple bodies in a single batch.
    ///
    /// More efficient than calling [`add_body`] in a loop: the energy
    /// baseline is invalidated once. Each body's name is resolved
    /// against the running set so auto-numbered placeholders stay
    /// monotonic and explicit names are preserved.
    pub fn add_bodies(&mut self, new_bodies: Vec<Body>) {
        for mut body in new_bodies {
            body.sync_physical_properties();
            self.total_mass += body.mass;
            let existing = self.existing_names();
            body.name = Some(resolved_name(body.name.take(), DEFAULT_NAME_PREFIX, &existing));
            self.bodies.push(body);
        }

        self.initial_energy = None;
        self.template_source = None;
    }

    /// Removes a body from the simulation.
    ///
    /// Uses `swap_remove` for O(1) removal. The trail buffer is reset
    /// because body indices change.
    pub fn remove_body(&mut self, index: usize) {
        if index < self.bodies.len() {
            let removed = self.bodies.swap_remove(index);
            self.total_mass -= removed.mass;

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
    /// All previous state is cleared, the trail buffer is reset, and
    /// the system is normalised to its COM rest frame. Each body's
    /// name is resolved against the running set during the load —
    /// explicit names (e.g. preset display names from a template via
    /// [`Body::with_name`]) ride through; `None` slots get
    /// `"Body N"` placeholders.
    pub fn load_bodies(&mut self, bodies: Vec<Body>) {
        self.bodies.clear();
        self.scratch_acc.clear();
        self.total_mass = 0.0;

        for mut body in bodies {
            body.sync_physical_properties();
            self.total_mass += body.mass;
            let existing = self.existing_names();
            body.name = Some(resolved_name(body.name.take(), DEFAULT_NAME_PREFIX, &existing));
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

    /// All body names, in registration order. Always `Some` because
    /// [`Self::add_body`] auto-fills `body.name` with a `"Body N"`
    /// placeholder when the caller didn't supply one.
    pub fn names(&self) -> Vec<&str> {
        self.bodies.iter().map(|b| b.name.as_deref().unwrap_or("")).collect()
    }

    /// Rename body `idx`. Silently ignores out-of-range indices.
    pub fn set_name(&mut self, idx: usize, name: String) {
        if let Some(slot) = self.bodies.get_mut(idx) {
            slot.name = Some(name);
        }
    }

    /// Helper: snapshot the current set of body names as owned
    /// `String`s. Used by the auto-name policy at insert time so
    /// `resolved_name` can compare against the running set without
    /// borrowing self twice.
    fn existing_names(&self) -> Vec<String> {
        self.bodies.iter().filter_map(|b| b.name.clone()).collect()
    }
}
