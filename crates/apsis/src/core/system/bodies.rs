//! Body management: add, remove, update, load, naming, COM calibration.

use crate::core::calibration;
use crate::core::diagnostics::{DiagnosticsComputer, SimulationDiagnostics};
use crate::core::system::System;
use crate::core::system::helpers::{
    auto_name, compute_closeness, l_sun, mass_to_solar, radius_to_solar, resolved_name,
};
use crate::domain::body::{Body, NamedBody};

impl System {
    /// Adds a new body to the simulation.
    ///
    /// The trail buffer is reset to accommodate the new body count; trail
    /// history is lost.  Energy baseline is reset because the system
    /// topology has changed. Also invalidates any template-rebuild source
    /// remembered by [`from_template`](System::from_template): a later
    /// [`with_seed`](System::with_seed) will no longer overwrite this
    /// manual addition.
    pub fn add_body(&mut self, mut body: Body) {
        use crate::domain::body::default_softening;
        body.sync_physical_properties();
        if (self.softening_scale - 1.0).abs() > 1e-15 {
            body.softening = default_softening(body.mass) * self.softening_scale;
        }
        self.total_mass += body.mass;
        self.names.push(auto_name(body.material, &self.names));
        body.update_luminosity(mass_to_solar(), radius_to_solar(), l_sun());
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
        use crate::domain::body::default_softening;
        for mut body in new_bodies {
            body.sync_physical_properties();
            if (self.softening_scale - 1.0).abs() > 1e-15 {
                body.softening = default_softening(body.mass) * self.softening_scale;
            }
            self.total_mass += body.mass;
            self.names.push(auto_name(body.material, &self.names));
            body.update_luminosity(mass_to_solar(), radius_to_solar(), l_sun());
            self.bodies.push(body);
        }

        self.initial_energy = None;
        self.template_source = None;
    }

    /// Add multiple bodies in a single batch while preserving explicit names.
    pub fn add_named_bodies(&mut self, new_bodies: Vec<NamedBody>) {
        use crate::domain::body::default_softening;
        for mut named_body in new_bodies {
            let mut body = named_body.body;
            body.sync_physical_properties();
            if (self.softening_scale - 1.0).abs() > 1e-15 {
                body.softening = default_softening(body.mass) * self.softening_scale;
            }
            self.total_mass += body.mass;
            let name = resolved_name(named_body.name.take(), body.material, &self.names);
            body.update_luminosity(mass_to_solar(), radius_to_solar(), l_sun());
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
            self.rel_energy_error = 0.0;
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

            body.update_luminosity(mass_to_solar(), radius_to_solar(), l_sun());
            *slot = body;

            if mass_changed {
                self.initial_energy = None;
                self.rel_energy_error = 0.0;
            }
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
            b.update_luminosity(mass_to_solar(), radius_to_solar(), l_sun());
            self.total_mass += b.mass;
            self.names.push(auto_name(b.material, &self.names));
            self.bodies.push(b);
        }

        self.initial_energy = None;
        self.rel_energy_error = 0.0;
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

        let (r_min, softening_max) = compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }

    /// Removes the centre-of-mass velocity so the system is in its rest frame.
    pub fn zero_com_velocity(&mut self) {
        calibration::zero_com_velocity(&mut self.bodies, self.total_mass);
    }

    /// Recenters the system so that the centre of mass is at the origin.
    pub fn recenter_com(&mut self) {
        if let Some((dx, dy)) = calibration::com_offset(&self.bodies, self.total_mass) {
            calibration::apply_body_shift(&mut self.bodies, dx, dy);
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
