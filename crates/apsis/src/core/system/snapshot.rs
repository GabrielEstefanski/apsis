//! Save/load via [`SimSnapshot`].

use crate::core::system::System;
use crate::core::system::helpers::{auto_name, compute_closeness};
use crate::domain::body::Body;
use crate::physics::integrator::make_integrator;

impl System {
    /// Capture the minimal state required for deterministic reproduction.
    pub fn to_snapshot(&self) -> crate::io::snapshot::SimSnapshot {
        use crate::io::snapshot::{BodyRecord, SimSnapshot};
        SimSnapshot {
            save_id: 0,
            t: self.t,
            steps: self.steps,
            dt: self.current_dt,
            theta: self.force_model.theta(),
            softening_scale: self.softening_scale,
            g_factor: self.g_factor,
            integrator_kind: self.integrator.kind(),
            // trail_every is now owned by TrailRecorder; persisted at app level.
            trail_every: 1,
            sim_name: String::new(),
            seed: self.seed,
            trail: None,
            bodies: self.bodies.iter().map(BodyRecord::from_body).collect(),
            names: self.names.clone(),
        }
    }

    /// Replace the current simulation state with a saved snapshot.
    ///
    /// The trail buffer is cleared (cosmetic — cannot be restored).
    /// Energy and angular-momentum references are reset so the first
    /// post-load step establishes new baselines.
    pub fn restore_from_snapshot(&mut self, snap: &crate::io::snapshot::SimSnapshot) {
        let bodies: Vec<Body> = snap.bodies.iter().map(|r| r.into_body()).collect();

        use crate::core::system::helpers::DEFAULT_NAME_PREFIX;
        self.names = if snap.names.len() == bodies.len() {
            snap.names.clone()
        } else {
            let mut acc: Vec<String> = Vec::with_capacity(bodies.len());
            for _ in &bodies {
                acc.push(auto_name(DEFAULT_NAME_PREFIX, &acc));
            }
            acc
        };

        let n = bodies.len();
        self.bodies = bodies;
        self.total_mass = self.bodies.iter().map(|b| b.mass).sum();
        self.scratch_acc.clear();
        // Trail buffer is now owned by the render-side TrailRecorder.
        // The physics snapshot carries TrailSnapshot for the app to forward
        // to TrailRecorder.restore_from_snapshot() after this call.
        let _ = n;

        self.t = snap.t;
        self.steps = snap.steps;
        self.current_dt = snap.dt;
        self.user_dt = snap.dt;
        self.force_model.set_theta(snap.theta);
        self.softening_scale = snap.softening_scale;
        self.g_factor = snap.g_factor;
        self.integrator = make_integrator(snap.integrator_kind);
        // trail_every is now managed by TrailRecorder at app level.
        self.seed = snap.seed;
        // Bodies came from disk, not a preset builder — invalidate the
        // template-rebuild path so a later `.with_seed(...)` won't
        // overwrite the restored state.
        self.template_source = None;
        self.initial_energy = None;
        self.initial_angular_momentum = None;
        self.rel_energy_error = 0.0;
        self.rel_angular_momentum_error = 0.0;
        self.abs_angular_momentum_error = 0.0;
        self.last_kinetic = 0.0;
        self.last_potential = 0.0;
        self.diagnostics = crate::core::diagnostics::DiagnosticsComputer::new();
        self.last_diag = crate::core::diagnostics::SimulationDiagnostics::default();
        self.last_step_degraded = false;
        self.orbital_cache.clear();
        self.dt_ctrl.reset();
        self.theta_ctrl.set(snap.theta);

        let (r_min, softening_max) = compute_closeness(&self.bodies);
        self.r_min = r_min;
        self.softening_max = softening_max;
    }
}
