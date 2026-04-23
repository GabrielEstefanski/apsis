//! Metrics assembly and recommended-dt computation.

use crate::core::metrics::Metrics;
use crate::core::system::System;
use crate::physics::energy::{angular_momentum_z, center_of_mass_state, kinetic_energy, total_energy};

impl System {
    /// Assemble a [`Metrics`] snapshot of the current simulation state.
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

            integrator_kind: self.integrator.kind(),
            g_factor: self.g_factor,
            theta: self.force_model.theta(),
            force_is_direct: self.force_model.is_deterministic(),
            dt: self.current_dt,
            user_dt: self.user_dt,
            dt_mode: self.dt_mode,
            adaptive_theta: self.adaptive_theta,

            max_acc: self.last_diag.max_acc,
            jerk: self.last_diag.jerk,
            max_vel: self.last_diag.max_vel,
            last_step_degraded: self.last_step_degraded,

            r_min: self.r_min,
            softening_max: self.softening_max,

            recommended_dt: self.recommended_dt(),
            adaptive_stats: self.integrator.adaptive_stats(),
        }
    }

    /// Accelerations computed during the last integration step.
    pub fn last_accelerations(&self) -> &[(f64, f64)] {
        &self.scratch_acc
    }

    /// Lightweight accessor for adaptive-integrator counters. Avoids
    /// the full [`Metrics`] assembly when the caller only needs the
    /// adaptive state — useful inside the physics thread's hot loop
    /// (e.g. Precision Run telemetry updates) where rebuilding every
    /// observable every tick is wasteful.
    pub fn adaptive_stats(
        &self,
    ) -> Option<crate::physics::integrator::traits::AdaptiveStats> {
        self.integrator.adaptive_stats()
    }

    /// Current relative energy error `δE/E₀` (signed). The same value
    /// is available via [`Metrics::rel_energy_error`] but this
    /// accessor avoids the allocation-cost of a full metrics build.
    pub fn rel_energy_error(&self) -> f64 {
        self.rel_energy_error
    }

    // ── Direct scalar queries ────────────────────────────────────────────────
    //
    // Each returns the single value most scripts reach for, computed or
    // cached with no DTO allocation. Use [`metrics`](Self::metrics) when you
    // need everything at once (e.g. periodic CSV dumps).

    /// Total energy `E = K + U` at the last completed step.
    pub fn energy(&self) -> f64 {
        total_energy(self.last_kinetic, self.last_potential)
    }

    /// Relative energy drift `δE = (E − E₀) / |E₀|` at the last step.
    ///
    /// Alias for [`rel_energy_error`](Self::rel_energy_error); named to
    /// match `energy()` for script ergonomics (`sys.energy()` /
    /// `sys.energy_delta()`).
    #[inline]
    pub fn energy_delta(&self) -> f64 {
        self.rel_energy_error
    }

    /// Kinetic energy `K = Σ ½ mᵢ vᵢ²` at the last completed step.
    #[inline]
    pub fn kinetic_energy(&self) -> f64 {
        self.last_kinetic
    }

    /// Potential energy `U = Σᵢ Σⱼ<ᵢ −G mᵢ mⱼ / |rᵢ−rⱼ|` at the last step.
    #[inline]
    pub fn potential_energy(&self) -> f64 {
        self.last_potential
    }

    /// Total angular momentum (z-component) at the current body state.
    pub fn lz(&self) -> f64 {
        angular_momentum_z(&self.bodies)
    }

    /// Relative angular-momentum drift `δLz = (Lz − Lz₀) / |Lz₀|`.
    #[inline]
    pub fn lz_delta(&self) -> f64 {
        self.rel_angular_momentum_error
    }

    /// Centre-of-mass state `(x, y, vx, vy)` at the current body state.
    pub fn center_of_mass(&self) -> (f64, f64, f64, f64) {
        center_of_mass_state(&self.bodies)
    }

    /// Physics-justified recommended timestep from the current system state.
    ///
    /// Uses two criteria (Power et al. 2003 acceleration + Aarseth jerk).
    /// Returns `None` before the first force evaluation or when N = 0.
    pub fn recommended_dt(&self) -> Option<f64> {
        if self.bodies.is_empty() || self.last_diag.max_acc <= 1e-30 {
            return None;
        }

        let eps_min = self.bodies.iter().map(|b| b.softening).fold(f64::MAX, f64::min);
        if eps_min >= f64::MAX || eps_min <= 0.0 {
            return None;
        }

        const ETA: f64 = 0.05;

        let dt_acc = ETA * (eps_min / self.last_diag.max_acc).sqrt();

        let dt_jerk = if self.last_diag.jerk > 1e-30 {
            ETA * (self.last_diag.max_acc / self.last_diag.jerk).sqrt()
        } else {
            f64::MAX
        };

        Some(dt_acc.min(dt_jerk).clamp(1e-9, 1e6))
    }
}
