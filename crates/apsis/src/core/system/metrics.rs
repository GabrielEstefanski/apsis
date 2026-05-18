//! Metrics assembly and recommended-dt computation.

use crate::core::metrics::Metrics;
use crate::core::system::System;
use crate::math::Vec3;
use crate::physics::energy::{angular_momentum_z, center_of_mass_state, total_energy};

impl System {
    /// Assemble a [`Metrics`] snapshot of the current simulation state.
    pub fn metrics(&self) -> Metrics {
        let kinetic = self.last_kinetic;
        let potential = self.last_potential;
        let total = total_energy(kinetic, potential);

        let lz = angular_momentum_z(&self.bodies);
        let (com_pos, com_vel) = center_of_mass_state(&self.bodies);

        Metrics {
            kinetic,
            potential,
            total_energy: total,
            initial_energy: self.initial_energy.unwrap_or(total),
            abs_energy_error: self.abs_energy_error,
            rel_energy_error: self.rel_energy_error,

            angular_momentum_z: lz,
            initial_angular_momentum_z: self.initial_angular_momentum.unwrap_or(lz),
            rel_angular_momentum_error: self.rel_angular_momentum_error,
            abs_angular_momentum_error: self.abs_angular_momentum_error,

            com_x: com_pos.x,
            com_y: com_pos.y,
            com_z: com_pos.z,
            com_vx: com_vel.x,
            com_vy: com_vel.y,
            com_vz: com_vel.z,

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
            kernel_epsilon_squared: self.force_model.kernel().epsilon_squared(),

            recommended_dt: self.recommended_dt(),
            adaptive_stats: self.integrator.adaptive_stats(),
        }
    }

    /// Accelerations computed during the last integration step.
    pub fn last_accelerations(&self) -> &[Vec3] {
        &self.scratch_acc
    }

    /// Lightweight accessor for adaptive-integrator counters. Avoids
    /// the full [`Metrics`] assembly when the caller only needs the
    /// adaptive state — useful inside the physics thread's hot loop
    /// (e.g. Precision Run telemetry updates) where rebuilding every
    /// observable every tick is wasteful.
    pub fn adaptive_stats(&self) -> Option<crate::physics::integrator::traits::AdaptiveStats> {
        self.integrator.adaptive_stats()
    }

    /// Relative energy drift `(E − E₀) / |E₀|`, or `None` when
    /// `|E₀| < MIN_RELATIVE_DENOMINATOR` (precision-limited regime).
    pub fn rel_energy_error(&self) -> Option<f64> {
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

    /// Initial total energy. Falls back to the current energy if the
    /// first force evaluation has not run yet.
    #[inline]
    pub fn initial_energy(&self) -> f64 {
        self.initial_energy.unwrap_or_else(|| self.energy())
    }

    /// Absolute energy drift `E − E₀` (signed).
    #[inline]
    pub fn abs_energy_drift(&self) -> f64 {
        self.abs_energy_error
    }

    /// Relative energy drift `(E − E₀) / |E₀|`, or `None` when
    /// `|E₀| < MIN_RELATIVE_DENOMINATOR`. Alias for
    /// [`rel_energy_error`](Self::rel_energy_error), named to match
    /// `energy()`.
    #[inline]
    pub fn energy_delta(&self) -> Option<f64> {
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

    /// Initial Lz. Falls back to the current Lz if the first state
    /// evaluation has not run yet.
    #[inline]
    pub fn initial_lz(&self) -> f64 {
        self.initial_angular_momentum.unwrap_or_else(|| self.lz())
    }

    /// Absolute Lz drift `|Lz − Lz₀|`.
    #[inline]
    pub fn abs_lz_drift(&self) -> f64 {
        self.abs_angular_momentum_error
    }

    /// Relative Lz drift `(Lz − Lz₀) / |Lz₀|`, or `None` when
    /// `|Lz₀| < MIN_RELATIVE_DENOMINATOR`.
    #[inline]
    pub fn lz_delta(&self) -> Option<f64> {
        self.rel_angular_momentum_error
    }

    /// Planar projection of the centre-of-mass state: `(x, y, vx, vy)`.
    ///
    /// The `xy`-projection of [`center_of_mass_3d`](Self::center_of_mass_3d).
    /// Useful for 2D plots and overlay code; callers operating in
    /// three dimensions read the full state through `_3d` directly.
    pub fn center_of_mass(&self) -> (f64, f64, f64, f64) {
        let (pos, vel) = center_of_mass_state(&self.bodies);
        (pos.x, pos.y, vel.x, vel.y)
    }

    /// Centre-of-mass position and velocity in the inertial frame.
    pub fn center_of_mass_3d(&self) -> (Vec3, Vec3) {
        center_of_mass_state(&self.bodies)
    }

    /// Physics-justified recommended timestep from the current system state.
    ///
    /// Two regimes, selected by whether the active gravity kernel uses
    /// softening:
    ///
    /// - **Softened (η = 0.05):** Power-style acceleration criterion
    ///   `η · √(ε / a_max)` (Power et al. 2003) combined with Aarseth's
    ///   jerk criterion `η · √(a_max / |jerk|)` (Aarseth 2003 §2). Uses
    ///   the kernel's softening length as the resolution scale.
    ///
    /// - **Unsoftened fallback (η = 0.01):** shortest pairwise Kepler
    ///   period `T_ij = 2π · √(r_ij³ / μ_ij)` (Aarseth 2003 §2). Closed-
    ///   form, scale-aware, naturally tightens during close encounters.
    ///
    /// Returns `None` before the first force evaluation, when N = 0, or
    /// when every pair degenerates.
    pub fn recommended_dt(&self) -> Option<f64> {
        if self.bodies.is_empty() || self.last_diag.max_acc <= 1e-30 {
            return None;
        }

        let eps_sq = self.force_model.kernel().epsilon_squared();
        let eps = eps_sq.sqrt();

        if eps > 0.0 {
            const ETA: f64 = 0.05;

            let dt_acc = ETA * (eps / self.last_diag.max_acc).sqrt();

            let dt_jerk = if self.last_diag.jerk > 1e-30 {
                ETA * (self.last_diag.max_acc / self.last_diag.jerk).sqrt()
            } else {
                f64::MAX
            };

            return Some(dt_acc.min(dt_jerk).clamp(1e-9, 1e6));
        }

        // Unsoftened fallback: shortest pairwise Kepler period.
        if self.bodies.len() < 2 {
            return None;
        }
        const ETA_PAIR: f64 = 0.01;
        let g = self.g_factor;
        let mut min_period = f64::MAX;
        for i in 0..self.bodies.len() {
            let bi = &self.bodies[i];
            for j in (i + 1)..self.bodies.len() {
                let bj = &self.bodies[j];
                let dx = bi.pos_x - bj.pos_x;
                let dy = bi.pos_y - bj.pos_y;
                let dz = bi.pos_z - bj.pos_z;
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < 1e-60 {
                    continue;
                }
                let mu = g * (bi.mass + bj.mass);
                if mu < 1e-30 {
                    continue;
                }
                // T = 2π · √(r³ / μ); compute as 2π · √(r² · r / μ).
                let r = r2.sqrt();
                let t_pair = std::f64::consts::TAU * (r2 * r / mu).sqrt();
                if t_pair < min_period {
                    min_period = t_pair;
                }
            }
        }
        if min_period >= f64::MAX || !min_period.is_finite() {
            return None;
        }
        Some((ETA_PAIR * min_period).clamp(1e-9, 1e6))
    }
}

#[cfg(test)]
mod tests {
    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::physics::integrator::IntegratorKind;
    use crate::units::UnitSystem;

    /// Reconstruct the expected `recommended_dt` from the current state of a
    /// fully-unsoftened system: η · min over all pairs of `2π · √(r³/μ)` with
    /// `μ = G · (m_i + m_j)`. Computed from `sys.bodies()` directly so the
    /// assertion is on the formula, not on a static a-priori IC value (the
    /// integrator drifts r slightly during the warm-up step that populates
    /// `last_diag`).
    fn expected_dt_unsoftened(sys: &System) -> f64 {
        const ETA_PAIR: f64 = 0.01;
        let bodies = sys.bodies();
        let g = sys.metrics().g_factor;
        let mut min_period = f64::MAX;
        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                let dx = bodies[i].pos_x - bodies[j].pos_x;
                let dy = bodies[i].pos_y - bodies[j].pos_y;
                let dz = bodies[i].pos_z - bodies[j].pos_z;
                let r2 = dx * dx + dy * dy + dz * dz;
                let r = r2.sqrt();
                let mu = g * (bodies[i].mass + bodies[j].mass);
                let t_pair = std::f64::consts::TAU * (r2 * r / mu).sqrt();
                if t_pair < min_period {
                    min_period = t_pair;
                }
            }
        }
        ETA_PAIR * min_period
    }

    /// Two equal-mass unsoftened bodies. The pair-Kepler fallback returns
    /// `η · 2π · √(r³/μ)` evaluated at the current geometry. Asserts the
    /// formula by recomputing the expected value from `sys.bodies()` after
    /// the warm-up step — verifies that the implementation iterates pairs,
    /// uses `r³/μ` correctly, and applies `η = 0.01`.
    #[test]
    fn unsoftened_two_body_returns_pair_kepler_period() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_integrator(IntegratorKind::VelocityVerlet)
            .with_dt(0.01);
        sys.step();

        let dt = sys.recommended_dt().expect("unsoftened two-body should yield Some");
        let expected = expected_dt_unsoftened(&sys);
        let rel_err = (dt - expected).abs() / expected;
        assert!(
            rel_err < 1e-14,
            "two-body pair-Kepler formula mismatch: expected {expected:.6e}, got {dt:.6e}",
        );
    }

    /// Three unsoftened bodies at distinct separations. The closest pair
    /// (separation ≈ 1, reduced mass 2) must dominate; the spectator at
    /// `x = 10` produces longer-period pairs that must not be selected.
    #[test]
    fn unsoftened_three_body_closest_pair_dominates() {
        let bodies = vec![
            Body::rocky(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.5, 0.5),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 0.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_integrator(IntegratorKind::VelocityVerlet)
            .with_dt(0.01);
        sys.step();

        let dt = sys.recommended_dt().expect("unsoftened three-body should yield Some");
        let expected = expected_dt_unsoftened(&sys);
        let rel_err = (dt - expected).abs() / expected;
        assert!(
            rel_err < 1e-14,
            "three-body closest-pair selection mismatch: expected {expected:.6e}, got {dt:.6e}",
        );

        // Sanity: the dt selected must correspond to the closest pair —
        // not the {0, body_at_10} or {1, body_at_10} pairs whose periods
        // are an order of magnitude longer.
        let bodies_now = sys.bodies();
        let r_close = (bodies_now[1].pos_x - bodies_now[0].pos_x)
            .hypot(bodies_now[1].pos_y - bodies_now[0].pos_y);
        assert!(r_close < 2.0, "closest pair must remain close after one step");
    }

    /// Default Newton kernel routes through the pair-Kepler fallback.
    #[test]
    fn newton_kernel_falls_back_to_pair_kepler() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_integrator(IntegratorKind::VelocityVerlet)
            .with_dt(0.01);
        sys.step();

        let dt = sys.recommended_dt().expect("Newton kernel should yield Some via pair-Kepler");
        let expected = expected_dt_unsoftened(&sys);
        let rel_err = (dt - expected).abs() / expected;
        assert!(
            rel_err < 1e-14,
            "Newton kernel did not route through pair-Kepler fallback: \
             expected {expected:.6e}, got {dt:.6e}",
        );
    }
}
