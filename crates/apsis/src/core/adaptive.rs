/// Timestep management policy for the N-body integrator.
///
/// # Symplecticity and why this distinction matters
///
/// Symplectic integrators (Velocity Verlet, Yoshida 4th-order, Wisdom–Holman)
/// are derived from a generating function that exactly preserves the symplectic
/// 2-form on phase space. This guarantees that the numerical flow is conjugate to
/// the exact flow of a *modified Hamiltonian* H̃, which differs from the true H
/// by terms of order O(dtᵖ⁺¹) (p = integrator order). As a consequence, energy
/// error is **bounded and oscillatory** for all time — it cannot drift secularly.
///
/// This guarantee holds **only when dt is constant**. When dt varies between
/// steps, the per-step maps Φ_{dt₁}, Φ_{dt₂}, … are each individually symplectic
/// but their composition Φ_{dt₁} ∘ Φ_{dt₂} ∘ … is not a single symplectic map.
/// The modified Hamiltonian changes with each step; there is no single conserved
/// quantity bounding the error. Energy error may then drift monotonically, at a
/// rate proportional to the magnitude and frequency of dt changes.
///
/// For a result intended for publication or long-term stability analysis this is
/// not an acceptable trade-off. The only rigorous options are:
///
/// 1. **Fixed dt** — the standard approach; choose dt small enough for the system.
/// 2. **Time-transformed Hamiltonians** (Mikkola & Tanikawa 1999; Preto & Tremaine
///    1999) — reformulate in a fictitious time variable so that dt_fictitious is
///    constant while the physical step size adapts; still symplectic. Not
///    implemented here.
///
/// `Adaptive` mode is provided as a **convenience for exploration and
/// dt tuning only** — finding a stable dt, qualitative demonstrations.
/// It must not be used for any run whose results are cited.
///
/// # Selection guide
///
/// | Mode | Energy error | Suitable for |
/// |------|-------------|--------------|
/// | [`Fixed`] | Bounded, oscillatory — O(dtᵖ) amplitude | Long-term integration, publication |
/// | [`Adaptive`] | Potentially secular drift | Scenario exploration, dt tuning |
///
/// # References
///
/// - Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical
///   Integration: Structure-Preserving Algorithms for Ordinary Differential
///   Equations* (2nd ed.). Springer. §VI.
/// - Leimkuhler, B., & Reich, S. (2004). *Simulating Hamiltonian Dynamics*.
///   Cambridge University Press. §4.2.
/// - Yoshida, H. (1990). Construction of higher order symplectic integrators.
///   *Phys. Lett. A* 150, 262–268.
/// - Mikkola, S., & Tanikawa, K. (1999). Explicit symplectic algorithms for
///   time-transformed Hamiltonians. *Cel. Mech. Dyn. Astron.* 74, 287–295.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DtMode {
    /// **Constant timestep.** The integration is symplectic — energy error is
    /// bounded and oscillatory for all time. The magnitude of oscillation scales
    /// as O(dtᵖ), where p is the integrator order (2 for Velocity Verlet and
    /// Wisdom–Holman; 4 for Yoshida 4th-order).
    ///
    /// This is the only mode compatible with publication-quality results and
    /// long-term stability studies. The [`DtController`] is not consulted when
    /// this mode is active; `current_dt` equals `user_dt` exactly.
    Fixed,

    /// **Adaptive timestep.** [`DtController`] modulates `dt` each step based
    /// on the relative energy error and an acceleration-based CFL criterion.
    ///
    /// # ⚠ Symplecticity is broken in this mode
    ///
    /// Varying dt breaks the symplectic structure of the integrator. Energy
    /// error is no longer guaranteed to remain bounded; it may drift secularly.
    /// The severity depends on the rate and magnitude of dt changes, but even
    /// small variations accumulate over long integrations.
    ///
    /// **Use only for scenario exploration and initial condition setup.**
    /// Switch to [`Fixed`] before any run whose output will be analysed or cited.
    Adaptive,
}

impl DtMode {
    /// Stable short label for serialisation / display. Used by
    /// `apsis::records` to write `integrator.dt_mode` in the record
    /// header without coupling to `Debug` formatting.
    pub const fn label(self) -> &'static str {
        match self {
            Self::Fixed => "Fixed",
            Self::Adaptive => "Adaptive",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct AccelerationStats {
    pub max_acc: f64,
    pub jerk: f64,
}

impl AccelerationStats {
    pub fn new(max_acc: f64, jerk: f64) -> Self {
        Self { max_acc, jerk }
    }
}

#[derive(Debug, Clone)]
pub struct DtAdaptationConfig {
    pub enabled: bool,
    pub min_dt: f64,
    pub max_dt: f64,
    pub target_rel_energy_error: f64,
    pub accel_epsilon: f64,
    pub grow_limit: f64,
    pub shrink_limit: f64,
    pub dt_slew_fraction: f64,
}

#[derive(Debug, Clone)]
pub struct ThetaTuning {
    pub ema_time_constant: f64,
    pub response_gain: f64,
    pub tighten_gain: f64,
    pub relax_gain: f64,
    pub min_step_fraction: f64,
    pub max_step_fraction: f64,
}

#[derive(Debug, Clone)]
pub struct ThetaController {
    pub target_error: f64,
    pub min_theta: f64,
    pub max_theta: f64,
    theta: f64,
    error_ema: f64,
    tuning: ThetaTuning,
}

impl ThetaController {
    pub fn new(target_error: f64, min_theta: f64, max_theta: f64) -> Self {
        let min_theta = min_theta.max(1e-6);

        Self {
            target_error: target_error.max(1e-12),
            min_theta,
            max_theta: max_theta.max(min_theta),
            theta: min_theta,
            error_ema: 1.0,
            tuning: ThetaTuning {
                ema_time_constant: 5.0,
                response_gain: 0.5,
                tighten_gain: 0.6,
                relax_gain: 0.25,
                min_step_fraction: 0.01,
                max_step_fraction: 0.08,
            },
        }
    }

    pub fn with_initial_theta(mut self, theta: f64) -> Self {
        self.theta = theta.clamp(self.min_theta, self.max_theta);
        self
    }

    pub fn current(&self) -> f64 {
        self.theta
    }

    pub fn set(&mut self, theta: f64) {
        self.theta = theta.clamp(self.min_theta, self.max_theta);
    }

    pub fn update(&mut self, e_theta: f64, _dt: f64) -> f64 {
        let current_theta = self.theta;

        let e_norm = e_theta / self.target_error;

        let alpha = (1.0 / (self.tuning.ema_time_constant + 1.0)).clamp(0.05, 0.4);

        self.error_ema = alpha * e_norm + (1.0 - alpha) * self.error_ema;
        let e = self.error_ema;

        let factor = if e > 1.0 {
            1.0 / (1.0 + self.tuning.tighten_gain * (e - 1.0))
        } else {
            1.0 + self.tuning.relax_gain * (1.0 - e)
        };

        let desired = (current_theta * factor).clamp(self.min_theta, self.max_theta);

        let responsiveness = (e - 1.0).abs().clamp(0.0, 1.0);
        let response = self.tuning.response_gain * (0.3 + 0.7 * responsiveness);

        let blended = current_theta + response * (desired - current_theta);

        let step_scale = (e - 1.0).abs().clamp(0.1, 2.0);

        let max_step = current_theta.abs().max(self.min_theta)
            * (self.tuning.min_step_fraction
                + (self.tuning.max_step_fraction - self.tuning.min_step_fraction) * step_scale);

        let delta = (blended - current_theta).clamp(-max_step, max_step);

        let next = (current_theta + delta).clamp(self.min_theta, self.max_theta);

        self.theta = next;

        next
    }

    pub fn error(&self) -> f64 {
        self.error_ema
    }
}

/// Whether the dt controller has a well-conditioned relative-error
/// signal to feed back on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeedbackMode {
    /// Relative error is defined; controller adjusts dt against
    /// `target_rel_energy_error`.
    Active,
    /// `|E_initial|` is below the conditioning threshold; controller
    /// disables feedback and returns `user_dt` unchanged.
    DisabledPrecisionLimited,
}

#[derive(Debug, Clone)]
pub struct DtController {
    pub config: DtAdaptationConfig,
    last_dt: f64,
    feedback_mode: FeedbackMode,
}

impl DtController {
    pub fn new(config: DtAdaptationConfig) -> Self {
        Self {
            config: DtAdaptationConfig {
                enabled: config.enabled,
                min_dt: config.min_dt.max(1e-12),
                max_dt: config.max_dt.max(config.min_dt.max(1e-12)),
                target_rel_energy_error: config.target_rel_energy_error.max(1e-12),
                accel_epsilon: config.accel_epsilon.max(1e-18),
                grow_limit: config.grow_limit.max(1.0),
                shrink_limit: config.shrink_limit.clamp(1e-6, 1.0),
                dt_slew_fraction: config.dt_slew_fraction.clamp(0.02, 1.0),
            },
            last_dt: 0.0,
            feedback_mode: FeedbackMode::Active,
        }
    }

    pub fn reset(&mut self) {
        self.last_dt = 0.0;
        self.feedback_mode = FeedbackMode::Active;
    }

    pub fn feedback_mode(&self) -> FeedbackMode {
        self.feedback_mode
    }

    pub fn update(
        &mut self,
        proposed_dt: f64,
        rel_energy_error: Option<f64>,
        stats: AccelerationStats,
    ) -> f64 {
        let cfg = &self.config;

        let clamp = |dt: f64| dt.clamp(cfg.min_dt, cfg.max_dt);

        if !cfg.enabled {
            self.feedback_mode = FeedbackMode::Active;
            let dt = clamp(proposed_dt);
            self.last_dt = dt;
            return dt;
        }

        let Some(rel) = rel_energy_error else {
            self.feedback_mode = FeedbackMode::DisabledPrecisionLimited;
            let dt = clamp(proposed_dt);
            self.last_dt = dt;
            return dt;
        };

        self.feedback_mode = FeedbackMode::Active;

        let prev = if self.last_dt > 0.0 { self.last_dt } else { clamp(proposed_dt) };

        let mut dt = clamp(proposed_dt).min(prev);

        let ratio = (rel.abs() / cfg.target_rel_energy_error).max(1e-12);

        let energy_scale =
            if ratio > 1.0 { 1.0 / (1.0 + 0.5 * (ratio - 1.0)) } else { 1.0 + 0.2 * (1.0 - ratio) };

        dt *= energy_scale;

        let effective_acc = stats.max_acc + stats.jerk * dt;

        let accel_dt = if effective_acc > 0.0 {
            (cfg.accel_epsilon / effective_acc).sqrt()
        } else {
            cfg.max_dt
        };

        let candidate = dt.min(accel_dt);
        let lo = prev * cfg.shrink_limit;
        let hi = prev * cfg.grow_limit;
        let smoothed = candidate.clamp(lo.max(cfg.min_dt), hi.min(cfg.max_dt));

        let out = clamp(smoothed);

        self.last_dt = out;

        out
    }

    pub fn last_dt(&self) -> f64 {
        self.last_dt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> DtAdaptationConfig {
        DtAdaptationConfig {
            enabled: true,
            min_dt: 1e-9,
            max_dt: 1e6,
            target_rel_energy_error: 1e-6,
            accel_epsilon: 0.1,
            grow_limit: 1.2,
            shrink_limit: 0.5,
            dt_slew_fraction: 0.1,
        }
    }

    fn stats() -> AccelerationStats {
        AccelerationStats { max_acc: 1.0, jerk: 0.0 }
    }

    #[test]
    fn feedback_mode_active_when_rel_error_is_some() {
        let mut ctrl = DtController::new(config());
        ctrl.update(1e-3, Some(1e-9), stats());
        assert_eq!(ctrl.feedback_mode(), FeedbackMode::Active);
    }

    #[test]
    fn feedback_mode_disabled_when_rel_error_is_none() {
        let mut ctrl = DtController::new(config());
        let out = ctrl.update(1e-3, None, stats());
        assert_eq!(ctrl.feedback_mode(), FeedbackMode::DisabledPrecisionLimited);
        assert_eq!(out, 1e-3, "proposed_dt returned unchanged in precision-limited regime");
    }

    #[test]
    fn feedback_mode_transitions_when_signal_recovers() {
        let mut ctrl = DtController::new(config());
        ctrl.update(1e-3, None, stats());
        assert_eq!(ctrl.feedback_mode(), FeedbackMode::DisabledPrecisionLimited);
        ctrl.update(1e-3, Some(1e-9), stats());
        assert_eq!(ctrl.feedback_mode(), FeedbackMode::Active);
    }
}
