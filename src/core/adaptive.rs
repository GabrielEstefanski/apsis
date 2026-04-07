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

    pub fn update(&mut self, e_theta: f64, dt: f64) -> f64 {
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

#[derive(Debug, Clone)]
pub struct DtController {
    pub config: DtAdaptationConfig,
    last_dt: f64,
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
        }
    }

    pub fn reset(&mut self) {
        self.last_dt = 0.0;
    }

    pub fn update(
        &mut self,
        proposed_dt: f64,
        rel_energy_error: f64,
        stats: AccelerationStats,
    ) -> f64 {
        let cfg = &self.config;

        let clamp = |dt: f64| dt.clamp(cfg.min_dt, cfg.max_dt);

        if !cfg.enabled {
            let dt = clamp(proposed_dt);
            self.last_dt = dt;
            return dt;
        }

        let prev = if self.last_dt > 0.0 {
            self.last_dt
        } else {
            clamp(proposed_dt)
        };

        let mut dt = clamp(proposed_dt).min(prev);

        let ratio = (rel_energy_error.abs() / cfg.target_rel_energy_error).max(1e-12);

        let energy_scale = if ratio > 1.0 {
            1.0 / (1.0 + 0.5 * (ratio - 1.0))
        } else {
            1.0 + 0.2 * (1.0 - ratio)
        };

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
