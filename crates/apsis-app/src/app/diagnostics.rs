//! Render-loop frame diagnostics (FPS, frame time).
//!
//! `tick(animating)` must be called once per frame; the `animating`
//! flag freezes the EMA whenever the loop is paced by egui's idle
//! throttle instead of by real workload, so the readout never reports
//! the ~20 Hz throttle as the application's actual throughput.

use std::time::Instant;

const EMA_TAU_S: f64 = 1.0;
const WARMUP_S: f64 = 1.0;
const DT_MIN_S: f64 = 0.001;
const DT_MAX_S: f64 = 0.250;

#[derive(Debug, Default)]
pub struct Diagnostics {
    last_instant: Option<Instant>,
    start_instant: Option<Instant>,
    ema_dt: f64,
    idle: bool,
}

impl Diagnostics {
    pub fn new() -> Self {
        Self::default()
    }

    /// `animating` should be `true` only when the application is
    /// driving continuous repaints (sim running, camera animating,
    /// drag in progress). Resuming from idle restarts warmup so the
    /// first post-resume sample does not poison the EMA.
    pub fn tick(&mut self, animating: bool) {
        let now = Instant::now();
        let prev = self.last_instant.replace(now);

        if !animating {
            self.idle = true;
            return;
        }

        if self.idle {
            self.idle = false;
            self.start_instant = Some(now);
            self.ema_dt = 0.0;
            return;
        }

        let prev = match prev {
            Some(t) => t,
            None => {
                self.start_instant = Some(now);
                return;
            },
        };

        let dt = (now - prev).as_secs_f64().clamp(DT_MIN_S, DT_MAX_S);

        if self.ema_dt == 0.0 {
            self.ema_dt = dt;
        } else {
            let alpha = 1.0 - (-dt / EMA_TAU_S).exp();
            self.ema_dt += alpha * (dt - self.ema_dt);
        }
    }

    pub fn is_idle(&self) -> bool {
        self.idle
    }

    pub fn warming(&self) -> bool {
        if self.idle {
            return false;
        }
        match self.start_instant {
            Some(t) => t.elapsed().as_secs_f64() < WARMUP_S,
            None => true,
        }
    }

    pub fn fps(&self) -> f64 {
        if self.ema_dt > 0.0 { 1.0 / self.ema_dt } else { 0.0 }
    }

    pub fn frame_ms(&self) -> f64 {
        self.ema_dt * 1000.0
    }
}
