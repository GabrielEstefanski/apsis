//! Observable snapshot of an in-flight precision run.
//!
//! Populated by the physics thread on a regular cadence (typically
//! after each accepted sub-step or every wall-second, whichever is
//! less frequent); read by the UI to render the Precision Run panel.
//!
//! # Why a snapshot, not a stream
//!
//! UI consumers do not need every intermediate value — they render
//! at ~60 Hz and would discard 99% of a high-frequency feed. A
//! snapshot that reflects "the state at time T" matches the render
//! cadence with zero synchronisation overhead beyond the single
//! lock on [`PrecisionRunController`](super::PrecisionRunController).
//!
//! # Design of the throughput window
//!
//! Raw instantaneous throughput (substeps / wall-second-of-last-batch)
//! is noisy in adaptive integrators — one close encounter can tank
//! the rate for a substep, then recover. The fields below record a
//! **moving-window** view: [`substeps_per_second_window`] and
//! [`sim_time_per_second_window`] are rolling averages kept by the
//! physics thread, so the UI can render a stable number without
//! re-implementing the smoothing itself.

use std::time::Duration;

/// One-shot snapshot of a precision run's observable state.
///
/// All fields are plain owned values so the snapshot can be cloned
/// cheaply from the physics thread into the controller without
/// holding locks.
#[derive(Debug, Clone)]
pub struct Telemetry {
    // ── Controller counters (cumulative) ──────────────────────────────────────
    /// Accepted sub-steps in this run (zeroed at
    /// [`start`](super::PrecisionRunController::start)).
    pub substeps: u64,
    /// Rejections caused by Picard non-convergence in this run.
    pub rejections_picard: u64,
    /// Rejections caused by truncation-error breach in this run.
    pub rejections_truncation: u64,
    /// Picard iterations totalled across all attempts (accepted +
    /// rejected) in this run.
    pub picard_iters: u64,
    /// Degraded accepts in this run (any cause — floor or deadline).
    pub degraded: u64,
    /// Subset of [`degraded`](Self::degraded) attributable to the
    /// `DT_MIN` floor specifically (scenario-stiffness signal).
    pub floor_hits: u64,

    // ── Current step size ─────────────────────────────────────────────────────
    /// Size of the most recent accepted sub-step.
    pub current_dt: f64,
    /// 50th-percentile sub-step size in the rolling window.
    pub dt_p50: f64,
    /// 95th-percentile sub-step size in the rolling window.
    pub dt_p95: f64,

    // ── Energy diagnostics ────────────────────────────────────────────────────
    /// Peak absolute relative energy error since run start
    /// (`max_t |δE/E₀|`).
    pub peak_energy_err: f64,
    /// Current relative energy error (`δE/E₀`), signed.
    pub current_energy_err: f64,

    // ── Throughput (rolling-window) ───────────────────────────────────────────
    /// Accepted sub-steps per wall-second, smoothed over a moving
    /// window (typically 500 ms to 1 s). Zero until the first
    /// window has elapsed.
    pub substeps_per_second_window: f64,
    /// Simulation-time units advanced per wall-second, smoothed.
    /// Drives the primary progress label ("1.4 yr/s"). Zero until
    /// the first window has elapsed.
    pub sim_time_per_second_window: f64,

    // ── Progress cache ────────────────────────────────────────────────────────
    /// Latest `(t_sim - t_start) / (t_target - t_start)` — cached
    /// so the UI can render a frozen bar after a run ends with
    /// partial progress (abort). Range `[0.0, 1.0]`.
    pub last_progress_fraction: f32,
    /// Wall time the run has spent *running* (pause time excluded).
    pub wall_elapsed: Duration,
}

impl Default for Telemetry {
    fn default() -> Self {
        Self {
            substeps: 0,
            rejections_picard: 0,
            rejections_truncation: 0,
            picard_iters: 0,
            degraded: 0,
            floor_hits: 0,
            current_dt: 0.0,
            dt_p50: 0.0,
            dt_p95: 0.0,
            peak_energy_err: 0.0,
            current_energy_err: 0.0,
            substeps_per_second_window: 0.0,
            sim_time_per_second_window: 0.0,
            last_progress_fraction: 0.0,
            wall_elapsed: Duration::ZERO,
        }
    }
}

impl Telemetry {
    /// Total rejections this run (Picard + truncation).
    pub fn rejections_total(&self) -> u64 {
        self.rejections_picard.saturating_add(self.rejections_truncation)
    }

    /// Acceptance rate `[0.0, 1.0]`. Returns `1.0` when no substeps
    /// have been attempted yet (division-by-zero guard).
    pub fn acceptance_rate(&self) -> f32 {
        let attempts = self.substeps + self.rejections_total();
        if attempts == 0 { 1.0 } else { (self.substeps as f32) / (attempts as f32) }
    }
}

/// Builder that accumulates readings between two snapshots.
///
/// The physics thread constructs one, feeds it per-substep
/// observations via the setter methods, then calls [`finish`](Self::finish)
/// to produce a [`Telemetry`] value it hands to
/// [`PrecisionRunController::update_telemetry`](super::PrecisionRunController::update_telemetry).
///
/// Today this is a thin wrapper; its purpose is to give the physics
/// thread a single place to compose a snapshot without the UI side
/// seeing intermediate states through a half-updated `Telemetry`
/// field.
#[derive(Debug, Default, Clone)]
pub struct TelemetryBuilder {
    snapshot: Telemetry,
}

impl TelemetryBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_substeps(mut self, n: u64) -> Self {
        self.snapshot.substeps = n;
        self
    }
    pub fn with_rejections_picard(mut self, n: u64) -> Self {
        self.snapshot.rejections_picard = n;
        self
    }
    pub fn with_rejections_truncation(mut self, n: u64) -> Self {
        self.snapshot.rejections_truncation = n;
        self
    }
    pub fn with_picard_iters(mut self, n: u64) -> Self {
        self.snapshot.picard_iters = n;
        self
    }
    pub fn with_degraded(mut self, n: u64) -> Self {
        self.snapshot.degraded = n;
        self
    }
    pub fn with_floor_hits(mut self, n: u64) -> Self {
        self.snapshot.floor_hits = n;
        self
    }

    pub fn with_current_dt(mut self, dt: f64) -> Self {
        self.snapshot.current_dt = dt;
        self
    }
    pub fn with_dt_p50(mut self, dt: f64) -> Self {
        self.snapshot.dt_p50 = dt;
        self
    }
    pub fn with_dt_p95(mut self, dt: f64) -> Self {
        self.snapshot.dt_p95 = dt;
        self
    }

    pub fn with_peak_energy_err(mut self, err: f64) -> Self {
        self.snapshot.peak_energy_err = err;
        self
    }
    pub fn with_current_energy_err(mut self, err: f64) -> Self {
        self.snapshot.current_energy_err = err;
        self
    }

    pub fn with_substeps_per_second(mut self, rate: f64) -> Self {
        self.snapshot.substeps_per_second_window = rate;
        self
    }
    pub fn with_sim_time_per_second(mut self, rate: f64) -> Self {
        self.snapshot.sim_time_per_second_window = rate;
        self
    }

    pub fn with_progress_fraction(mut self, f: f32) -> Self {
        self.snapshot.last_progress_fraction = f;
        self
    }
    pub fn with_wall_elapsed(mut self, d: Duration) -> Self {
        self.snapshot.wall_elapsed = d;
        self
    }

    pub fn finish(self) -> Telemetry {
        self.snapshot
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_all_zero() {
        let t = Telemetry::default();
        assert_eq!(t.substeps, 0);
        assert_eq!(t.rejections_total(), 0);
        assert_eq!(t.acceptance_rate(), 1.0);
        assert_eq!(t.wall_elapsed, Duration::ZERO);
    }

    #[test]
    fn acceptance_rate_with_mixed_outcomes() {
        let t = TelemetryBuilder::new()
            .with_substeps(80)
            .with_rejections_picard(5)
            .with_rejections_truncation(15)
            .finish();
        assert_eq!(t.rejections_total(), 20);
        assert!((t.acceptance_rate() - 0.8).abs() < 1e-6);
    }

    #[test]
    fn builder_sets_every_field() {
        let t = TelemetryBuilder::new()
            .with_substeps(10)
            .with_current_dt(1.5e-3)
            .with_dt_p50(1.0e-3)
            .with_dt_p95(3.0e-3)
            .with_peak_energy_err(2e-12)
            .with_substeps_per_second(400.0)
            .with_sim_time_per_second(0.5)
            .with_progress_fraction(0.25)
            .with_wall_elapsed(Duration::from_millis(250))
            .finish();
        assert_eq!(t.substeps, 10);
        assert_eq!(t.current_dt, 1.5e-3);
        assert_eq!(t.dt_p50, 1.0e-3);
        assert_eq!(t.dt_p95, 3.0e-3);
        assert_eq!(t.peak_energy_err, 2e-12);
        assert_eq!(t.substeps_per_second_window, 400.0);
        assert_eq!(t.sim_time_per_second_window, 0.5);
        assert_eq!(t.last_progress_fraction, 0.25);
        assert_eq!(t.wall_elapsed, Duration::from_millis(250));
    }
}
