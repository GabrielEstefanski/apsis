//! [`TrailRecorder`] — render-side owner of the trail ring buffer.
//!
//! Decouples trail sampling entirely from the physics domain: the physics
//! thread publishes body positions at ~60 Hz via `RenderState`; this
//! recorder runs on the UI thread, calls its [`TrailSampler`] strategy to
//! decide when to push a sample, and owns the [`TrailBuffer`] that the GPU
//! renderer consumes.
//!
//! # Responsibilities
//!
//! | Concern | Method |
//! |---------|--------|
//! | Sampling cadence | [`tick`] — called once per rendered frame |
//! | COM recenter correction | [`apply_com_shift`] |
//! | Topology changes (add/remove body) | auto-detected in [`tick`] |
//! | Colour updates | auto-detected in [`tick`] |
//! | Snapshot save/restore | [`to_snapshot`] / [`restore_from_snapshot`] |

use crate::domain::body::Body;
use crate::io::snapshot::TrailSnapshot;
use crate::render::trail::sampler::{AdaptiveSampler, SampleDecision, TrailSampler};
use crate::render::trail_buffer::{TrailBuffer, adaptive_capacity};

/// Render-side trail recorder.
///
/// Owns the ring buffer and the sampling strategy. Lives in
/// [`crate::app::ui::SimulationApp`] and is ticked once per render frame.
pub struct TrailRecorder {
    buffer: TrailBuffer,
    sampler: AdaptiveSampler,
    last_sim_t: f64,
    /// Baseline interval multiplier corresponding to the `trail_every` UI
    /// control (kept for snapshot round-trip parity).
    interval_multiplier: usize,
    // Track whether colours need refreshing (happens when body properties change).
    last_n_bodies: usize,
}

impl TrailRecorder {
    const BASE_INTERVAL: f64 = 0.01;
    const MAX_PER_FRAME: u32 = 64;

    pub fn new() -> Self {
        Self {
            buffer: TrailBuffer::new(0),
            sampler: AdaptiveSampler::new(Self::BASE_INTERVAL, Self::MAX_PER_FRAME),
            last_sim_t: 0.0,
            interval_multiplier: 1,
            last_n_bodies: 0,
        }
    }

    // ── Frame tick ─────────────────────────────────────────────────────────────

    /// Called once per rendered frame with the latest published body state.
    ///
    /// Handles topology changes (reset + recolour), colour sync, and sampling.
    /// `sim_t` is the cumulative simulation time; `steps_per_frame` is the SPF
    /// hint forwarded to the adaptive sampler.
    pub fn tick(&mut self, bodies: &[Body], sim_t: f64, steps_per_frame: u32) {
        let n = bodies.len();

        // ── Topology change ───────────────────────────────────────────────────
        if n != self.last_n_bodies {
            let cap = adaptive_capacity(n.max(1));
            self.buffer.reset(n, cap);
            self.buffer.update_colors(bodies);
            self.last_n_bodies = n;
            self.last_sim_t = sim_t;
            self.sampler.reset();
            return;
        }

        // ── Colour refresh (body properties may change each frame) ────────────
        self.buffer.update_colors(bodies);

        // ── Sampling ─────────────────────────────────────────────────────────
        if n > 0 {
            let dt = (sim_t - self.last_sim_t).max(0.0);
            self.last_sim_t = sim_t;

            if dt > 0.0 && self.sampler.decide(dt, steps_per_frame) == SampleDecision::Record {
                self.buffer.push(bodies);
            }
        }

        self.sampler.tick_frame();
    }

    // ── COM shift correction ───────────────────────────────────────────────────

    /// Applies an accumulated COM translation to all stored positions.
    ///
    /// Called each frame with the value published from the physics thread.
    /// `(0.0, 0.0)` is a no-op.
    pub fn apply_com_shift(&mut self, dx: f32, dy: f32) {
        if dx != 0.0 || dy != 0.0 {
            self.buffer.translate(dx, dy);
        }
    }

    // ── Sampler configuration ─────────────────────────────────────────────────

    /// Sets the interval multiplier (legacy `trail_every` control).
    ///
    /// The effective base interval becomes `BASE_INTERVAL * n`. Larger values
    /// produce sparser trails, matching the old behaviour.
    pub fn set_interval_multiplier(&mut self, n: usize) {
        let n = n.max(1);
        self.interval_multiplier = n;
        let interval = Self::BASE_INTERVAL * n as f64;
        self.sampler = AdaptiveSampler::new(interval, Self::MAX_PER_FRAME);
    }

    pub fn interval_multiplier(&self) -> usize {
        self.interval_multiplier
    }

    // ── Buffer access ─────────────────────────────────────────────────────────

    pub fn buffer(&self) -> &TrailBuffer {
        &self.buffer
    }

    // ── Snapshot ──────────────────────────────────────────────────────────────

    pub fn to_snapshot(&self) -> TrailSnapshot {
        self.buffer.to_snapshot()
    }

    /// Restores trail positions from a snapshot.
    ///
    /// `n_bodies` must match the current simulation body count.
    /// Silently skips restoration if the snapshot is incompatible (topology
    /// change between save and load).
    pub fn restore_from_snapshot(&mut self, snap: &TrailSnapshot, bodies: &[Body]) {
        let n = bodies.len();
        let cap = adaptive_capacity(n.max(1));
        self.buffer.reset(n, cap);
        self.buffer.update_colors(bodies);
        self.last_n_bodies = n;

        if snap.n_bodies == n as u32
            && snap.positions.len() == (snap.n_bodies * snap.capacity) as usize
        {
            self.buffer.restore_from_snapshot(snap);
        }
        self.sampler.reset();
    }

    /// Resets to an empty state (used when the simulation is fully reloaded).
    pub fn clear(&mut self) {
        self.buffer.reset(0, 1);
        self.last_n_bodies = 0;
        self.last_sim_t = 0.0;
        self.sampler.reset();
    }
}

impl Default for TrailRecorder {
    fn default() -> Self {
        Self::new()
    }
}
