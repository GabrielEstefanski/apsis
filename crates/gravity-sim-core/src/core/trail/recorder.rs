//! [`TrailRecorder`] — owner of the trail ring buffer.
//!
//! The recorder is a **pure consumer** of sample columns produced by the
//! physics thread. Sampling cadence decisions live inside the physics
//! thread's [`TrailSampler`](crate::core::trail::sampler::TrailSampler) —
//! they must happen at physics-step granularity so that fast orbits never
//! alias to polygons and high-SPF runs never starve the trail.
//!
//! # Responsibilities
//!
//! | Concern                           | Method                      |
//! |-----------------------------------|-----------------------------|
//! | Ingest new samples this frame     | [`ingest`]                  |
//! | COM recenter correction           | [`apply_com_shift`]         |
//! | Topology changes (add/remove body)| auto-detected in [`ingest`] |
//! | Colour updates                    | auto-detected in [`ingest`] |
//! | Snapshot save/restore             | [`to_snapshot`] / [`restore_from_snapshot`] |
//! | Sampler config (for physics)      | [`sampler_kind`]            |

use crate::core::trail::buffer::{TrailBuffer, adaptive_capacity};
use crate::core::trail::sampler::TrailSamplerKind;
use crate::domain::body::Body;
use crate::io::snapshot::TrailSnapshot;

/// Baseline arc-length trigger: a body must traverse ~2 % of the scene scale
/// before a new sample is recorded. Yields ~314 samples per orbital period
/// for a circular orbit — above the eye's polygonal-aliasing threshold.
const BASE_THRESHOLD_FRAC: f32 = 0.02;

/// Core-side trail recorder.
///
/// Owns the ring buffer. The sampling strategy itself lives in the physics
/// thread; this struct only exposes the *config* that drives it so the UI
/// can tune density (via [`set_interval_multiplier`]) and the snapshot
/// layer can round-trip the multiplier.
pub struct TrailRecorder {
    buffer: TrailBuffer,
    /// `trail_every` semantic preserved from the legacy UI control:
    /// multiplies the arc-length threshold (larger → sparser samples).
    interval_multiplier: usize,
    last_n_bodies: usize,
}

impl TrailRecorder {
    pub fn new() -> Self {
        Self {
            buffer: TrailBuffer::new(0),
            interval_multiplier: 1,
            last_n_bodies: 0,
        }
    }

    // ── Frame ingest ──────────────────────────────────────────────────────────

    /// Consume the samples produced by the physics thread this frame.
    ///
    /// `samples` is the list of position columns as delivered by
    /// [`PhysicsHandle::take_trail_samples`](crate::core::physics_thread::PhysicsHandle::take_trail_samples).
    /// `bodies` is the current body list — used for topology detection and
    /// colour sync.
    ///
    /// On topology change the buffer is reset and any stale samples are
    /// silently dropped (their column widths no longer match).
    pub fn ingest(&mut self, samples: Vec<Vec<[f32; 2]>>, bodies: &[Body]) {
        self.ingest_with_colors(samples, bodies, None);
    }

    /// Same as [`ingest`] but takes an optional per-body RGB override.
    ///
    /// `colors_override = Some(rgb)` is the entry point for the data-driven
    /// colour pipeline: trails pick up the exact colours produced by the
    /// active [`ColorView`](crate::render::color::ColorViewSelection) that
    /// frame, so the trail gradient matches the body rendering. When
    /// `colors_override = None` the recorder falls back to each body's
    /// material colour — the pre-existing behaviour.
    pub fn ingest_with_colors(
        &mut self,
        samples: Vec<Vec<[f32; 2]>>,
        bodies: &[Body],
        colors_override: Option<&[[u8; 3]]>,
    ) {
        let n = bodies.len();

        if n != self.last_n_bodies {
            let cap = adaptive_capacity(n.max(1));
            self.buffer.reset(n, cap);
            self.last_n_bodies = n;
        }

        if n == 0 {
            return;
        }

        match colors_override {
            Some(rgb) if rgb.len() == n => self.buffer.set_colors_rgb(rgb),
            _ => self.buffer.update_colors(bodies),
        }

        // Any column whose width doesn't match the current body count was
        // produced before a topology change and is silently discarded.
        for col in samples {
            if col.len() == n {
                self.buffer.push_column(&col);
            }
        }
    }

    // ── COM shift correction ───────────────────────────────────────────────────

    /// Applies an accumulated COM translation to all stored positions.
    pub fn apply_com_shift(&mut self, dx: f32, dy: f32) {
        if dx != 0.0 || dy != 0.0 {
            self.buffer.translate(dx, dy);
        }
    }

    // ── Sampler configuration ─────────────────────────────────────────────────

    /// Sets the legacy `trail_every` multiplier. The effective arc-length
    /// threshold is `BASE_THRESHOLD_FRAC × multiplier` — larger values
    /// produce sparser trails, matching the old step-based semantic.
    pub fn set_interval_multiplier(&mut self, n: usize) {
        self.interval_multiplier = n.max(1);
    }

    pub fn interval_multiplier(&self) -> usize {
        self.interval_multiplier
    }

    /// The sampler configuration that should be sent to the physics thread
    /// whenever this recorder's density settings change.
    pub fn sampler_kind(&self) -> TrailSamplerKind {
        TrailSamplerKind::ArcLength {
            threshold_frac: BASE_THRESHOLD_FRAC * self.interval_multiplier as f32,
        }
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
    pub fn restore_from_snapshot(&mut self, snap: &TrailSnapshot, bodies: &[Body]) {
        let n = bodies.len();
        let cap = adaptive_capacity(n.max(1));
        self.buffer.reset(n, cap);
        if n > 0 {
            self.buffer.update_colors(bodies);
        }
        self.last_n_bodies = n;

        if snap.n_bodies == n as u32
            && snap.positions.len() == (snap.n_bodies * snap.capacity) as usize
        {
            self.buffer.restore_from_snapshot(snap);
        }
    }

    /// Resets to an empty state (used when the simulation is fully reloaded).
    pub fn clear(&mut self) {
        self.buffer.reset(0, 1);
        self.last_n_bodies = 0;
    }
}

impl Default for TrailRecorder {
    fn default() -> Self {
        Self::new()
    }
}
