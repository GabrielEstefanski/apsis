//! [`TrailSampler`] — strategy deciding *when* to record a position sample.
//!
//! Sampling cadence is a **visualization concern**, not a physics one: the
//! physics step runs at whatever rate gives numerical accuracy, while the
//! trail wants temporal/spatial density appropriate for human perception.
//! Decoupling the two lets us solve two pathologies independently:
//!
//! - **Low steps-per-frame + fast orbits** → sparse samples → polygonal aliasing
//! - **High steps-per-frame** → thousands of samples/frame → opaque stacking
//!
//! # Strategies
//!
//! | Strategy           | Signal used                  | Best for                 |
//! |--------------------|------------------------------|--------------------------|
//! | [`TimeSampler`]    | elapsed sim-time since last  | deterministic replay     |
//! | [`AdaptiveSampler`]| sim-time + per-frame cap     | live interactive viewing |
//! | *(future)* `ArcLengthSampler` | world-space distance   | zoom-independent detail  |
//!
//! All strategies are pure — they hold counters, never bodies or buffers.

/// Result of asking a sampler whether the current simulation state should
/// be recorded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleDecision {
    /// Do not record this step.
    Skip,
    /// Record this step into the trail buffer.
    Record,
}

/// Decides when a simulation state warrants a trail sample.
///
/// Implementors are called once per physics step (or once per UI tick in
/// Phase 2); they maintain whatever internal state they need. The method is
/// `&mut self` so counters can be updated without interior mutability.
pub trait TrailSampler: Send {
    /// Should the current step be recorded?
    ///
    /// - `dt`: simulation time elapsed since the previous call (sim units).
    /// - `steps_per_frame`: hint about render cadence; strategies may use
    ///   this to cap samples-per-frame. Pass `1` if unknown.
    fn decide(&mut self, dt: f64, steps_per_frame: u32) -> SampleDecision;

    /// Reset any internal state (e.g. after a seek or snapshot restore).
    fn reset(&mut self);
}

// ── TimeSampler ───────────────────────────────────────────────────────────────

/// Records whenever accumulated sim-time crosses a fixed `interval`.
///
/// Deterministic: replaying the same physics trace yields identical samples.
/// Suffers from the known pathologies when `interval` is fixed — use
/// [`AdaptiveSampler`] for interactive viewing.
#[derive(Clone, Debug)]
pub struct TimeSampler {
    interval: f64,
    acc: f64,
}

impl TimeSampler {
    pub fn new(interval: f64) -> Self {
        Self { interval: interval.max(1e-9), acc: 0.0 }
    }
}

impl TrailSampler for TimeSampler {
    fn decide(&mut self, dt: f64, _spf: u32) -> SampleDecision {
        self.acc += dt;
        if self.acc >= self.interval {
            self.acc -= self.interval;
            SampleDecision::Record
        } else {
            SampleDecision::Skip
        }
    }

    fn reset(&mut self) {
        self.acc = 0.0;
    }
}

// ── AdaptiveSampler ───────────────────────────────────────────────────────────

/// Time-based sampler that additionally caps **samples per frame**.
///
/// Solves the high-SPF stacking pathology: at `steps_per_frame = 10⁶` a pure
/// time sampler would record hundreds of thousands of columns per render
/// frame, all of which end up drawn on top of each other and saturate the
/// alpha blend to opaque.
///
/// The cap is enforced by auto-scaling the effective interval so that no
/// more than `max_per_frame` decisions come back `Record` per frame window.
#[derive(Clone, Debug)]
pub struct AdaptiveSampler {
    base_interval: f64,
    effective_interval: f64,
    acc: f64,
    frame_acc: u32,
    steps_this_frame: u32,
    max_per_frame: u32,
}

impl AdaptiveSampler {
    pub fn new(base_interval: f64, max_per_frame: u32) -> Self {
        Self {
            base_interval: base_interval.max(1e-9),
            effective_interval: base_interval.max(1e-9),
            acc: 0.0,
            frame_acc: 0,
            steps_this_frame: 0,
            max_per_frame: max_per_frame.max(1),
        }
    }

    fn rescale(&mut self, spf: u32) {
        // When spf is large, widen the interval so we expect roughly
        // `max_per_frame` records per frame. Formula:
        //   expected_records_per_frame ≈ spf / (interval / dt_avg)
        // We don't know dt here; approximate by assuming the caller passes
        // honest `dt`s and let `acc` drive selection. The widening factor
        // therefore depends on *observed* record density, not raw spf.
        let widen = (spf as f64 / self.max_per_frame as f64).max(1.0);
        self.effective_interval = self.base_interval * widen;
    }
}

impl TrailSampler for AdaptiveSampler {
    fn decide(&mut self, dt: f64, spf: u32) -> SampleDecision {
        if spf != self.steps_this_frame {
            self.rescale(spf);
            self.steps_this_frame = spf;
        }

        self.acc += dt;
        if self.acc >= self.effective_interval && self.frame_acc < self.max_per_frame {
            self.acc -= self.effective_interval;
            self.frame_acc += 1;
            SampleDecision::Record
        } else {
            SampleDecision::Skip
        }
    }

    fn reset(&mut self) {
        self.acc = 0.0;
        self.frame_acc = 0;
        self.effective_interval = self.base_interval;
    }
}

impl AdaptiveSampler {
    /// Called by the driver once per rendered frame to reset the per-frame
    /// record counter. Separate from [`TrailSampler::reset`] which wipes
    /// long-lived state.
    pub fn tick_frame(&mut self) {
        self.frame_acc = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_sampler_fires_on_interval() {
        let mut s = TimeSampler::new(0.1);
        assert_eq!(s.decide(0.05, 1), SampleDecision::Skip);
        assert_eq!(s.decide(0.05, 1), SampleDecision::Record);
        assert_eq!(s.decide(0.05, 1), SampleDecision::Skip);
    }

    #[test]
    fn adaptive_caps_per_frame() {
        let mut s = AdaptiveSampler::new(0.01, 3);
        // A single frame with 1_000_000 tiny dts; without the cap this would
        // fire ~10_000 times.
        let mut fired = 0;
        for _ in 0..1_000_000 {
            if s.decide(0.001, 1_000_000) == SampleDecision::Record {
                fired += 1;
            }
        }
        assert!(fired <= 3, "fired {fired}");
    }
}
