//! [`TrailSampler`] — strategy deciding *when* to record a position sample.
//!
//! Trail sampling is a **visualization concern** driven by geometry, not
//! wall-clock. The pathologies we need to avoid:
//!
//! - **Low steps-per-frame + fast orbits** → if we sample per-frame, short
//!   orbital periods get only 2-3 samples per loop and render as polygons.
//! - **High steps-per-frame** → if we sample per-step with a fixed time
//!   interval, thousands of samples stack into opaque bands, then starve
//!   when the auto-widened interval overshoots simulation time.
//!
//! The fix is to run sampling **at physics-step granularity** using an
//! **arc-length criterion** — fire whenever a body has moved enough (as a
//! fraction of the scene scale) to warrant a new trail vertex. This makes
//! trail density invariant to both SPF and orbital period.
//!
//! # Strategies
//!
//! | Strategy             | Signal used                 | Use case               |
//! |----------------------|-----------------------------|------------------------|
//! | [`ArcLengthSampler`] | max body displacement ratio | live interactive view  |
//! | [`StepSampler`]      | physics-step count          | deterministic replay   |

use crate::domain::body::Body;

/// Decides — per physics step — whether the current body state should be
/// appended to the trail buffer as a new sample column.
///
/// Implementations own the state needed to make that decision (e.g. anchor
/// positions). They never own the buffer itself.
pub trait TrailSampler: Send {
    /// Should the current body positions be recorded as a trail sample?
    ///
    /// Called once per physics step. When the implementation returns `true`
    /// it must have updated its internal state so that subsequent calls with
    /// the same inputs return `false` (i.e. the decision is idempotent for
    /// a given displacement).
    fn should_sample(&mut self, bodies: &[Body]) -> bool;

    /// Discard internal state. Called after snapshot restore or a reset that
    /// invalidates anchor positions.
    fn reset(&mut self);
}

// ── ArcLengthSampler ──────────────────────────────────────────────────────────

const EPS_SCALE_SQ: f32 = 1e-12;

/// Samples whenever the largest body displacement exceeds a configurable
/// **fraction of the scene scale**.
///
/// The scene scale is the RMS distance of all bodies from the origin,
/// recomputed each step. Using a ratio rather than an absolute distance makes
/// the sampler scale-invariant: a solar system at 1 AU and a TRAPPIST-1
/// system at 0.05 AU both produce visually identical trail densities.
///
/// # Density in orbital terms
///
/// For a circular orbit, `displacement / radius ≈ Δθ`. A threshold of
/// `0.02` therefore yields roughly `2π / 0.02 ≈ 314` samples per revolution —
/// well above the polygonal-aliasing threshold.
pub struct ArcLengthSampler {
    /// Last recorded world-space position of each body (anchor).
    /// Resized on topology change.
    anchors: Vec<(f32, f32)>,
    /// Squared trigger threshold, in units of `(displacement / scene_scale)²`.
    threshold_sq: f32,
}

impl ArcLengthSampler {
    pub fn new(threshold_frac: f32) -> Self {
        let t = threshold_frac.max(1e-6);
        Self { anchors: Vec::new(), threshold_sq: t * t }
    }
}

impl TrailSampler for ArcLengthSampler {
    fn should_sample(&mut self, bodies: &[Body]) -> bool {
        let n = bodies.len();

        // Topology change — reseat anchors and force a sample so the new
        // body's trail starts immediately.
        if self.anchors.len() != n {
            self.anchors.clear();
            self.anchors.extend(bodies.iter().map(|b| (b.x as f32, b.y as f32)));
            return n > 0;
        }
        if n == 0 {
            return false;
        }

        // Scene scale = RMS distance of bodies from origin. Clamped to a
        // small epsilon so a fully-centred cluster (all at origin) doesn't
        // blow up the ratio test.
        let mut scale_sq_acc = 0.0_f32;
        for b in bodies {
            let x = b.x as f32;
            let y = b.y as f32;
            scale_sq_acc += x * x + y * y;
        }
        let scene_scale_sq = (scale_sq_acc / n as f32).max(EPS_SCALE_SQ);

        // Maximum squared displacement since the last recorded sample.
        let mut max_disp_sq = 0.0_f32;
        for (i, b) in bodies.iter().enumerate() {
            let dx = b.x as f32 - self.anchors[i].0;
            let dy = b.y as f32 - self.anchors[i].1;
            let d = dx * dx + dy * dy;
            if d > max_disp_sq {
                max_disp_sq = d;
            }
        }

        let ratio_sq = max_disp_sq / scene_scale_sq;
        if ratio_sq >= self.threshold_sq {
            for (i, b) in bodies.iter().enumerate() {
                self.anchors[i] = (b.x as f32, b.y as f32);
            }
            true
        } else {
            false
        }
    }

    fn reset(&mut self) {
        self.anchors.clear();
    }
}

// ── StepSampler ───────────────────────────────────────────────────────────────

/// Fires every `interval` physics steps. Deterministic, replay-stable.
pub struct StepSampler {
    interval: u32,
    counter: u32,
}

impl StepSampler {
    pub fn new(interval: u32) -> Self {
        Self { interval: interval.max(1), counter: 0 }
    }
}

impl TrailSampler for StepSampler {
    fn should_sample(&mut self, _bodies: &[Body]) -> bool {
        self.counter += 1;
        if self.counter >= self.interval {
            self.counter = 0;
            true
        } else {
            false
        }
    }

    fn reset(&mut self) {
        self.counter = 0;
    }
}

// ── Config / factory ──────────────────────────────────────────────────────────

/// Serializable configuration describing *which* sampler to build.
///
/// Lives in the render layer (samplers are a visualization concern) but is
/// sent to the physics thread via [`crate::core::physics_thread::PhysicsCmd::SetTrailSampler`]
/// so sampling can happen at physics-step granularity.
#[derive(Clone, Copy, Debug)]
pub enum TrailSamplerKind {
    /// Fire when max-body displacement exceeds `threshold_frac` × scene scale.
    ArcLength { threshold_frac: f32 },
    /// Fire every `interval` physics steps.
    Step { interval: u32 },
}

impl TrailSamplerKind {
    /// Constructs the concrete sampler for this configuration.
    pub fn build(&self) -> Box<dyn TrailSampler> {
        match *self {
            Self::ArcLength { threshold_frac } => Box::new(ArcLengthSampler::new(threshold_frac)),
            Self::Step { interval } => Box::new(StepSampler::new(interval)),
        }
    }
}

impl Default for TrailSamplerKind {
    /// Arc-length with a threshold that yields ~314 samples per revolution.
    fn default() -> Self {
        Self::ArcLength { threshold_frac: 0.02 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;

    fn body_at(x: f64, y: f64) -> Body {
        Body::rocky(1.0).at(x, y).with_velocity(0.0, 0.0)
    }

    #[test]
    fn step_sampler_fires_on_interval() {
        let mut s = StepSampler::new(3);
        let bodies = vec![body_at(0.0, 0.0)];
        assert!(!s.should_sample(&bodies));
        assert!(!s.should_sample(&bodies));
        assert!(s.should_sample(&bodies));
        assert!(!s.should_sample(&bodies));
    }

    #[test]
    fn arc_length_fires_on_topology_change() {
        let mut s = ArcLengthSampler::new(0.1);
        let bodies = vec![body_at(1.0, 0.0), body_at(-1.0, 0.0)];
        assert!(s.should_sample(&bodies));
        // Same positions → no new sample.
        assert!(!s.should_sample(&bodies));
    }

    #[test]
    fn arc_length_fires_after_displacement() {
        let mut s = ArcLengthSampler::new(0.05);
        let mut bodies = vec![body_at(1.0, 0.0), body_at(-1.0, 0.0)];
        let _ = s.should_sample(&bodies); // seat anchors
        // Scene scale sqrt = 1.0. Move body 0 by 0.1 → ratio 0.1 > 0.05.
        bodies[0].x = 1.1;
        assert!(s.should_sample(&bodies));
        // Anchors updated — no further fire without more motion.
        assert!(!s.should_sample(&bodies));
    }
}
