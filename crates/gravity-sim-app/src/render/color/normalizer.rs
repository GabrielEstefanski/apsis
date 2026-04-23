//! [`Normalizer`] — maps a physical scalar value to `[0, 1]` given a data
//! range.
//!
//! Separated from [`Colormap`](super::colormap::Colormap) because the
//! "which transform before the colour ramp" choice (linear / log / …) is
//! orthogonal to the colour ramp itself. This is the SPLASH / yt pattern:
//! `field → normalize(value, range) → colormap.sample(t)`.

/// Normalizes a physical scalar into the unit interval given a data range.
pub trait Normalizer: Send + Sync {
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;
    /// Map `v` into `[0, 1]` using `(min, max)` from the data. Implementations
    /// must handle degenerate ranges (`min == max`) and non-finite inputs by
    /// returning a finite value in `[0, 1]`.
    fn normalize(&self, v: f64, range: (f64, f64)) -> f32;
}

/// Linear normalizer: `t = (v - min) / (max - min)`.
pub struct LinearNormalizer;

impl Normalizer for LinearNormalizer {
    fn id(&self) -> &'static str {
        "linear"
    }
    fn name(&self) -> &'static str {
        "Linear"
    }
    fn normalize(&self, v: f64, (lo, hi): (f64, f64)) -> f32 {
        if !v.is_finite() {
            return 0.0;
        }
        let span = hi - lo;
        if span.abs() < 1e-300 {
            return 0.5;
        }
        (((v - lo) / span) as f32).clamp(0.0, 1.0)
    }
}

/// Log normalizer: `t = (log(v) - log(min)) / (log(max) - log(min))`.
///
/// Non-positive values and range endpoints are clamped to a small floor so
/// the function is total — SPLASH takes the same approach for surface
/// density plots where a small fraction of SPH particles can legitimately
/// have zero contribution.
pub struct LogNormalizer;

impl Normalizer for LogNormalizer {
    fn id(&self) -> &'static str {
        "log"
    }
    fn name(&self) -> &'static str {
        "Log"
    }
    fn normalize(&self, v: f64, (lo, hi): (f64, f64)) -> f32 {
        const FLOOR: f64 = 1e-30;
        let lv = v.max(FLOOR).ln();
        let llo = lo.max(FLOOR).ln();
        let lhi = hi.max(FLOOR).ln();
        let span = lhi - llo;
        if !lv.is_finite() || span.abs() < 1e-300 {
            return 0.5;
        }
        (((lv - llo) / span) as f32).clamp(0.0, 1.0)
    }
}
