//! [`Colormap`] — maps a normalized scalar `t ∈ [0, 1]` to an RGB triple.
//!
//! Implementations are typically defined by a fixed set of *stops* (as used
//! in the matplotlib colormap source) and linearly interpolated between
//! them; [`sample_stops`] is the helper they share.

/// A colormap: `t → RGB`. Must accept any finite `t` and clamp internally.
pub trait Colormap: Send + Sync {
    fn id(&self) -> &'static str;
    fn name(&self) -> &'static str;
    /// Samples the colormap at `t`. Implementations clamp `t` to `[0, 1]`.
    fn sample(&self, t: f32) -> [u8; 3];
}

/// Piecewise-linear interpolation between a fixed stop table.
///
/// `stops[i] = (position, rgb)` with positions non-decreasing over
/// `[0, 1]`. Used by every built-in colormap; callers hand their const
/// stop table in directly.
pub fn sample_stops(stops: &[(f32, [u8; 3])], t: f32) -> [u8; 3] {
    if stops.is_empty() {
        return [0, 0, 0];
    }
    let t = t.clamp(0.0, 1.0);

    // Below first / above last stop.
    if t <= stops[0].0 {
        return stops[0].1;
    }
    if t >= stops[stops.len() - 1].0 {
        return stops[stops.len() - 1].1;
    }

    // Linear search is fine — a typical table has 5–16 entries.
    for w in stops.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t >= t0 && t <= t1 {
            let span = (t1 - t0).max(1e-6);
            let f = (t - t0) / span;
            return [lerp_u8(c0[0], c1[0], f), lerp_u8(c0[1], c1[1], f), lerp_u8(c0[2], c1[2], f)];
        }
    }
    stops[stops.len() - 1].1
}

#[inline]
fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    let af = a as f32;
    let bf = b as f32;
    (af + (bf - af) * t).round().clamp(0.0, 255.0) as u8
}
