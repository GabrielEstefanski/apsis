//! Auto-exposure domain state — the scalar bookkeeping that sits between
//! the GPU luminance reducer (see [`crate::render::luminance_reducer`])
//! and the tonemap pipeline.
//!
//! Keeping this file pure (no `wgpu` types) means the EMA math, clamps,
//! and parameter interpretation are unit-testable without a GPU, and the
//! `wgpu_backend` wiring has one obvious adapter point: take a measured
//! soft-max luminance in, get back an exposure scale to feed into the
//! existing [`crate::render::tonemap::TonemapPipeline::set_exposure`] API.
//!
//! # Why soft-max (power-mean), not average or raw max
//!
//! The simulation's canvas is dominated by empty black. A classical
//! log-average exposure meter would "open up" the aperture in a zoomed-
//! out frame and pop as soon as a star enters view — the opposite of
//! what the eye wants. A raw max is the other extreme: one hot pixel
//! (a star centre in HDR) dominates the reading and the camera chases
//! that pixel forever.
//!
//! The GPU reducer computes a **Minkowski p-norm** of luminance:
//!
//! ```text
//! L_soft = ( mean(L^p) )^(1/p)
//! ```
//!
//! For `p → 1` this becomes the average; for `p → ∞` it becomes the max.
//! `p = 4` sits in the sweet spot: strongly weighted toward bright
//! regions (so empty space doesn't sway the meter), but still averages
//! over those bright regions so a single outlier pixel doesn't drive
//! the whole exposure. This module's job is to take that scalar,
//! temporally smooth it, and clamp to usable bounds.
//!
//! # EMA in log space
//!
//! The exposure scale is smoothed in log space so a 2× brightness jump
//! takes the same adaptation time as a 0.5× dimming — matching how the
//! human eye adapts (log-luminance response). Linear EMA would make
//! up-adaptation feel faster than down-adaptation and is what produces
//! the classic "game engine pops when you turn a corner" artefact.

/// Luminance contribution weights for Rec. 709 / sRGB primaries.
///
/// Mirrored in the luminance reducer's WGSL init shader. Kept here so
/// the coefficient is documented in one place; any future colour-space
/// change (e.g. Rec. 2020) updates both call sites together.
pub const LUMA_WEIGHTS: [f32; 3] = [0.2126, 0.7152, 0.0722];

/// Power used in the Minkowski p-norm soft-max on the GPU. Kept in
/// domain code because the CPU EMA undoes the exponent (`L = m^(1/p)`)
/// after reading back the reduced texel.
///
/// `p = 4` biases strongly toward bright regions without letting a
/// single hot pixel dominate — smaller values drift toward mean-metering
/// (bad for sparse scenes), larger values toward raw-max (brittle).
pub const SOFT_MAX_P: f32 = 4.0;

/// Auto-exposure state + tuning knobs.
///
/// Owned by the renderer; a UI layer or save file can override the
/// tunable fields. The `current_scale` field is runtime state — the
/// smoothed exposure value the tonemap actually consumes — and gets
/// updated by [`ExposureState::tick`] each frame when a new GPU
/// measurement is available.
#[derive(Clone, Debug)]
pub struct ExposureState {
    /// Master on/off. When `false`, [`ExposureState::current_scale`]
    /// returns [`ExposureState::manual_scale`] and the GPU reducer can
    /// be skipped entirely.
    pub enabled: bool,

    /// What we want the soft-max luminance to *look like* after
    /// exposure is applied. The tonemap then runs ACES over that,
    /// which knees down anything above ~1.
    ///
    /// We aim at **0.18** (photography's "middle-grey" key) rather than
    /// at the top of ACES' linear region: with a p-norm soft-max the
    /// metric is biased toward bright regions already, so targeting
    /// middle-grey keeps *peaks* in ACES' linear range and reserves the
    /// upper knee for genuine highlights (a star's hot core, a nearby
    /// flyby) instead of parking the average bright body at the knee.
    /// 0.8 produced visible blow-out on planets once the scale opened
    /// up for a sparse zoomed-out Solar System frame.
    pub target_luminance: f32,

    /// Time for the EMA to cover half of a change, in seconds. Smaller
    /// values = snappier adaptation; larger = more film-like lag.
    /// 0.5s is a good default — noticeable but not laggy.
    pub ema_half_life_sec: f32,

    /// Lower bound on the exposure scale. Prevents the meter from
    /// collapsing to zero in a completely black frame (all bodies
    /// off-screen) and producing a totally black tonemap output.
    pub min_scale: f32,

    /// Upper bound. Prevents an isolated dim frame (e.g. between
    /// template loads) from blowing the gain up to the point where
    /// the first bright pixel on the next frame clips white.
    pub max_scale: f32,

    /// Used when `enabled = false` — passed straight through.
    /// Also used as the initial value of `current_scale`.
    pub manual_scale: f32,

    /// Smoothed exposure scale currently fed to the tonemap. Updated
    /// by [`tick`](Self::tick); read by the backend.
    pub current_scale: f32,
}

impl Default for ExposureState {
    fn default() -> Self {
        Self {
            enabled: true,
            target_luminance: 0.18,
            ema_half_life_sec: 0.5,
            // Asymmetric bounds: darkening is cheap (HDR peaks take care
            // of themselves once compressed by ACES) so we leave 3 stops
            // below unity; brightening is where over-exposure lives, so
            // we cap at 2 stops above unity. A wider `max_scale` (32×)
            // produced blown-out grids and planets on sparse zoomed-out
            // frames — the meter reads a small bright region, divides
            // into `target`, and inflates the gain unrealistically.
            min_scale: 1.0 / 8.0,
            max_scale: 4.0,
            manual_scale: 1.0,
            current_scale: 1.0,
        }
    }
}

impl ExposureState {
    /// Advances the EMA toward the target computed from a measured
    /// soft-max luminance, returns the new `current_scale`.
    ///
    /// * `measured_soft_max` — the power-mean already raised to `1/p`
    ///   on the CPU side (see [`decode_reduced_texel`]). Must be
    ///   non-negative; values ≤ 0 are treated as "no valid reading"
    ///   and leave `current_scale` untouched.
    /// * `dt_sec` — wall-clock seconds since the previous tick. Used
    ///   to convert the half-life into a per-frame blend factor so
    ///   adaptation speed is framerate-independent.
    ///
    /// When `enabled` is false this is a no-op that returns
    /// `manual_scale` — the backend can still call it every frame
    /// without branching.
    pub fn tick(&mut self, measured_soft_max: f32, dt_sec: f32) -> f32 {
        if !self.enabled {
            self.current_scale = self.manual_scale.max(0.0);
            return self.current_scale;
        }
        if !(measured_soft_max > 0.0) || !(dt_sec > 0.0) {
            return self.current_scale;
        }

        // Target scale this frame — what exposure would make
        // soft_max land exactly at target_luminance.
        let raw_target =
            (self.target_luminance / measured_soft_max).clamp(self.min_scale, self.max_scale);

        // EMA in log space. Half-life h → per-step alpha such that
        // applying (1-alpha)^n to the error halves it when n·dt = h:
        //     (1-alpha) = 0.5^(dt / h)
        //     alpha     = 1 - 0.5^(dt / h)
        //
        // Working on log(scale) so a 2× ramp-up and 0.5× ramp-down
        // take the same perceptual time. Linear EMA in the scale
        // itself is asymmetric and produces the "engine pops" feel.
        let alpha = 1.0 - 0.5_f32.powf(dt_sec / self.ema_half_life_sec.max(1e-4));
        let log_curr = self.current_scale.max(1e-6).ln();
        let log_target = raw_target.max(1e-6).ln();
        let log_new = log_curr + alpha * (log_target - log_curr);
        self.current_scale = log_new.exp().clamp(self.min_scale, self.max_scale);

        self.current_scale
    }
}

/// Converts the single f16 texel read back from the reducer's final
/// mip (the `mean(L^p)` over the whole HDR frame, stored at 1×1) into
/// a soft-max luminance ready for [`ExposureState::tick`].
///
/// The GPU stores `L^p` in `R16Float` and reduces by bilinear averaging
/// through the mip chain, so the final texel holds `mean(L^p)`. The
/// 1/p root happens here on the CPU — doing it on the GPU would cost a
/// dedicated single-pixel pass for one `pow`, not worth it.
#[inline]
pub fn decode_reduced_texel(mean_l_to_p: f32) -> f32 {
    if mean_l_to_p <= 0.0 {
        return 0.0;
    }
    mean_l_to_p.powf(1.0 / SOFT_MAX_P)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_passes_manual_scale() {
        let mut e = ExposureState { enabled: false, manual_scale: 2.5, ..Default::default() };
        let s = e.tick(999.0, 1.0 / 60.0);
        assert_eq!(s, 2.5);
        assert_eq!(e.current_scale, 2.5);
    }

    #[test]
    fn invalid_measurement_leaves_state_untouched() {
        let mut e = ExposureState { current_scale: 1.2, ..Default::default() };
        let s = e.tick(0.0, 1.0 / 60.0);
        assert!((s - 1.2).abs() < 1e-6);
        let s = e.tick(-5.0, 1.0 / 60.0);
        assert!((s - 1.2).abs() < 1e-6);
        // Zero dt is also a no-op (prevents division-by-zero paths).
        let s = e.tick(0.5, 0.0);
        assert!((s - 1.2).abs() < 1e-6);
    }

    #[test]
    fn half_life_takes_roughly_one_half_life_to_halve_error() {
        // Starting at scale=1.0 and walking toward target_scale=4.0.
        // After one half-life's worth of ticks, the log-space error
        // should be halved: log(2.0) ≈ 0.5 · log(4.0).
        let mut e = ExposureState {
            target_luminance: 0.8,
            ema_half_life_sec: 0.5,
            min_scale: 1e-3,
            max_scale: 1e3,
            ..Default::default()
        };
        // target = 0.8 / measured → measured = 0.8 / 4 = 0.2 for target_scale = 4.0
        let measured = 0.8 / 4.0;
        let dt = 0.5; // exactly one half-life in one tick

        e.tick(measured, dt);
        // After one half-life the remaining log-error should be ≈ half.
        // log(4) = 1.386; half = 0.693 ≈ log(2), so scale ≈ 2.
        assert!(
            (e.current_scale - 2.0).abs() < 0.05,
            "expected ~2.0 after one half-life, got {}",
            e.current_scale
        );
    }

    #[test]
    fn clamps_target_to_min_max() {
        let mut e = ExposureState {
            target_luminance: 0.8,
            min_scale: 0.25,
            max_scale: 4.0,
            ema_half_life_sec: 0.0001, // snap to target
            ..Default::default()
        };

        // Wildly bright frame → target_scale would be tiny; clamps to min.
        e.current_scale = 1.0;
        e.tick(100.0, 1.0);
        assert!(
            (e.current_scale - 0.25).abs() < 1e-3,
            "expected clamp to min_scale, got {}",
            e.current_scale
        );

        // Wildly dim frame → clamps to max.
        e.current_scale = 1.0;
        e.tick(1e-6, 1.0);
        assert!(
            (e.current_scale - 4.0).abs() < 1e-3,
            "expected clamp to max_scale, got {}",
            e.current_scale
        );
    }

    #[test]
    fn log_space_adaptation_is_symmetric() {
        // A factor-of-4 brighten and a factor-of-4 dim should take the
        // same number of ticks to halve the log-error.
        let cfg = ExposureState {
            target_luminance: 0.8,
            ema_half_life_sec: 0.2,
            min_scale: 1e-4,
            max_scale: 1e4,
            ..Default::default()
        };

        let mut up = cfg.clone();
        up.current_scale = 1.0;
        up.tick(0.2, 0.2); // target_scale = 4

        let mut down = cfg.clone();
        down.current_scale = 1.0;
        down.tick(3.2, 0.2); // target_scale = 0.25

        // In log space, up went from log(1)=0 toward log(4), down from 0 toward log(0.25).
        // Both should be at ±log(2) after one half-life.
        assert!(
            (up.current_scale.ln() + down.current_scale.ln()).abs() < 0.01,
            "log-space adaptation not symmetric: up={} down={}",
            up.current_scale,
            down.current_scale
        );
    }

    #[test]
    fn decode_root_matches_power() {
        // If L_actual = 0.6, then mean(L^p) in a uniform frame is 0.6^p,
        // and decode_reduced_texel should invert that cleanly.
        let l_actual = 0.6_f32;
        let stored = l_actual.powf(SOFT_MAX_P);
        let decoded = decode_reduced_texel(stored);
        assert!(
            (decoded - l_actual).abs() < 1e-4,
            "round-trip lost precision: actual={}, decoded={}",
            l_actual,
            decoded
        );
    }

    #[test]
    fn decode_handles_zero_and_negatives() {
        assert_eq!(decode_reduced_texel(0.0), 0.0);
        assert_eq!(decode_reduced_texel(-1.0), 0.0);
    }
}
