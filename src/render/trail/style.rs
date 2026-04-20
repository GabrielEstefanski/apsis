//! [`TrailStyle`] — immutable value object describing trail appearance.
//!
//! Holds every parameter the renderer needs to turn a `TrailBuffer` into
//! pixels. Injected into [`crate::render::WgpuBackend`] and forwarded into
//! the shader as uniform state. Swappable at runtime via
//! [`TrailStylePreset`].

/// Parameters that drive the trail shader.
///
/// All fields are visually orthogonal; change one without surprising side
/// effects on the others. Units where relevant are noted inline.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrailStyle {
    /// Half-width of the trail quad at its brightest (pixels, at tip).
    pub width: f32,

    /// Exponential decay constant for tail alpha: `alpha ∝ exp(-k * (1-age))`.
    /// Higher → tail vanishes faster. Typical range 3.0 – 8.0.
    pub decay_k: f32,

    /// Amount of desaturation applied at the tail end (0 = original colour,
    /// 1 = fully grayscale).
    pub tail_desaturate: f32,

    /// Base opacity multiplier applied on top of the decay curve. Lower
    /// values prevent opaque stacking when many segments overlap at high
    /// steps-per-frame. Range 0.0 – 1.0.
    pub base_alpha: f32,

    /// Soft-feather width as a fraction of the quad's transverse half-width.
    /// 0 → hard edges (old behaviour); 0.4 → ~40 % of the width is a smooth
    /// SDF falloff on each side. Tunable aesthetic.
    pub feather: f32,

    /// Extra core brightness near the leading tip (age ≈ 1). Scales the
    /// RGB vector; 1.0 disables the effect.
    pub core_boost: f32,
}

impl TrailStyle {
    pub const fn new(width: f32) -> Self {
        Self {
            width,
            decay_k: 6.0,
            tail_desaturate: 0.5,
            base_alpha: 1.0,
            feather: 0.35,
            core_boost: 1.25,
        }
    }
}

impl Default for TrailStyle {
    fn default() -> Self {
        TrailStylePreset::UniverseSandbox.style(1.5)
    }
}

// ── Presets ───────────────────────────────────────────────────────────────────

/// Named visual presets. Encapsulate aesthetic decisions so callers pick
/// an intent rather than twiddling six floats.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TrailStylePreset {
    /// Softer, glowing look. Lower base alpha and strong feathering so trails
    /// do not obscure the body at the leading tip and overlaps remain
    /// translucent.
    UniverseSandbox,
    /// Crisper, higher-contrast trail suitable for scientific presentation
    /// where the exact path matters more than atmosphere.
    Scientific,
    /// Thin, single-hue trail with no core boost — minimal visual noise for
    /// dense simulations.
    Minimal,
}

impl TrailStylePreset {
    /// Returns the concrete [`TrailStyle`] for this preset, using `width` as
    /// the user-facing size control (the one thing the UI exposes today).
    pub fn style(self, width: f32) -> TrailStyle {
        match self {
            Self::UniverseSandbox => TrailStyle {
                width,
                decay_k: 5.5,
                tail_desaturate: 0.55,
                base_alpha: 0.55,
                feather: 0.45,
                core_boost: 1.3,
            },
            Self::Scientific => TrailStyle {
                width,
                decay_k: 7.0,
                tail_desaturate: 0.25,
                base_alpha: 0.85,
                feather: 0.15,
                core_boost: 1.1,
            },
            Self::Minimal => TrailStyle {
                width: width * 0.75,
                decay_k: 8.0,
                tail_desaturate: 0.7,
                base_alpha: 0.45,
                feather: 0.3,
                core_boost: 1.0,
            },
        }
    }
}
