//! Scene lighting — domain types + GPU uniform packing.
//!
//! The renderer splits its lighting responsibilities in two:
//!
//! 1. **High-level state** — [`LightSpec`] and [`SceneLighting`] describe
//!    what the *scene* currently looks like: which bodies are luminous,
//!    where they sit in world space, how the shading should behave
//!    (ambient floor, terminator softening, distance falloff). These are
//!    the types the app layer deals with; they know nothing about bytes,
//!    binding slots, or alignment.
//!
//! 2. **GPU-packed state** — [`LightingUniform`] is the exact byte layout
//!    the WGSL shader expects. It's produced from a [`SceneLighting`] via
//!    [`LightingUniform::pack`], which handles sorting/clipping of
//!    lights and pre-squaring of distances so the fragment shader can run
//!    a tight loop without per-sample normalisation.
//!
//! Keeping the two separate means the app layer can evolve the shading
//! inputs (e.g. tint per light, coloured shadows) without touching
//! backend plumbing, and the GPU packing can change (16-byte stride
//! adjustments, field reordering) without rippling into callers.
//!
//! # 3D readiness
//!
//! All positions are [`[f32; 3]`]. 2D callers set `z = 0`; a future 3D
//! camera populates the full vector. Nothing in the shader cares which
//! it is.

use bytemuck::{Pod, Zeroable};

/// Maximum concurrent light sources the renderer supports.
///
/// Raising this requires a matching bump to the WGSL `array<PackedLight, N>`
/// dimension in [`crate::render::wgpu_backend`]. Four covers binary and
/// triple-star systems with headroom; going higher is cheap on the GPU
/// (still a few dot-products per fragment) but wastes uniform bandwidth
/// for the common 1-light case.
pub const MAX_LIGHTS: usize = 4;

// ── High-level domain ────────────────────────────────────────────────────────

/// A single active light source in world coordinates.
#[derive(Clone, Copy, Debug)]
pub struct LightSpec {
    /// World-space position. 2D callers use `z = 0`.
    pub world_pos: [f32; 3],
    /// Relative luminosity — weights this source's contribution against
    /// its siblings. A reference star sets `intensity = 1.0`; a binary
    /// companion at 30% of the primary's output sets `0.3`. The absolute
    /// scale is free (the renderer clamps per-fragment output), so
    /// callers can use "luminosity / max_luminosity" as a cheap default.
    pub intensity: f32,
}

/// Per-frame scene lighting configuration. Built by the app layer and
/// handed to the renderer via
/// [`WgpuBackend::set_scene_lighting`](crate::render::WgpuBackend::set_scene_lighting).
#[derive(Clone, Debug)]
pub struct SceneLighting {
    /// Active sources. More than [`MAX_LIGHTS`] are truncated after
    /// sorting by intensity (brightest wins) — rare in practice and
    /// physically acceptable: the faintest contributors to Lambert sum
    /// are also the ones the eye can't distinguish from noise.
    pub lights: Vec<LightSpec>,

    /// Multiplicative ambient floor in [0, 1]. The body shader computes
    ///
    /// ```text
    /// lit_factor = mix(ambient_floor, 1.0, saturate(diffuse_total))
    /// ```
    ///
    /// so `0.0` is pure Lambert (black terminator), `1.0` fully flat
    /// (no shading). Small values (~0.05) prevent distant / back-facing
    /// bodies from disappearing into the HDR black without washing out
    /// the day side. Multiplicative, not additive — we want a *floor* on
    /// the albedo, not extra brightness piled on top of lit surfaces.
    pub ambient_floor: f32,

    /// Reference distance at which a light contributes attenuation = 1.
    ///
    /// The falloff in the shader is
    ///
    /// ```text
    /// attenuation = r_ref² / (r² + falloff_bias²)
    /// ```
    ///
    /// so picking `r_ref = characteristic_distance` of the scene
    /// (e.g. 1 AU for a Solar-System view) keeps the "typical" body
    /// rendered at natural brightness regardless of how the user
    /// scales the simulation.
    pub r_ref: f32,

    /// Soft bias added to the attenuation denominator. Prevents the
    /// `1/r²` singularity when a body overlaps its light source and
    /// gives the artist a knob to soften the falloff near the primary.
    /// Units: same as world coordinates.
    pub falloff_bias: f32,

    /// Terminator-softening knob in [0, 1].
    ///
    /// * `0.0` — pure Lambert: `max(dot(n, L), 0)`. Hard cutoff.
    /// * `1.0` — full half-Lambert wrap: `((dot(n, L) + 1) / 2)²`.
    ///           No dark side, but physically plausible for diffuse
    ///           interreflection in atmospheric / dusty bodies.
    ///
    /// Intermediate values mix the two, giving a gentle terminator
    /// without erasing the phase look entirely.
    pub wrap: f32,
}

impl Default for SceneLighting {
    fn default() -> Self {
        Self {
            lights: Vec::new(),
            ambient_floor: 0.05,
            r_ref: 1.0,
            falloff_bias: 0.05,
            wrap: 0.25,
        }
    }
}

// ── GPU packing ──────────────────────────────────────────────────────────────

/// One light source in its GPU-uniform layout.
///
/// The WGSL mirror (`struct PackedLight`) matches this byte-for-byte:
/// `vec3<f32>` for the position (size 12, align 16 — trailing 4 bytes
/// swallowed by the struct's own align-16) followed by a scalar f32.
/// Total 16 bytes per entry; the array stride in the uniform matches.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct PackedLight {
    world_pos: [f32; 3],
    intensity: f32,
}

/// GPU-uniform layout consumed by the body fragment shader.
///
/// Field order and padding are deliberate — see the module-level
/// comment for the WGSL alignment constraints that shape them. In
/// particular the `_pad*` scalars bring the struct up to a 16-byte
/// multiple so `array<PackedLight, N>` at offset 0 keeps its natural
/// alignment when the buffer is reallocated.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct LightingUniform {
    lights: [PackedLight; MAX_LIGHTS],
    num_lights: u32,
    ambient_floor: f32,
    r_ref_sq: f32,
    bias_sq: f32,
    wrap: f32,
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
}

impl Default for LightingUniform {
    fn default() -> Self {
        Self {
            lights: [PackedLight { world_pos: [0.0; 3], intensity: 0.0 }; MAX_LIGHTS],
            num_lights: 0,
            ambient_floor: 0.0,
            r_ref_sq: 1.0,
            bias_sq: 1e-4,
            wrap: 0.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
        }
    }
}

impl LightingUniform {
    /// Packs high-level [`SceneLighting`] into the GPU-uniform shape.
    ///
    /// Lights are sorted by intensity descending; any beyond
    /// [`MAX_LIGHTS`] are dropped. Distances are pre-squared so the
    /// fragment shader runs a scalar divide per light rather than a
    /// `sqrt` + `mul`.
    pub fn pack(scene: &SceneLighting) -> Self {
        let mut out = Self::default();
        out.ambient_floor = scene.ambient_floor.clamp(0.0, 1.0);
        out.r_ref_sq = scene.r_ref.max(1e-6).powi(2);
        out.bias_sq = scene.falloff_bias.max(0.0).powi(2);
        out.wrap = scene.wrap.clamp(0.0, 1.0);

        // Sort a local copy so the caller's slice keeps insertion order
        // (callers may rely on it for UI display, logging, etc.).
        let mut sorted: Vec<LightSpec> = scene.lights.clone();
        sorted.sort_by(|a, b| {
            b.intensity
                .partial_cmp(&a.intensity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let n = sorted.len().min(MAX_LIGHTS);
        for (i, src) in sorted.iter().take(n).enumerate() {
            out.lights[i] = PackedLight {
                world_pos: src.world_pos,
                intensity: src.intensity.max(0.0),
            };
        }
        out.num_lights = n as u32;
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_size_is_multiple_of_16() {
        // WGSL requires uniform buffers to be a multiple of 16 bytes so
        // the trailing array element strides remain aligned. If this
        // assertion ever fails, re-check the _pad* fields at the end of
        // LightingUniform.
        assert_eq!(std::mem::size_of::<LightingUniform>() % 16, 0);
    }

    #[test]
    fn pack_sorts_lights_by_intensity_descending() {
        let scene = SceneLighting {
            lights: vec![
                LightSpec { world_pos: [1.0, 0.0, 0.0], intensity: 0.3 },
                LightSpec { world_pos: [2.0, 0.0, 0.0], intensity: 1.0 },
                LightSpec { world_pos: [3.0, 0.0, 0.0], intensity: 0.6 },
            ],
            ..Default::default()
        };
        let u = LightingUniform::pack(&scene);
        assert_eq!(u.num_lights, 3);
        // Brightest (intensity 1.0, pos x=2) must occupy slot 0.
        assert_eq!(u.lights[0].world_pos[0], 2.0);
        assert_eq!(u.lights[1].world_pos[0], 3.0);
        assert_eq!(u.lights[2].world_pos[0], 1.0);
    }

    #[test]
    fn pack_truncates_past_max_lights() {
        let mut lights = Vec::new();
        for i in 0..(MAX_LIGHTS + 3) {
            lights.push(LightSpec {
                world_pos: [i as f32, 0.0, 0.0],
                intensity: (i + 1) as f32, // increasing → last wins ordering
            });
        }
        let scene = SceneLighting { lights, ..Default::default() };
        let u = LightingUniform::pack(&scene);
        assert_eq!(u.num_lights as usize, MAX_LIGHTS);
        // Brightest (intensity = MAX_LIGHTS+3) must be slot 0.
        assert_eq!(u.lights[0].intensity, (MAX_LIGHTS + 3) as f32);
    }

    #[test]
    fn pack_pre_squares_distances() {
        let scene = SceneLighting {
            r_ref: 4.0,
            falloff_bias: 0.5,
            ..Default::default()
        };
        let u = LightingUniform::pack(&scene);
        assert!((u.r_ref_sq - 16.0).abs() < 1e-6);
        assert!((u.bias_sq - 0.25).abs() < 1e-6);
    }

    #[test]
    fn pack_clamps_ambient_floor_and_wrap() {
        let scene = SceneLighting {
            ambient_floor: 1.5,
            wrap: -0.3,
            ..Default::default()
        };
        let u = LightingUniform::pack(&scene);
        assert_eq!(u.ambient_floor, 1.0);
        assert_eq!(u.wrap, 0.0);
    }
}
