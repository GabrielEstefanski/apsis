//! Camera-relative coordinates for the render path.
//!
//! Each frame the renderer picks a `render_origin` (the camera's eye)
//! and every world position handed to the GPU is shifted by it before
//! the `f64 → f32` cast. The cast then operates on a small magnitude
//! near zero, so the `f32` mantissa retains full precision regardless
//! of where the camera is in the absolute world — the same technique
//! KSP, Outerra, and the NASA Eyes simulators use to keep AU-scale
//! scenes free of `f32` jitter.
//!
//! [`RenderRelativeVec3`] is a newtype around `glam::Vec3` so the
//! type system can distinguish render-frame positions from raw
//! direction vectors. Helper paths that upload geometry to the GPU
//! consume this type; passing an absolute `Vec3` is a compile error.

use glam::{DVec3, Vec3};

/// Position in the render frame.
///
/// Constructed from an absolute world position (`DVec3`) by
/// subtracting the current `render_origin` in `f64` and casting the
/// small difference to `f32`. The wrapped value is `world - origin`
/// — magnitude is distance from the camera, not from the absolute
/// origin.
#[derive(Clone, Copy, Debug, Default)]
pub struct RenderRelativeVec3(pub Vec3);

impl RenderRelativeVec3 {
    /// `f64` subtraction → `f32` cast. The cast happens on a small
    /// magnitude, so the result keeps full `f32` precision down to
    /// sub-metre scale at any solar-system distance.
    #[inline]
    pub fn from_world(world: DVec3, origin: DVec3) -> Self {
        let rel = world - origin;
        Self(Vec3::new(rel.x as f32, rel.y as f32, rel.z as f32))
    }

    #[inline]
    pub fn as_vec3(self) -> Vec3 {
        self.0
    }

    #[inline]
    pub fn as_array(self) -> [f32; 3] {
        self.0.to_array()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_at_origin() {
        let rel = RenderRelativeVec3::from_world(DVec3::new(1.0, 2.0, 3.0), DVec3::ZERO);
        assert_eq!(rel.as_vec3(), Vec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn precision_preserved_at_au_scale() {
        // Body 1 km from the camera at 9.4 AU absolute world position.
        // Without the f64 subtraction first, this collapses to noise
        // in the f32 cast.
        let au = 1.495_978_707e11_f64; // metres per AU, scaled here in metres for clarity
        let camera = DVec3::new(9.4 * au, 0.0, 0.0);
        let body_offset_metres = 1000.0_f64;
        let body = DVec3::new(9.4 * au + body_offset_metres, 0.0, 0.0);

        let rel = RenderRelativeVec3::from_world(body, camera);
        let recovered_metres = rel.as_vec3().x as f64;

        let err = (recovered_metres - body_offset_metres).abs();
        assert!(err < 1e-3, "expected sub-millimetre precision at AU scale, got error {err} m");
    }

    #[test]
    fn naive_f32_cast_loses_precision_at_au_scale() {
        // Sanity check that the problem we're solving is real: doing
        // the subtraction in `f32` after the cast destroys the signal.
        let au = 1.495_978_707e11_f32;
        let camera = au * 9.4;
        let body = au * 9.4 + 1000.0;

        let naive_diff = body - camera;
        let err = (naive_diff as f64 - 1000.0).abs();
        assert!(
            err > 100.0,
            "expected catastrophic precision loss in naive f32 path; got error {err} m"
        );
    }
}
