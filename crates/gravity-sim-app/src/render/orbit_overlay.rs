//! Predicted Keplerian orbit rendering.
//!
//! Projects a sampled orbit (world coordinates, 3D in general) onto the
//! current viewport and submits it as a sequence of line segments through
//! [`WgpuBackend::draw_line_segment`]. This module is a pure consumer of
//! the sampling produced by
//! [`OrbitalElements::sample_orbit`](gravity_sim_core::physics::orbital::OrbitalElements::sample_orbit);
//! it owns no physics state and holds no buffers of its own.
//!
//! # Why a separate module
//!
//! The overlay sits between `physics/orbital` (pure math) and
//! `render/wgpu_backend` (pure GPU submission). Putting the projection +
//! style logic in its own module means:
//!
//! * `physics/orbital` stays math-only — it never learns about pixels.
//! * `render/wgpu_backend` stays primitive — it never learns about orbits.
//! * The canvas layer only has to resolve elements + call one function.
//!
//! # 3D readiness
//!
//! The projection is passed as a closure so the same code path serves
//! the current 2D viewport and a future 3D camera. 2D callers ignore
//! the z-axis; a later 3D camera will project (x, y, z) through a view
//! matrix. Nothing else in the overlay changes.

use crate::render::wgpu_backend::WgpuBackend;
use gravity_sim_core::physics::orbital::OrbitalElements;

/// Visual style for one predicted-orbit polyline.
#[derive(Debug, Clone, Copy)]
pub struct OrbitOverlayStyle {
    /// RGBA (0-255). Alpha controls how prominent the prediction is
    /// relative to the live trail.
    pub color: [u8; 4],

    /// Line width in screen pixels.
    pub width_px: f32,
}

impl OrbitOverlayStyle {
    /// Default for the *selected-body* overlay: dim cyan, thin.
    ///
    /// The colour is deliberately neutral (not a material colour) so the
    /// predicted track reads as an annotation — a hint about where the
    /// body is heading — rather than competing with the body itself.
    pub const fn selected_default() -> Self {
        Self { color: [140, 210, 255, 205], width_px: 1.2 }
    }

    /// Default for the *all-bodies* overlay: same hue, fainter alpha.
    ///
    /// Used when the user enables "Orbit ellipses" globally and every
    /// bound body draws its own predicted track. The lower alpha lets
    /// many overlapping orbits coexist without overwhelming the scene
    /// or drowning out the brighter selected-body overlay on top.
    pub const fn background_default() -> Self {
        Self { color: [140, 210, 255, 95], width_px: 1.0 }
    }
}

impl Default for OrbitOverlayStyle {
    fn default() -> Self {
        Self::selected_default()
    }
}

/// Projects `points` (world coordinates, 3D) through `world_to_screen`
/// and submits them to `backend` as a polyline of line segments.
///
/// * `points` is the output of
///   [`sample_orbit`](gravity_sim_core::physics::orbital::OrbitalElements::sample_orbit).
/// * `world_to_screen` is caller-supplied: 2D mode drops the z-axis,
///   a future 3D mode applies its camera view matrix.
///
/// No-op when fewer than two points are supplied.
/// Draws the periapsis and (when defined) apoapsis markers for one orbit.
///
/// Periapsis is rendered as a small filled disk, apoapsis as a hollow
/// ring — both in the same hue as `style` so they read as annotations
/// attached to the polyline, not independent objects. Hyperbolic and
/// parabolic orbits have no apoapsis; only the periapsis marker is
/// emitted for them.
///
/// Marker sizes scale with `style.width_px` so the background overlay's
/// markers stay subordinate to the selected-orbit overlay's — same
/// visual hierarchy the polylines already follow.
pub fn draw_orbit_apsides<F>(
    backend: &mut WgpuBackend,
    el: &OrbitalElements,
    primary_pos: [f64; 3],
    mut world_to_screen: F,
    style: &OrbitOverlayStyle,
) where
    F: FnMut([f64; 3]) -> [f32; 2],
{
    // Periapsis: small filled disk. `draw_circle_stroke` with
    // width = 2·radius yields inner = 0, i.e. a solid disk.
    let r_peri_px = (style.width_px * 3.0).max(3.5);
    if let Some(peri_world) = el.periapsis_world(primary_pos) {
        let p = world_to_screen(peri_world);
        backend.draw_circle_stroke(p, r_peri_px, r_peri_px * 2.0, style.color);
    }
    // Apoapsis: hollow ring, slightly larger than periapsis so the two
    // read as distinct even when the orbit is nearly edge-on.
    let r_apo_px = (style.width_px * 3.8).max(4.5);
    if let Some(apo_world) = el.apoapsis_world(primary_pos) {
        let p = world_to_screen(apo_world);
        backend.draw_circle_stroke(p, r_apo_px, 1.5_f32.max(style.width_px), style.color);
    }
}

pub fn draw_orbit_polyline<F>(
    backend: &mut WgpuBackend,
    points: &[[f64; 3]],
    mut world_to_screen: F,
    style: &OrbitOverlayStyle,
) where
    F: FnMut([f64; 3]) -> [f32; 2],
{
    if points.len() < 2 {
        return;
    }
    let mut prev = world_to_screen(points[0]);
    for pt in &points[1..] {
        let cur = world_to_screen(*pt);
        backend.draw_line_segment(prev, cur, style.width_px, style.color);
        prev = cur;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn style_default_matches_selected_default() {
        let a = OrbitOverlayStyle::default();
        let b = OrbitOverlayStyle::selected_default();
        assert_eq!(a.color, b.color);
        assert_eq!(a.width_px, b.width_px);
    }

    #[test]
    fn background_is_fainter_than_selected() {
        // Invariant: the global "all bodies" overlay must not outshine
        // the focused selected-body overlay, otherwise the user loses
        // the annotation they're actively looking at.
        let sel = OrbitOverlayStyle::selected_default();
        let bg = OrbitOverlayStyle::background_default();
        assert!(bg.color[3] < sel.color[3], "background alpha must be lower");
        assert!(bg.width_px <= sel.width_px, "background width must not exceed selected");
        // Same hue family — only alpha/width differ.
        assert_eq!(bg.color[0..3], sel.color[0..3]);
    }
}
