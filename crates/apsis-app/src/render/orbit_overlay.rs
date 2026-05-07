//! Predicted Keplerian orbit rendering.
//!
//! Pure consumer of `OrbitalElements::sample_orbit`. Caller supplies a
//! 3D-camera projection closure that returns `(screen, w)` per world
//! point, with `w` the clip-space depth (positive in front of the
//! camera). Segments crossing the near plane are linearly interpolated
//! in world space until `w ≈ NEAR_BIAS` so the polyline lands exactly
//! on the camera horizon instead of leaving a visible gap.

use crate::render::wgpu_backend::WgpuBackend;
use apsis::physics::orbital::OrbitalElements;

/// Minimum positive clip-space `w` we accept for projection. Slightly
/// past the camera's near plane so the resulting screen position stays
/// numerically stable after the perspective divide.
const NEAR_BIAS: f32 = 1.05e-3;

#[derive(Debug, Clone, Copy)]
pub struct OrbitOverlayStyle {
    pub color: [u8; 4],
    pub width_px: f32,
}

impl OrbitOverlayStyle {
    pub const fn selected_default() -> Self {
        Self { color: [180, 220, 255, 235], width_px: 1.6 }
    }

    pub const fn background_default() -> Self {
        Self { color: [180, 220, 255, 110], width_px: 1.0 }
    }
}

impl Default for OrbitOverlayStyle {
    fn default() -> Self {
        Self::selected_default()
    }
}

#[inline]
fn lerp_world(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

/// Submit a polyline. `tint` overrides the style's RGB; `directional_anchor`,
/// when present, biases alpha toward that sample (60–100% of style alpha).
pub fn draw_orbit_polyline<F>(
    backend: &mut WgpuBackend,
    points: &[[f64; 3]],
    mut world_project: F,
    style: &OrbitOverlayStyle,
    tint: Option<[u8; 3]>,
    directional_anchor: Option<usize>,
) where
    F: FnMut([f64; 3]) -> ([f32; 2], f32),
{
    if points.len() < 2 {
        return;
    }

    let rgb = tint.unwrap_or_else(|| [style.color[0], style.color[1], style.color[2]]);
    let base_alpha = style.color[3] as f32;

    // Floor at 60% so even the dimmest part stays legible against HDR
    // bloom; delta of 40% is above the alpha-perception threshold without
    // reading as a particle effect.
    const ALPHA_FLOOR: f32 = 0.60;

    let n = points.len();
    let alpha_at = |i: usize| -> u8 {
        let scale = match directional_anchor {
            Some(anchor) => {
                let raw = i.abs_diff(anchor);
                let dist = raw.min(n - raw);
                let half = (n / 2).max(1) as f32;
                let t = (dist as f32 / half).clamp(0.0, 1.0);
                ALPHA_FLOOR + (1.0 - ALPHA_FLOOR) * (1.0 - t)
            },
            None => 1.0,
        };
        (base_alpha * scale).round().clamp(0.0, 255.0) as u8
    };

    let (mut s_prev, mut w_prev) = world_project(points[0]);
    for i in 1..points.len() {
        let p_prev = points[i - 1];
        let p_cur = points[i];
        let (s_cur, w_cur) = world_project(p_cur);

        let in_prev = w_prev > NEAR_BIAS;
        let in_cur = w_cur > NEAR_BIAS;
        let a = ((alpha_at(i - 1) as u16 + alpha_at(i) as u16) / 2) as u8;
        let color = [rgb[0], rgb[1], rgb[2], a];

        if in_prev && in_cur {
            backend.draw_line_segment(s_prev, s_cur, style.width_px, color);
        } else if in_prev || in_cur {
            // One endpoint sits behind the near plane. Find the world
            // point where w == NEAR_BIAS and draw the partial segment
            // up to it instead of dropping the whole edge.
            let denom = w_cur - w_prev;
            if denom.abs() > 1e-12 {
                let t = ((NEAR_BIAS - w_prev) / denom) as f64;
                let t = t.clamp(0.0, 1.0);
                let p_clip = lerp_world(p_prev, p_cur, t);
                let (s_clip, _) = world_project(p_clip);
                if in_prev {
                    backend.draw_line_segment(s_prev, s_clip, style.width_px, color);
                } else {
                    backend.draw_line_segment(s_clip, s_cur, style.width_px, color);
                }
            }
        }

        s_prev = s_cur;
        w_prev = w_cur;
    }
}

/// Polyline with a wider, faint underlay. Halo stays inside the style's
/// alpha budget so it doesn't push extra energy into the bloom pass.
pub fn draw_orbit_polyline_with_halo<F>(
    backend: &mut WgpuBackend,
    points: &[[f64; 3]],
    mut world_project: F,
    style: &OrbitOverlayStyle,
    directional_anchor: Option<usize>,
) where
    F: FnMut([f64; 3]) -> ([f32; 2], f32),
{
    let halo = OrbitOverlayStyle {
        color: [
            style.color[0],
            style.color[1],
            style.color[2],
            ((style.color[3] as f32) * 0.30).round().clamp(0.0, 255.0) as u8,
        ],
        width_px: style.width_px * 2.6,
    };
    draw_orbit_polyline(backend, points, &mut world_project, &halo, None, None);
    draw_orbit_polyline(backend, points, world_project, style, None, directional_anchor);
}

pub fn draw_orbit_apsides<F>(
    backend: &mut WgpuBackend,
    el: &OrbitalElements,
    primary_pos: [f64; 3],
    mut world_project: F,
    style: &OrbitOverlayStyle,
) where
    F: FnMut([f64; 3]) -> ([f32; 2], f32),
{
    let r_peri_px = (style.width_px * 3.0).max(3.5);
    if let Some(peri_world) = el.periapsis_world(primary_pos) {
        let (p, w) = world_project(peri_world);
        if w > NEAR_BIAS {
            backend.draw_circle_stroke(p, r_peri_px, r_peri_px * 2.0, style.color);
        }
    }
    let r_apo_px = (style.width_px * 3.8).max(4.5);
    if let Some(apo_world) = el.apoapsis_world(primary_pos) {
        let (p, w) = world_project(apo_world);
        if w > NEAR_BIAS {
            backend.draw_circle_stroke(p, r_apo_px, 1.5_f32.max(style.width_px), style.color);
        }
    }
}

/// Index of the sample in `points` closest to `body_pos` in 3D world coords.
pub fn closest_sample_index(points: &[[f64; 3]], body_pos: [f64; 3]) -> Option<usize> {
    if points.is_empty() {
        return None;
    }
    let mut best_idx = 0usize;
    let mut best_d2 = f64::INFINITY;
    for (i, p) in points.iter().enumerate() {
        let dx = p[0] - body_pos[0];
        let dy = p[1] - body_pos[1];
        let dz = p[2] - body_pos[2];
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < best_d2 {
            best_d2 = d2;
            best_idx = i;
        }
    }
    Some(best_idx)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn background_is_fainter_than_selected() {
        let sel = OrbitOverlayStyle::selected_default();
        let bg = OrbitOverlayStyle::background_default();
        assert!(bg.color[3] < sel.color[3]);
        assert!(bg.width_px <= sel.width_px);
    }

    #[test]
    fn closest_sample_index_picks_nearest() {
        let points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        assert_eq!(closest_sample_index(&points, [1.4, 0.0, 0.0]), Some(1));
        assert_eq!(closest_sample_index(&points, [2.6, 0.0, 0.0]), Some(3));
    }

    #[test]
    fn closest_sample_index_uses_3d_distance() {
        let points = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!(closest_sample_index(&points, [0.0, 0.0, 0.9]), Some(2));
    }

    #[test]
    fn closest_sample_index_empty_returns_none() {
        assert_eq!(closest_sample_index(&[], [0.0, 0.0, 0.0]), None);
    }
}
