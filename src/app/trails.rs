use eframe::egui::{Color32, Painter, Pos2, Rect, Stroke};
use std::collections::VecDeque;

const MAX_TOTAL_SEGMENTS: usize = 60_000;

pub fn draw_trails(
    painter: &Painter,
    trails: &[VecDeque<(f64, f64)>],
    body_colors: &[Color32],
    center: Pos2,
    scale: f32,
    rect: Rect,
) {
    let n = trails.len();
    if n == 0 {
        return;
    }

    let segs_per_body = (MAX_TOTAL_SEGMENTS / n).clamp(20, 2000);

    for i in 0..n {
        let trail = &trails[i];
        let len = trail.len();
        if len < 2 {
            continue;
        }

        let base = if i < body_colors.len() {
            body_colors[i]
        } else {
            Color32::WHITE
        };

        let step = (len / segs_per_body).max(1);
        let sampled_len = (len / step).max(1);

        let mut prev: Option<(Pos2, f32)> = None;

        for (j, (tx, ty)) in trail.iter().step_by(step).enumerate() {
            let t = j as f32 / sampled_len as f32;
            let p = Pos2::new(center.x + *tx as f32 * scale, center.y + *ty as f32 * scale);

            if let Some((prev_p, prev_t)) = prev {
                if rect.contains(prev_p) || rect.contains(p) {
                    let mid_t = (prev_t + t) * 0.5;
                    // Fade in over the first 15% of the trail, fully opaque after that.
                    let alpha = ((mid_t / 0.15).min(1.0) * 160.0) as u8;
                    let width = 0.4 + mid_t * 0.8;
                    painter.line_segment(
                        [prev_p, p],
                        Stroke::new(
                            width,
                            Color32::from_rgba_premultiplied(base.r(), base.g(), base.b(), alpha),
                        ),
                    );
                }
            }

            prev = Some((p, t));
        }
    }
}
