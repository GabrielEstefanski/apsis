use crate::app::theme::{ACCENT, BG, TEXT_DIM, body_radius, fmt_world, nice_grid_world};
use crate::app::trails::draw_trails;
use crate::app::ui::{SelectionForm, SimulationApp};
use crate::domain::body::{
    Body, default_moment_inertia, default_softening, radius_from_density_mass,
};
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Stroke};

/// Visual radius in pixels: R ∝ m^(1/3), minimum 2 px.
/// Independent of the physical body radius (which is calibration-dependent);
/// this purely encodes mass so that a body twice as massive appears visually
/// larger by a factor of ∛2 ≈ 1.26.
fn visual_r_adaptive(
    physical_radius: f64,
    min_r: f64,
    max_r: f64,
    min_px: f32,
    max_px: f32,
) -> f32 {
    if (max_r - min_r).abs() < 1e-12 {
        return (min_px + max_px) * 0.5;
    }

    let log_min = min_r.max(1e-30).ln();
    let log_max = max_r.max(1e-30).ln();
    let log_range = (log_max - log_min).max(1e-6);

    let t = ((physical_radius.max(1e-30).ln() - log_min) / log_range).clamp(0.0, 1.0);

    let t = t.powf(0.6);

    min_px + (max_px - min_px) * t as f32
}

/// Dim a Color32 by a linear factor in [0, 1].
fn dim(c: Color32, f: f32) -> Color32 {
    Color32::from_rgba_premultiplied(
        (c.r() as f32 * f) as u8,
        (c.g() as f32 * f) as u8,
        (c.b() as f32 * f) as u8,
        c.a(),
    )
}

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default()
            .frame(egui::Frame::none().fill(BG))
            .show(ctx, |ui| {
                let rect = ui.max_rect();
                let center = rect.center() + self.offset;

                let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());

                if self.place_mode {
                    if response.drag_started() {
                        self.place_drag_start = ctx.input(|i| i.pointer.press_origin());
                    }
                    if response.drag_stopped() {
                        if let Some(start) = self.place_drag_start.take() {
                            let end = ctx.input(|i| i.pointer.interact_pos()).unwrap_or(start);
                            let wx = (start.x - center.x) as f64 / self.scale as f64;
                            let wy = (start.y - center.y) as f64 / self.scale as f64;
                            let vx = (end.x - start.x) as f64 / self.scale as f64;
                            let vy = (end.y - start.y) as f64 / self.scale as f64;

                            let mut b = Body::new(
                                wx,
                                wy,
                                vx,
                                vy,
                                self.place_mass,
                                crate::domain::materials::Material::Rocky,
                            );

                            b.density = self.place_density;
                            let physical_radius = radius_from_density_mass(b.density, b.mass);
                            b.physical_radius = physical_radius;
                            b.radius = physical_radius;
                            b.softening = default_softening(b.mass).max(physical_radius * 2.0);

                            b.moment_inertia = default_moment_inertia(b.mass, physical_radius);
                            self.system.add_body(b);
                        }
                    }
                } else {
                    if response.dragged() {
                        self.offset += response.drag_delta();
                    }
                    if response.clicked() {
                        let cursor = ctx.input(|i| i.pointer.interact_pos());
                        if let Some(cursor) = cursor {
                            let hit = self.find_body_at(cursor, center);
                            match hit {
                                Some(idx) if self.selected_body != Some(idx) => {
                                    let body = self.system.bodies()[idx];
                                    self.selection_form = Some(SelectionForm::from_body(&body));
                                    self.selected_body = Some(idx);
                                }
                                Some(_) => {}
                                None => {
                                    self.selected_body = None;
                                    self.selection_form = None;
                                }
                            }
                        }
                    }
                }

                let scroll = ctx.input(|i| i.raw_scroll_delta.y);
                if scroll != 0.0 {
                    self.scale = (self.scale * (1.0 + scroll * 0.001)).clamp(0.5, 500.0);
                }

                let painter = ui.painter();

                // ── Grid ─────────────────────────────────────────────────── //
                if self.show_grid {
                    let grid_world = nice_grid_world(self.scale);
                    let grid_px = grid_world * self.scale;
                    let line_col = Color32::from_rgba_premultiplied(32, 32, 42, 130);
                    let axis_col = Color32::from_rgba_premultiplied(55, 55, 70, 200);

                    let first_x = ((rect.left() - center.x) / grid_px).ceil() * grid_px + center.x;
                    let mut gx = first_x;
                    while gx <= rect.right() + grid_px {
                        let is_axis = (gx - center.x).abs() < 1.0;
                        painter.line_segment(
                            [Pos2::new(gx, rect.top()), Pos2::new(gx, rect.bottom())],
                            Stroke::new(0.5, if is_axis { axis_col } else { line_col }),
                        );
                        gx += grid_px;
                    }
                    let first_y = ((rect.top() - center.y) / grid_px).ceil() * grid_px + center.y;
                    let mut gy = first_y;
                    while gy <= rect.bottom() + grid_px {
                        let is_axis = (gy - center.y).abs() < 1.0;
                        painter.line_segment(
                            [Pos2::new(rect.left(), gy), Pos2::new(rect.right(), gy)],
                            Stroke::new(0.5, if is_axis { axis_col } else { line_col }),
                        );
                        gy += grid_px;
                    }
                    painter.text(
                        Pos2::new(rect.left() + 10.0, rect.bottom() - 10.0),
                        Align2::LEFT_BOTTOM,
                        format!("grid {}", fmt_world(grid_world)),
                        FontId::proportional(9.0),
                        TEXT_DIM,
                    );
                }

                // ── COM crosshair ────────────────────────────────────────── //
                {
                    let m = self.system.metrics();
                    let cx = center.x + m.com_x as f32 * self.scale;
                    let cy = center.y + m.com_y as f32 * self.scale;
                    let s = 4.0;
                    let c = Color32::from_rgba_premultiplied(80, 160, 110, 140);
                    painter.line_segment(
                        [Pos2::new(cx - s, cy), Pos2::new(cx + s, cy)],
                        Stroke::new(1.0, c),
                    );
                    painter.line_segment(
                        [Pos2::new(cx, cy - s), Pos2::new(cx, cy + s)],
                        Stroke::new(1.0, c),
                    );
                    painter.text(
                        Pos2::new(cx + 6.0, cy - 7.0),
                        Align2::LEFT_CENTER,
                        "COM",
                        FontId::proportional(8.5),
                        c,
                    );
                }

                // ── Trails ───────────────────────────────────────────────── //
                // Trails use palette colours (body identity), not velocity colours,
                // so the orbital paths are always distinguishable.
                if self.show_trails {
                    let colors: Vec<Color32> = self
                        .system
                        .bodies()
                        .iter()
                        .enumerate()
                        .map(|(_i, b)| {
                            let [r, g, b_] = b.color;
                            Color32::from_rgb(r, g, b_)
                        })
                        .collect();
                    draw_trails(
                        painter,
                        self.system.trails(),
                        &colors,
                        center,
                        self.scale,
                        rect,
                    );
                }

                // ── Velocity → colour normalisation ──────────────────────── //
                let v_max = self
                    .system
                    .bodies()
                    .iter()
                    .map(|b| (b.vx * b.vx + b.vy * b.vy).sqrt())
                    .fold(0.0_f64, f64::max)
                    .max(1e-30);

                let bodies = self.system.bodies();

                let (min_r, max_r) =
                    bodies
                        .iter()
                        .fold((f64::INFINITY, 0.0_f64), |(min_r, max_r), b| {
                            (min_r.min(b.physical_radius), max_r.max(b.physical_radius))
                        });

                // Clone accelerations to avoid simultaneous borrow of self.system.
                let accs: Vec<(f64, f64)> = self.system.last_accelerations().to_vec();

                // ── Bodies ───────────────────────────────────────────────── //
                for (i, b) in self.system.bodies().iter().enumerate() {
                    let px = center.x + b.x as f32 * self.scale;
                    let py = center.y + b.y as f32 * self.scale;
                    let pos = Pos2::new(px, py);

                    let r = visual_r_adaptive(
                        b.physical_radius,
                        min_r,
                        max_r,
                        3.0,  // mínimo visível
                        40.0, // máximo aceitável
                    );

                    // let speed = (b.vx * b.vx + b.vy * b.vy).sqrt();
                    let [cr, cg, cb] = b.color;
                    let col = Color32::from_rgb(cr, cg, cb);

                    painter.circle_filled(pos, r, col);
                    painter.circle_stroke(
                        pos,
                        r,
                        Stroke::new(0.5, Color32::from_rgba_premultiplied(255, 255, 255, 18)),
                    );

                    // ── Rotation spoke ────────────────────────────────────── //
                    // A line from centre to the surface at the current rotation
                    // angle makes spin visible.  Only drawn when the body is
                    // large enough and actually spinning.
                    if r > 3.5 && b.omega_z.abs() > 1e-6 {
                        let angle = self.body_angles.get(i).copied().unwrap_or(0.0) as f32;
                        let spoke = Pos2::new(px + angle.cos() * r, py + angle.sin() * r);
                        painter.line_segment([pos, spoke], Stroke::new(1.5, dim(col, 0.55)));
                        // Dot at the surface point for visibility
                        painter.circle_filled(spoke, 1.5, dim(col, 0.8));
                    }

                    // ── Velocity vectors ──────────────────────────────────── //
                    if self.show_vectors {
                        let vscale = self.scale * 0.3;
                        let tip = Pos2::new(px + b.vx as f32 * vscale, py + b.vy as f32 * vscale);
                        let vcol = dim(col, 0.75);
                        painter.line_segment([pos, tip], Stroke::new(1.0, vcol));
                        // Arrow head
                        let dx = tip.x - px;
                        let dy = tip.y - py;
                        let len = (dx * dx + dy * dy).sqrt().max(1e-5);
                        if len > 4.0 {
                            let (nx, ny) = (dx / len, dy / len);
                            let hs = 4.0_f32;
                            painter.line_segment(
                                [
                                    tip,
                                    Pos2::new(
                                        tip.x - nx * hs - ny * hs * 0.5,
                                        tip.y - ny * hs + nx * hs * 0.5,
                                    ),
                                ],
                                Stroke::new(0.8, vcol),
                            );
                            painter.line_segment(
                                [
                                    tip,
                                    Pos2::new(
                                        tip.x - nx * hs + ny * hs * 0.5,
                                        tip.y - ny * hs - nx * hs * 0.5,
                                    ),
                                ],
                                Stroke::new(0.8, vcol),
                            );
                        }
                    }

                    // ── Force / acceleration vectors ──────────────────────── //
                    if self.show_force_vectors {
                        if let Some(&(ax, ay)) = accs.get(i) {
                            let fscale = self.scale * 0.5;
                            let ftip = Pos2::new(px + ax as f32 * fscale, py + ay as f32 * fscale);
                            let fcol = Color32::from_rgba_premultiplied(220, 100, 40, 210);
                            painter.line_segment([pos, ftip], Stroke::new(1.0, fcol));
                            // Arrow head
                            let dx = ftip.x - px;
                            let dy = ftip.y - py;
                            let len = (dx * dx + dy * dy).sqrt().max(1e-5);
                            if len > 3.0 {
                                let (nx, ny) = (dx / len, dy / len);
                                let hs = 4.0_f32;
                                painter.line_segment(
                                    [
                                        ftip,
                                        Pos2::new(
                                            ftip.x - nx * hs - ny * hs * 0.5,
                                            ftip.y - ny * hs + nx * hs * 0.5,
                                        ),
                                    ],
                                    Stroke::new(0.8, fcol),
                                );
                                painter.line_segment(
                                    [
                                        ftip,
                                        Pos2::new(
                                            ftip.x - nx * hs + ny * hs * 0.5,
                                            ftip.y - ny * hs - nx * hs * 0.5,
                                        ),
                                    ],
                                    Stroke::new(0.8, fcol),
                                );
                            }
                        }
                    }

                    // ── Mass label ────────────────────────────────────────── //
                    if r > 6.0 {
                        painter.text(
                            Pos2::new(px, py + r + 8.0),
                            Align2::CENTER_CENTER,
                            format!("{:.1}", b.mass),
                            FontId::proportional(9.0),
                            Color32::from_rgba_premultiplied(col.r(), col.g(), col.b(), 120),
                        );
                    }
                }

                // ── Selection ring ───────────────────────────────────────── //
                if let Some(sel) = self.selected_body {
                    if sel < self.system.bodies().len() {
                        let b = self.system.bodies()[sel];
                        let px = center.x + b.x as f32 * self.scale;
                        let py = center.y + b.y as f32 * self.scale;
                        let r = visual_r_adaptive(b.physical_radius, min_r, max_r, 3.0, 40.0);
                        painter.circle_stroke(Pos2::new(px, py), r + 4.0, Stroke::new(1.0, ACCENT));
                    } else {
                        self.selected_body = None;
                        self.selection_form = None;
                    }
                }

                // ── Impact visual effects ─────────────────────────────────── //
                for effect in &self.impact_effects {
                    let sx = center.x + effect.world_x as f32 * self.scale;
                    let sy = center.y + effect.world_y as f32 * self.scale;
                    let t = effect.age;
                    let fade = (1.0 - t).max(0.0);
                    let alpha = (fade.powf(0.5) * 255.0) as u8;

                    // Expanding flash ring
                    let ring_r = t * 28.0;
                    if ring_r > 0.5 {
                        painter.circle_stroke(
                            Pos2::new(sx, sy),
                            ring_r,
                            Stroke::new(
                                (1.8 * fade).max(0.3),
                                Color32::from_rgba_premultiplied(255, 210, 80, alpha),
                            ),
                        );
                    }

                    // Inner bright flash (fades quickly in first 30% of lifetime)
                    if t < 0.3 {
                        let i_fade = (0.3 - t) / 0.3;
                        painter.circle_filled(
                            Pos2::new(sx, sy),
                            9.0 * i_fade,
                            Color32::from_rgba_premultiplied(255, 240, 180, (i_fade * 220.0) as u8),
                        );
                    }

                    // Collision normal line (shown when "nrm" toggle is on)
                    if self.show_impact_normals {
                        let len = 32.0 * fade;
                        let na = Pos2::new(sx + effect.nx * len, sy + effect.ny * len);
                        let nb = Pos2::new(sx - effect.nx * len, sy - effect.ny * len);
                        painter.line_segment(
                            [na, nb],
                            Stroke::new(
                                1.0,
                                Color32::from_rgba_premultiplied(100, 200, 255, alpha),
                            ),
                        );
                    }

                    // Burst particles
                    for p in &effect.particles {
                        let px_p = center.x + p[0] as f32 * self.scale;
                        let py_p = center.y + p[1] as f32 * self.scale;
                        let p_alpha = (fade.powf(1.5) * 220.0) as u8;
                        painter.circle_filled(
                            Pos2::new(px_p, py_p),
                            1.5,
                            Color32::from_rgba_premultiplied(255, 175, 55, p_alpha),
                        );
                    }
                }

                // ── Place-mode overlay ───────────────────────────────────── //
                if self.place_mode {
                    if let Some(start) = self.place_drag_start {
                        let current = ctx.input(|i| i.pointer.hover_pos()).unwrap_or(start);
                        let r = body_radius(self.place_mass);
                        painter.circle_stroke(start, r, Stroke::new(1.0, ACCENT));
                        let delta = current - start;
                        if delta.length() > 4.0 {
                            painter.line_segment([start, current], Stroke::new(1.0, ACCENT));
                            painter.circle_filled(current, 2.0, ACCENT);
                        }
                    } else if let Some(cursor) = ctx.input(|i| i.pointer.hover_pos()) {
                        if rect.contains(cursor) {
                            painter.circle_stroke(
                                cursor,
                                body_radius(self.place_mass),
                                Stroke::new(
                                    0.5,
                                    Color32::from_rgba_premultiplied(200, 200, 210, 60),
                                ),
                            );
                        }
                    }
                    painter.text(
                        Pos2::new(rect.right() - 10.0, rect.bottom() - 10.0),
                        Align2::RIGHT_BOTTOM,
                        "click+drag = place body",
                        FontId::proportional(9.5),
                        TEXT_DIM,
                    );
                }
            });
    }

    fn find_body_at(&self, cursor: Pos2, center: Pos2) -> Option<usize> {
        let bodies = self.system.bodies();
        let (min_r, max_r) = bodies
            .iter()
            .fold((f64::INFINITY, 0.0_f64), |(min_r, max_r), b| {
                (min_r.min(b.physical_radius), max_r.max(b.physical_radius))
            });

        for i in (0..bodies.len()).rev() {
            let b = &bodies[i];
            let px = center.x + b.x as f32 * self.scale;
            let py = center.y + b.y as f32 * self.scale;
            let r = visual_r_adaptive(b.physical_radius, min_r, max_r, 3.0, 40.0).max(6.0);
            let dx = cursor.x - px;
            let dy = cursor.y - py;
            if dx * dx + dy * dy <= r * r {
                return Some(i);
            }
        }
        None
    }
}
