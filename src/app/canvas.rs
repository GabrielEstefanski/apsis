use crate::app::theme::{
    body_color, body_radius, fmt_world, nice_grid_world, BG, TEXT_DIM, ACCENT,
};
use crate::app::trails::draw_trails;
use crate::app::ui::SimulationApp;
use crate::domain::body::Body;
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Stroke};

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default()
            .frame(egui::Frame::none().fill(BG))
            .show(ctx, |ui| {
                let rect = ui.max_rect();
                let center = rect.center() + self.offset;

                let response = ui.allocate_rect(rect, egui::Sense::drag());

                if self.place_mode {
                    if response.drag_started() {
                        self.place_drag_start = ctx.input(|i| i.pointer.press_origin());
                    }
                    if response.drag_released() {
                        if let Some(start) = self.place_drag_start.take() {
                            let end =
                                ctx.input(|i| i.pointer.interact_pos()).unwrap_or(start);
                            let wx = (start.x - center.x) as f64 / self.scale as f64;
                            let wy = (start.y - center.y) as f64 / self.scale as f64;
                            let vx = (end.x - start.x) as f64 / self.scale as f64;
                            let vy = (end.y - start.y) as f64 / self.scale as f64;
                            self.system.add_body(Body {
                                x: wx,
                                y: wy,
                                vx,
                                vy,
                                mass: self.place_mass,
                            });
                        }
                    }
                } else if response.dragged() {
                    self.offset += response.drag_delta();
                }

                let scroll = ctx.input(|i| i.raw_scroll_delta.y);
                if scroll != 0.0 {
                    self.scale =
                        (self.scale * (1.0 + scroll * 0.001)).clamp(0.5, 500.0);
                }

                let painter = ui.painter();

                if self.show_grid {
                    let grid_world = nice_grid_world(self.scale);
                    let grid_px = grid_world * self.scale;
                    let line_col = Color32::from_rgba_premultiplied(32, 32, 42, 130);
                    let axis_col = Color32::from_rgba_premultiplied(55, 55, 70, 200);

                    let first_x =
                        ((rect.left() - center.x) / grid_px).ceil() * grid_px + center.x;
                    let mut gx = first_x;
                    while gx <= rect.right() + grid_px {
                        let is_axis = (gx - center.x).abs() < 1.0;
                        painter.line_segment(
                            [Pos2::new(gx, rect.top()), Pos2::new(gx, rect.bottom())],
                            Stroke::new(0.5, if is_axis { axis_col } else { line_col }),
                        );
                        gx += grid_px;
                    }

                    let first_y =
                        ((rect.top() - center.y) / grid_px).ceil() * grid_px + center.y;
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

                if self.show_trails {
                    let colors: Vec<Color32> = self
                        .system
                        .bodies()
                        .iter()
                        .enumerate()
                        .map(|(i, b)| body_color(i, b.mass))
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

                for (i, b) in self.system.bodies().iter().enumerate() {
                    let px = center.x + b.x as f32 * self.scale;
                    let py = center.y + b.y as f32 * self.scale;
                    let pos = Pos2::new(px, py);
                    let r = body_radius(b.mass);
                    let col = body_color(i, b.mass);

                    painter.circle_filled(pos, r, col);
                    painter.circle_stroke(
                        pos,
                        r,
                        Stroke::new(0.5, Color32::from_rgba_premultiplied(255, 255, 255, 20)),
                    );

                    if self.show_vectors {
                        let vscale = self.scale * 0.3;
                        let tip = Pos2::new(
                            px + b.vx as f32 * vscale,
                            py + b.vy as f32 * vscale,
                        );
                        painter.line_segment(
                            [pos, tip],
                            Stroke::new(
                                0.8,
                                Color32::from_rgba_premultiplied(
                                    col.r(),
                                    col.g(),
                                    col.b(),
                                    150,
                                ),
                            ),
                        );
                        painter.circle_filled(tip, 1.5, col);
                    }

                    if r > 6.0 {
                        painter.text(
                            Pos2::new(px, py + r + 8.0),
                            Align2::CENTER_CENTER,
                            format!("{:.1}", b.mass),
                            FontId::proportional(9.0),
                            Color32::from_rgba_premultiplied(col.r(), col.g(), col.b(), 140),
                        );
                    }
                }

                if self.place_mode {
                    if let Some(start) = self.place_drag_start {
                        let current =
                            ctx.input(|i| i.pointer.hover_pos()).unwrap_or(start);
                        let r = body_radius(self.place_mass);
                        painter.circle_stroke(start, r, Stroke::new(1.0, ACCENT));
                        let delta = current - start;
                        if delta.length() > 4.0 {
                            painter.line_segment(
                                [start, current],
                                Stroke::new(1.0, ACCENT),
                            );
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
}
