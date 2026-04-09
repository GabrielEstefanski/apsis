use crate::app::render_params::{RenderParams, compute_render_radius};
use crate::app::theme::{ACCENT, BG, TEXT_DIM, body_radius, fmt_world, nice_grid_world};
use crate::app::ui::{SelectionForm, SemanticScaleMode, SimulationApp};
use crate::domain::body::Body;
use crate::render::{RenderBackend, TrailRenderer, WgpuBackend};
use crate::templates::instantiate_at;
use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Stroke, Vec2};

fn dim(c: Color32, f: f32) -> Color32 {
    Color32::from_rgba_premultiplied(
        (c.r() as f32 * f) as u8,
        (c.g() as f32 * f) as u8,
        (c.b() as f32 * f) as u8,
        c.a(),
    )
}

fn alpha(c: Color32, a: u8) -> Color32 {
    Color32::from_rgba_premultiplied(c.r(), c.g(), c.b(), a)
}

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default()
            .frame(egui::Frame::NONE.fill(BG))
            .show(ctx, |ui| {
                let rect = ui.max_rect();
                let center = rect.center() + self.offset;

                let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());

                if self.template_drag.is_some() {
                    self.handle_template_drag(ctx, rect, center);
                } else if self.place_mode {
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
                            b.sync_physical_properties();
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
                    // Drag selected body to reposition
                    if let Some(sel_idx) = self.selected_body {
                        if response.drag_started() {
                            self.dragging_body = Some(sel_idx);
                            let cursor_pos =
                                ctx.input(|i| i.pointer.hover_pos()).unwrap_or(Pos2::ZERO);
                            let wx = (cursor_pos.x - center.x) as f64 / self.scale as f64;
                            let wy = (cursor_pos.y - center.y) as f64 / self.scale as f64;
                            self.drag_start_world = Some((wx, wy));
                        }
                        if response.dragged() && self.dragging_body.is_some() {
                            if let Some(start) = self.drag_start_world {
                                let cur =
                                    ctx.input(|i| i.pointer.hover_pos()).unwrap_or(Pos2::ZERO);
                                let cur_wx = (cur.x - center.x) as f64 / self.scale as f64;
                                let cur_wy = (cur.y - center.y) as f64 / self.scale as f64;
                                let dx = cur_wx - start.0;
                                let dy = cur_wy - start.1;
                                let mut body = self.system.bodies()[sel_idx];
                                body.x += dx;
                                body.y += dy;
                                self.system.update_body(sel_idx, body);
                                self.drag_start_world = Some((cur_wx, cur_wy));
                            }
                        }
                        if response.drag_stopped() {
                            self.dragging_body = None;
                            self.drag_start_world = None;
                        }
                    }
                }

                let scroll = ctx.input(|i| i.smooth_scroll_delta.y);

                if scroll.abs() > 0.0 {
                    let zoom_factor = (1.0 + scroll * 0.001).clamp(0.9, 1.1);
                    self.scale = (self.scale * zoom_factor).clamp(0.5, 5000.0);
                }

                let painter = ui.painter();
                let mut backend = WgpuBackend::new(ui, rect);

                backend.begin();

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
                if self.show_trails {
                    TrailRenderer::submit(
                        ui,
                        rect,
                        self.system.trail_buf_mut(),
                        [center.x, center.y],
                        self.scale,
                    );
                }

                let bodies = self.system.bodies();
                let render_params = RenderParams {
                    world_scale: self.scale,
                    mode: self.semantic_scale_mode,
                    min_px: match self.semantic_scale_mode {
                        SemanticScaleMode::Physical => 0.0,
                        SemanticScaleMode::Comparative => 3.0,
                        SemanticScaleMode::Illustrative => 5.0,
                    },
                };

                // ── Bodies ───────────────────────────────────────────────── //
                let mut indices: Vec<usize> = (0..bodies.len()).collect();

                indices.sort_by(|&a, &b| {
                    bodies[a]
                        .physical_radius
                        .partial_cmp(&bodies[b].physical_radius)
                        .unwrap()
                });

                for &i in &indices {
                    let b = &bodies[i];

                    let [cr, cg, cb] = b.color;

                    let render_r = compute_render_radius(b.physical_radius, render_params);

                    let px = center.x + b.x as f32 * self.scale;
                    let py = center.y + b.y as f32 * self.scale;

                    if px < rect.left() - 50.0
                        || px > rect.right() + 50.0
                        || py < rect.top() - 50.0
                        || py > rect.bottom() + 50.0
                    {
                        continue;
                    }

                    backend.draw_circle([px, py], render_r, [cr, cg, cb]);

                    //painter.circle_stroke(pos, render_r, Stroke::new(outline, Color32::BLACK));

                    //painter.circle_stroke(pos, render_r + 2.0, Stroke::new(1.0, alpha(col, 60)));

                    // ── VELOCITY COLOR HINT (subtle) ─────────────────────────────── //
                    let speed = (b.vx * b.vx + b.vy * b.vy).sqrt();
                    if speed > 0.0 && render_r > 4.0 {
                        let t = (speed * 0.15).clamp(0.0, 1.0);
                        let alpha = (t * 120.0) as u8;

                        backend.draw_circle_stroke(
                            [px, py],
                            render_r + 1.5,
                            1.0,
                            [255, 120, 40, alpha],
                        );
                    }
                }

                // ── Selection ring ───────────────────────────────────────── //
                if let Some(sel) = self.selected_body {
                    if sel < bodies.len() {
                        let b = bodies[sel];

                        let px = center.x + b.x as f32 * self.scale;
                        let py = center.y + b.y as f32 * self.scale;

                        let mut r = compute_render_radius(b.physical_radius, render_params);

                        let r = r.max(6.0);

                        painter.circle_stroke(Pos2::new(px, py), r + 4.0, Stroke::new(1.0, ACCENT));
                    } else {
                        self.selected_body = None;
                        self.selection_form = None;
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

                backend.end();
            });
    }

    /// Called every frame while `template_drag` is Some.
    /// Renders a ghost preview at the cursor and commits the drop on mouse release.
    fn handle_template_drag(&mut self, ctx: &egui::Context, rect: Rect, center: Pos2) {
        ctx.set_cursor_icon(egui::CursorIcon::Grabbing);

        let hover = ctx.input(|i| i.pointer.hover_pos());

        // Ghost: draw template bodies semi-transparently under the cursor
        if let Some(cursor) = hover {
            if rect.contains(cursor) {
                let build_fn = self.template_drag.as_ref().unwrap();
                let template = build_fn();

                let total_mass: f64 = template.bodies.iter().map(|t| t.mass).sum();
                let (com_x, com_y) = template.bodies.iter().fold((0.0, 0.0), |(ax, ay), t| {
                    let [px, py] = t.position.unwrap_or([0.0, 0.0]);
                    (ax + t.mass * px / total_mass, ay + t.mass * py / total_mass)
                });

                let painter = ctx.layer_painter(egui::LayerId::new(
                    egui::Order::Tooltip,
                    egui::Id::new("template_ghost"),
                ));

                for body in &template.bodies {
                    let [bx, by] = body.position.unwrap_or([0.0, 0.0]);
                    let rel_x = (bx - com_x) * self.scale as f64;
                    let rel_y = (by - com_y) * self.scale as f64;
                    let screen = Pos2::new(cursor.x + rel_x as f32, cursor.y + rel_y as f32);
                    let r = (body.radius * self.scale as f64).max(4.0) as f32;
                    let [cr, cg, cb] = body.material.props().base_color;

                    painter.circle_filled(
                        screen,
                        r,
                        Color32::from_rgba_premultiplied(cr, cg, cb, 90),
                    );

                    painter.circle_stroke(
                        screen,
                        r,
                        Stroke::new(1.0, Color32::from_rgba_premultiplied(cr, cg, cb, 180)),
                    );
                }

                // Drop hint label
                painter.text(
                    cursor + Vec2::new(12.0, 12.0),
                    Align2::LEFT_TOP,
                    "release to place",
                    FontId::proportional(9.0),
                    Color32::from_rgba_premultiplied(200, 200, 210, 140),
                );
            }
        }

        // Drop: commit on mouse release
        if ctx.input(|i| i.pointer.primary_released()) {
            let build_fn = self.template_drag.take().unwrap();
            if let Some(cursor) = hover {
                if rect.contains(cursor) {
                    let wx = (cursor.x - center.x) as f64 / self.scale as f64;
                    let wy = (cursor.y - center.y) as f64 / self.scale as f64;
                    let template = build_fn();
                    for b in instantiate_at(&template, wx, wy) {
                        self.system.add_body(b);
                    }
                }
                // Released outside canvas → cancel (template_drag is already None)
            } else {
                // Pointer left the window → cancel
            }
        }
    }

    fn find_body_at(&self, cursor: Pos2, center: Pos2) -> Option<usize> {
        let bodies = self.system.bodies();

        let render_params = RenderParams {
            world_scale: self.scale,
            mode: self.semantic_scale_mode,
            min_px: match self.semantic_scale_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            },
        };

        for i in (0..bodies.len()).rev() {
            let b = &bodies[i];

            let px = center.x + b.x as f32 * self.scale;
            let py = center.y + b.y as f32 * self.scale;

            let mut r = compute_render_radius(b.physical_radius, render_params);

            let r = r.max(6.0);

            let dx = cursor.x - px;
            let dy = cursor.y - py;

            if dx * dx + dy * dy <= r * r {
                return Some(i);
            }
        }

        None
    }
}
