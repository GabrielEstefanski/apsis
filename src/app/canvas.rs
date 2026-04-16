use crate::app::render_params::{RenderParams, compute_render_radius};
use crate::app::ui::{SelectionForm, SemanticScaleMode, SimulationApp, UndoRecord};
use crate::render::CallbackFn;
use crate::render::wgpu_backend::LightSource;
use crate::templates::instantiate_at;
use eframe::egui::{self, Color32, FontId, Pos2, Stroke};
use eframe::egui_wgpu;

// ── Tunables ──────────────────────────────────────────────────────────────────

/// Minimum hit-test radius in pixels. Makes small/distant bodies easier to click.
const MIN_HIT_PX: f32 = 10.0;

/// Maximum distance in pixels from the body centre to place its name label.
/// Prevents labels drifting far from large bodies at high zoom.
const MAX_LABEL_OFFSET_PX: f32 = 48.0;

/// Label font size.
const LABEL_FONT_SIZE: f32 = 10.5;

/// Selection ring: base gap outside the body disc.
const RING_GAP: f32 = 5.0;

/// Camera pan animation: fraction of remaining distance applied each frame.
const CAM_LERP: f32 = 0.16;
const FOLLOW_LERP: f32 = 0.32;

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ui: &mut egui::Ui) {
        let rect = ui.max_rect();
        let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
        let ctx = ui.ctx();

        let (scroll_y, hover_pos, dt, time) = ctx.input(|i| {
            (
                i.smooth_scroll_delta.y,
                i.pointer.hover_pos(),
                i.stable_dt.min(0.05_f32),
                i.time as f32,
            )
        });

        // ── Camera pan animation (spring toward target) ───────────────────────
        if let Some(target) = self.camera_anim_target {
            let delta = target - self.offset;
            if delta.length_sq() < 0.1 {
                self.offset = target;
                self.camera_anim_target = None;
            } else {
                self.offset += delta * CAM_LERP;
                ctx.request_repaint();
            }
        }

        if self.follow_selected_body {
            if let Some(idx) = self.selected_body {
                if let Some(body) = self.system.bodies().get(idx) {
                    let target =
                        egui::vec2(-body.x as f32 * self.scale, -body.y as f32 * self.scale);
                    let delta = target - self.offset;
                    if delta.length_sq() > 0.0001 {
                        self.offset += delta * FOLLOW_LERP;
                        ctx.request_repaint();
                    }
                } else {
                    self.follow_selected_body = false;
                    self.selected_body = None;
                    self.selection_form = None;
                }
            } else {
                self.follow_selected_body = false;
            }
        }

        // centre is derived AFTER any animation update
        let center = rect.center() + self.offset;

        // ── Zoom with inertia ─────────────────────────────────────────────────
        if scroll_y.abs() > 0.0 {
            self.zoom_vel += scroll_y * 0.0004;
            self.zoom_vel = self.zoom_vel.clamp(-0.15, 0.15);
        }

        if self.zoom_vel.abs() > 0.00005 {
            let old_scale = self.scale;
            self.scale = (self.scale * (1.0 + self.zoom_vel)).clamp(0.001, 50_000.0);

            if let Some(mouse) = hover_pos {
                let ratio = self.scale / old_scale;
                self.offset += (mouse - center) * (1.0 - ratio);
                // keep animation target in sync with the scale-adjusted offset
                if let Some(t) = self.camera_anim_target.as_mut() {
                    *t *= ratio;
                }
            }

            self.zoom_vel *= 0.80;
            if self.zoom_vel.abs() < 0.00005 {
                self.zoom_vel = 0.0;
            } else {
                ctx.request_repaint();
            }
        }

        // ── Pan with inertia ──────────────────────────────────────────────────
        if response.dragged() {
            let delta = response.drag_delta();
            self.offset += delta;
            // Cancel smooth-pan if the user takes manual control
            self.camera_anim_target = None;
            self.follow_selected_body = false;
            let frame_vel = delta / dt.max(1.0 / 120.0);
            self.pan_vel = self.pan_vel * 0.4 + frame_vel * 0.6;
            ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        } else {
            if self.pan_vel.length_sq() > 1.0 {
                self.offset += self.pan_vel * dt;
                self.pan_vel *= (-dt * 2.8_f32).exp();
                if self.pan_vel.length_sq() < 1.0 {
                    self.pan_vel = egui::Vec2::ZERO;
                }
                ctx.request_repaint();
            }
        }

        // ── Render params ─────────────────────────────────────────────────────
        let render_params = RenderParams {
            world_scale: self.scale,
            mode: self.semantic_scale_mode,
            min_px: match self.semantic_scale_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            },
        };

        // ── Hover detection (before click, drives cursor + ring) ──────────────
        let center_after_pan = rect.center() + self.offset;
        self.hovered_body = hover_pos
            .filter(|p| rect.contains(*p))
            .and_then(|p| self.find_body_at(p, center_after_pan, render_params));

        // Cursor: pointer over bodies, grab over empty space, grabbing while dragging
        if response.dragged() {
            ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        } else if self.hovered_body.is_some() {
            ctx.set_cursor_icon(egui::CursorIcon::PointingHand);
        } else if hover_pos.map_or(false, |p| rect.contains(p)) {
            ctx.set_cursor_icon(egui::CursorIcon::Grab);
        }

        // ── GPU body rendering ────────────────────────────────────────────────
        {
            let mut backend = self.backend.lock().unwrap();
            backend.begin();
            backend.show_grid = self.show_grid;

            let bodies = self.system.bodies();

            let mut screen_positions = Vec::with_capacity(bodies.len());

            for b in bodies {
                let px = center_after_pan.x + b.x as f32 * self.scale;
                let py = center_after_pan.y + b.y as f32 * self.scale;

                screen_positions.push(([px, py], b));
            }

            for (sp, b) in &screen_positions {
                if b.is_luminous() {
                    backend.add_light_source(LightSource {
                        screen_pos: *sp,
                        luminosity: b.luminosity as f32,
                    });
                }
            }

            for (sp, b) in &screen_positions {
                let [cr, cg, cb] = b.color;
                let r = compute_render_radius(b.physical_radius, render_params);

                backend.draw_circle(*sp, r, [cr, cg, cb]);
            }

            backend.set_lighting_params(0.55, 0.7);
            backend.center = [center_after_pan.x, center_after_pan.y];
            backend.scale = self.scale;
            backend.trail_width = self.trail_width;
            backend.trail_buffer = if self.show_trails {
                let bodies = self.system.bodies();
                let dom_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);
                backend.trail_visibility = Some(
                    bodies
                        .iter()
                        .map(|b| dom_mass == 0.0 || b.mass / dom_mass >= self.trail_min_mass_ratio)
                        .collect(),
                );
                Some(self.system.clone_trail_buf())
            } else {
                backend.trail_visibility = None;
                None
            };
        }

        // ── Place-mode: click or drag-release to spawn a body ────────────────
        if self.place_mode {
            let pointer = ctx.input(|i| i.pointer.clone());

            // Track drag start (only when pressing on empty space)
            if response.drag_started() {
                if let Some(pos) = pointer.press_origin() {
                    // Only start a place-drag if not clicking an existing body
                    if self.find_body_at(pos, center_after_pan, render_params).is_none() {
                        self.place_drag_start = Some(pos);
                    }
                }
            }

            // Draw velocity arrow while dragging
            if let Some(start) = self.place_drag_start {
                if let Some(cur) = pointer.hover_pos() {
                    let painter = ui.painter();
                    let delta = cur - start;
                    let len = delta.length();
                    if len > 4.0 {
                        // Line
                        painter.line_segment(
                            [start, cur],
                            Stroke::new(1.5, Color32::from_rgba_premultiplied(120, 200, 255, 200)),
                        );
                        // Arrowhead
                        let dir = delta / len;
                        let perp = egui::vec2(-dir.y, dir.x);
                        let tip = cur;
                        let base = tip - dir * 8.0;
                        painter.add(egui::Shape::convex_polygon(
                            vec![tip, base + perp * 4.0, base - perp * 4.0],
                            Color32::from_rgba_premultiplied(120, 200, 255, 200),
                            Stroke::NONE,
                        ));
                        // Velocity label
                        let v_scale = 0.5 / self.scale as f64;
                        let vx = delta.x as f64 * v_scale;
                        let vy = delta.y as f64 * v_scale;
                        let speed = (vx * vx + vy * vy).sqrt();
                        painter.text(
                            cur + egui::vec2(8.0, -12.0),
                            egui::Align2::LEFT_CENTER,
                            format!("v={:.3}", speed),
                            egui::FontId::monospace(9.5),
                            Color32::from_rgba_premultiplied(140, 210, 255, 200),
                        );
                    }
                }
            }

            // On release: spawn body
            if response.drag_stopped() || response.clicked() {
                let start = self.place_drag_start.take();
                if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                    // Don't spawn if clicking an existing body
                    if self.find_body_at(cursor, center_after_pan, render_params).is_none() {
                        let spawn_pos = start.unwrap_or(cursor);
                        let wx = (spawn_pos.x - center_after_pan.x) as f64 / self.scale as f64;
                        let wy = (spawn_pos.y - center_after_pan.y) as f64 / self.scale as f64;

                        // Velocity from drag delta
                        let v_scale = 0.5 / self.scale as f64;
                        let (vx, vy) = if let Some(s) = start {
                            let d = cursor - s;
                            (d.x as f64 * v_scale, d.y as f64 * v_scale)
                        } else {
                            (0.0, 0.0)
                        };

                        use crate::domain::materials::density as mat_density;
                        let mut body = crate::domain::body::Body::new(
                            wx,
                            wy,
                            vx,
                            vy,
                            self.place_mass,
                            self.place_material,
                        );
                        body.density = mat_density(self.place_material, self.place_mass);
                        body.sync_physical_properties();

                        self.push_undo(UndoRecord::AddedBodies(1));
                        self.system.add_body(body);
                    }
                }
            }

            // Override cursor to crosshair while in place-mode (unless dragging a body)
            if self.dragging_body.is_none() {
                ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
            }

            ctx.request_repaint();
        } else {
            // ── Normal click: select + pan to center ──────────────────────────
            self.place_drag_start = None;
            if response.clicked() {
                if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                    match self.find_body_at(cursor, center_after_pan, render_params) {
                        Some(idx) => {
                            let body = self.system.bodies()[idx];
                            self.selected_body = Some(idx);
                            self.follow_selected_body = true;
                            let name = self.system.name(idx).to_owned();
                            self.selection_form = Some(SelectionForm::from_body(&body, &name));

                            self.pan_vel = egui::Vec2::ZERO;
                            self.zoom_vel = 0.0;

                            let screen_r =
                                compute_render_radius(body.physical_radius, render_params);
                            if screen_r < 6.0 && body.physical_radius > 1e-30 {
                                let desired_px = 24.0_f32;
                                let new_scale = (desired_px / body.physical_radius as f32)
                                    .clamp(self.scale * 2.0, self.scale * 500.0)
                                    .min(50_000.0);
                                self.scale = new_scale;
                            }

                            self.camera_anim_target = Some(egui::vec2(
                                -body.x as f32 * self.scale,
                                -body.y as f32 * self.scale,
                            ));
                        },
                        None => {
                            self.selected_body = None;
                            self.follow_selected_body = false;
                            self.selection_form = None;
                        },
                    }
                }
            }
        }

        // ── Template drag-drop ────────────────────────────────────────────────
        // Check if a drag from the template panel was released over this canvas.
        if self.template_drag.is_some() {
            let released = ctx.input(|i| i.pointer.any_released());
            let drop_pos = ctx.input(|i| i.pointer.interact_pos());

            if released {
                if let (Some(build_fn), Some(screen_pos)) = (self.template_drag.take(), drop_pos) {
                    if rect.contains(screen_pos) {
                        // Convert screen pos → world pos
                        let wx = (screen_pos.x - center_after_pan.x) as f64 / self.scale as f64;
                        let wy = (screen_pos.y - center_after_pan.y) as f64 / self.scale as f64;
                        let template = (build_fn)();
                        self.active_units = template.units;
                        let bodies = instantiate_at(&template, wx, wy);
                        self.push_undo(UndoRecord::AddedBodies(bodies.len()));
                        self.system.add_named_bodies(bodies);
                        self.pending_fit = false; // dropped at explicit position — no auto-fit
                    } else {
                        // Released outside canvas — discard
                        self.template_drag = None;
                    }
                } else {
                    self.template_drag = None;
                }
            } else {
                // Still dragging — show a ghost cursor
                if hover_pos.map_or(false, |p| rect.contains(p)) {
                    ctx.set_cursor_icon(egui::CursorIcon::Crosshair);
                }
            }
        }

        // Keep animating while a body is selected (pulsing ring needs repaints)
        if self.selected_body.is_some() {
            ctx.request_repaint();
        }

        // ── GPU paint callback ────────────────────────────────────────────────
        let device = self.device.as_ref().unwrap().clone();
        let queue = self.queue.as_ref().unwrap().clone();
        let format = self.format.unwrap();

        ui.painter().add(egui_wgpu::Callback::new_paint_callback(
            rect,
            CallbackFn {
                backend: self.backend.clone(),
                device,
                queue,
                format,
                // Canvas dimensions (not full window) so `to_ndc` maps correctly
                // into the canvas-rect viewport that egui_wgpu sets for callbacks.
                screen: [rect.width(), rect.height()],
                viewport_min: [rect.min.x, rect.min.y],
            },
        ));

        // ── Overlay: rings + labels (on top of GPU layer) ─────────────────────
        self.draw_overlay(ui, center_after_pan, time);

        // ── Loading overlay ───────────────────────────────────────────────────
        if self.system.is_loading() {
            self.draw_loading_overlay(ui, rect, time);
            ctx.request_repaint();
        }

        // ── Playbar ───────────────────────────────────────────────────────────
        self.draw_playbar(ctx, rect, time);
    }

    // ── Overlay ───────────────────────────────────────────────────────────────

    fn draw_overlay(&self, ui: &egui::Ui, center: Pos2, time: f32) {
        let bodies = self.system.bodies();
        let names = self.system.names();

        if bodies.is_empty() {
            return;
        }

        let max_mass = bodies.iter().map(|b| b.mass).fold(0.0_f64, f64::max);

        // !! Use the SAME min_px as the GPU renderer so rings/labels are anchored
        //    to the visible disc, not to the (possibly sub-pixel) physical radius.
        let render_params = RenderParams {
            world_scale: self.scale,
            mode: self.semantic_scale_mode,
            min_px: match self.semantic_scale_mode {
                SemanticScaleMode::Physical => 0.0,
                SemanticScaleMode::Comparative => 3.0,
                SemanticScaleMode::Illustrative => 5.0,
            },
        };

        // Label visibility: show when body is large enough on screen OR important
        // enough given current zoom. threshold ↑ at low zoom → fewer labels.
        let importance_threshold = (2.0_f64 / self.scale as f64).clamp(0.001, 1.0);

        let painter = ui.painter();
        let font = FontId::proportional(LABEL_FONT_SIZE);

        // Pulse for selection ring: ±1.5 px at ~3.5 Hz
        let pulse = (time * 3.5).sin() * 1.5_f32;

        for (i, (body, name)) in bodies.iter().zip(names.iter()).enumerate() {
            // visual_r = how big the body actually appears on screen (matches GPU)
            let visual_r = compute_render_radius(body.physical_radius, render_params);

            let px = center.x + body.x as f32 * self.scale;
            let py = center.y + body.y as f32 * self.scale;
            let body_pos = egui::pos2(px, py);

            let is_selected = self.selected_body == Some(i);
            let is_hovered = self.hovered_body == Some(i) && !is_selected;

            // ── Hover ring ───────────────────────────────────────────────
            if is_hovered {
                painter.circle_stroke(
                    body_pos,
                    visual_r + RING_GAP - 1.0,
                    Stroke::new(1.0, Color32::from_rgba_premultiplied(160, 160, 200, 90)),
                );
            }

            // ── Selection ring (pulsing) ─────────────────────────────────
            if is_selected {
                let r1 = visual_r + RING_GAP + pulse;
                let r2 = r1 + 3.5;

                // Outer dim halo
                painter.circle_stroke(
                    body_pos,
                    r2,
                    Stroke::new(0.8, Color32::from_rgba_premultiplied(130, 130, 200, 55)),
                );
                // Main ring
                painter.circle_stroke(
                    body_pos,
                    r1,
                    Stroke::new(1.5, Color32::from_rgba_premultiplied(200, 200, 255, 210)),
                );
                // Inner tick marks at compass points for clarity
                for angle in [0.0_f32, 90.0, 180.0, 270.0] {
                    let rad = angle.to_radians();
                    let dir = egui::vec2(rad.cos(), rad.sin());
                    let inner = body_pos + dir * (r1 - 2.5);
                    let outer = body_pos + dir * (r1 + 2.5);
                    painter.line_segment(
                        [inner, outer],
                        Stroke::new(1.5, Color32::from_rgba_premultiplied(200, 200, 255, 180)),
                    );
                }
            }

            // ── Name label ───────────────────────────────────────────────
            let importance = if max_mass > 0.0 { body.mass / max_mass } else { 0.0 };
            let show_label = visual_r >= 5.0 || importance >= importance_threshold || is_selected;

            if show_label {
                // Cap offset so the label never drifts far from the body
                let offset_y = (visual_r + 4.0).min(MAX_LABEL_OFFSET_PX);
                let label_pos = egui::pos2(px, py + offset_y);

                let color = if is_selected {
                    Color32::from_rgb(220, 220, 255)
                } else {
                    Color32::from_rgba_premultiplied(175, 175, 195, 195)
                };

                // Shadow
                painter.text(
                    label_pos + egui::vec2(1.0, 1.0),
                    egui::Align2::CENTER_TOP,
                    name.as_str(),
                    font.clone(),
                    Color32::from_black_alpha(110),
                );
                painter.text(
                    label_pos,
                    egui::Align2::CENTER_TOP,
                    name.as_str(),
                    font.clone(),
                    color,
                );
            }
        }
    }

    // ── Loading overlay ───────────────────────────────────────────────────────

    fn draw_loading_overlay(&self, ui: &egui::Ui, rect: egui::Rect, time: f32) {
        let painter = ui.painter();

        // Dim backdrop — subtle, doesn't obliterate the scene.
        painter.rect_filled(rect, 0.0, Color32::from_black_alpha(120));

        let cx = rect.center();

        // ── Spinner: three arcs rotating at different phases ──────────────────
        // Outer ring
        draw_spinner_arc(
            painter,
            cx,
            22.0,
            2.5,
            time,
            1.0,
            Color32::from_rgba_premultiplied(160, 160, 255, 220),
        );
        // Middle ring (counter-rotate, dimmer)
        draw_spinner_arc(
            painter,
            cx,
            14.0,
            2.0,
            -time * 1.4,
            0.75,
            Color32::from_rgba_premultiplied(120, 120, 200, 150),
        );
        // Inner dot
        painter.circle_filled(cx, 3.5, Color32::from_rgba_premultiplied(180, 180, 255, 200));

        // ── "LOADING" label ────────────────────────────────────────────────────
        let label_pos = egui::pos2(cx.x, cx.y + 36.0);
        // Shadow
        painter.text(
            label_pos + egui::vec2(1.0, 1.0),
            egui::Align2::CENTER_TOP,
            "LOADING",
            eframe::egui::FontId::proportional(11.0),
            Color32::from_black_alpha(120),
        );
        painter.text(
            label_pos,
            egui::Align2::CENTER_TOP,
            "LOADING",
            eframe::egui::FontId::proportional(11.0),
            Color32::from_rgba_premultiplied(170, 170, 220, 210),
        );
    }

    // ── Hit-test ──────────────────────────────────────────────────────────────

    // ── Playbar ───────────────────────────────────────────────────────────────

    fn draw_playbar(&mut self, ctx: &egui::Context, canvas_rect: egui::Rect, time: f32) {
        use crate::app::theme::{ACCENT, ACCENT_DIM, SUCCESS, TEXT_DIM, TEXT_SEC};

        let bar_w = 400.0_f32;
        let bar_h = 44.0_f32;
        let anchor =
            egui::pos2(canvas_rect.center().x - bar_w * 0.5, canvas_rect.max.y - bar_h - 16.0);

        egui::Area::new(egui::Id::new("playbar"))
            .fixed_pos(anchor)
            .order(egui::Order::Foreground)
            .show(ctx, |ui| {
                egui::Frame::NONE
                    .fill(egui::Color32::from_rgba_unmultiplied(14, 14, 20, 220))
                    .stroke(egui::Stroke::new(
                        1.0,
                        egui::Color32::from_rgba_unmultiplied(50, 50, 70, 180),
                    ))
                    .corner_radius(10.0)
                    .inner_margin(egui::Margin::symmetric(12, 0))
                    .show(ui, |ui| {
                        ui.set_width(bar_w - 24.0);
                        ui.set_height(bar_h);

                        ui.horizontal_centered(|ui| {
                            ui.spacing_mut().item_spacing.x = 5.0;

                            // ── Simulation time ───────────────────────────────
                            let m = self.system.metrics();
                            let t_str = fmt_sim_time(m.t);
                            ui.label(
                                egui::RichText::new(t_str).monospace().size(10.5).color(TEXT_SEC),
                            );

                            ui.add(egui::Separator::default().vertical().spacing(4.0));

                            // ── Play / Pause ──────────────────────────────────
                            let (icon, icon_col) =
                                if self.paused { ("▶", SUCCESS) } else { ("⏸", ACCENT) };

                            // Pulsing glow ring when running
                            let btn_pos = ui.next_widget_position() + egui::vec2(18.0, 18.0);
                            if !self.paused {
                                let pulse = ((time * 2.0).sin() * 0.5 + 0.5) * 0.35 + 0.1;
                                ui.painter().circle_stroke(
                                    btn_pos,
                                    22.0,
                                    egui::Stroke::new(
                                        1.5,
                                        egui::Color32::from_rgba_unmultiplied(
                                            ACCENT.r(),
                                            ACCENT.g(),
                                            ACCENT.b(),
                                            (pulse * 180.0) as u8,
                                        ),
                                    ),
                                );
                            }

                            let play_btn = ui
                                .add(
                                    egui::Button::new(
                                        egui::RichText::new(icon).size(18.0).color(icon_col),
                                    )
                                    .fill(if self.paused {
                                        ACCENT_DIM
                                    } else {
                                        egui::Color32::from_rgba_unmultiplied(30, 50, 35, 180)
                                    })
                                    .stroke(egui::Stroke::new(1.0, icon_col.gamma_multiply(0.5)))
                                    .min_size(egui::vec2(36.0, 36.0)),
                                )
                                .on_hover_text(if self.paused {
                                    "Play (Space)"
                                } else {
                                    "Pause (Space)"
                                });
                            if play_btn.clicked() {
                                self.paused = !self.paused;
                            }

                            ui.add(egui::Separator::default().vertical().spacing(4.0));

                            // ── Speed slider (logarithmic steps/frame) ────────
                            // Label shows current multiplier; slider gives fine control.
                            let spf_col = if self.steps_per_frame > 1 { ACCENT } else { TEXT_DIM };
                            ui.label(
                                egui::RichText::new(format!("×{}", self.steps_per_frame))
                                    .monospace()
                                    .size(10.0)
                                    .color(spf_col),
                            );

                            let mut spf_f = self.steps_per_frame as f32;
                            let slider_r = ui.add_sized(
                                [110.0, 20.0],
                                egui::Slider::new(&mut spf_f, 1.0..=10000.0)
                                    .logarithmic(true)
                                    .show_value(false),
                            );
                            if slider_r.changed() {
                                self.steps_per_frame = spf_f.round().max(1.0) as u32;
                            }

                            ui.add(egui::Separator::default().vertical().spacing(4.0));

                            // ── dt ────────────────────────────────────────────
                            ui.label(egui::RichText::new("dt").size(9.5).color(TEXT_DIM));
                            let mut dt = self.system.dt();
                            let dt_speed = (dt * 0.05).max(1e-7);
                            let dt_r = ui
                                .add(
                                    egui::DragValue::new(&mut dt)
                                        .speed(dt_speed)
                                        .range(1e-7_f64..=10.0)
                                        .max_decimals(6)
                                        .min_decimals(1),
                                )
                                .on_hover_text("Integration timestep — smaller = more accurate");
                            if dt_r.changed() {
                                self.system.set_dt(dt);
                            }
                        });
                    });
            });
    }

    fn find_body_at(
        &self,
        cursor: Pos2,
        center: Pos2,
        render_params: RenderParams,
    ) -> Option<usize> {
        let bodies = self.system.bodies();

        // Iterate in reverse so top-rendered (last) body wins ties
        for i in (0..bodies.len()).rev() {
            let b = &bodies[i];

            let px = center.x + b.x as f32 * self.scale;
            let py = center.y + b.y as f32 * self.scale;

            // Hit radius: visual size + generous minimum for easy clicking
            let r = compute_render_radius(b.physical_radius, render_params).max(MIN_HIT_PX);

            let dx = cursor.x - px;
            let dy = cursor.y - py;

            if dx * dx + dy * dy <= r * r {
                return Some(i);
            }
        }

        None
    }
}

// ── Sim-time formatter ────────────────────────────────────────────────────────

/// Compact display of simulated time: uses natural units when small, sci notation when large.
fn fmt_sim_time(t: f64) -> String {
    if t == 0.0 {
        return "t=0".into();
    }
    let a = t.abs();
    if a < 1e-3 {
        format!("t={:+.2e}", t)
    } else if a < 1_000.0 {
        format!("t={:.4}", t)
    } else {
        format!("t={:.3e}", t)
    }
}

// ── Spinner helpers ───────────────────────────────────────────────────────────

/// Draw a rotating arc at `center` with the given `radius` and line `width`.
/// `time` drives rotation speed; `arc_frac` (0..1) is how much of the circle
/// is filled (e.g. 0.75 = 270°). Alpha fades from full at the head to 0 at
/// the tail for a "comet tail" effect.
fn draw_spinner_arc(
    painter: &egui::Painter,
    center: Pos2,
    radius: f32,
    width: f32,
    time: f32,
    arc_frac: f32,
    color: Color32,
) {
    let segments = 32usize;
    let arc_radians = arc_frac * std::f32::consts::TAU;
    let angle_step = arc_radians / segments as f32;
    let base_angle = time * 2.2; // rotation speed

    let [r, g, b, _] = color.to_array();

    for i in 0..segments {
        let t_frac = i as f32 / segments as f32; // 0 = tail, 1 = head
        let alpha = (t_frac * t_frac * 255.0) as u8; // quadratic fade

        let a0 = base_angle + i as f32 * angle_step;
        let a1 = a0 + angle_step * 1.1; // slight overlap to avoid gaps

        let p0 = center + egui::vec2(a0.cos(), a0.sin()) * radius;
        let p1 = center + egui::vec2(a1.cos(), a1.sin()) * radius;

        painter.line_segment(
            [p0, p1],
            Stroke::new(width, Color32::from_rgba_premultiplied(r, g, b, alpha)),
        );
    }
}
