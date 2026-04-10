use crate::app::render_params::{RenderParams, compute_render_radius};
use crate::app::theme::{BG, nice_grid_world};
use crate::app::ui::{SelectionForm, SemanticScaleMode, SimulationApp};
use crate::domain::body::Body;
use crate::render::CallbackFn;
use eframe::egui::{self, Pos2};
use eframe::egui_wgpu;

impl SimulationApp {
    pub(super) fn draw_canvas(&mut self, ui: &mut egui::Ui) {
        // 🔴 REMOVIDO CentralPanel daqui

        let rect = ui.max_rect();

        let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());

        let center = rect.center() + self.offset;

        let ctx = ui.ctx();

        let (scroll_y, hover_pos, dt) = ctx.input(|i| (
            i.smooth_scroll_delta.y,
            i.pointer.hover_pos(),
            i.stable_dt.min(0.05_f32),
        ));

        // ── Zoom com inércia ─────────────────────────────────────────────────
        if scroll_y.abs() > 0.0 {
            // Acumula velocidade — clamp para evitar zoom explosivo em trackpads
            self.zoom_vel += scroll_y * 0.0004;
            self.zoom_vel = self.zoom_vel.clamp(-0.15, 0.15);
        }

        if self.zoom_vel.abs() > 0.00005 {
            let old_scale = self.scale;
            self.scale = (self.scale * (1.0 + self.zoom_vel)).clamp(0.001, 50_000.0);

            // Mantém o ponto do mundo sob o cursor fixo na tela
            if let Some(mouse) = hover_pos {
                let ratio = self.scale / old_scale;
                self.offset += (mouse - center) * (1.0 - ratio);
            }

            // Fricção: decai ~20% por frame → suave em ~15 frames (0.25s)
            self.zoom_vel *= 0.80;
            if self.zoom_vel.abs() < 0.00005 {
                self.zoom_vel = 0.0;
            } else {
                ctx.request_repaint();
            }
        }

        // ── Pan com inércia ──────────────────────────────────────────────────
        if response.dragged() {
            let delta = response.drag_delta();
            self.offset += delta;
            // Atualiza velocidade como média exponencial do delta instantâneo
            let frame_vel = delta / dt.max(1.0 / 120.0);
            self.pan_vel = self.pan_vel * 0.4 + frame_vel * 0.6;
            ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        } else {
            if self.pan_vel.length_sq() > 1.0 {
                self.offset += self.pan_vel * dt;
                // Fricção contínua: T½ ≈ 0.25s
                self.pan_vel *= (-dt * 2.8_f32).exp();
                if self.pan_vel.length_sq() < 1.0 {
                    self.pan_vel = egui::Vec2::ZERO;
                }
                ctx.request_repaint();
            }
            if hover_pos.map_or(false, |p| rect.contains(p)) {
                ctx.set_cursor_icon(egui::CursorIcon::Grab);
            }
        }

        let mut backend = self.backend.lock().unwrap();
        backend.begin();
        backend.show_grid = self.show_grid;

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

        for b in bodies {
            let [cr, cg, cb] = b.color;

            let r = compute_render_radius(b.physical_radius, render_params);

            let px = center.x + b.x as f32 * self.scale;
            let py = center.y + b.y as f32 * self.scale;

            backend.draw_circle([px, py], r, [cr, cg, cb]);
        }

        // Passa estado de câmera e trail buffer pro callback de GPU
        backend.center = [center.x, center.y];
        backend.scale = self.scale;
        backend.trail_width = self.trail_width;
        if self.show_trails {
            backend.trail_buffer = Some(self.system.trail_buf().clone());
        } else {
            backend.trail_buffer = None;
        }

        drop(backend);

        // ── Click to select / deselect body ──────────────────────────────────
        if response.clicked() {
            if let Some(cursor) = ctx.input(|i| i.pointer.interact_pos()) {
                let hit = self.find_body_at(cursor, center);
                match hit {
                    Some(idx) => {
                        let body = self.system.bodies()[idx];
                        self.selected_body = Some(idx);
                        self.selection_form = Some(SelectionForm::from_body(&body));
                    }
                    None => {
                        self.selected_body = None;
                        self.selection_form = None;
                    }
                }
            }
        }

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
                screen: {
                    let r = ctx.input(|i| i.content_rect());
                    [r.width(), r.height()]
                },
            },
        ));
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

            let r = compute_render_radius(b.physical_radius, render_params).max(6.0);

            let dx = cursor.x - px;
            let dy = cursor.y - py;

            if dx * dx + dy * dy <= r * r {
                return Some(i);
            }
        }

        None
    }
}
