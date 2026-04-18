use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{SemanticScaleMode, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

// Frame-counter for the rec dot pulse (toolbar has no `time` param)
// We borrow ctx.input time instead.

impl SimulationApp {
    pub(super) fn toolbar_content(&mut self, ui: &mut egui::Ui) {
        let time = ui.ctx().input(|i| i.time as f32);

        ui.horizontal(|ui| {
            ui.label(RichText::new("GRAVITY SIM").size(11.0).color(TEXT_PRI).strong());

            // ── Simulation name (editable inline) ───────────────────────────── //
            ui.separator();
            let name_display =
                if self.sim_name.is_empty() { "Unnamed".to_owned() } else { self.sim_name.clone() };
            let name_resp = ui.add(
                egui::TextEdit::singleline(&mut self.sim_name)
                    .desired_width(120.0)
                    .hint_text(name_display)
                    .font(egui::FontId::proportional(10.0)),
            );
            if name_resp.hovered() {
                name_resp.on_hover_text("Simulation name — used in save files and recordings");
            }

            ui.separator();

            // ── Camera / view controls ──────────────────────────────── //
            ui.label(RichText::new("zoom").size(10.0).color(TEXT_SEC));
            let zoom_speed = self.scale * 0.01;
            ui.add(
                egui::DragValue::new(&mut self.scale)
                    .speed(zoom_speed)
                    .range(0.001..=50_000.0_f32)
                    .max_decimals(2),
            );

            ui.label(RichText::new("body sz").size(10.0).color(TEXT_SEC));
            ui.add(
                egui::DragValue::new(&mut self.body_size_boost)
                    .speed(0.5)
                    .range(1.0..=500.0_f32)
                    .max_decimals(1),
            )
            .on_hover_text("Visual size multiplier for bodies (does not affect physics)");

            ui.label(RichText::new("scale mode").size(10.0).color(TEXT_SEC));
            egui::ComboBox::from_id_salt("semantic_scale_mode")
                .selected_text(self.semantic_scale_mode.label())
                .width(96.0)
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.semantic_scale_mode,
                        SemanticScaleMode::Physical,
                        SemanticScaleMode::Physical.label(),
                    );
                    ui.selectable_value(
                        &mut self.semantic_scale_mode,
                        SemanticScaleMode::Comparative,
                        SemanticScaleMode::Comparative.label(),
                    );
                    ui.selectable_value(
                        &mut self.semantic_scale_mode,
                        SemanticScaleMode::Illustrative,
                        SemanticScaleMode::Illustrative.label(),
                    );
                });

            ui.separator();

            // ── Display toggles ─────────────────────────────────────── //
            toggle(ui, &mut self.show_grid, "grid");
            toggle(ui, &mut self.show_trails, "trails");
            if self.show_trails {
                ui.add(
                    egui::DragValue::new(&mut self.trail_width)
                        .speed(0.1)
                        .range(0.5_f32..=20.0)
                        .max_decimals(1)
                        .prefix("w:"),
                )
                .on_hover_text("Trail ribbon width in pixels");
            }
            toggle(ui, &mut self.show_orbit_ellipses, "orbits");
            toggle(ui, &mut self.show_vectors, "vel");
            toggle(ui, &mut self.show_force_vectors, "force");
            toggle(ui, &mut self.show_belts, "structure");

            ui.separator();

            // ── Right-side actions ──────────────────────────────────── //
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                // Body count
                let n = self.system.bodies().len();
                ui.label(RichText::new(format!("{n} bodies")).size(10.0).color(TEXT_DIM));
                ui.separator();

                // ⚙ Settings
                if ui
                    .add(
                        egui::Button::new(RichText::new("⚙").size(13.0).color(TEXT_DIM))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(24.0, 20.0)),
                    )
                    .on_hover_text("Settings — unit labels, recording")
                    .clicked()
                {
                    self.show_settings_modal = !self.show_settings_modal;
                }

                // ? Shortcuts
                if ui
                    .add(
                        egui::Button::new(RichText::new("?").size(11.0).color(TEXT_DIM))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(22.0, 20.0)),
                    )
                    .on_hover_text("Keyboard shortcuts (H)")
                    .clicked()
                {
                    self.show_shortcuts_modal = !self.show_shortcuts_modal;
                }

                // ● Record dot
                let is_rec = self.recorder.is_some();
                let pulse_alpha =
                    if is_rec { (((time * 2.5).sin() * 0.4 + 0.6) * 255.0) as u8 } else { 100 };
                let rec_col = Color32::from_rgba_unmultiplied(
                    DANGER.r(),
                    DANGER.g(),
                    DANGER.b(),
                    pulse_alpha,
                );
                let rec_btn = ui
                    .add(
                        egui::Button::new(RichText::new("●").size(13.0).color(rec_col))
                            .fill(if is_rec {
                                Color32::from_rgba_unmultiplied(50, 10, 10, 80)
                            } else {
                                Color32::TRANSPARENT
                            })
                            .stroke(Stroke::new(
                                0.5,
                                if is_rec { DANGER.gamma_multiply(0.5) } else { BORDER },
                            ))
                            .min_size(egui::vec2(24.0, 20.0)),
                    )
                    .on_hover_text(if is_rec {
                        "Recording — click to stop"
                    } else {
                        "Start recording (open Settings)"
                    });
                if rec_btn.clicked() {
                    if is_rec {
                        if let Some(mut rec) = self.recorder.take() {
                            let _ = rec.flush();
                        }
                    } else {
                        self.show_settings_modal = true;
                    }
                }

                ui.separator();

                // Saves
                if ui
                    .add(
                        egui::Button::new(RichText::new("Saves").size(10.0).color(ACCENT))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, ACCENT))
                            .min_size(egui::vec2(48.0, 20.0)),
                    )
                    .on_hover_text("Browse and load saved states")
                    .clicked()
                {
                    self.open_save_modal();
                }

                if ui
                    .add(
                        egui::Button::new(RichText::new("Save").size(10.0).color(SUCCESS))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, SUCCESS))
                            .min_size(egui::vec2(40.0, 20.0)),
                    )
                    .on_hover_text("Quick-save current state")
                    .clicked()
                {
                    let _ = self.do_save();
                }

                ui.separator();

                if ui
                    .add(
                        egui::Button::new(RichText::new("Clear").size(10.0).color(DANGER))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, DANGER))
                            .min_size(egui::vec2(46.0, 20.0)),
                    )
                    .clicked()
                {
                    self.system.load_bodies(vec![]);
                    self.paused = true;
                    self.reset_drift_peaks();
                    self.sim_name = String::new();
                }

                if ui
                    .add(
                        egui::Button::new(RichText::new("Zero COM").size(10.0).color(TEXT_SEC))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(60.0, 20.0)),
                    )
                    .on_hover_text("Zero centre-of-mass velocity")
                    .clicked()
                {
                    self.system.zero_com_velocity();
                }

                if ui
                    .add(
                        egui::Button::new(RichText::new("Fit view").size(10.0).color(TEXT_SEC))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(52.0, 20.0)),
                    )
                    .clicked()
                {
                    self.fit_to_view();
                }
            });
        });
    }
}

fn toggle(ui: &mut egui::Ui, value: &mut bool, label: &str) {
    let col = if *value { TEXT_PRI } else { TEXT_SEC };
    if ui
        .add(
            egui::Button::new(RichText::new(label).size(10.5).color(col))
                .fill(Color32::TRANSPARENT)
                .stroke(Stroke::new(0.5, if *value { BORDER } else { Color32::TRANSPARENT }))
                .min_size(egui::vec2(0.0, 20.0)),
        )
        .clicked()
    {
        *value = !*value;
    }
}
