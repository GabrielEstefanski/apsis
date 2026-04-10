use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{SemanticScaleMode, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn toolbar_content(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(
                RichText::new("GRAVITY SIM")
                    .size(11.0)
                    .color(TEXT_PRI)
                    .strong(),
            );

            ui.separator();

            // ── Play / Pause ────────────────────────────────────────── //
            let (lbl, col) = if self.paused {
                ("▶  Run", SUCCESS)
            } else {
                ("⏸  Pause", ACCENT)
            };
            if ui
                .add(
                    egui::Button::new(RichText::new(lbl).size(11.0).color(col))
                        .fill(Color32::TRANSPARENT)
                        .stroke(Stroke::new(1.0, col))
                        .min_size(egui::vec2(72.0, 22.0)),
                )
                .clicked()
            {
                self.paused = !self.paused;
            }

            ui.separator();

            // ── Integration timestep ────────────────────────────────── //
            ui.label(RichText::new("dt").size(10.0).color(TEXT_SEC));
            let mut dt = self.system.dt();
            let speed = (dt * 0.05).max(1e-7);
            let r = ui.add(
                egui::DragValue::new(&mut dt)
                    .speed(speed)
                    .range(1e-6..=1.0)
                    .max_decimals(6),
            );
            if r.changed() {
                self.system.set_dt(dt);
            }
            if r.hovered() {
                egui::show_tooltip_text(ui.ctx(), ui.layer_id(), egui::Id::new("dt_tip"),
                    "Integration timestep. Smaller = more accurate but slower.");
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
                let n = self.system.bodies().len();
                ui.label(
                    RichText::new(format!("{n} bodies"))
                        .size(10.0)
                        .color(TEXT_DIM),
                );
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
                }

                if ui
                    .add(
                        egui::Button::new(
                            RichText::new("Zero COM").size(10.0).color(TEXT_SEC),
                        )
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
                        egui::Button::new(
                            RichText::new("Fit view").size(10.0).color(TEXT_SEC),
                        )
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
