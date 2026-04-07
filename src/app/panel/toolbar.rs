use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{SemanticScaleMode, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn toolbar_content(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(
                RichText::new("GRAVITY SIM")
                    .size(12.0)
                    .color(TEXT_PRI)
                    .strong(),
            );

            ui.separator();

            // Play / Pause
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

            // Simulation parameters
            let dt_speed = self.proposed_dt * 0.05;

            ui.label(RichText::new("dt").size(10.0).color(TEXT_SEC));
            ui.add(
                egui::DragValue::new(&mut self.proposed_dt)
                    .speed(dt_speed)
                    .clamp_range(1e-5..=0.5)
                    .min_decimals(3)
                    .max_decimals(5),
            );

            ui.label(RichText::new("zoom").size(10.0).color(TEXT_SEC));
            ui.add(
                egui::DragValue::new(&mut self.scale)
                    .speed(0.5)
                    .clamp_range(1.0..=5000.0f32)
                    .max_decimals(1),
            );

            ui.label(RichText::new("size").size(10.0).color(TEXT_SEC));
            ui.add(
                egui::DragValue::new(&mut self.body_size_boost)
                    .speed(0.5)
                    .clamp_range(1.0..=500.0f32)
                    .max_decimals(1),
            );

            ui.label(RichText::new("camera").size(10.0).color(TEXT_SEC));
            egui::ComboBox::from_id_source("semantic_scale_mode")
                .selected_text(self.semantic_scale_mode.label())
                .width(86.0)
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

            // Display toggles
            ui.checkbox(
                &mut self.show_grid,
                RichText::new("grid").size(11.0).color(TEXT_SEC),
            );
            ui.checkbox(
                &mut self.show_trails,
                RichText::new("trails").size(11.0).color(TEXT_SEC),
            );
            ui.checkbox(
                &mut self.show_vectors,
                RichText::new("vel").size(11.0).color(TEXT_SEC),
            );
            ui.checkbox(
                &mut self.show_force_vectors,
                RichText::new("force").size(11.0).color(TEXT_SEC),
            );
            ui.checkbox(
                &mut self.show_impact_normals,
                RichText::new("nrm").size(11.0).color(TEXT_SEC),
            );

            ui.separator();

            // Right-aligned system buttons
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let n = self.system.bodies().len();
                ui.label(
                    RichText::new(format!("{} bodies", n))
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
                            RichText::new("Reset E₀").size(10.0).color(TEXT_SEC),
                        )
                        .fill(Color32::TRANSPARENT)
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(60.0, 20.0)),
                    )
                    .clicked()
                {
                    self.system.reset_energy_baseline();
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
                    .clicked()
                {
                    self.system.zero_com_velocity();
                    self.system.reset_energy_baseline();
                }
            });
        });
    }
}
