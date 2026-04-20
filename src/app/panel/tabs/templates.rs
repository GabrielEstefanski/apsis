use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{SimulationApp, UndoRecord};
use crate::templates::{TEMPLATES, TemplateCategory, instantiate_at};
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn panel_tab_templates(&mut self, ui: &mut egui::Ui) {
        ui.add_space(2.0);
        ui.label(
            RichText::new("click to load at center  ·  drag to canvas")
                .size(9.5)
                .color(TEXT_DIM),
        );
        ui.add_space(8.0);

        for cat in [
            TemplateCategory::Bodies,
            TemplateCategory::Systems,
            TemplateCategory::ThreeBodyProblems,
        ] {
            let entries: Vec<_> = TEMPLATES.iter().filter(|e| e.category == cat).collect();
            if entries.is_empty() {
                continue;
            }

            // Category header with hairline rule
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new(cat.label())
                        .size(9.0)
                        .color(TEXT_SEC)
                        .strong(),
                );
                ui.add_space(4.0);
                let r = ui.available_rect_before_wrap();
                ui.painter().line_segment(
                    [
                        egui::pos2(r.left(), r.center().y),
                        egui::pos2(r.right(), r.center().y),
                    ],
                    Stroke::new(0.5, BORDER),
                );
            });
            ui.add_space(4.0);

            // 2-column button grid
            egui::Grid::new(cat.grid_id())
                .num_columns(2)
                .spacing([4.0, 4.0])
                .show(ui, |ui| {
                    for (i, entry) in entries.iter().enumerate() {
                        let response = ui.add(
                            egui::Button::new(
                                RichText::new(entry.name).size(10.0).color(TEXT_PRI),
                            )
                            .fill(Color32::from_rgb(20, 20, 26))
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(88.0, 22.0))
                            .sense(egui::Sense::click_and_drag()),
                        );

                        if response.drag_started() {
                            let seed = self.system.seed();
                            let build = entry.build;
                            self.template_drag = Some(Box::new(move || (build)(seed)));
                        }

                        if response.clicked() {
                            let template = (entry.build)(self.system.seed());
                            let bodies = instantiate_at(&template, 0.0, 0.0);
                            self.push_undo(UndoRecord::AddedBodies(bodies.len()));
                            self.system.add_named_bodies(bodies);
                            self.pending_fit = true;
                            self.reset_drift_peaks();
                            if self.sim_name.is_empty() {
                                self.sim_name = entry.name.to_owned();
                            }
                        }

                        if response.hovered() && self.template_drag.is_none() {
                            response.on_hover_text_at_pointer(
                                RichText::new("click: add at center\ndrag: place on canvas")
                                    .size(10.0)
                                    .color(TEXT_SEC),
                            );
                        }

                        if i % 2 == 1 {
                            ui.end_row();
                        }
                    }
                    if entries.len() % 2 == 1 {
                        ui.end_row();
                    }
                });

            ui.add_space(10.0);
        }

        let _ = (ACCENT, ACCENT_DIM); // palette reserved for future active state
    }
}
