use crate::app::theme::{ACCENT_DIM, BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{SimulationApp, UndoRecord};
use crate::templates::{TEMPLATES, TemplateCategory, instantiate_at};
use eframe::egui::{self, RichText, Stroke};

impl SimulationApp {
    pub(super) fn panel_tab_templates(&mut self, ui: &mut egui::Ui) {
        ui.add_space(2.0);
        ui.label(
            RichText::new("drag → canvas  or  click → viewport center")
                .size(9.0)
                .color(TEXT_DIM)
                .italics(),
        );
        ui.add_space(6.0);

        for cat in [
            TemplateCategory::Bodies,
            TemplateCategory::Systems,
            TemplateCategory::ThreeBodyProblems,
        ] {
            let entries: Vec<_> = TEMPLATES.iter().filter(|e| e.category == cat).collect();
            if entries.is_empty() {
                continue;
            }

            ui.label(RichText::new(cat.label()).size(9.5).color(TEXT_DIM).strong());
            ui.add_space(1.0);
            ui.add(egui::Separator::default().spacing(4.0));
            ui.add_space(3.0);

            egui::Grid::new(cat.grid_id()).num_columns(2).spacing([4.0, 4.0]).show(ui, |ui| {
                for (i, entry) in entries.iter().enumerate() {
                    let response = ui.add(
                        egui::Button::new(RichText::new(entry.name).size(10.0).color(TEXT_PRI))
                            .fill(ACCENT_DIM)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(88.0, 22.0))
                            .sense(egui::Sense::click_and_drag()),
                    );

                    // Drag → initiate cross-widget drag to canvas
                    if response.drag_started() {
                        let seed = self.system.seed();
                        let build = entry.build;
                        self.template_drag = Some(Box::new(move || (build)(seed)));
                    }

                    // Click (no drag) → spawn at viewport center (world origin of
                    // current offset, i.e. the screen centre in world space)
                    if response.clicked() {
                        let template = (entry.build)(self.system.seed());
                        let bodies = instantiate_at(&template, 0.0, 0.0);
                        self.push_undo(UndoRecord::AddedBodies(bodies.len()));
                        self.system.add_named_bodies(bodies);
                        self.pending_fit = true;
                        self.reset_drift_peaks();
                        // Seed the simulation name from the template (user can rename later)
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

            ui.add_space(8.0);
        }
    }
}
