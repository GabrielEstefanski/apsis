use crate::app::theme::{ACCENT, TEXT_DIM};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn panel_time_speed(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        ui.horizontal(|ui| {
            ui.label(
                RichText::new("TIME SPEED")
                    .size(9.5)
                    .color(TEXT_DIM)
                    .strong(),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let col = if self.steps_per_frame > 1 {
                    ACCENT
                } else {
                    TEXT_DIM
                };
                ui.label(
                    RichText::new(format!("×{}", self.steps_per_frame))
                        .monospace()
                        .size(10.0)
                        .color(col),
                );
            });
        });

        ui.add_space(2.0);
        ui.add(
            egui::Slider::new(&mut self.steps_per_frame, 1..=1000u32)
                .logarithmic(true)
                .show_value(false),
        );
    }
}
