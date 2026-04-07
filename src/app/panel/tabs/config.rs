use crate::app::config::PhysicsConfig;
use crate::app::theme::{TEXT_DIM, TEXT_SEC};
use crate::app::theme::secondary_btn;
use crate::app::ui::SimulationApp;
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn panel_tab_config(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);

        ui.label(RichText::new("G multiplier").size(10.0).color(TEXT_SEC));
        let g_resp = ui.add(
            egui::Slider::new(&mut self.physics_cfg.g_factor, 0.01..=10.0)
                .logarithmic(true)
                .show_value(true)
                .text(""),
        );
        if g_resp.changed() {
            self.system.set_g_factor(self.physics_cfg.g_factor);
            self.system.reset_energy_baseline();
        }
        ui.label(
            RichText::new(format!("G_eff = {:.4}", self.physics_cfg.g_factor))
                .size(9.0)
                .color(TEXT_DIM),
        );

        ui.add_space(8.0);
        ui.label(
            RichText::new("restitution  (0 = merge)")
                .size(10.0)
                .color(TEXT_SEC),
        );
        if ui
            .add(
                egui::Slider::new(&mut self.collision_cor, 0.0..=1.0)
                    .show_value(true)
                    .text(""),
            )
            .changed()
        {
            self.system.set_cor(self.collision_cor);
        }

        ui.add_space(8.0);
        ui.label(
            RichText::new("UNIT LABELS")
                .size(9.5)
                .color(TEXT_DIM)
                .strong(),
        );
        ui.add_space(3.0);
        egui::Grid::new("unit_labels")
            .num_columns(2)
            .spacing([4.0, 2.0])
            .show(ui, |ui| {
                ui.label(RichText::new("mass").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.mass_label)
                        .desired_width(48.0),
                );
                ui.end_row();
                ui.label(RichText::new("dist").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.dist_label)
                        .desired_width(48.0),
                );
                ui.end_row();
                ui.label(RichText::new("time").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.time_label)
                        .desired_width(48.0),
                );
                ui.end_row();
            });

        ui.add_space(8.0);
        if secondary_btn(ui, "Reset to defaults") {
            self.physics_cfg = PhysicsConfig::default();
            self.system.set_g_factor(1.0);
            self.system.set_cor(0.0);
            self.collision_cor = 0.0;
            self.system.reset_energy_baseline();
        }
    }
}
