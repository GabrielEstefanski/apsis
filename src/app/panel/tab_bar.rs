use crate::app::theme::{ACCENT_DIM, BORDER, TEXT_PRI, TEXT_SEC};
use crate::app::ui::{PanelTab, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn panel_tab_bar(&mut self, ui: &mut egui::Ui) {
        const TABS: &[(PanelTab, &str)] = &[
            (PanelTab::Add, "Add"),
            (PanelTab::Templates, "Library"),
            (PanelTab::Config, "Config"),
        ];
        let w = (196.0 - 4.0 * (TABS.len() as f32 - 1.0)) / TABS.len() as f32;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 4.0;
            for (tab, label) in TABS {
                let active = self.panel_tab == *tab;
                let col = if active { TEXT_PRI } else { TEXT_SEC };
                let fill = if active { ACCENT_DIM } else { Color32::TRANSPARENT };
                if ui
                    .add(
                        egui::Button::new(RichText::new(*label).size(10.5).color(col))
                            .fill(fill)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(w, 22.0)),
                    )
                    .clicked()
                {
                    self.panel_tab = *tab;
                }
            }
        });
    }
}
