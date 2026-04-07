mod add;
mod config;
mod templates;

use crate::app::ui::{PanelTab, SimulationApp};
use eframe::egui;

impl SimulationApp {
    /// Dispatches to the active tab's content. Called from the scroll area in
    /// `panel/mod.rs`. Adding a new tab = add a variant to `PanelTab` + a new
    /// file here.
    pub(super) fn panel_tab_dispatch(&mut self, ui: &mut egui::Ui) {
        match self.panel_tab {
            PanelTab::Add => self.panel_tab_add(ui),
            PanelTab::Templates => self.panel_tab_templates(ui),
            PanelTab::Config => self.panel_tab_config(ui),
        }
    }
}
