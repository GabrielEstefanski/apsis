mod add;
mod camera;
mod config;
mod overview;
mod templates;
mod view;

use crate::app::ui::{PanelTab, SimulationApp};
use eframe::egui;

impl SimulationApp {
    /// Dispatches to the active tool's contextual panel content.
    /// Adding a new tool = add a variant to `PanelTab` + a new file here +
    /// a match arm below.
    pub(super) fn panel_tab_dispatch(&mut self, ui: &mut egui::Ui) {
        match self.panel_tab {
            PanelTab::Overview => self.panel_tab_overview(ui),
            PanelTab::Add => self.panel_tab_add(ui),
            PanelTab::Templates => self.panel_tab_templates(ui),
            PanelTab::View => self.panel_tab_view(ui),
            PanelTab::Camera => self.panel_tab_camera(ui),
            PanelTab::Config => self.panel_tab_config(ui),
        }
    }
}
