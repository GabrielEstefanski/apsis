mod inspector;
mod metrics;
mod tab_bar;
mod tabs;
mod time_speed;
mod toolbar;

use crate::app::ui::SimulationApp;
use eframe::egui;

impl SimulationApp {
    // ── PUBLIC ENTRY POINTS ───────────────────────────────────────────────── //
    // These methods are called from the main update loop.
    // They own only the egui frame/shell; all content is delegated to sub-modules.

    pub(super) fn draw_toolbar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("toolbar")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(12, 5)),
            )
            .show(ctx, |ui| self.toolbar_content(ui));
    }

    pub(super) fn draw_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("controls")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(12, 10)),
            )
            .min_width(220.0)
            .max_width(220.0)
            .show(ctx, |ui| {
                ui.set_width(196.0);

                self.panel_metrics_compact(ui);
                self.panel_time_speed(ui);
                self.panel_tab_bar(ui);

                ui.add_space(4.0);

                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.set_width(196.0);
                    self.panel_tab_dispatch(ui);
                });
            });
    }

    pub(super) fn draw_inspector(&mut self, ctx: &egui::Context) {
        let idx = match self.selected_body {
            Some(i) => i,
            None => return,
        };

        if idx >= self.system.bodies().len() {
            self.selected_body = None;
            self.selection_form = None;
            return;
        }

        egui::SidePanel::right("inspector")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14, 14)),
            )
            .min_width(200.0)
            .max_width(200.0)
            .show(ctx, |ui| {
                ui.set_width(172.0);
                self.inspector_content(ui, idx);
            });
    }

    // Visible to sub-modules that need to reset the viewport after loading a
    // template or any other operation that changes the body set.
    pub(in crate::app::panel) fn fit_to_view(&mut self) {
        let bodies = self.system.bodies();

        if bodies.is_empty() {
            return;
        }

        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for b in bodies {
            min_x = min_x.min(b.x);
            max_x = max_x.max(b.x);
            min_y = min_y.min(b.y);
            max_y = max_y.max(b.y);
        }

        let width = (max_x - min_x) as f32;
        let height = (max_y - min_y) as f32;
        let size = width.max(height).max(1e-9);

        self.scale = 400.0 / (size * 1.2);

        let center_x = (min_x + max_x) as f32 * 0.5;
        let center_y = (min_y + max_y) as f32 * 0.5;
        self.offset = egui::vec2(-center_x * self.scale, -center_y * self.scale);
    }
}
