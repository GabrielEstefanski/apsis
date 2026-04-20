//! Zone orchestration for the main UI.
//!
//! Layout (see `project_ui_redesign.md` in memory for the full spec):
//!
//! ```text
//! ┌ top bar (session) ──────────────────────────────────────────────┐
//! ├──┬──────────────┬──────────────────────────────┬────────────────┤
//! │🔍│              │                              │                │
//! │➕│  contextual  │         canvas (wgpu)        │   inspector    │
//! │★ │  panel       │                              │   (auto-show)  │
//! │👁│  (280 px)    │                              │   (320 px)     │
//! │⎆ │              │                              │                │
//! │⚙ │              │                              │                │
//! ├──┴──────────────┴──────────────────────────────┴────────────────┤
//! │ playbar (play/pause · dt · SPF · integrator · reset)            │
//! └─────────────────────────────────────────────────────────────────┘
//!   48 px   280 px                                   320 px
//! ```
//!
//! Panel registration order matters in egui: each panel carves space out of
//! the remaining rect. Top → bottom → left-most → next-left → right →
//! central. Keep this order stable so the canvas always occupies the
//! residual rectangle.

mod inspector;
mod metrics;
mod playbar;
mod save_modal;
mod settings_modal;
mod shortcuts_modal;
mod tabs;
mod toolbar;

use crate::app::theme::{BORDER, PANEL_BG};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Stroke};

const CONTEXTUAL_MIN: f32 = 240.0;
const CONTEXTUAL_DEFAULT: f32 = 300.0;
const INSPECTOR_WIDTH: f32 = 320.0;

impl SimulationApp {
    // ── Top bar (session) ──────────────────────────────────────────────────
    pub(super) fn draw_toolbar(&mut self, ctx: &egui::Context) {
        egui::Panel::top("toolbar")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14, 10))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .show(ctx, |ui| self.toolbar_content(ui));
    }

    // ── Contextual panel (driven by tool rail selection) ───────────────────
    pub(super) fn draw_panel(&mut self, ctx: &egui::Context) {
        egui::Panel::left("contextual")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14, 12))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .default_size(CONTEXTUAL_DEFAULT)
            .min_size(CONTEXTUAL_MIN)
            .resizable(true)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    self.panel_tab_dispatch(ui);
                });
            });
    }

    // ── Inspector (right, auto-show when a body is selected) ───────────────
    pub(super) fn draw_inspector(&mut self, ctx: &egui::Context) {
        let idx = match self.selected_body {
            Some(i) => i,
            None => return,
        };

        if idx >= self.system.bodies().len() {
            self.selected_body = None;
            self.follow_selected_body = false;
            self.selection_form = None;
            return;
        }

        egui::Panel::right("inspector")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(16, 14))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .default_size(INSPECTOR_WIDTH)
            .min_size(INSPECTOR_WIDTH)
            .max_size(INSPECTOR_WIDTH)
            .resizable(false)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    self.inspector_content(ui, idx);
                });
            });
    }

    // ── Fit view (shared by Camera tool and F shortcut) ────────────────────
    pub(in crate::app) fn fit_to_view(&mut self) {
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
        self.follow_selected_body = false;
    }
}
