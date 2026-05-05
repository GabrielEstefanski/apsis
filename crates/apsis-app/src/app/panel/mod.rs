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
mod notifications_panel;
mod playbar;
mod precision_confirmation_modal;
mod precision_panel;
mod save_modal;
mod settings_modal;
mod shortcuts_modal;
mod tabs;
mod tool_rail;
mod toolbar;

use crate::app::theme::{BORDER, PANEL_BG};
use crate::app::ui::{BodySelection, SimulationApp};
use eframe::egui::{self, Stroke};
use std::collections::BTreeSet;

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
        if matches!(self.selection, BodySelection::None) {
            self.system.set_orbital_elements_needed(false);
            return;
        }

        // Stale-index guard: prune selected indices that have gone out of range.
        let n_bodies = self.system.bodies().len();
        let stale = match &self.selection {
            BodySelection::Single(i) if *i >= n_bodies => Some(BodySelection::default()),
            BodySelection::Multi(set) => {
                let valid: BTreeSet<usize> =
                    set.iter().copied().filter(|&i| i < n_bodies).collect();
                (valid.len() < set.len()).then(|| match valid.len() {
                    0 => BodySelection::default(),
                    1 => BodySelection::Single(*valid.iter().next().unwrap()),
                    _ => BodySelection::Multi(valid),
                })
            },
            _ => None,
        };
        if let Some(sel) = stale {
            if matches!(sel, BodySelection::None) {
                self.follow_selected_body = false;
                self.selection_form = None;
                self.system.set_orbital_elements_needed(false);
            }
            self.selection = sel;
            if matches!(self.selection, BodySelection::None) {
                return;
            }
        }

        // Extract dispatch info before the closure to avoid borrow conflicts.
        let single_idx = self.selection.single();
        let multi_set: Option<BTreeSet<usize>> = match &self.selection {
            BodySelection::Multi(s) => Some(s.clone()),
            _ => None,
        };

        self.system.set_orbital_elements_needed(single_idx.is_some());

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
                    if let Some(idx) = single_idx {
                        self.inspector_content(ui, idx);
                    } else if let Some(ref indices) = multi_set {
                        self.aggregate_content(ui, indices);
                    }
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
