//! Vertical tool rail — 48 px strip on the left that selects which
//! contextual panel is visible.
//!
//! Owns the tool-selection interaction for the whole app:
//!   * click an inactive tool → switch and open the sidebar
//!   * click the active tool   → collapse the sidebar (VS Code pattern)
//!   * press `1`..`6`          → switch to that tool and open the sidebar
//!                               (keyboard always opens — predictable)
//!
//! Shortcut wiring lives in [`crate::app::ui::SimulationApp::handle_keys`];
//! this file owns the visual feedback (hover, active, accent stripe).

use crate::app::icons;
use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, PANEL_BG, SURFACE_STRIP, TEXT_DIM, TEXT_SEC};
use crate::app::ui::{PanelTab, SimulationApp};
use eframe::egui::{self, Color32, Frame, Pos2, Rect, RichText, Sense, Stroke, Vec2};

/// Fixed rail width. 48 px is enough for a 40×40 hit target with 4 px
/// breathing room on each side; narrower than that and the active
/// accent stripe on the left edge starts to feel cramped.
const RAIL_W: f32 = 48.0;
/// Per-button hit area. Square, centred horizontally in the rail.
const BTN_SIZE: f32 = 40.0;
/// Gap after Templates and Camera — visually groups the 6 tools into
/// three clusters (context · view · config). Small enough to stay a
/// grouping cue, not a hard separator.
const GROUP_GAP: f32 = 8.0;
/// Vertical spacing between adjacent buttons inside a group.
const BTN_GAP: f32 = 2.0;
/// Width of the left-edge accent stripe shown on the active tool.
/// VS Code uses ~2 px and it reads clearly without stealing attention.
const ACTIVE_STRIPE_W: f32 = 2.0;

impl SimulationApp {
    /// Draw the tool rail. Always rendered (even when the sidebar is
    /// collapsed) — the rail itself is how the user reopens the sidebar.
    pub(in crate::app) fn draw_tool_rail(&mut self, ctx: &egui::Context) {
        egui::Panel::left("tool_rail")
            .frame(
                Frame::NONE
                    .fill(PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(0, 6))
                    .stroke(Stroke::new(0.5, BORDER)),
            )
            .default_size(RAIL_W)
            .min_size(RAIL_W)
            .max_size(RAIL_W)
            .resizable(false)
            .show(ctx, |ui| {
                ui.vertical_centered(|ui| {
                    ui.spacing_mut().item_spacing = Vec2::new(0.0, BTN_GAP);

                    // Context group — what the user looks at or builds with.
                    self.rail_button(ui, PanelTab::Overview);
                    self.rail_button(ui, PanelTab::Add);
                    self.rail_button(ui, PanelTab::Templates);

                    ui.add_space(GROUP_GAP);

                    // View group — how the scene is presented.
                    self.rail_button(ui, PanelTab::View);
                    self.rail_button(ui, PanelTab::Camera);

                    ui.add_space(GROUP_GAP);

                    // Config group — numerical knobs. Slated to move into the
                    // Settings modal in F5, but lives here for now so nothing
                    // disappears during the skeleton rewrite.
                    self.rail_button(ui, PanelTab::Config);
                });
            });
    }

    /// Allocate + render + handle one rail button.
    ///
    /// We bypass `egui::Button` and paint manually because the three
    /// visual states (idle / hover / active + accent stripe) don't map
    /// cleanly onto `Button`'s fill/stroke API — the stripe in
    /// particular needs to be painted *inside* the button rect, flush
    /// to the left edge.
    fn rail_button(&mut self, ui: &mut egui::Ui, tab: PanelTab) {
        let is_active = self.panel_tab == tab && !self.sidebar_collapsed;

        let (rect, resp) = ui.allocate_exact_size(Vec2::splat(BTN_SIZE), Sense::click());
        let hovered = resp.hovered();
        if hovered {
            ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
        }

        // Background fill — three-state.
        let fill = if is_active {
            ACCENT_DIM
        } else if hovered {
            SURFACE_STRIP
        } else {
            Color32::TRANSPARENT
        };
        let painter = ui.painter();
        painter.rect_filled(rect, 4.0, fill);

        // Active indicator — 2 px accent stripe on the left edge,
        // matching VS Code's activity bar. Renders *after* the fill so
        // it always sits on top.
        if is_active {
            let stripe = Rect::from_min_max(
                Pos2::new(rect.left(), rect.top() + 4.0),
                Pos2::new(rect.left() + ACTIVE_STRIPE_W, rect.bottom() - 4.0),
            );
            painter.rect_filled(stripe, 1.0, ACCENT);
        }

        // Icon — colour follows state.
        let icon_col = if is_active {
            ACCENT
        } else if hovered {
            TEXT_SEC
        } else {
            TEXT_DIM
        };
        let icon = tool_icon(tab, is_active);
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            icon,
            egui::FontId::proportional(17.0),
            icon_col,
        );

        // Rich two-line tooltip: name + shortcut on line 1, contextual
        // description on line 2. Plain `.on_hover_text` would collapse
        // both into one visual beat; the UI variant lets us colour the
        // description dim so the hierarchy reads at a glance.
        let (title, desc) = tooltip_text(tab);
        let shortcut = tab_shortcut(tab);
        resp.clone().on_hover_ui(|ui| {
            ui.set_max_width(260.0);
            ui.label(RichText::new(format!("{title}  [{shortcut}]")).strong().size(12.0));
            ui.label(RichText::new(desc).size(11.0).color(TEXT_SEC));
        });

        if resp.clicked() {
            // VS Code-style: clicking the active tool collapses; any
            // other click (inactive tool, or collapsed sidebar) opens
            // onto that tool.
            if self.sidebar_collapsed {
                self.sidebar_collapsed = false;
                self.panel_tab = tab;
            } else if self.panel_tab == tab {
                self.sidebar_collapsed = true;
            } else {
                self.panel_tab = tab;
            }
        }
    }

    /// Activate a tool via keyboard. Always opens the sidebar —
    /// keyboard shortcuts should be predictable, not toggles.
    pub(in crate::app) fn activate_tool(&mut self, tab: PanelTab) {
        self.panel_tab = tab;
        self.sidebar_collapsed = false;
    }
}

// ── Per-tab metadata ──────────────────────────────────────────────────────────

/// Icon glyph (active / inactive variants defined in `crate::app::icons`).
fn tool_icon(tab: PanelTab, active: bool) -> &'static str {
    match (tab, active) {
        (PanelTab::Overview, false) => icons::TOOL_OVERVIEW,
        (PanelTab::Overview, true) => icons::TOOL_OVERVIEW_ON,
        (PanelTab::Add, false) => icons::TOOL_ADD,
        (PanelTab::Add, true) => icons::TOOL_ADD_ON,
        (PanelTab::Templates, false) => icons::TOOL_TEMPLATES,
        (PanelTab::Templates, true) => icons::TOOL_TEMPLATES_ON,
        (PanelTab::View, false) => icons::TOOL_VIEW,
        (PanelTab::View, true) => icons::TOOL_VIEW_ON,
        (PanelTab::Camera, false) => icons::TOOL_CAMERA,
        (PanelTab::Camera, true) => icons::TOOL_CAMERA_ON,
        (PanelTab::Config, false) => icons::TOOL_CONFIG,
        (PanelTab::Config, true) => icons::TOOL_CONFIG_ON,
    }
}

/// Two-line tooltip content: `(title, description)`. The shortcut is
/// appended to the title at paint time so we don't duplicate it here.
fn tooltip_text(tab: PanelTab) -> (&'static str, &'static str) {
    match tab {
        PanelTab::Overview => ("Overview", "Scene summary, masses, energy diagnostics"),
        PanelTab::Add => ("Add", "Create body — form, template, or click-place"),
        PanelTab::Templates => ("Templates", "Load preset systems (Solar, TRAPPIST, …)"),
        PanelTab::View => ("View", "Toggle grid, trails, orbits, labels"),
        PanelTab::Camera => ("Camera", "Fit to view, follow selected, re-center"),
        PanelTab::Config => ("Config", "Integrator, θ, G (moves to Settings later)"),
    }
}

fn tab_shortcut(tab: PanelTab) -> &'static str {
    match tab {
        PanelTab::Overview => "1",
        PanelTab::Add => "2",
        PanelTab::Templates => "3",
        PanelTab::View => "4",
        PanelTab::Camera => "5",
        PanelTab::Config => "6",
    }
}
