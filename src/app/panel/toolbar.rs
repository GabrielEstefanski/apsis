//! Top bar — brand, hamburger menu, tool tabs, scene identity, session actions.
//!
//! Design rule: no parameter controls here. Everything that changes sim
//! behaviour belongs to the contextual panel or the bottom playbar.
//!
//! Layout zones (left → right):
//!   [brand] [☰] [tool tabs] ............ [actions | scene | rec | utils]
//!
//! The left cluster carries navigation (what panel is active) and global
//! actions (behind the hamburger). The right cluster carries session state
//! (scene name, body count) and its immediate actions (Clear / Save / Load).
//! Vertical rhythm is 24px for every interactive element.

use crate::app::icons;
use crate::app::theme::{
    ACCENT, ACCENT_DIM, BORDER, DANGER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC,
};
use crate::app::ui::{PanelTab, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

const BTN_H: f32 = 24.0;
const BTN_BG: Color32 = Color32::from_rgb(20, 20, 26);

impl SimulationApp {
    pub(super) fn toolbar_content(&mut self, ui: &mut egui::Ui) {
        let time = ui.ctx().input(|i| i.time as f32);

        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 6.0;

            // ── Zone 1: Brand ─────────────────────────────────────────────── //
            let (logo_rect, _) =
                ui.allocate_exact_size(egui::vec2(16.0, 16.0), egui::Sense::hover());
            ui.painter().circle_filled(logo_rect.center(), 6.5, ACCENT);
            ui.painter().circle_stroke(
                logo_rect.center(),
                6.5,
                Stroke::new(1.0, Color32::from_rgba_unmultiplied(255, 244, 200, 70)),
            );

            ui.label(
                RichText::new("Gravity Sim")
                    .size(13.0)
                    .color(TEXT_PRI)
                    .strong(),
            );
            ui.label(
                RichText::new("v0.4.2")
                    .size(9.5)
                    .color(TEXT_DIM)
                    .monospace(),
            );

            ui.add_space(4.0);
            ui.add(egui::Separator::default().vertical().spacing(6.0));

            // ── Zone 2: Hamburger (global actions) ────────────────────────── //
            self.menu_hamburger(ui);

            ui.add(egui::Separator::default().vertical().spacing(6.0));

            // ── Zone 3: Tool tabs ─────────────────────────────────────────── //
            ui.spacing_mut().item_spacing.x = 2.0;
            for tab in PanelTab::ALL {
                self.tool_tab_btn(ui, tab);
            }

            // Sidebar collapse toggle — last tab-adjacent control.
            ui.add_space(2.0);
            let (chevron, tip) = if self.sidebar_collapsed {
                (icons::SIDEBAR_OPEN, "Show sidebar  [B]")
            } else {
                (icons::SIDEBAR_CLOSE, "Hide sidebar  [B]")
            };
            if tb_icon_btn(ui, chevron, tip).clicked() {
                self.sidebar_collapsed = !self.sidebar_collapsed;
            }

            // ── Right side (RTL) — UNCHANGED ──────────────────────────────── //
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.spacing_mut().item_spacing.x = 6.0;

                if tb_icon_btn(ui, icons::SETTINGS, "Settings").clicked() {
                    self.show_settings_modal = !self.show_settings_modal;
                }
                if tb_icon_btn(ui, icons::HELP, "Keyboard shortcuts (H)").clicked() {
                    self.show_shortcuts_modal = !self.show_shortcuts_modal;
                }

                ui.add(egui::Separator::default().vertical().spacing(6.0));

                let is_rec = self.recorder.is_some();
                let pulse_alpha = if is_rec {
                    (((time * 2.5).sin() * 0.4 + 0.6) * 255.0) as u8
                } else {
                    80
                };
                let rec_col = Color32::from_rgba_unmultiplied(
                    DANGER.r(),
                    DANGER.g(),
                    DANGER.b(),
                    pulse_alpha,
                );
                let rec_btn = ui
                    .add(
                        egui::Button::new(RichText::new(icons::RECORD).size(13.0).color(rec_col))
                            .fill(BTN_BG)
                            .stroke(Stroke::new(
                                0.5,
                                if is_rec {
                                    DANGER.gamma_multiply(0.55)
                                } else {
                                    BORDER
                                },
                            ))
                            .min_size(egui::vec2(26.0, BTN_H))
                            .corner_radius(4.0),
                    )
                    .on_hover_text(if is_rec {
                        "Recording CSV — click to stop"
                    } else {
                        "Start CSV recording (Export menu)"
                    });
                if rec_btn.clicked() {
                    if is_rec {
                        if let Some(mut rec) = self.recorder.take() {
                            let _ = rec.flush();
                        }
                    } else {
                        self.show_settings_modal = true;
                    }
                }

                ui.add(egui::Separator::default().vertical().spacing(6.0));

                if tb_action_btn(ui, icons::LOAD, "Load", ACCENT, "Browse saved states").clicked() {
                    self.open_save_modal();
                }
                if tb_action_btn(ui, icons::SAVE, "Save", SUCCESS, "Quick-save current state")
                    .clicked()
                {
                    let _ = self.do_save();
                }
                if tb_action_btn(
                    ui,
                    icons::CLEAR,
                    "Clear",
                    DANGER,
                    "Remove all bodies and reset simulation",
                )
                .clicked()
                {
                    self.system.load_bodies(vec![]);
                    self.paused = true;
                    self.reset_drift_peaks();
                    self.sim_name = String::new();
                }

                ui.add(egui::Separator::default().vertical().spacing(6.0));

                let n = self.system.bodies().len();
                ui.label(
                    RichText::new(format!("• {n} bodies"))
                        .size(10.5)
                        .monospace()
                        .color(TEXT_DIM),
                );
                ui.add(
                    egui::TextEdit::singleline(&mut self.sim_name)
                        .desired_width(150.0)
                        .hint_text("Unnamed")
                        .font(egui::FontId::proportional(12.0))
                        .text_color(TEXT_PRI),
                );
                ui.label(
                    RichText::new("scene:")
                        .size(10.5)
                        .color(TEXT_DIM),
                );
            });
        });

        let _ = TEXT_SEC; // suppress unused warning
        let _ = PANEL_BG;
    }

    // ── Tool tab button ──────────────────────────────────────────────────── //

    fn tool_tab_btn(&mut self, ui: &mut egui::Ui, tab: PanelTab) {
        // "Active" means both currently selected AND sidebar visible.
        // When sidebar is collapsed, no tab appears visually active — keeps
        // the collapsed state unambiguous.
        let is_active = self.panel_tab == tab && !self.sidebar_collapsed;
        let (icon, label) = (tool_icon(tab, is_active), tab.label());

        let fill = if is_active { ACCENT_DIM } else { BTN_BG };
        let stroke_col = if is_active {
            ACCENT.gamma_multiply(0.6)
        } else {
            BORDER
        };
        let text_col = if is_active { TEXT_PRI } else { TEXT_DIM };

        let text = RichText::new(format!("{icon}  {label}"))
            .size(11.0)
            .color(text_col);

        let hover = if is_active {
            format!("{}  [{}]   ·   click again to hide sidebar", label, tab_shortcut(tab))
        } else if self.sidebar_collapsed {
            format!("{}  [{}]   ·   opens sidebar", label, tab_shortcut(tab))
        } else {
            format!("{}  [{}]", label, tab_shortcut(tab))
        };

        let resp = ui
            .add(
                egui::Button::new(text)
                    .fill(fill)
                    .stroke(Stroke::new(if is_active { 1.0 } else { 0.5 }, stroke_col))
                    .min_size(egui::vec2(0.0, BTN_H))
                    .corner_radius(4.0),
            )
            .on_hover_text(hover);

        if resp.clicked() {
            // Click rules:
            //   sidebar hidden            → open sidebar + switch to tab
            //   sidebar visible, other tab → switch tab
            //   sidebar visible, same tab → hide sidebar (toggle)
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

    // ── Hamburger: single entry point for all global actions ────────────── //

    fn menu_hamburger(&mut self, ui: &mut egui::Ui) {
        let trigger = RichText::new(icons::MENU).size(14.0).color(TEXT_SEC);

        ui.menu_button(trigger, |ui| {
            ui.set_min_width(190.0);

            // File
            ui.label(RichText::new("FILE").size(9.0).color(TEXT_DIM));
            if ui
                .add(egui::Button::new("Save").shortcut_text("Ctrl+S"))
                .clicked()
            {
                let _ = self.do_save();
                ui.close_menu();
            }
            if ui.button("Load…").clicked() {
                self.open_save_modal();
                ui.close_menu();
            }
            if ui.button("Clear all bodies").clicked() {
                self.system.load_bodies(vec![]);
                self.paused = true;
                self.reset_drift_peaks();
                self.sim_name = String::new();
                ui.close_menu();
            }

            ui.separator();

            // Edit
            ui.label(RichText::new("EDIT").size(9.0).color(TEXT_DIM));
            if ui
                .add(egui::Button::new("Undo").shortcut_text("Ctrl+Z"))
                .clicked()
            {
                self.perform_undo();
                ui.close_menu();
            }

            ui.separator();

            // Simulation
            ui.label(RichText::new("SIMULATION").size(9.0).color(TEXT_DIM));
            let play_lbl = if self.paused { "Play" } else { "Pause" };
            if ui
                .add(egui::Button::new(play_lbl).shortcut_text("Space"))
                .clicked()
            {
                self.paused = !self.paused;
                ui.close_menu();
            }
            if ui
                .add(egui::Button::new("Step one frame").shortcut_text("→"))
                .clicked()
            {
                self.paused = false;
                self.step_pending = true;
                ui.close_menu();
            }
            if ui
                .add(egui::Button::new("Fit to view").shortcut_text("F"))
                .clicked()
            {
                self.fit_to_view();
                ui.close_menu();
            }
            if ui.button("Reset drift peaks").clicked() {
                self.reset_drift_peaks();
                ui.close_menu();
            }

            ui.separator();

            // Export
            ui.label(RichText::new("EXPORT").size(9.0).color(TEXT_DIM));
            let is_rec = self.recorder.is_some();
            let lbl = if is_rec {
                "Stop recording CSV"
            } else {
                "Record CSV…"
            };
            if ui.button(lbl).clicked() {
                if is_rec {
                    if let Some(mut rec) = self.recorder.take() {
                        let _ = rec.flush();
                    }
                } else {
                    self.show_settings_modal = true;
                }
                ui.close_menu();
            }

            ui.separator();

            // Help
            ui.label(RichText::new("HELP").size(9.0).color(TEXT_DIM));
            if ui
                .add(egui::Button::new("Keyboard shortcuts").shortcut_text("H"))
                .clicked()
            {
                self.show_shortcuts_modal = !self.show_shortcuts_modal;
                ui.close_menu();
            }
            if ui.button("Settings").clicked() {
                self.show_settings_modal = !self.show_settings_modal;
                ui.close_menu();
            }
        });
    }
}

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

fn tb_icon_btn(ui: &mut egui::Ui, icon: &str, hover: &str) -> egui::Response {
    ui.add(
        egui::Button::new(RichText::new(icon).size(14.0).color(TEXT_DIM))
            .fill(BTN_BG)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(26.0, BTN_H))
            .corner_radius(4.0),
    )
    .on_hover_text(hover)
}

fn tb_action_btn(
    ui: &mut egui::Ui,
    icon: &str,
    label: &str,
    col: Color32,
    hover: &str,
) -> egui::Response {
    let text = RichText::new(format!("{icon}  {label}"))
        .size(11.0)
        .color(col);
    ui.add(
        egui::Button::new(text)
            .fill(BTN_BG)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(66.0, BTN_H))
            .corner_radius(4.0),
    )
    .on_hover_text(hover)
}
