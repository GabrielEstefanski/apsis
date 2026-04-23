//! Top bar — brand, scene identity, session actions.
//!
//! Three logical zones (left → center → right):
//!   [brand │ ☰]   [scene name · N bodies]   [Load Save Clear │ Rec │ Help Sett]
//!
//! The center zone fills remaining space between the two fixed zones,
//! giving the scene name a prominent, publication-style position.
//!
//! Tool selection lives in [`crate::app::panel::tool_rail`] (vertical
//! rail on the left); this bar is session-only.

use crate::app::icons;
use crate::app::theme::{ACCENT, BORDER, DANGER, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use eframe::egui::containers::menu::MenuButton;
use eframe::egui::{self, Color32, RichText, Stroke};

/// Full toolbar height. Drives icon button heights.
const TOOLBAR_ROW_H: f32 = 28.0;
const BTN_BG: Color32 = Color32::from_rgb(20, 20, 26);
/// Approximate pixel width kept for the right action zone.
/// Update if buttons are added/removed from that zone.
const RIGHT_RESERVE: f32 = 192.0;

impl SimulationApp {
    pub(super) fn toolbar_content(&mut self, ui: &mut egui::Ui) {
        let time = ui.ctx().input(|i| i.time as f32);

        ui.horizontal(|ui| {
            ui.set_min_height(TOOLBAR_ROW_H);
            ui.spacing_mut().item_spacing.x = 6.0;

            // ── Brand ──────────────────────────────────────────────────── //
            let (logo_rect, _) =
                ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
            ui.painter().circle_filled(logo_rect.center(), 5.0, ACCENT);
            ui.painter().circle_stroke(
                logo_rect.center(),
                5.0,
                Stroke::new(0.7, Color32::from_rgba_unmultiplied(255, 244, 200, 55)),
            );

            ui.horizontal(|ui| {
                ui.spacing_mut().item_spacing.x = 4.0;
                ui.label(RichText::new("Gravity Sim").size(12.0).color(TEXT_PRI).strong());
                ui.label(RichText::new("v0.4.2").size(8.5).color(TEXT_DIM).monospace());
            });

            ui.add_space(2.0);
            vsep(ui);

            // ── Menu ────────────────────────────────────────────────────── //
            self.menu_hamburger(ui);

            // ── Center: scene identity ─────────────────────────────────── //
            let n = self.system.bodies().len();
            let center_w = (ui.available_width() - RIGHT_RESERVE).max(80.0);

            ui.allocate_ui_with_layout(
                egui::vec2(center_w, TOOLBAR_ROW_H),
                egui::Layout::left_to_right(egui::Align::Center),
                |ui| {
                    ui.spacing_mut().item_spacing.x = 5.0;

                    // Reserve space for the body count on the right of the center zone.
                    let count_w: f32 = if n > 0 { 72.0 } else { 0.0 };
                    let edit_w = (center_w - count_w - 10.0).max(60.0);

                    ui.add_sized(
                        egui::vec2(edit_w, TOOLBAR_ROW_H - 6.0),
                        egui::TextEdit::singleline(&mut self.sim_name)
                            .desired_width(edit_w)
                            .hint_text("Unnamed simulation")
                            .font(egui::FontId::proportional(11.0))
                            .text_color(TEXT_PRI),
                    );

                    if n > 0 {
                        ui.label(RichText::new("·").size(9.0).color(TEXT_DIM));
                        ui.label(
                            RichText::new(format!("{n}"))
                                .size(10.5)
                                .monospace()
                                .color(TEXT_SEC),
                        );
                        ui.label(RichText::new("bodies").size(9.0).color(TEXT_DIM));
                    }
                },
            );

            // ── Right: actions + record + utils ────────────────────────── //
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.spacing_mut().item_spacing.x = 4.0;

                if tb_icon_btn(ui, icons::SETTINGS, "Settings").clicked() {
                    self.show_settings_modal = !self.show_settings_modal;
                }
                if tb_icon_btn(ui, icons::HELP, "Keyboard shortcuts  [H]").clicked() {
                    self.show_shortcuts_modal = !self.show_shortcuts_modal;
                }

                // Notification bell with unread badge.
                let (unread, total) = {
                    let store = self.notifications.lock().unwrap();
                    (store.unread_count(), store.len())
                };
                let bell_icon = if unread > 0 { icons::BELL_ON } else { icons::BELL };
                let bell_color = if unread > 0 { ACCENT } else { TEXT_DIM };
                let bell_btn = ui
                    .add(
                        egui::Button::new(RichText::new(bell_icon).size(13.0).color(bell_color))
                            .fill(BTN_BG)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(24.0, TOOLBAR_ROW_H))
                            .corner_radius(3.0),
                    )
                    .on_hover_text(if unread > 0 {
                        format!("Notifications — {} unread / {} total", unread, total)
                    } else if total > 0 {
                        format!("Notifications — {} total", total)
                    } else {
                        "Notifications".into()
                    });
                if bell_btn.clicked() {
                    self.show_notifications_panel = !self.show_notifications_panel;
                    if self.show_notifications_panel {
                        self.notifications.lock().unwrap().mark_all_read();
                    }
                }
                if unread > 0 {
                    // Small badge glyph overlaid on the bell corner
                    // via an inline count label next to it. Kept as
                    // a sibling so egui layout handles wrapping; no
                    // absolute positioning.
                    ui.label(
                        RichText::new(if unread > 99 {
                            "99+".to_string()
                        } else {
                            format!("{}", unread)
                        })
                        .size(9.0)
                        .monospace()
                        .color(ACCENT),
                    );
                }

                vsep(ui);

                // CSV record indicator
                let is_rec = self.recorder.is_some();
                let pulse_alpha = if is_rec {
                    (((time * 2.2).sin() * 0.25 + 0.55) * 255.0) as u8
                } else {
                    75
                };
                let rec_col = Color32::from_rgba_unmultiplied(
                    DANGER.r(),
                    DANGER.g(),
                    DANGER.b(),
                    pulse_alpha,
                );
                let rec_btn = ui
                    .add(
                        egui::Button::new(
                            RichText::new(icons::RECORD).size(12.0).color(rec_col),
                        )
                        .fill(BTN_BG)
                        .stroke(Stroke::new(
                            0.5,
                            if is_rec { DANGER.gamma_multiply(0.45) } else { BORDER },
                        ))
                        .min_size(egui::vec2(24.0, TOOLBAR_ROW_H))
                        .corner_radius(3.0),
                    )
                    .on_hover_text(if is_rec {
                        "Recording CSV — click to stop"
                    } else {
                        "Record CSV data"
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

                vsep(ui);

                // Mutation-capable file actions are disabled during a
                // Precision Run. Save stays enabled — snapshotting the
                // in-flight state is a read-only operation.
                let edit_locked = self.is_editing_locked();
                let clear_hint = if edit_locked {
                    self.editing_lock_hint()
                } else {
                    "Clear all bodies"
                };
                let clear_btn = ui.add_enabled(
                    !edit_locked,
                    egui::Button::new(RichText::new(icons::CLEAR).size(13.0).color(TEXT_DIM))
                        .fill(BTN_BG)
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(24.0, TOOLBAR_ROW_H))
                        .corner_radius(3.0),
                )
                .on_hover_text(clear_hint);
                if clear_btn.clicked() {
                    self.system.load_bodies(vec![]);
                    self.paused = true;
                    self.reset_drift_peaks();
                    self.sim_name = String::new();
                }
                if tb_icon_btn(ui, icons::SAVE, "Save  Ctrl+S").clicked() {
                    let _ = self.do_save();
                }
                let load_hint = if edit_locked {
                    self.editing_lock_hint()
                } else {
                    "Load saved state"
                };
                let load_btn = ui.add_enabled(
                    !edit_locked,
                    egui::Button::new(RichText::new(icons::LOAD).size(13.0).color(TEXT_DIM))
                        .fill(BTN_BG)
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(24.0, TOOLBAR_ROW_H))
                        .corner_radius(3.0),
                )
                .on_hover_text(load_hint);
                if load_btn.clicked() {
                    self.open_save_modal();
                }
            });
        });
    }

    // ── Hamburger: global actions ────────────────────────────────────────── //

    fn menu_hamburger(&mut self, ui: &mut egui::Ui) {
        let btn = egui::Button::new(RichText::new(icons::MENU).size(13.0).color(TEXT_SEC))
            .fill(BTN_BG)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(26.0, TOOLBAR_ROW_H))
            .corner_radius(4.0);

        let _ = MenuButton::from_button(btn).ui(ui, |ui: &mut egui::Ui| {
            ui.set_min_width(190.0);

            ui.label(RichText::new("FILE").size(9.0).color(TEXT_DIM));
            if ui.add(egui::Button::new("Save").shortcut_text("Ctrl+S")).clicked() {
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

            ui.label(RichText::new("EDIT").size(9.0).color(TEXT_DIM));
            if ui.add(egui::Button::new("Undo").shortcut_text("Ctrl+Z")).clicked() {
                self.perform_undo();
                ui.close_menu();
            }

            ui.separator();

            ui.label(RichText::new("SIMULATION").size(9.0).color(TEXT_DIM));
            let play_lbl = if self.paused { "Play" } else { "Pause" };
            if ui.add(egui::Button::new(play_lbl).shortcut_text("Space")).clicked() {
                self.paused = !self.paused;
                ui.close_menu();
            }
            if ui.add(egui::Button::new("Step one frame").shortcut_text("→")).clicked() {
                self.paused = false;
                self.step_pending = true;
                ui.close_menu();
            }
            if ui.add(egui::Button::new("Fit to view").shortcut_text("F")).clicked() {
                self.fit_to_view();
                ui.close_menu();
            }
            if ui.button("Reset drift peaks").clicked() {
                self.reset_drift_peaks();
                ui.close_menu();
            }

            ui.separator();

            ui.label(RichText::new("EXPORT").size(9.0).color(TEXT_DIM));
            let is_rec = self.recorder.is_some();
            let lbl = if is_rec { "Stop recording CSV" } else { "Record CSV…" };
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

            ui.label(RichText::new("HELP").size(9.0).color(TEXT_DIM));
            if ui.add(egui::Button::new("Keyboard shortcuts").shortcut_text("H")).clicked() {
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

// ── Free helpers ──────────────────────────────────────────────────────────────

fn vsep(ui: &mut egui::Ui) {
    ui.add(egui::Separator::default().vertical());
}

fn tb_icon_btn(ui: &mut egui::Ui, icon: &str, hover: &str) -> egui::Response {
    ui.add(
        egui::Button::new(RichText::new(icon).size(13.0).color(TEXT_DIM))
            .fill(BTN_BG)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(24.0, TOOLBAR_ROW_H))
            .corner_radius(3.0),
    )
    .on_hover_text(hover)
}
