use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::snapshot::SimSnapshot;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(in crate::app) fn draw_save_modal(&mut self, ctx: &egui::Context) {
        if !self.show_save_modal {
            return;
        }

        // ── Backdrop ────────────────────────────────────────────────────────
        egui::Area::new(egui::Id::new("save_modal_backdrop"))
            .fixed_pos(egui::pos2(0.0, 0.0))
            .order(egui::Order::Background)
            .show(ctx, |ui| {
                let screen = ctx.screen_rect();
                ui.painter().rect_filled(screen, 0.0, Color32::from_black_alpha(160));
            });

        // ── Modal window ────────────────────────────────────────────────────
        let mut open = true;

        egui::Window::new("Saves")
            .id(egui::Id::new("save_modal"))
            .collapsible(false)
            .resizable(false)
            .min_width(420.0)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .open(&mut open)
            .show(ctx, |ui| {
                ui.set_width(420.0);

                // ── Header row ──────────────────────────────────────────────
                ui.horizontal(|ui| {
                    ui.label(RichText::new("SAVE / LOAD").size(11.0).color(TEXT_PRI).strong());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(
                            RichText::new(format!("dir: {}", self.save_dir))
                                .size(9.0)
                                .color(TEXT_DIM),
                        );
                    });
                });

                ui.separator();

                // ── Save-dir & autosave controls ────────────────────────────
                ui.horizontal(|ui| {
                    ui.label(RichText::new("save dir").size(10.0).color(TEXT_SEC));
                    if ui
                        .add(
                            egui::TextEdit::singleline(&mut self.save_dir)
                                .desired_width(180.0),
                        )
                        .changed()
                    {
                        // Refresh listing when path changes
                        self.save_modal_entries = crate::core::snapshot::list_saves(
                            std::path::Path::new(&self.save_dir),
                        );
                    }

                    ui.label(RichText::new("auto (s)").size(10.0).color(TEXT_SEC));
                    ui.add(
                        egui::DragValue::new(&mut self.autosave_interval_secs)
                            .speed(10.0)
                            .range(0.0_f64..=3600.0)
                            .max_decimals(0),
                    )
                    .on_hover_text("Auto-save interval in real seconds. 0 = disabled.");
                });

                ui.add_space(4.0);

                // ── Manual save button ──────────────────────────────────────
                if ui
                    .add(
                        egui::Button::new(RichText::new("  Save now  ").size(11.0).color(SUCCESS))
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(1.0, SUCCESS))
                            .min_size(egui::vec2(100.0, 24.0)),
                    )
                    .clicked()
                {
                    match self.do_save() {
                        Ok(_) => {
                            self.save_modal_error = None;
                            // Refresh listing
                            self.save_modal_entries = crate::core::snapshot::list_saves(
                                std::path::Path::new(&self.save_dir),
                            );
                        }
                        Err(e) => self.save_modal_error = Some(e),
                    }
                }

                if let Some(err) = &self.save_modal_error.clone() {
                    ui.label(RichText::new(err).size(10.0).color(DANGER));
                }

                ui.add_space(6.0);
                ui.separator();

                // ── Pending load confirmation ────────────────────────────────
                if let Some(snap) = &self.pending_load {
                    let snap = snap.clone();
                    ui.add_space(4.0);
                    ui.label(
                        RichText::new(format!(
                            "Load \"{}\" — t={:.4e}, {} bodies?",
                            if snap.sim_name.is_empty() { "Unnamed" } else { &snap.sim_name },
                            snap.t,
                            snap.bodies.len()
                        ))
                        .size(10.5)
                        .color(ACCENT),
                    );
                    ui.horizontal(|ui| {
                        if ui
                            .add(
                                egui::Button::new(RichText::new("Confirm load").size(10.5).color(SUCCESS))
                                    .fill(Color32::TRANSPARENT)
                                    .stroke(Stroke::new(1.0, SUCCESS)),
                            )
                            .clicked()
                        {
                            self.system.restore_from_snapshot(&snap);
                            // Sync app UI state back from snapshot
                            self.physics_cfg.integrator      = snap.integrator;
                            self.physics_cfg.theta           = snap.theta;
                            self.physics_cfg.softening_scale = snap.softening_scale;
                            self.physics_cfg.g_factor        = snap.g_factor;
                            self.physics_cfg.trail_every     = snap.trail_every;
                            self.sim_name = snap.sim_name.clone();
                            // Restore seed (0 = old save without seed → generate fresh)
                            self.sim_seed = if snap.seed != 0 { snap.seed } else { crate::core::snapshot::SimSnapshot::new_seed() };
                            self.paused = true;
                            self.selected_body = None;
                            self.follow_selected_body = false;
                            self.selection_form = None;
                            self.pending_load = None;
                            self.pending_fit = true;
                            self.show_save_modal = false;
                            self.reset_drift_peaks();
                        }
                        if ui
                            .add(
                                egui::Button::new(RichText::new("Cancel").size(10.5).color(TEXT_DIM))
                                    .fill(Color32::TRANSPARENT)
                                    .stroke(Stroke::new(0.5, BORDER)),
                            )
                            .clicked()
                        {
                            self.pending_load = None;
                        }
                    });
                    ui.separator();
                }

                // ── Save list ────────────────────────────────────────────────
                ui.add_space(4.0);
                ui.label(RichText::new("SAVED STATES").size(9.5).color(TEXT_DIM).strong());
                ui.add_space(2.0);

                if self.save_modal_entries.is_empty() {
                    ui.label(
                        RichText::new("  No saves found in this directory.")
                            .size(10.0)
                            .color(TEXT_DIM)
                            .italics(),
                    );
                } else {
                    egui::ScrollArea::vertical()
                        .max_height(280.0)
                        .show(ui, |ui| {
                            ui.set_width(400.0);

                            // Collect paths for delete first to avoid borrow issues
                            let mut delete_idx: Option<usize> = None;
                            let mut load_idx: Option<usize> = None;

                            for (i, entry) in self.save_modal_entries.iter().enumerate() {
                                ui.horizontal(|ui| {
                                    // Name + date + seed column
                                    ui.vertical(|ui| {
                                        ui.label(
                                            RichText::new(entry.display_name())
                                                .size(10.5)
                                                .color(TEXT_PRI)
                                                .strong(),
                                        );
                                        ui.label(
                                            RichText::new(entry.display_date())
                                                .size(9.0)
                                                .color(TEXT_DIM)
                                                .monospace(),
                                        );
                                        if entry.seed != 0 {
                                            ui.label(
                                                RichText::new(format!("seed {}", entry.seed))
                                                    .size(9.0)
                                                    .color(TEXT_DIM)
                                                    .monospace(),
                                            ).on_hover_text("Reproducibility seed — share this to reproduce the initial state");
                                        }
                                    });

                                    // Stats
                                    ui.label(
                                        RichText::new(format!(
                                            "t={:.3e}  {} steps  {} bodies",
                                            entry.t, entry.steps, entry.n_bodies
                                        ))
                                        .size(9.5)
                                        .color(TEXT_SEC),
                                    );

                                    ui.with_layout(
                                        egui::Layout::right_to_left(egui::Align::Center),
                                        |ui| {
                                            if ui
                                                .add(
                                                    egui::Button::new(
                                                        RichText::new("✕").size(9.5).color(DANGER),
                                                    )
                                                    .fill(Color32::TRANSPARENT)
                                                    .stroke(Stroke::new(0.5, DANGER))
                                                    .min_size(egui::vec2(20.0, 18.0)),
                                                )
                                                .on_hover_text("Delete this save")
                                                .clicked()
                                            {
                                                delete_idx = Some(i);
                                            }

                                            if ui
                                                .add(
                                                    egui::Button::new(
                                                        RichText::new("Load").size(9.5).color(ACCENT),
                                                    )
                                                    .fill(Color32::TRANSPARENT)
                                                    .stroke(Stroke::new(0.5, ACCENT))
                                                    .min_size(egui::vec2(38.0, 18.0)),
                                                )
                                                .clicked()
                                            {
                                                load_idx = Some(i);
                                            }
                                        },
                                    );
                                });

                                ui.add(egui::Separator::default().spacing(1.0));
                            }

                            // Process load
                            if let Some(i) = load_idx {
                                let path = self.save_modal_entries[i].path.clone();
                                match SimSnapshot::load_from(&path) {
                                    Ok(snap) => {
                                        self.pending_load = Some(snap);
                                        self.save_modal_error = None;
                                    }
                                    Err(e) => {
                                        self.save_modal_error =
                                            Some(format!("Load failed: {e}"));
                                    }
                                }
                            }

                            // Process delete
                            if let Some(i) = delete_idx {
                                let path = self.save_modal_entries[i].path.clone();
                                if std::fs::remove_file(&path).is_ok() {
                                    self.save_modal_entries.remove(i);
                                }
                            }
                        });
                }
            });

        if !open {
            self.show_save_modal = false;
            self.pending_load = None;
        }
    }
}
