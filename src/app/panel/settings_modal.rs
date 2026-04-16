use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::recorder::{RecordMetadata, SimRecorder};
use eframe::egui::{self, Color32, RichText, Stroke};
use std::path::Path;

impl SimulationApp {
    pub(in crate::app) fn draw_settings_modal(&mut self, ctx: &egui::Context) {
        if !self.show_settings_modal {
            return;
        }

        let mut open = true;

        egui::Window::new("Settings")
            .id(egui::Id::new("settings_modal"))
            .collapsible(false)
            .resizable(false)
            .min_width(340.0)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .open(&mut open)
            .show(ctx, |ui| {
                ui.set_width(340.0);

                ui.label(RichText::new("SETTINGS").size(11.0).color(TEXT_PRI).strong());
                ui.separator();
                ui.add_space(4.0);

                // ── UNIT LABELS ───────────────────────────────────────────────
                ui.label(RichText::new("UNIT LABELS").size(9.5).color(TEXT_DIM).strong());
                ui.add_space(2.0);
                ui.label(
                    RichText::new("Cosmetic only — do not affect physics or recorded data.")
                        .size(9.0)
                        .color(TEXT_DIM),
                );
                ui.add_space(4.0);

                egui::Grid::new("settings_unit_labels").num_columns(2).spacing([8.0, 4.0]).show(
                    ui,
                    |ui| {
                        let lw = 40.0_f32;
                        let fw = 80.0_f32;

                        ui.add_sized(
                            egui::vec2(lw, 18.0),
                            egui::Label::new(RichText::new("mass").size(10.0).color(TEXT_SEC)),
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.physics_cfg.mass_label)
                                .desired_width(fw),
                        );
                        ui.end_row();

                        ui.add_sized(
                            egui::vec2(lw, 18.0),
                            egui::Label::new(RichText::new("dist").size(10.0).color(TEXT_SEC)),
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.physics_cfg.dist_label)
                                .desired_width(fw),
                        );
                        ui.end_row();

                        ui.add_sized(
                            egui::vec2(lw, 18.0),
                            egui::Label::new(RichText::new("time").size(10.0).color(TEXT_SEC)),
                        );
                        ui.add(
                            egui::TextEdit::singleline(&mut self.physics_cfg.time_label)
                                .desired_width(fw),
                        );
                        ui.end_row();
                    },
                );

                ui.add_space(10.0);
                ui.separator();
                ui.add_space(6.0);

                // ── RECORDING ────────────────────────────────────────────────
                ui.label(RichText::new("RECORDING").size(9.5).color(TEXT_DIM).strong());
                ui.add_space(2.0);
                ui.label(
                    RichText::new("Exports two CSV files: _bodies.csv and _system.csv")
                        .size(9.0)
                        .color(TEXT_DIM),
                );
                ui.add_space(6.0);

                let is_recording = self.recorder.is_some();

                egui::Grid::new("settings_recording").num_columns(2).spacing([8.0, 5.0]).show(
                    ui,
                    |ui| {
                        ui.label(RichText::new("output path").size(10.0).color(TEXT_SEC));
                        ui.add_enabled(
                            !is_recording,
                            egui::TextEdit::singleline(&mut self.record_base_path)
                                .desired_width(200.0),
                        );
                        ui.end_row();

                        ui.label(RichText::new("Δt_record").size(10.0).color(TEXT_SEC))
                            .on_hover_text(
                                "Simulated time gap between CSV rows.\nRule of thumb: 10–100 × dt.",
                            );
                        let speed = (self.record_interval * 0.05).max(1e-6);
                        ui.add_enabled(
                            !is_recording,
                            egui::DragValue::new(&mut self.record_interval)
                                .speed(speed)
                                .range(1e-6_f64..=1e4)
                                .max_decimals(6),
                        );
                        ui.end_row();
                    },
                );

                ui.add_space(6.0);

                if is_recording {
                    let records = self.recorder.as_ref().map(|r| r.records_written).unwrap_or(0);
                    let path = self
                        .recorder
                        .as_ref()
                        .map(|r| r.base_path.display().to_string())
                        .unwrap_or_default();

                    // Recording status badge
                    ui.horizontal(|ui| {
                        let (dot, _) =
                            ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                        ui.painter().circle_filled(dot.center(), 4.0, DANGER);
                        ui.label(RichText::new("Recording").size(10.5).color(DANGER).strong());
                        ui.label(
                            RichText::new(format!("— {records} rows")).size(10.0).color(TEXT_DIM),
                        );
                    });
                    ui.label(
                        RichText::new(format!("{path}_bodies.csv"))
                            .monospace()
                            .size(9.0)
                            .color(TEXT_DIM),
                    );

                    ui.add_space(6.0);
                    if ui
                        .add(
                            egui::Button::new(
                                RichText::new("■  Stop recording").size(10.5).color(DANGER),
                            )
                            .fill(Color32::from_rgba_unmultiplied(60, 20, 20, 80))
                            .stroke(Stroke::new(1.0, DANGER.gamma_multiply(0.5)))
                            .min_size(egui::vec2(ui.available_width(), 26.0)),
                        )
                        .clicked()
                    {
                        if let Some(mut rec) = self.recorder.take() {
                            let _ = rec.flush();
                        }
                    }
                } else {
                    if let Some(err) = &self.record_error.clone() {
                        ui.label(RichText::new(err).size(9.0).color(DANGER));
                        ui.add_space(4.0);
                    }

                    if ui
                        .add(
                            egui::Button::new(
                                RichText::new("●  Start recording").size(10.5).color(SUCCESS),
                            )
                            .fill(Color32::from_rgba_unmultiplied(20, 60, 35, 80))
                            .stroke(Stroke::new(1.0, SUCCESS.gamma_multiply(0.5)))
                            .min_size(egui::vec2(ui.available_width(), 26.0)),
                        )
                        .clicked()
                    {
                        // Auto-generate a versioned path when still using the default.
                        if self.record_base_path == "./records/sim_export"
                            || self.record_base_path.is_empty()
                        {
                            let name_slug = if self.sim_name.is_empty() {
                                "unnamed".to_owned()
                            } else {
                                self.sim_name
                                    .to_lowercase()
                                    .chars()
                                    .map(|c| if c.is_alphanumeric() { c } else { '_' })
                                    .collect::<String>()
                            };
                            let ts = std::time::SystemTime::now()
                                .duration_since(std::time::UNIX_EPOCH)
                                .map(|d| d.as_secs())
                                .unwrap_or(0);
                            self.record_base_path = format!("./records/{name_slug}_{ts}");
                        }
                        let meta = RecordMetadata {
                            n_bodies: self.system.bodies().len(),
                            integrator_label: self.physics_cfg.integrator.label(),
                            integrator_order: self.physics_cfg.integrator.order(),
                            dt: self.system.dt(),
                            theta: self.physics_cfg.theta,
                            softening_scale: self.physics_cfg.softening_scale,
                            g_factor: self.physics_cfg.g_factor,
                            record_interval: self.record_interval,
                            units: self.active_units,
                        };
                        match SimRecorder::create(
                            Path::new(&self.record_base_path),
                            self.record_interval,
                            &meta,
                        ) {
                            Ok(rec) => {
                                self.recorder = Some(rec);
                                self.record_error = None;
                            },
                            Err(e) => {
                                self.record_error = Some(format!("Failed: {e}"));
                            },
                        }
                    }
                }
            });

        if !open {
            self.show_settings_modal = false;
        }
    }
}
