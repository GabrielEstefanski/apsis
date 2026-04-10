use crate::app::config::PhysicsConfig;
use crate::app::theme::{secondary_btn, ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::recorder::{RecordMetadata, SimRecorder};
use crate::physics::integrator::Integrator;
use eframe::egui::{self, RichText, Stroke};
use std::path::Path;

// ── Section heading helper ────────────────────────────────────────────────────

fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.label(RichText::new(title).size(9.5).color(TEXT_DIM).strong());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let rect = ui.max_rect();
            let y = rect.center().y;
            ui.painter().line_segment(
                [egui::pos2(ui.min_rect().left(), y), egui::pos2(rect.right(), y)],
                Stroke::new(0.5, BORDER),
            );
        });
    });
    ui.add_space(3.0);
}

fn param_row<R>(ui: &mut egui::Ui, label: &str, tip: &str, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let r = ui.horizontal(|ui| {
        ui.add_sized(
            egui::vec2(72.0, 0.0),
            egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)),
        );
        add(ui)
    });
    r.response.on_hover_text(tip);
    ui.add_space(1.0);
    r.inner
}

impl SimulationApp {
    pub(super) fn panel_tab_config(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);

        // ── FORCE ACCURACY ──────────────────────────────────────────────────
        section(ui, "FORCE ACCURACY");

        // θ — Barnes–Hut opening angle
        let theta_tip = "Barnes–Hut opening angle θ.\n\
            θ → 0   exact O(N²), maximum accuracy\n\
            θ = 0.5 balanced (default)\n\
            θ → 1.5 fast O(N log N), less accurate\n\
            Rule of thumb: θ < 0.3 for publication-quality runs.";

        let mut theta = self.physics_cfg.theta;
        let changed = param_row(ui, "θ (Barnes–Hut)", theta_tip, |ui| {
            ui.add(
                egui::Slider::new(&mut theta, 0.05_f64..=1.5)
                    .step_by(0.05)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.2}"))
                    .text(""),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.theta = theta;
            self.system.set_theta(theta);
        }

        // Accuracy hint
        let hint = if theta < 0.2 {
            ("≈ O(N²)", ACCENT)
        } else if theta < 0.5 {
            ("high accuracy", ACCENT)
        } else if theta < 0.8 {
            ("balanced", TEXT_SEC)
        } else {
            ("approximate", TEXT_DIM)
        };
        ui.label(
            RichText::new(format!("  → {}", hint.0))
                .size(9.0)
                .color(hint.1),
        );

        ui.add_space(4.0);

        // ε — Plummer softening scale
        let eps_tip = "Global Plummer softening scale.\n\
            Per-body default: ε = 0.02 · m^(1/3)\n\
            This multiplier scales that value for all bodies.\n\
            1.0 = default (physically motivated)\n\
            > 1  suppresses close-encounter singularities\n\
            < 1  sharper forces, physically closer to point-mass\n\
            Changes take effect immediately on all bodies.";

        let mut eps = self.physics_cfg.softening_scale;
        let changed = param_row(ui, "ε scale", eps_tip, |ui| {
            ui.add(
                egui::Slider::new(&mut eps, 0.01_f64..=10.0)
                    .logarithmic(true)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.3}"))
                    .text(""),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.softening_scale = eps;
            self.system.set_softening_scale(eps);
        }

        ui.label(
            RichText::new(format!(
                "  ε_eff = 0.02 · m^(1/3) · {:.3}",
                self.physics_cfg.softening_scale
            ))
            .size(9.0)
            .color(TEXT_DIM),
        );

        // ── GRAVITY ─────────────────────────────────────────────────────────
        section(ui, "GRAVITY");

        let g_tip = "Effective gravitational constant G_eff = G₀ · factor.\n\
            1.0 = natural simulation units (default)\n\
            Scales all pairwise forces simultaneously.\n\
            Useful for studying virial scaling or non-standard cosmologies.";

        let mut g = self.physics_cfg.g_factor;
        let changed = param_row(ui, "G multiplier", g_tip, |ui| {
            ui.add(
                egui::Slider::new(&mut g, 0.01_f64..=100.0)
                    .logarithmic(true)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.4}"))
                    .text(""),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.g_factor = g;
            self.system.set_g_factor(g);
        }

        ui.label(
            RichText::new(format!("  G_eff = {:.6}", self.physics_cfg.g_factor))
                .size(9.0)
                .color(TEXT_DIM),
        );

        // ── INTEGRATION ──────────────────────────────────────────────────────
        section(ui, "INTEGRATION");

        // Integrator selector
        let integ_tip = "Integration algorithm.\n\
            Velocity Verlet (2nd-order): standard leapfrog, 2 force evals/step.\n\
            Yoshida 4th-order: Forest–Ruth composition, 3 force evals/step but\n\
            allows 5–10× larger Δt for equal energy conservation.";

        param_row(ui, "algorithm", integ_tip, |ui| {
            egui::ComboBox::from_id_salt("integrator_sel")
                .selected_text(self.physics_cfg.integrator.label())
                .width(130.0)
                .show_ui(ui, |ui| {
                    for variant in [Integrator::VelocityVerlet, Integrator::Yoshida4] {
                        let r = ui.selectable_value(
                            &mut self.physics_cfg.integrator,
                            variant,
                            variant.label(),
                        );
                        if r.clicked() {
                            self.system.set_integrator(variant);
                        }
                    }
                });
        });

        ui.label(
            RichText::new(format!(
                "  order {} · {} F-eval/step",
                self.physics_cfg.integrator.order(),
                self.physics_cfg.integrator.force_evals_per_step(),
            ))
            .size(9.0)
            .color(TEXT_DIM),
        );
        ui.label(
            RichText::new(format!("  {}", self.physics_cfg.integrator.description()))
                .size(8.5)
                .color(TEXT_DIM)
                .italics(),
        );

        ui.add_space(4.0);

        let dt_tip = "Fixed time step Δt.\n\
            Smaller Δt → better energy conservation, slower simulation.\n\
            Larger Δt → faster but energy drift increases.\n\
            With Yoshida-4 you can typically use 3–5× larger Δt than VV\n\
            for the same dE/E₀. Monitor the metrics panel.";

        param_row(ui, "Δt", dt_tip, |ui| {
            let mut dt = self.system.dt();
            let speed = (dt * 0.05).max(1e-7);
            let r = ui.add(
                egui::DragValue::new(&mut dt)
                    .speed(speed)
                    .range(1e-7_f64..=10.0)
                    .max_decimals(7),
            );
            if r.changed() {
                self.system.set_dt(dt);
            }
        });

        let spf_tip = "Physics steps computed per rendered frame.\n\
            Increase to speed up simulated time without changing Δt.\n\
            Very high values may cause UI lag.";

        param_row(ui, "steps / frame", spf_tip, |ui| {
            ui.add(
                egui::DragValue::new(&mut self.steps_per_frame)
                    .speed(1)
                    .range(1..=10000u32),
            );
        });

        // ── TRAILS ──────────────────────────────────────────────────────────
        section(ui, "TRAILS");

        let te_tip = "Trail sampling interval: record a trail point every N frames.\n\
            1 = maximum trail density (one sample per rendered frame)\n\
            Higher values produce sparser but longer-lived trails.\n\
            Useful when steps_per_frame is high — prevents trail aliasing.";

        let mut trail_every = self.system.trail_every();
        let changed = param_row(ui, "sample every", te_tip, |ui| {
            ui.add(
                egui::DragValue::new(&mut trail_every)
                    .speed(1)
                    .range(1..=256usize),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.trail_every = trail_every;
            self.system.set_trail_every(trail_every);
        }

        ui.label(
            RichText::new(format!("  record every {} frame(s)", trail_every))
                .size(9.0)
                .color(TEXT_DIM),
        );

        // ── UNIT LABELS ──────────────────────────────────────────────────────
        section(ui, "UNIT LABELS");

        ui.label(
            RichText::new("Cosmetic only — do not affect physics.")
                .size(9.0)
                .color(TEXT_DIM),
        );
        ui.add_space(3.0);

        egui::Grid::new("unit_labels")
            .num_columns(2)
            .spacing([4.0, 3.0])
            .show(ui, |ui| {
                ui.label(RichText::new("mass").size(10.0).color(TEXT_SEC));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.mass_label)
                        .desired_width(48.0),
                );
                ui.end_row();

                ui.label(RichText::new("dist").size(10.0).color(TEXT_SEC));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.dist_label)
                        .desired_width(48.0),
                );
                ui.end_row();

                ui.label(RichText::new("time").size(10.0).color(TEXT_SEC));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.time_label)
                        .desired_width(48.0),
                );
                ui.end_row();
            });

        // ── RECORDING ───────────────────────────────────────────────────────
        section(ui, "RECORDING");

        let is_recording = self.recorder.is_some();

        // Output path
        param_row(ui, "output path", "Base file path (no extension).\nTwo files will be created:\n  <path>_bodies.csv\n  <path>_system.csv", |ui| {
            ui.add_enabled(
                !is_recording,
                egui::TextEdit::singleline(&mut self.record_base_path).desired_width(140.0),
            );
        });

        // Sample interval
        let interval_tip = "Simulated time between successive CSV records.\n\
            Smaller → denser data, larger files.\n\
            Larger → sparser data, suitable for long runs.\n\
            Rule of thumb: set to ~10–100 × Δt for smooth plots.";
        param_row(ui, "Δt_record", interval_tip, |ui| {
            let speed = (self.record_interval * 0.05).max(1e-6);
            ui.add_enabled(
                !is_recording,
                egui::DragValue::new(&mut self.record_interval)
                    .speed(speed)
                    .range(1e-6_f64..=1e4)
                    .max_decimals(6),
            );
        });

        ui.add_space(4.0);

        if is_recording {
            let records = self.recorder.as_ref().map(|r| r.records_written).unwrap_or(0);
            let path = self.recorder.as_ref().map(|r| r.base_path.display().to_string()).unwrap_or_default();
            ui.label(RichText::new(format!("  recording… {records} records")).size(9.0).color(SUCCESS));
            ui.label(RichText::new(format!("  {path}_bodies.csv")).size(8.5).color(TEXT_DIM));

            if secondary_btn(ui, "Stop recording") {
                if let Some(mut rec) = self.recorder.take() {
                    let _ = rec.flush();
                }
            }
        } else {
            if let Some(err) = &self.record_error.clone() {
                ui.label(RichText::new(format!("  {err}")).size(9.0).color(DANGER));
            }

            if secondary_btn(ui, "Start recording") {
                let meta = RecordMetadata {
                    n_bodies: self.system.bodies().len(),
                    integrator_label: self.physics_cfg.integrator.label(),
                    integrator_order: self.physics_cfg.integrator.order(),
                    dt: self.system.dt(),
                    theta: self.physics_cfg.theta,
                    softening_scale: self.physics_cfg.softening_scale,
                    g_factor: self.physics_cfg.g_factor,
                    record_interval: self.record_interval,
                };
                match SimRecorder::create(Path::new(&self.record_base_path), self.record_interval, &meta) {
                    Ok(rec) => {
                        self.recorder = Some(rec);
                        self.record_error = None;
                    }
                    Err(e) => {
                        self.record_error = Some(format!("Failed to create files: {e}"));
                    }
                }
            }
        }

        // ── RESET ────────────────────────────────────────────────────────────
        ui.add_space(10.0);

        if secondary_btn(ui, "Reset to defaults") {
            let defaults = PhysicsConfig::default();
            self.system.set_theta(defaults.theta);
            self.system.set_softening_scale(defaults.softening_scale);
            self.system.set_g_factor(defaults.g_factor);
            self.system.set_trail_every(defaults.trail_every);
            self.physics_cfg = defaults;
        }
    }
}
