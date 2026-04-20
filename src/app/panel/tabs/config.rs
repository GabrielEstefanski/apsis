use crate::app::config::PhysicsConfig;
use crate::app::theme::secondary_btn;
use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::physics::integrator::IntegratorKind;
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke};

// ── Layout constants ──────────────────────────────────────────────────────────

/// Width of all `DragValue` number inputs (px).
const DV_W: f32 = 72.0;

/// Width of all `Slider` widgets (px).
const SL_W: f32 = 120.0;

/// Width of the label column in `param_row` (px).
const LBL_W: f32 = 80.0;

// ── Section / row helpers ─────────────────────────────────────────────────────

/// Renders a labelled section divider with a hairline rule.
fn section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.label(RichText::new(title).size(9.5).color(TEXT_DIM).strong());
        ui.add_space(4.0);
        let r = ui.available_rect_before_wrap();
        ui.painter().line_segment(
            [egui::pos2(r.left(), r.center().y), egui::pos2(r.right(), r.center().y)],
            Stroke::new(0.5, BORDER),
        );
    });
    ui.add_space(3.0);
}

/// Renders a two-column parameter row: a fixed-width label on the left and a
/// right-aligned widget on the right.
///
/// `label_w` sets the label column width in pixels; the widget receives the
/// remaining horizontal space.  Returns whatever `add` returns.
fn param_row<R>(
    ui: &mut egui::Ui,
    label: &str,
    tip: &str,
    label_w: f32,
    add: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    ui.horizontal(|ui| {
        ui.add_sized(
            egui::vec2(label_w, 18.0),
            egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)),
        )
        .on_hover_text(tip);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| add(ui)).inner
    })
    .inner
}

// ── Config panel ──────────────────────────────────────────────────────────────

impl SimulationApp {
    /// Renders the **Config** side-panel tab.
    ///
    /// Exposes controls for force accuracy (Barnes–Hut θ, Plummer softening),
    /// gravity scaling, integration algorithm and time step, trail sampling,
    /// and a one-click reset to factory defaults.
    pub(super) fn panel_tab_config(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);

        // ── DIAGNOSTICS (collapsible) ─────────────────────────────────────────

        egui::CollapsingHeader::new(
            egui::RichText::new("DIAGNOSTICS").size(9.5).color(TEXT_DIM).strong(),
        )
        .default_open(true)
        .show(ui, |ui| {
            ui.add_space(2.0);
            self.panel_diagnostics_detail(ui);
            ui.add_space(4.0);
        });

        ui.add_space(4.0);

        // ── FORCE ACCURACY ────────────────────────────────────────────────────

        section(ui, "FORCE ACCURACY");

        let theta_tip = "Barnes–Hut opening angle θ.\n\
            θ → 0   exact O(N²), maximum accuracy\n\
            θ = 0.5 balanced default\n\
            θ → 1.5 fast O(N log N), less accurate\n\
            Rule of thumb: θ < 0.3 for publication-quality runs.";

        let mut theta = self.physics_cfg.theta;
        let changed = param_row(ui, "θ (Barnes–Hut)", theta_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(SL_W, 18.0),
                egui::Slider::new(&mut theta, 0.05_f64..=1.5)
                    .step_by(0.05)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.2}")),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.theta = theta;
            self.system.set_theta(theta);
        }

        let thr_tip = "Direct O(N²) threshold.\n\
            N ≤ this value → exact pairwise sum (always used for benchmarks).\n\
            N > this value → Barnes-Hut tree approximation.\n\
            Default 64.  Set to 0 (= 1 after clamp) to force BH at all N.\n\
            Set high (e.g. 10000) to force exact evaluation at all N.";

        let mut thr = self.physics_cfg.exact_threshold;
        let changed = param_row(ui, "direct N ≤", thr_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut thr).speed(1).range(1..=10_000usize),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.exact_threshold = thr;
            self.system.set_exact_threshold(thr);
        }

        let eps_tip = "Global Plummer softening scale.\n\
            Per-body default: ε = 0.02 · m^(1/3)\n\
            1.0 = default  |  > 1 suppresses singularities  |  < 1 sharper forces";

        let mut eps = self.physics_cfg.softening_scale;
        let changed = param_row(ui, "ε scale", eps_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(SL_W, 18.0),
                egui::Slider::new(&mut eps, 0.01_f64..=10.0)
                    .logarithmic(true)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.3}")),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.softening_scale = eps;
            self.system.set_softening_scale(eps);
        }

        ui.label(
            RichText::new(format!("  ε_eff = 0.02·m^⅓·{:.3}", self.physics_cfg.softening_scale))
                .size(9.0)
                .color(TEXT_DIM),
        );

        // Softening validity indicator: estimates the fractional force error
        // from the ratio ε_max / r_min and flags critical close encounters.
        {
            let m = self.system.metrics();
            let r_min = m.r_min;
            let soft_max = m.softening_max;
            let has_data = r_min < f64::MAX && r_min > 1e-30 && soft_max > 0.0;

            if has_data {
                let ratio = soft_max / r_min;

                // Fractional force error ≈ (3/2)(ε/r)² for ε/r ≪ 1.
                let force_err_pct = (1.5 * ratio * ratio * 100.0).min(9999.0);

                let (dot, color, sev_label) = if ratio > 0.3 {
                    ("▲", DANGER, "critical")
                } else if ratio > 0.1 {
                    ("▲", ACCENT, "warning")
                } else {
                    ("●", SUCCESS, "ok")
                };

                ui.add_space(2.0);
                ui.horizontal(|ui| {
                    ui.label(RichText::new(dot).size(9.0).color(color));
                    ui.label(RichText::new("softening").size(9.5).color(TEXT_DIM));
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        ui.label(RichText::new(sev_label).size(9.5).color(color).strong());
                    });
                });

                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new(format!("  ε/r_min = {ratio:.3e}")).size(9.0).color(TEXT_DIM),
                    );
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        ui.label(
                            RichText::new(format!("~{force_err_pct:.1}% err"))
                                .size(9.0)
                                .color(color),
                        );
                    });
                });

                if ratio > 0.1 {
                    egui::Frame::NONE
                        .fill(Color32::from_rgba_unmultiplied(color.r(), color.g(), color.b(), 18))
                        .stroke(Stroke::new(0.5, color.gamma_multiply(0.4)))
                        .corner_radius(3.0)
                        .inner_margin(egui::Margin::symmetric(6, 3))
                        .show(ui, |ui| {
                            ui.set_width(ui.available_width());
                            let hint = if ratio > 0.3 {
                                "force accuracy severely compromised — \
                                 reduce ε scale or increase body separations"
                            } else {
                                "close encounter detected — \
                                 reduce ε scale or switch to Yoshida 4th-order"
                            };
                            ui.add(
                                egui::Label::new(RichText::new(hint).size(9.0).color(color)).wrap(),
                            );
                        });
                } else {
                    ui.add_space(18.0);
                }
            }
        }

        // ── GRAVITY ───────────────────────────────────────────────────────────

        section(ui, "GRAVITY");

        let g_tip = "Effective gravitational constant G_eff = G₀ · factor.\n\
            1.0 = natural simulation units (default)\n\
            Scales all pairwise forces simultaneously.";

        let mut g = self.physics_cfg.g_factor;
        let changed = param_row(ui, "G multiplier", g_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(SL_W, 18.0),
                egui::Slider::new(&mut g, 0.01_f64..=100.0)
                    .logarithmic(true)
                    .show_value(true)
                    .custom_formatter(|v, _| format!("{v:.4}")),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.g_factor = g;
            self.system.set_g_factor(g);
        }

        // ── REPRODUCIBILITY ───────────────────────────────────────────────────

        section(ui, "REPRODUCIBILITY");

        let seed_tip = "Reproducibility seed.\n\
            Presets with random elements (solar system, trojans, etc.) use this\n\
            seed so the same initial conditions can be regenerated.\n\
            0 = randomised each time a template is loaded.\n\
            Any nonzero value → fully deterministic preset.";

        // Keep local copy in sync with system seed
        self.physics_cfg.seed = self.system.seed();
        let mut seed = self.physics_cfg.seed;
        let changed = param_row(ui, "seed", seed_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut seed).speed(1).range(0..=u64::MAX),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.seed = seed;
            self.system.set_seed(seed);
        }

        // ── INTEGRATION ───────────────────────────────────────────────────────

        section(ui, "INTEGRATION");

        let integ_tip = "Integration algorithm.\n\
            Velocity Verlet (2nd): standard leapfrog, 2 force evals/step.\n\
            Yoshida 4th: Forest–Ruth composition, 3 evals/step but\n\
            allows 5–10× larger Δt for equal energy conservation.\n\
            Wisdom–Holman (2nd): exact Keplerian drift + perturbation kicks;\n\
            optimal for nearly-Keplerian systems (star + planets).";

        param_row(ui, "algorithm", integ_tip, LBL_W, |ui| {
            egui::ComboBox::from_id_salt("integrator_sel")
                .selected_text(
                    RichText::new(self.physics_cfg.integrator.label()).size(10.0).color(TEXT_PRI),
                )
                .width(SL_W)
                .show_ui(ui, |ui| {
                    for variant in IntegratorKind::ALL {
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
                "  O({}) · {}F/step — {}",
                self.physics_cfg.integrator.order(),
                self.physics_cfg.integrator.force_evals_per_step(),
                self.physics_cfg.integrator.description(),
            ))
            .size(9.0)
            .color(TEXT_DIM),
        );

        // Wisdom–Holman applicability warning.
        if self.physics_cfg.integrator == IntegratorKind::WisdomHolman {
            egui::Frame::NONE
                .fill(Color32::from_rgba_unmultiplied(ACCENT.r(), ACCENT.g(), ACCENT.b(), 18))
                .stroke(Stroke::new(0.5, ACCENT.gamma_multiply(0.4)))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(6, 3))
                .show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.add(
                        egui::Label::new(
                            RichText::new(
                                "bodies[0] must be the dominant central mass. \
                             Accuracy degrades near planet–planet close encounters.",
                            )
                            .size(9.0)
                            .color(ACCENT),
                        )
                        .wrap(),
                    );
                });
        }

        ui.add_space(4.0);

        let dt_tip = "Fixed time step Δt.\n\
            Smaller Δt → better energy conservation, slower simulation.\n\
            Yoshida-4 can use 3–5× larger Δt than VV for same dE/E₀.";

        param_row(ui, "Δt", dt_tip, LBL_W, |ui| {
            let mut dt = self.system.dt();
            let speed = (dt * 0.05).max(1e-7);
            let r = ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut dt).speed(speed).range(1e-7_f64..=10.0).max_decimals(7),
            );
            if r.changed() {
                self.system.set_dt(dt);
            }
        });

        // Recommended Δt hint — Physics-justified suggestion derived from the
        // Power et al. (2003) acceleration criterion and Aarseth jerk criterion.
        // Clicking "apply" sets dt directly; the run stays fully symplectic.
        {
            let current_dt = self.system.dt();
            let rec = self.system.metrics().recommended_dt;

            if let Some(rec) = rec {
                // Ratio > 1 means the user's dt is larger (coarser) than recommended.
                let ratio = current_dt / rec;
                let (status_color, status_label) = if ratio <= 2.0 {
                    (SUCCESS, "ok")
                } else if ratio <= 10.0 {
                    (ACCENT, "coarse")
                } else {
                    (DANGER, "too large")
                };

                ui.horizontal(|ui| {
                    ui.add_space(4.0);
                    ui.label(
                        RichText::new(format!("suggested  {:.2e}", rec))
                            .size(9.0)
                            .color(TEXT_DIM),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .add_sized(
                                egui::vec2(42.0, 14.0),
                                egui::Button::new(RichText::new("apply").size(9.0)),
                            )
                            .on_hover_text(
                                "Set Δt to the recommended value (Power et al. 2003 + Aarseth criterion).\n\
                                 Keeps DtMode::Fixed — integration remains fully symplectic.",
                            )
                            .clicked()
                        {
                            self.system.set_dt(rec);
                        }
                        ui.add_space(4.0);
                        ui.label(
                            RichText::new(status_label)
                                .size(9.0)
                                .color(status_color)
                                .strong(),
                        );
                    });
                });
            } else {
                // Reserve height so layout is stable before first step.
                ui.add_space(16.0);
            }
        }

        let spf_tip = "Physics steps computed per rendered frame.\n\
            Increase to speed up simulated time without changing Δt.\n\
            Also controllable with the speed slider in the playbar.";

        param_row(ui, "steps / frame", spf_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut self.steps_per_frame).speed(1).range(1..=10_000u32),
            );
        });

        // ── TRAILS ────────────────────────────────────────────────────────────

        section(ui, "TRAILS");

        let te_tip = "Trail sampling: record one trail point every N frames.\n\
            1 = max density  |  Higher = sparser, longer-lived trails\n\
            Useful at high steps/frame to prevent trail aliasing.";

        let mut trail_every = self.system.trail_every();
        let changed = param_row(ui, "sample every", te_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut trail_every).speed(1).range(1..=256usize),
            )
            .changed()
        });
        if changed {
            self.physics_cfg.trail_every = trail_every;
            self.system.set_trail_every(trail_every);
        }

        let mr_tip = "Minimum mass ratio (body / dominant body) for a trail to be shown.\n\
            Raise to hide more bodies (e.g. 1e-4 hides asteroid-mass objects).\n\
            Lower to show trails for smaller bodies (0 = all bodies).";

        param_row(ui, "min mass ratio", mr_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut self.trail_min_mass_ratio)
                    .speed(1e-8)
                    .range(0.0_f64..=1.0)
                    .custom_formatter(|v, _| format!("{v:.1e}"))
                    .custom_parser(|s| s.parse::<f64>().ok()),
            )
        });

        // ── RESET ─────────────────────────────────────────────────────────────

        ui.add_space(10.0);
        if secondary_btn(ui, "Reset to defaults") {
            let defaults = PhysicsConfig::default();
            self.system.set_exact_threshold(defaults.exact_threshold);
            self.system.set_seed(defaults.seed);
            self.system.set_theta(defaults.theta);
            self.system.set_softening_scale(defaults.softening_scale);
            self.system.set_g_factor(defaults.g_factor);
            self.system.set_trail_every(defaults.trail_every);
            self.physics_cfg = defaults;
        }
    }
}
