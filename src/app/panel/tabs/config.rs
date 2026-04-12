use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::theme::secondary_btn;
use crate::app::config::PhysicsConfig;
use crate::app::ui::SimulationApp;
use crate::physics::integrator::Integrator;
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke};

// ── Helpers ───────────────────────────────────────────────────────────────────

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

/// A label + fixed-width right-aligned widget row.
/// `label_w` is the label column width; the widget gets the rest.
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
        ).on_hover_text(tip);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            add(ui)
        }).inner
    }).inner
}

// The DragValue width we want all number inputs to be (px).
const DV_W: f32 = 72.0;
// The Slider width.
const SL_W: f32 = 120.0;
// The label column.
const LBL_W: f32 = 80.0;

// ── Config tab ────────────────────────────────────────────────────────────────

impl SimulationApp {
    pub(super) fn panel_tab_config(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);

        // ── FORCE ACCURACY ────────────────────────────────────────────────────
        section(ui, "FORCE ACCURACY");

        // θ — Barnes-Hut opening angle
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
            ).changed()
        });
        if changed {
            self.physics_cfg.theta = theta;
            self.system.set_theta(theta);
        }

        // ε — Plummer softening scale
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
            ).changed()
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

        // Softening validity indicator
        {
            let m = self.system.metrics();
            let r_min = m.r_min;
            let soft_max = m.softening_max;
            let has_data = r_min < f64::MAX && r_min > 1e-30 && soft_max > 0.0;

            if has_data {
                let ratio = soft_max / r_min;
                // Fractional force error ≈ (3/2)(ε/r)² for ε/r << 1
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
                        ui.label(
                            RichText::new(sev_label)
                                .size(9.5)
                                .color(color)
                                .strong(),
                        );
                    });
                });

                ui.horizontal(|ui| {
                    ui.label(RichText::new(format!("  ε/r_min = {ratio:.3e}")).size(9.0).color(TEXT_DIM));
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
                        .fill(Color32::from_rgba_unmultiplied(
                            color.r(), color.g(), color.b(), 18,
                        ))
                        .stroke(Stroke::new(0.5, color.gamma_multiply(0.4)))
                        .corner_radius(3.0)
                        .inner_margin(egui::Margin::symmetric(6, 3))
                        .show(ui, |ui| {
                            ui.set_width(ui.available_width());
                            let hint = if ratio > 0.3 {
                                "force accuracy severely compromised — reduce ε scale or increase body separations"
                            } else {
                                "close encounter detected — reduce ε scale or switch to Yoshida 4th-order"
                            };
                            ui.add(egui::Label::new(
                                RichText::new(hint).size(9.0).color(color),
                            ).wrap());
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
            ).changed()
        });
        if changed {
            self.physics_cfg.g_factor = g;
            self.system.set_g_factor(g);
        }

        // ── INTEGRATION ───────────────────────────────────────────────────────
        section(ui, "INTEGRATION");

        let integ_tip = "Integration algorithm.\n\
            Velocity Verlet (2nd): standard leapfrog, 2 force evals/step.\n\
            Yoshida 4th: Forest–Ruth composition, 3 evals/step but\n\
            allows 5–10× larger Δt for equal energy conservation.";

        param_row(ui, "algorithm", integ_tip, LBL_W, |ui| {
            egui::ComboBox::from_id_salt("integrator_sel")
                .selected_text(
                    RichText::new(self.physics_cfg.integrator.label())
                        .size(10.0)
                        .color(TEXT_PRI),
                )
                .width(SL_W)
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
                "  O({}) · {}F/step — {}",
                self.physics_cfg.integrator.order(),
                self.physics_cfg.integrator.force_evals_per_step(),
                self.physics_cfg.integrator.description(),
            ))
            .size(9.0)
            .color(TEXT_DIM),
        );

        ui.add_space(4.0);

        let dt_tip = "Fixed time step Δt.\n\
            Smaller Δt → better energy conservation, slower simulation.\n\
            Yoshida-4 can use 3–5× larger Δt than VV for same dE/E₀.";

        param_row(ui, "Δt", dt_tip, LBL_W, |ui| {
            let mut dt = self.system.dt();
            let speed = (dt * 0.05).max(1e-7);
            let r = ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut dt)
                    .speed(speed)
                    .range(1e-7_f64..=10.0)
                    .max_decimals(7),
            );
            if r.changed() { self.system.set_dt(dt); }
        });

        let spf_tip = "Physics steps computed per rendered frame.\n\
            Increase to speed up simulated time without changing Δt.\n\
            Also controllable with the speed slider in the playbar.";

        param_row(ui, "steps / frame", spf_tip, LBL_W, |ui| {
            ui.add_sized(
                egui::vec2(DV_W, 18.0),
                egui::DragValue::new(&mut self.steps_per_frame)
                    .speed(1)
                    .range(1..=10_000u32),
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
                egui::DragValue::new(&mut trail_every)
                    .speed(1)
                    .range(1..=256usize),
            ).changed()
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
            self.system.set_theta(defaults.theta);
            self.system.set_softening_scale(defaults.softening_scale);
            self.system.set_g_factor(defaults.g_factor);
            self.system.set_trail_every(defaults.trail_every);
            self.physics_cfg = defaults;
        }
    }
}
