use crate::app::config::PhysicsConfig;
use crate::app::theme::secondary_btn;
use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
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
        ui.add_space(2.0);
        ui.label(RichText::new("Physics").size(13.0).color(TEXT_PRI).strong());
        ui.label(
            RichText::new("Force model, gravity & reproducibility — affects results.")
                .size(10.0)
                .color(TEXT_DIM),
        );

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

        // ── RESET ─────────────────────────────────────────────────────────────

        ui.add_space(10.0);
        if secondary_btn(ui, "Reset to defaults") {
            let defaults = PhysicsConfig::default();
            self.system.set_exact_threshold(defaults.exact_threshold);
            self.system.set_seed(defaults.seed);
            self.system.set_theta(defaults.theta);
            self.system.set_softening_scale(defaults.softening_scale);
            self.system.set_g_factor(defaults.g_factor);
            self.physics_cfg = defaults;
        }
    }
}
