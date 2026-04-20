//! Bottom playbar — compact transport + temporal readouts.
//!
//! Groups (left → right):
//!   [▶ ⏭ ⟲]  |  T value · yr  |  DT  |  SPEED [slider] ×N  |  INTEGRATOR  |  Δ ENERGY  |  [trails]

use crate::app::icons;
use crate::app::theme::{
    ACCENT, ACCENT_DIM, BORDER, DANGER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC,
};
use crate::app::ui::SimulationApp;
use crate::physics::integrator::IntegratorKind;
use eframe::egui::{self, Color32, RichText, Stroke};
use std::f64::consts::PI;

pub const PLAYBAR_HEIGHT: f32 = 40.0;

const WARN: Color32 = Color32::from_rgb(210, 160, 50);

impl SimulationApp {
    pub(in crate::app) fn draw_playbar(&mut self, ctx: &egui::Context) {
        let time = ctx.input(|i| i.time as f32);

        egui::Panel::bottom("playbar")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .stroke(Stroke::new(0.5, BORDER))
                    .inner_margin(egui::Margin::symmetric(12, 0)),
            )
            .default_size(PLAYBAR_HEIGHT)
            .min_size(PLAYBAR_HEIGHT)
            .max_size(PLAYBAR_HEIGHT)
            .resizable(false)
            .show(ctx, |ui| {
                ui.horizontal_centered(|ui| {
                    ui.spacing_mut().item_spacing.x = 6.0;

                    // ── Transport ────────────────────────────────────────────
                    self.play_button(ui, time);
                    self.step_button(ui);
                    self.reset_btn(ui);

                    vsep(ui);

                    // ── Time ─────────────────────────────────────────────────
                    let t = self.system.t();
                    ui.label(
                        RichText::new("T").size(8.5).color(TEXT_DIM).strong(),
                    );
                    ui.label(
                        RichText::new(fmt_sci_signed(t, 5))
                            .size(10.5)
                            .monospace()
                            .color(TEXT_SEC),
                    );
                    ui.label(RichText::new("·").size(10.0).color(TEXT_DIM));
                    ui.label(
                        RichText::new(fmt_years(t))
                            .size(10.0)
                            .monospace()
                            .color(TEXT_DIM),
                    );

                    vsep(ui);

                    // ── dt ───────────────────────────────────────────────────
                    ui.label(RichText::new("DT").size(8.5).color(TEXT_DIM).strong());
                    let mut dt = self.system.dt();
                    let dt_speed = (dt * 0.05).max(1e-7);
                    if ui
                        .add(
                            egui::DragValue::new(&mut dt)
                                .speed(dt_speed)
                                .range(1e-7_f64..=10.0)
                                .max_decimals(6)
                                .min_decimals(1),
                        )
                        .on_hover_text(
                            "Integration timestep.\nSmaller = more accurate but slower.",
                        )
                        .changed()
                    {
                        self.system.set_dt(dt);
                    }

                    vsep(ui);

                    // ── Speed (steps-per-frame) ───────────────────────────────
                    ui.label(RichText::new("SPEED").size(8.5).color(TEXT_DIM).strong());
                    let mut spf_f = self.steps_per_frame as f32;
                    if ui
                        .add_sized(
                            [90.0, 14.0],
                            egui::Slider::new(&mut spf_f, 1.0..=10_000.0)
                                .logarithmic(true)
                                .show_value(false),
                        )
                        .changed()
                    {
                        self.steps_per_frame = spf_f.round().max(1.0) as u32;
                    }
                    let spf_col = if self.steps_per_frame > 1 { ACCENT } else { TEXT_DIM };
                    ui.label(
                        RichText::new(format!("×{}", self.steps_per_frame))
                            .monospace()
                            .size(10.0)
                            .color(spf_col),
                    )
                    .on_hover_text("Sub-steps rendered per frame");

                    vsep(ui);

                    // ── Integrator ────────────────────────────────────────────
                    ui.label(
                        RichText::new("INTEGRATOR").size(8.5).color(TEXT_DIM).strong(),
                    );
                    let current = self.system.integrator_kind();
                    egui::ComboBox::from_id_salt("playbar_integrator")
                        .selected_text(
                            RichText::new(short_label(current))
                                .size(10.5)
                                .color(TEXT_PRI),
                        )
                        .width(140.0)
                        .show_ui(ui, |ui| {
                            for kind in IntegratorKind::ALL {
                                let r = ui.selectable_label(current == kind, kind.label());
                                if r.clicked() && current != kind {
                                    self.system.set_integrator(kind);
                                    self.physics_cfg.integrator = kind;
                                }
                                r.on_hover_text(kind.description());
                            }
                        });

                    vsep(ui);

                    // ── Δ Energy ──────────────────────────────────────────────
                    let m = self.system.metrics();
                    let de = m.rel_energy_error;
                    ui.label(
                        RichText::new("Δ ENERGY").size(8.5).color(TEXT_DIM).strong(),
                    );
                    ui.label(
                        RichText::new(fmt_sci_signed(de, 3))
                            .size(10.5)
                            .monospace()
                            .color(energy_color(de)),
                    )
                    .on_hover_text(format!(
                        "Relative energy drift: ΔE / E₀\nPeak this run: {:.2e}",
                        self.energy_drift_peak
                    ));

                    // ── Right: trails toggle ──────────────────────────────────
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let col = if self.show_trails { ACCENT } else { TEXT_DIM };
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new(format!("{} TRAILS", icons::RESET))
                                        .size(9.5)
                                        .color(col),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, col.gamma_multiply(0.45)))
                                .min_size(egui::vec2(72.0, 22.0))
                                .corner_radius(3.0),
                            )
                            .on_hover_text("Toggle trail visibility")
                            .clicked()
                        {
                            self.show_trails = !self.show_trails;
                        }
                    });
                });
            });
    }

    fn play_button(&mut self, ui: &mut egui::Ui, time: f32) {
        let (icon, icon_col, fill_col) = if self.paused {
            (icons::PLAY, SUCCESS, ACCENT_DIM)
        } else {
            (
                icons::PAUSE,
                ACCENT,
                Color32::from_rgba_unmultiplied(30, 50, 35, 180),
            )
        };

        if !self.paused {
            let btn_pos = ui.next_widget_position() + egui::vec2(14.0, 14.0);
            let pulse = ((time * 2.0).sin() * 0.5 + 0.5) * 0.35 + 0.1;
            ui.painter().circle_stroke(
                btn_pos,
                18.0,
                Stroke::new(
                    1.5,
                    Color32::from_rgba_unmultiplied(
                        ACCENT.r(),
                        ACCENT.g(),
                        ACCENT.b(),
                        (pulse * 150.0) as u8,
                    ),
                ),
            );
        }

        if ui
            .add(
                egui::Button::new(RichText::new(icon).size(14.0).color(icon_col))
                    .fill(fill_col)
                    .stroke(Stroke::new(1.0, icon_col.gamma_multiply(0.5)))
                    .min_size(egui::vec2(28.0, 28.0))
                    .corner_radius(5.0),
            )
            .on_hover_text(if self.paused { "Play (Space)" } else { "Pause (Space)" })
            .clicked()
        {
            self.paused = !self.paused;
        }
    }

    fn step_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .add(
                egui::Button::new(RichText::new(icons::STEP).size(12.0).color(TEXT_DIM))
                    .fill(Color32::TRANSPARENT)
                    .stroke(Stroke::new(0.5, BORDER))
                    .min_size(egui::vec2(24.0, 24.0))
                    .corner_radius(3.0),
            )
            .on_hover_text(format!(
                "Step — advance {} physics step(s) then pause",
                self.steps_per_frame
            ))
            .clicked()
        {
            self.paused = false;
            self.step_pending = true;
        }
    }

    fn reset_btn(&mut self, ui: &mut egui::Ui) {
        if ui
            .add(
                egui::Button::new(RichText::new(icons::RESET).size(12.0).color(TEXT_DIM))
                    .fill(Color32::TRANSPARENT)
                    .stroke(Stroke::new(0.5, BORDER))
                    .min_size(egui::vec2(24.0, 24.0))
                    .corner_radius(3.0),
            )
            .on_hover_text("Reset energy/angular-momentum drift peaks")
            .clicked()
        {
            self.reset_drift_peaks();
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn vsep(ui: &mut egui::Ui) {
    ui.add(egui::Separator::default().vertical().spacing(2.0));
}

fn short_label(k: IntegratorKind) -> &'static str {
    match k {
        IntegratorKind::VelocityVerlet => "Velocity Verlet",
        IntegratorKind::Yoshida4 => "Yoshida 4",
        IntegratorKind::WisdomHolman => "Wisdom–Holman",
    }
}

/// Format a float as `±M.MMMxxx×10^N` with sign prefix.
fn fmt_sci_signed(v: f64, sig: usize) -> String {
    if v == 0.0 || v.abs() < f64::MIN_POSITIVE {
        return "+0".into();
    }
    let sign = if v >= 0.0 { '+' } else { '−' };
    let a = v.abs();
    let exp = a.log10().floor() as i32;
    let mantissa = a / 10f64.powi(exp);
    let prec = sig.saturating_sub(1);
    format!("{sign}{:.prec$}×10{}", mantissa, superscript(exp), prec = prec)
}

fn superscript(n: i32) -> String {
    n.to_string()
        .chars()
        .map(|c| match c {
            '-' => '⁻',
            '0' => '⁰',
            '1' => '¹',
            '2' => '²',
            '3' => '³',
            '4' => '⁴',
            '5' => '⁵',
            '6' => '⁶',
            '7' => '⁷',
            '8' => '⁸',
            '9' => '⁹',
            other => other,
        })
        .collect()
}

/// Convert sim-time to years (`t / 2π` convention used in natural-unit sims).
fn fmt_years(t: f64) -> String {
    let yr = t / (2.0 * PI);
    if yr.abs() < 0.01 {
        format!("{:.2e} yr", yr)
    } else if yr.abs() < 10_000.0 {
        format!("{:.2} yr", yr)
    } else {
        format!("{:.2e} yr", yr)
    }
}

fn energy_color(de: f64) -> Color32 {
    let a = de.abs();
    if a < 1e-6 {
        SUCCESS
    } else if a < 1e-3 {
        Color32::from_rgb(120, 200, 140)
    } else if a < 0.1 {
        WARN
    } else {
        DANGER
    }
}
