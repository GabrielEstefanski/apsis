//! Bottom playbar — transport controls + key simulation readouts.
//!
//! Layout (left → right):
//!   [▶ ⏭ ⟲]  |  T · yr  |  DT  |  SPEED [slider] ×N  |  ···  |  ● stability

use crate::app::icons;
use crate::app::panel::metrics::DriftSeverity;
use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::physics::integrator::IntegratorKind;
use eframe::egui::{self, Color32, RichText, Stroke};
use std::f64::consts::PI;

pub const PLAYBAR_HEIGHT: f32 = 36.0;

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

                    // ── Time readout ─────────────────────────────────────────
                    let t = self.system.t();
                    let yr = t / (2.0 * PI);
                    ui.label(RichText::new("T").size(8.5).color(TEXT_DIM).strong());
                    ui.label(
                        RichText::new(fmt_sci(t, 4))
                            .size(10.0)
                            .monospace()
                            .color(TEXT_SEC),
                    );
                    ui.label(RichText::new("·").size(9.0).color(TEXT_DIM));
                    ui.label(
                        RichText::new(fmt_years(yr))
                            .size(10.0)
                            .monospace()
                            .color(TEXT_DIM),
                    );

                    vsep(ui);

                    // ── DT / ACC ─────────────────────────────────────────────
                    // IAS15 controls accuracy via epsilon (ε), not a fixed DT.
                    // For all other integrators, show the DT drag-value as usual.
                    if self.physics_cfg.integrator == IntegratorKind::Ias15 {
                        ui.label(RichText::new("ACC").size(8.5).color(TEXT_DIM).strong());
                        let presets: &[(&str, f64, &str)] = &[
                            ("Fast",   1e-6,  "ε = 1×10⁻⁶  ·  faster, less accurate"),
                            ("Normal", 1e-9,  "ε = 1×10⁻⁹  ·  REBOUND default — recommended"),
                            ("Fine",   1e-12, "ε = 1×10⁻¹²  ·  high precision, slower"),
                        ];
                        let current_label = presets
                            .iter()
                            .min_by(|a, b| {
                                let da = (a.1.log10() - self.ias15_epsilon.log10()).abs();
                                let db = (b.1.log10() - self.ias15_epsilon.log10()).abs();
                                da.partial_cmp(&db).unwrap()
                            })
                            .map(|p| p.0)
                            .unwrap_or("Normal");
                        egui::ComboBox::from_id_salt("playbar_ias15_epsilon")
                            .selected_text(
                                RichText::new(current_label).size(10.0).color(TEXT_PRI),
                            )
                            .width(70.0)
                            .show_ui(ui, |ui| {
                                for (label, eps, hint) in presets {
                                    ui.selectable_value(
                                        &mut self.ias15_epsilon,
                                        *eps,
                                        RichText::new(*label).size(10.0),
                                    )
                                    .on_hover_text(*hint);
                                }
                            });
                    } else {
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
                            .on_hover_text("Integration timestep — smaller = more accurate but slower")
                            .changed()
                        {
                            self.system.set_dt(dt);
                        }
                    }

                    vsep(ui);

                    // ── Speed (sim-rate target) ───────────────────────────────
                    ui.label(RichText::new("SPEED").size(8.5).color(TEXT_DIM).strong());
                    // Slider operates in yr/s; convert to/from internal units (2π = 1 yr).
                    let tau = std::f64::consts::TAU;
                    let mut speed_yr = self.sim_rate_target / tau;
                    if ui
                        .add_sized(
                            [80.0, 14.0],
                            egui::Slider::new(&mut speed_yr, 0.01_f64..=100_000.0_f64)
                                .logarithmic(true)
                                .show_value(false),
                        )
                        .changed()
                    {
                        self.sim_rate_target = speed_yr.max(0.01) * tau;
                    }
                    let sim_rate = self.system.sim_rate();
                    let actual_yr = sim_rate / tau;
                    let shortfall = sim_rate > 0.0 && actual_yr < speed_yr * 0.8;
                    let speed_col = if shortfall { TEXT_DIM } else { ACCENT };
                    ui.label(
                        RichText::new(fmt_speed(speed_yr))
                            .monospace()
                            .size(10.0)
                            .color(speed_col),
                    )
                    .on_hover_text(
                        "Target simulation speed (yr/s).\n\
                         The physics thread advances this many simulated years per real second.\n\
                         If the sim can't keep up, the label dims and actual speed is shown below.",
                    );

                    vsep(ui);

                    // ── Integrator ────────────────────────────────────────────
                    ui.label(RichText::new("ALGO").size(8.5).color(TEXT_DIM).strong());
                    let current = self.physics_cfg.integrator;
                    let current_short = integrator_short_label(current);
                    egui::ComboBox::from_id_salt("playbar_integrator")
                        .selected_text(
                            RichText::new(current_short).size(10.0).color(TEXT_PRI),
                        )
                        .width(90.0)
                        .show_ui(ui, |ui| {
                            for variant in IntegratorKind::ALL {
                                let r = ui.selectable_value(
                                    &mut self.physics_cfg.integrator,
                                    variant,
                                    RichText::new(variant.label()).size(10.0),
                                )
                                .on_hover_text(format!(
                                    "O({}) · {}F/step\n{}",
                                    variant.order(),
                                    variant.force_evals_per_step(),
                                    variant.description(),
                                ));
                                if r.clicked() {
                                    self.system.set_integrator(variant);
                                }
                            }
                        });

                    vsep(ui);

                    // ── Sim throughput ────────────────────────────────────────
                    let sim_rate = self.system.sim_rate();
                    if sim_rate > 0.0 {
                        let yr_per_s = sim_rate / (2.0 * PI);
                        ui.label(
                            RichText::new(fmt_rate(yr_per_s))
                                .monospace()
                                .size(9.5)
                                .color(TEXT_DIM),
                        )
                        .on_hover_text(format!(
                            "Simulation throughput\n{:.3e} sim·units/s · {:.3e} yr/s",
                            sim_rate, yr_per_s
                        ));
                    }

                    // ── Right: stability badge ────────────────────────────────
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.spacing_mut().item_spacing.x = 6.0;

                        let m = self.system.metrics();
                        let e_sev = DriftSeverity::from_peak(self.energy_drift_peak);
                        let lz_sev = DriftSeverity::from_peak(self.lz_drift_peak);
                        let worst = e_sev.max(lz_sev);
                        let col = worst.color();

                        let (dot, label) = if m.steps <= 10 {
                            ("○", "warming up")
                        } else {
                            (worst.dot(), worst.label())
                        };

                        let hint = format!(
                            "Numerical stability\n\nPeak ΔE/E₀  = {:.2e}\nPeak ΔLz/Lz₀ = {:.2e}\n\n\
                             Excellent < 1×10⁻⁹  |  Good < 1×10⁻⁶  |  Acceptable < 1×10⁻³",
                            self.energy_drift_peak, self.lz_drift_peak
                        );

                        ui.add(
                            egui::Button::new(
                                RichText::new(format!("{dot}  {label}")).size(9.5).color(col),
                            )
                            .fill(Color32::TRANSPARENT)
                            .stroke(Stroke::new(0.5, col.gamma_multiply(0.35)))
                            .min_size(egui::vec2(96.0, 22.0))
                            .corner_radius(3.0),
                        )
                        .on_hover_text(hint);
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
            .on_hover_text(if self.paused { "Play  [Space]" } else { "Pause  [Space]" })
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
            .on_hover_text("Step — advance one physics batch then pause")
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
            .on_hover_text("Reset drift peak counters")
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

fn integrator_short_label(k: IntegratorKind) -> &'static str {
    match k {
        IntegratorKind::VelocityVerlet => "Verlet",
        IntegratorKind::Yoshida4 => "Yoshida-4",
        IntegratorKind::WisdomHolman => "W–H",
        IntegratorKind::Ias15 => "IAS15",
    }
}

fn fmt_speed(yr_per_s: f64) -> String {
    if yr_per_s < 1.0 / 365.25 {
        format!("{:.1}d/s", yr_per_s * 365.25)
    } else if yr_per_s < 1.0 {
        format!("{:.2}yr/s", yr_per_s)
    } else if yr_per_s < 1_000.0 {
        format!("{:.1}yr/s", yr_per_s)
    } else if yr_per_s < 1_000_000.0 {
        format!("{:.1}kyr/s", yr_per_s / 1_000.0)
    } else {
        format!("{:.1}Myr/s", yr_per_s / 1_000_000.0)
    }
}

fn fmt_rate(yr_per_s: f64) -> String {
    if yr_per_s < 1.0 / 365.25 {
        format!("{:.1} h/s", yr_per_s * 365.25 * 24.0)
    } else if yr_per_s < 1.0 {
        format!("{:.1} d/s", yr_per_s * 365.25)
    } else if yr_per_s < 1_000.0 {
        format!("{:.1} yr/s", yr_per_s)
    } else if yr_per_s < 1_000_000.0 {
        format!("{:.1} kyr/s", yr_per_s / 1_000.0)
    } else {
        format!("{:.1} Myr/s", yr_per_s / 1_000_000.0)
    }
}

fn fmt_sci(v: f64, sig: usize) -> String {
    if v == 0.0 || v.abs() < f64::MIN_POSITIVE {
        return "+0".into();
    }
    let sign = if v >= 0.0 { '+' } else { '−' };
    let a = v.abs();
    let exp = a.log10().floor() as i32;
    let mantissa = a / 10f64.powi(exp);
    let prec = sig.saturating_sub(1);
    format!("{sign}{:.prec$}e{exp}", mantissa, prec = prec)
}

fn fmt_years(yr: f64) -> String {
    if yr.abs() < 0.01 {
        format!("{:.2e} yr", yr)
    } else if yr.abs() < 10_000.0 {
        format!("{:.2} yr", yr)
    } else {
        format!("{:.2e} yr", yr)
    }
}
