//! Bottom playbar — transport controls + key simulation readouts.
//!
//! Layout (left → right):
//!   [▶ ⏭ ⟲]  |  T · yr  |  DT  |  SPEED [slider] ×N  |  ···  |  ● stability

use crate::app::icons;
use crate::app::panel::metrics::DriftSeverity;
use crate::app::theme::{
    ACCENT, ACCENT_DIM, BORDER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC,
};
use crate::app::ui::SimulationApp;
use apsis::physics::integrator::IntegratorKind;
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
                    ui.label(RichText::new(fmt_sci(t, 4)).size(10.0).monospace().color(TEXT_SEC));
                    ui.label(RichText::new("·").size(9.0).color(TEXT_DIM));
                    ui.label(RichText::new(fmt_years(yr)).size(10.0).monospace().color(TEXT_DIM));

                    vsep(ui);

                    // ── DT / ACC ─────────────────────────────────────────────
                    // IAS15 controls accuracy via epsilon (ε), not a fixed DT.
                    // For all other integrators, show the DT drag-value as usual.
                    if self.physics_cfg.integrator == IntegratorKind::Ias15 {
                        ui.label(RichText::new("ACC").size(8.5).color(TEXT_DIM).strong());
                        let presets: &[(&str, f64, &str)] = &[
                            ("Fast", 1e-6, "ε = 1×10⁻⁶  ·  faster, less accurate"),
                            ("Normal", 1e-9, "ε = 1×10⁻⁹  ·  REBOUND default — recommended"),
                            ("Fine", 1e-12, "ε = 1×10⁻¹²  ·  high precision, slower"),
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
                            .selected_text(RichText::new(current_label).size(10.0).color(TEXT_PRI))
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
                            .on_hover_text(
                                "Integration timestep — smaller = more accurate but slower",
                            )
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
                    let ratio = if speed_yr > 0.0 { actual_yr / speed_yr } else { 1.0 };
                    let shortfall =
                        sim_rate > 0.0 && shortfall_with_hysteresis(self.shortfall_active, ratio);
                    self.shortfall_active = shortfall;
                    let speed_col = if shortfall { TEXT_DIM } else { ACCENT };
                    // Shortfall surfaces the gap as a single percentage of
                    // delivered work, not as a second rate value: target
                    // and achieved cross unit boundaries (300 yr/s vs.
                    // 0.7 d/s) and reading them as the same thing requires
                    // mental conversion. Percentage is unit-free and
                    // tells the user directly how much of what they asked
                    // for the solver is keeping up with.
                    let speed_text = if shortfall {
                        format!("{} · {}", fmt_speed(speed_yr), fmt_percent(ratio))
                    } else {
                        fmt_speed(speed_yr)
                    };
                    let speed_tooltip = if shortfall {
                        format!(
                            "Target {} yr/s — solver delivering {}.\n\
                             Physics can't keep up; render slows to match.\n\
                             Lower the slider, switch to a faster integrator, or\n\
                             reduce body count to close the gap.",
                            fmt_speed(speed_yr),
                            fmt_rate(actual_yr),
                        )
                    } else {
                        "Target simulation speed (yr/s).\n\
                         The physics thread advances this many simulated years per real second."
                            .to_string()
                    };
                    ui.label(RichText::new(speed_text).monospace().size(10.0).color(speed_col))
                        .on_hover_text(speed_tooltip);

                    vsep(ui);

                    // ── Integrator ────────────────────────────────────────────
                    ui.label(RichText::new("ALGO").size(8.5).color(TEXT_DIM).strong());
                    let current = self.physics_cfg.integrator;
                    let current_short = integrator_short_label(current);
                    egui::ComboBox::from_id_salt("playbar_integrator")
                        .selected_text(RichText::new(current_short).size(10.0).color(TEXT_PRI))
                        .width(90.0)
                        .show_ui(ui, |ui| {
                            for variant in IntegratorKind::ALL {
                                // `selectable_label` renders selection state
                                // without binding `&mut physics_cfg.integrator` —
                                // the actual mutation goes through
                                // `request_integrator_change`, which is the
                                // single code path that can surface the
                                // precision-confirmation modal.
                                let selected = self.physics_cfg.integrator == variant;
                                let r = ui
                                    .selectable_label(
                                        selected,
                                        RichText::new(variant.label()).size(10.0),
                                    )
                                    .on_hover_text(format!(
                                        "O({}) · {}F/step\n{}",
                                        variant.order(),
                                        variant.force_evals_per_step(),
                                        variant.description(),
                                    ));
                                if r.clicked() {
                                    self.request_integrator_change(variant);
                                }
                            }
                        });

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
            (icons::PAUSE, ACCENT, Color32::from_rgba_unmultiplied(30, 50, 35, 180))
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
        IntegratorKind::WHFast => "WHFast",
        IntegratorKind::Ias15 => "IAS15",
        IntegratorKind::Mercurius => "Mercurius",
        IntegratorKind::ImplicitMidpoint => "IM2",
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

/// Format `ratio ∈ [0, 1]` as a percentage with adaptive precision.
///
/// The shortfall display can land anywhere from 79 % (just inside the
/// hysteresis band) to 0.0006 % (the WH-class collapse the user
/// reported in #63 review). Fixed `{:.0}%` reads "0%" for everything
/// below half a percent, hiding two orders of magnitude of severity.
/// Adaptive precision keeps the readout informative across the band.
fn fmt_percent(ratio: f64) -> String {
    let pct = ratio * 100.0;
    if pct >= 10.0 {
        format!("{:.0}%", pct)
    } else if pct >= 1.0 {
        format!("{:.1}%", pct)
    } else if pct >= 0.01 {
        format!("{:.2}%", pct)
    } else {
        // Below 0.01 % the integrator is essentially stalled; just
        // show that it's negligible without trailing zero noise.
        "<0.01%".to_string()
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

/// Hysteretic threshold for the playbar's "physics behind target" cue.
///
/// `ratio` is `achieved / target`. The cue activates at 80 % achieved
/// and only deactivates once achieved climbs back above 90 % — the
/// 10 pp gap kills the binary flicker that a single threshold
/// produces when the achieved rate hovers around the boundary.
///
/// Stateless (pass the previous value in, get the next one out) so
/// the policy is unit-testable without an `egui::Ui`.
pub(super) fn shortfall_with_hysteresis(currently_active: bool, ratio: f64) -> bool {
    const ENTER: f64 = 0.80;
    const EXIT: f64 = 0.90;
    if currently_active { ratio < EXIT } else { ratio < ENTER }
}

#[cfg(test)]
mod tests {
    use super::{fmt_percent, shortfall_with_hysteresis};

    // ── Percentage formatting ───────────────────────────────────────────────

    #[test]
    fn fmt_percent_uses_integer_above_ten() {
        assert_eq!(fmt_percent(0.50), "50%");
        assert_eq!(fmt_percent(0.123), "12%");
    }

    #[test]
    fn fmt_percent_keeps_one_decimal_in_single_digits() {
        assert_eq!(fmt_percent(0.05), "5.0%");
        assert_eq!(fmt_percent(0.087), "8.7%");
    }

    #[test]
    fn fmt_percent_keeps_two_decimals_below_one_percent() {
        // The WH-collapse case from #63 review (~0.06 %): a fixed
        // {:.0}% would read "0%" and erase the severity. Two
        // decimals show the actual order of magnitude.
        assert_eq!(fmt_percent(0.0006), "0.06%");
        assert_eq!(fmt_percent(0.001), "0.10%");
    }

    #[test]
    fn fmt_percent_collapses_negligible_to_threshold_marker() {
        // Below 0.01 % every digit is noise — physics is essentially
        // stalled. Show the qualitative state, not a tail of zeros.
        assert_eq!(fmt_percent(1e-5), "<0.01%");
        assert_eq!(fmt_percent(0.0), "<0.01%");
    }

    #[test]
    fn shortfall_activates_below_enter_threshold() {
        // Cold start: cue is off and achieved drops below 80 %.
        assert!(shortfall_with_hysteresis(false, 0.79));
        assert!(shortfall_with_hysteresis(false, 0.50));
        assert!(shortfall_with_hysteresis(false, 0.0));
    }

    #[test]
    fn shortfall_stays_off_in_hysteresis_band_when_starting_off() {
        // Off and ratio in [0.80, 0.90): doesn't trigger yet — needs
        // to drop below ENTER first.
        assert!(!shortfall_with_hysteresis(false, 0.80));
        assert!(!shortfall_with_hysteresis(false, 0.85));
        assert!(!shortfall_with_hysteresis(false, 0.89));
    }

    #[test]
    fn shortfall_stays_on_in_hysteresis_band_when_starting_on() {
        // On and ratio in [0.80, 0.90): stays on — needs to climb
        // above EXIT to clear.
        assert!(shortfall_with_hysteresis(true, 0.80));
        assert!(shortfall_with_hysteresis(true, 0.85));
        assert!(shortfall_with_hysteresis(true, 0.89));
    }

    #[test]
    fn shortfall_deactivates_above_exit_threshold() {
        // On and achieved climbs above 90 %: cue clears.
        assert!(!shortfall_with_hysteresis(true, 0.90));
        assert!(!shortfall_with_hysteresis(true, 0.95));
        assert!(!shortfall_with_hysteresis(true, 1.05));
    }

    #[test]
    fn oscillation_around_enter_threshold_latches_active() {
        // The single-threshold flicker scenario: achieved hovers
        // around 80 % of target. With one threshold, every tick that
        // crosses 0.80 flips the state. With hysteresis, the first
        // dip below 0.80 latches active, and any subsequent values
        // inside the [0.80, 0.90) band keep the state.
        let mut active = false;
        for ratio in [0.78, 0.82, 0.79, 0.83, 0.78, 0.85, 0.88] {
            active = shortfall_with_hysteresis(active, ratio);
        }
        assert!(active, "should latch active after dipping below ENTER and never reaching EXIT");
    }

    #[test]
    fn oscillation_inside_band_does_not_flip_state() {
        // Pure hysteresis-band test: ratio stays inside [0.80, 0.90)
        // the whole time. State must be preserved (whatever it was
        // when entering the band) across the entire sequence.
        let band = [0.80, 0.85, 0.82, 0.89, 0.81, 0.87];

        let mut active = true;
        for ratio in band {
            active = shortfall_with_hysteresis(active, ratio);
        }
        assert!(active, "starting active, band oscillation must not deactivate");

        let mut active = false;
        for ratio in band {
            active = shortfall_with_hysteresis(active, ratio);
        }
        assert!(!active, "starting inactive, band oscillation must not activate");
    }

    #[test]
    fn full_recovery_above_exit_clears_state() {
        // Sanity: hysteresis is not a permanent latch. Once achieved
        // climbs above EXIT, the state clears; a subsequent dip into
        // the band stays clear (consistent with the band-test above).
        let mut active = true;
        for ratio in [0.85, 0.92, 0.85] {
            active = shortfall_with_hysteresis(active, ratio);
        }
        assert!(
            !active,
            "single excursion above EXIT must clear state, even if next sample re-enters band"
        );
    }
}
