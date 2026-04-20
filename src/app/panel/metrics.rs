use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Color32, RichText, Stroke};
use std::f64::consts::PI;

// ── Drift severity ────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::app::panel) enum DriftSeverity {
    Excellent,
    Good,
    Acceptable,
    Warning,
    Critical,
}

impl DriftSeverity {
    pub(in crate::app::panel) fn from_peak(peak: f64) -> Self {
        let p = peak.abs();
        if p < 1e-9 { Self::Excellent }
        else if p < 1e-6 { Self::Good }
        else if p < 1e-3 { Self::Acceptable }
        else if p < 1e-1 { Self::Warning }
        else { Self::Critical }
    }

    pub(in crate::app::panel) fn color(self) -> Color32 {
        match self {
            Self::Excellent | Self::Good => SUCCESS,
            Self::Acceptable => TEXT_DIM,
            Self::Warning => ACCENT,
            Self::Critical => DANGER,
        }
    }

    pub(in crate::app::panel) fn dot(self) -> &'static str {
        match self {
            Self::Excellent | Self::Good | Self::Acceptable => "●",
            Self::Warning | Self::Critical => "▲",
        }
    }

    pub(in crate::app::panel) fn label(self) -> &'static str {
        match self {
            Self::Excellent => "excellent",
            Self::Good => "good",
            Self::Acceptable => "acceptable",
            Self::Warning => "warning",
            Self::Critical => "critical",
        }
    }

    fn hint(self) -> Option<&'static str> {
        match self {
            Self::Warning => Some("reduce dt or switch to Yoshida 4th-order"),
            Self::Critical => Some("simulation diverging — restart with smaller dt"),
            _ => None,
        }
    }
}

// ── Shared helper ─────────────────────────────────────────────────────────────

fn kv(ui: &mut egui::Ui, label: &str, value: &str, col: Color32) {
    ui.horizontal(|ui| {
        ui.add(egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)).truncate());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add(
                egui::Label::new(RichText::new(value).monospace().size(10.5).color(col))
                    .truncate(),
            );
        });
    });
}

fn sci(v: f64) -> String { format!("{:+.2e}", v) }

fn fmt_e(v: f64) -> String {
    let a = v.abs();
    if a == 0.0 { return "0.0000".into(); }
    if a < 1e-4 || a >= 1e5 { format!("{:+.4e}", v) } else { format!("{:+.4}", v) }
}

fn instant_color(v: f64) -> Color32 {
    let a = v.abs();
    if a < 1e-8 { SUCCESS }
    else if a < 1e-5 { TEXT_DIM }
    else if a < 1e-3 { ACCENT }
    else { DANGER }
}

// ── Compact status strip (shown at top of every panel tab) ────────────────────

impl SimulationApp {
    /// Minimal 3-row summary: stability + time + step count.
    /// Full diagnostics live in the Config tab (see `panel_diagnostics_detail`).
    pub(super) fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();
        let e_sev = DriftSeverity::from_peak(self.energy_drift_peak);
        let lz_sev = DriftSeverity::from_peak(self.lz_drift_peak);
        let worst = e_sev.max(lz_sev);
        let col = worst.color();

        // Row 1: stability indicator
        ui.horizontal(|ui| {
            ui.label(RichText::new(worst.dot()).size(9.0).color(col));
            ui.label(RichText::new("stability").size(10.0).color(TEXT_DIM));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(RichText::new(worst.label()).size(10.0).color(col).strong());
            });
        });

        // Row 2: sim time + years
        let t = m.t;
        let yr = t / (2.0 * PI);
        let yr_str = if yr.abs() < 0.01 { format!("{:.2e} yr", yr) }
                     else { format!("{:.2} yr", yr) };
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(format!("t  {:.4e}", t))
                    .monospace().size(10.0).color(TEXT_PRI),
            );
            ui.label(RichText::new("·").size(9.0).color(TEXT_DIM));
            ui.label(
                RichText::new(yr_str).monospace().size(10.0).color(TEXT_DIM),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(format!("{} steps", m.steps))
                        .monospace().size(9.5).color(TEXT_DIM),
                );
            });
        });

        // Row 3: peak drift (only after warm-up, only if non-trivial)
        if m.steps > 10 && self.energy_drift_peak > 0.0 {
            ui.horizontal(|ui| {
                ui.label(RichText::new("peak ΔE/E₀").size(9.0).color(TEXT_DIM));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        RichText::new(format!("{:.2e}", self.energy_drift_peak))
                            .monospace().size(9.0).color(e_sev.color()),
                    );
                });
            });
        } else if m.steps <= 10 {
            ui.label(RichText::new("warming up…").size(9.0).color(TEXT_DIM).italics());
        }

        // Actionable hint for Warning/Critical
        if let Some(hint) = worst.hint() {
            ui.add_space(2.0);
            egui::Frame::NONE
                .fill(Color32::from_rgba_unmultiplied(DANGER.r(), DANGER.g(), DANGER.b(), 18))
                .stroke(Stroke::new(0.5, col.gamma_multiply(0.4)))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(6, 3))
                .show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.add(egui::Label::new(RichText::new(hint).size(9.0).color(col)).wrap());
                });
        }

        let _ = TEXT_SEC;
        let _ = BORDER;
    }

    /// Full diagnostics table — called from the Config tab.
    pub(in crate::app::panel) fn panel_diagnostics_detail(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();
        let e_sev = DriftSeverity::from_peak(self.energy_drift_peak);
        let lz_sev = DriftSeverity::from_peak(self.lz_drift_peak);
        let worst = e_sev.max(lz_sev);

        // Conservation
        kv(ui, "E", &fmt_e(m.total_energy), TEXT_PRI);
        kv(ui, "ΔE / E₀", &sci(m.rel_energy_error), instant_color(m.rel_energy_error));
        let lz_triv = m.angular_momentum_z.abs() < 1e-10;
        let lz_err_str = if lz_triv { "—".to_owned() } else { sci(m.rel_angular_momentum_error) };
        let lz_err_col = if lz_triv { TEXT_DIM } else { instant_color(m.rel_angular_momentum_error) };
        kv(ui, "Lz", &fmt_e(m.angular_momentum_z), TEXT_PRI);
        kv(ui, "ΔLz / Lz₀", &lz_err_str, lz_err_col);

        ui.add_space(3.0);
        ui.add(egui::Separator::default().spacing(2.0));
        ui.add_space(2.0);

        // Energy breakdown
        kv(ui, "K (kinetic)", &fmt_e(m.kinetic), TEXT_DIM);
        kv(ui, "U (potential)", &fmt_e(m.potential), TEXT_DIM);

        ui.add_space(3.0);
        ui.add(egui::Separator::default().spacing(2.0));
        ui.add_space(2.0);

        // Solver diagnostics
        kv(ui, "dt", &format!("{:.3e}", m.dt), TEXT_DIM);
        kv(ui, "θ (opening)", &format!("{:.3}", m.theta), TEXT_DIM);
        kv(ui, "steps", &format!("{}", m.steps), TEXT_DIM);
        kv(ui, "vmax", &format!("{:.3e}", m.max_vel), TEXT_DIM);
        kv(ui, "amax", &format!("{:.3e}", m.max_acc), TEXT_DIM);

        if let Some(rec) = m.recommended_dt {
            let ratio = m.dt / rec;
            let col = if ratio <= 2.0 { TEXT_DIM } else if ratio <= 10.0 { ACCENT } else { DANGER };
            kv(ui, "suggested dt", &format!("{:.3e}", rec), col);
        }

        ui.add_space(3.0);
        ui.add(egui::Separator::default().spacing(2.0));
        ui.add_space(2.0);

        // Stability summary
        ui.horizontal(|ui| {
            ui.label(RichText::new("stability").size(10.0).color(TEXT_DIM));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(worst.label()).size(10.0).color(worst.color()).strong(),
                );
                ui.label(RichText::new(worst.dot()).size(9.0).color(worst.color()));
            });
        });
        if self.energy_drift_peak > 0.0 {
            kv(
                ui,
                "peak ΔE/E₀",
                &format!("{:.2e}", self.energy_drift_peak),
                e_sev.color(),
            );
        }
        if self.lz_drift_peak > 0.0 {
            kv(
                ui,
                "peak ΔLz/Lz₀",
                &format!("{:.2e}", self.lz_drift_peak),
                lz_sev.color(),
            );
        }
    }
}
