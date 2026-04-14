use crate::app::theme::{ACCENT, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Color32, RichText, Stroke};

// ── Layout constants ──────────────────────────────────────────────────────────

const LW: f32 = 34.0; // label column width
const VW: f32 = 74.0; // value column width

// ── Drift severity ────────────────────────────────────────────────────────────

/// Classification of numerical drift, based on the *peak* value seen so far.
/// Only ever moves toward worse severity mid-run; resets on scenario load.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DriftSeverity {
    Excellent,   // |peak| < 1e-9
    Good,        // |peak| < 1e-6
    Acceptable,  // |peak| < 1e-3
    Warning,     // |peak| < 1e-1
    Critical,    // |peak| >= 1e-1
}

impl DriftSeverity {
    fn from_peak(peak: f64) -> Self {
        let p = peak.abs();
        if p < 1e-9      { Self::Excellent }
        else if p < 1e-6 { Self::Good }
        else if p < 1e-3 { Self::Acceptable }
        else if p < 1e-1 { Self::Warning }
        else             { Self::Critical }
    }

    fn color(self) -> Color32 {
        match self {
            Self::Excellent  => SUCCESS,
            Self::Good       => SUCCESS,
            Self::Acceptable => TEXT_DIM,
            Self::Warning    => ACCENT,
            Self::Critical   => DANGER,
        }
    }

    fn dot(self) -> &'static str {
        match self {
            Self::Excellent  => "●",
            Self::Good       => "●",
            Self::Acceptable => "●",
            Self::Warning    => "▲",
            Self::Critical   => "▲",
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Excellent  => "excellent",
            Self::Good       => "good",
            Self::Acceptable => "acceptable",
            Self::Warning    => "warning",
            Self::Critical   => "critical",
        }
    }

    /// One-liner hint shown only at Warning/Critical.
    fn hint(self) -> Option<&'static str> {
        match self {
            Self::Warning  => Some("reduce dt or switch to Yoshida 4th-order"),
            Self::Critical => Some("simulation diverging — restart with smaller dt"),
            _              => None,
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn lbl(ui: &mut egui::Ui, text: &str) {
    ui.add_sized(
        egui::vec2(LW, 0.0),
        egui::Label::new(RichText::new(text).size(10.5).color(TEXT_SEC)),
    );
}

fn val(ui: &mut egui::Ui, text: &str, color: Color32, size: f32) {
    ui.add_sized(
        egui::vec2(VW, 0.0),
        egui::Label::new(RichText::new(text).monospace().size(size).color(color)).truncate(),
    );
}

fn sci(v: f64) -> String { format!("{:+.2e}", v) }

fn fmt_e(v: f64) -> String {
    let a = v.abs();
    if a == 0.0                { return "0.0000".into(); }
    if a < 1e-4 || a >= 1e5   { format!("{:+.4e}", v) }
    else                       { format!("{:+.4}", v) }
}

// ── Draw ──────────────────────────────────────────────────────────────────────

impl SimulationApp {
    pub(super) fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m   = self.system.metrics();
        let tl  = &self.physics_cfg.time_label;
        let dl  = &self.physics_cfg.dist_label;

        let e_sev  = DriftSeverity::from_peak(self.energy_drift_peak);
        let lz_sev = DriftSeverity::from_peak(self.lz_drift_peak);

        // Color for current instantaneous values (based on current, not peak)
        let instant_color = |v: f64| -> Color32 {
            let a = v.abs();
            if a < 1e-8      { SUCCESS }
            else if a < 1e-5 { TEXT_DIM }
            else if a < 1e-3 { ACCENT }
            else             { DANGER }
        };

        // ── Conserved quantities ──────────────────────────────────────────────
        egui::Grid::new("metrics_main")
            .num_columns(4)
            .spacing([4.0, 3.0])
            .show(ui, |ui| {
                // Row: E | dE/E₀ (instantaneous, colored by current value)
                lbl(ui, "E");
                val(ui, &fmt_e(m.total_energy), TEXT_PRI, 11.0);
                lbl(ui, "dE/E₀");
                val(ui, &sci(m.rel_energy_error), instant_color(m.rel_energy_error), 11.0);
                ui.end_row();

                // Row: Lz | dLz/Lz₀
                let lz_trivial = m.angular_momentum_z.abs() < 1e-10;
                let (lz_val_str, lz_err_str, lz_col) = if lz_trivial {
                    (fmt_e(m.angular_momentum_z), "—".into(), TEXT_DIM)
                } else {
                    (
                        fmt_e(m.angular_momentum_z),
                        sci(m.rel_angular_momentum_error),
                        instant_color(m.rel_angular_momentum_error),
                    )
                };
                lbl(ui, "Lz");
                val(ui, &lz_val_str, TEXT_PRI, 11.0);
                lbl(ui, "dLz/Lz₀");
                val(ui, &lz_err_str, lz_col, 11.0);
                ui.end_row();

                // Row: K | U
                lbl(ui, "K");
                val(ui, &fmt_e(m.kinetic), TEXT_DIM, 10.5);
                lbl(ui, "U");
                val(ui, &fmt_e(m.potential), TEXT_DIM, 10.5);
                ui.end_row();

                // Row: t | steps
                lbl(ui, &format!("t [{tl}]"));
                val(ui, &format!("{:.4e}", m.t), TEXT_PRI, 11.0);
                lbl(ui, "steps");
                val(ui, &format!("{}", m.steps), TEXT_DIM, 10.5);
                ui.end_row();

                // Row: dt | θ
                lbl(ui, "dt");
                val(ui, &format!("{:.2e}", m.dt), TEXT_SEC, 10.5);
                lbl(ui, "θ");
                val(ui, &format!("{:.3}", m.theta), TEXT_SEC, 10.5);
                ui.end_row();

                // Row: rec. dt (only when a physics-justified suggestion exists)
                if let Some(rec) = m.recommended_dt {
                    let ratio = m.dt / rec;
                    let color = if ratio <= 2.0 { TEXT_DIM }
                                else if ratio <= 10.0 { ACCENT }
                                else { DANGER };
                    lbl(ui, "rec. dt");
                    val(ui, &format!("{:.2e}", rec), color, 10.5);
                    // Intentionally leave the two right-side columns blank
                    // so the row still occupies exactly two columns.
                    lbl(ui, "");
                    val(ui, "", TEXT_DIM, 10.5);
                    ui.end_row();
                }
            });

        // ── Integrator ────────────────────────────────────────────────────────
        ui.add_space(2.0);
        ui.horizontal(|ui| {
            ui.label(RichText::new("integr.").size(10.5).color(TEXT_SEC));
            ui.add_space(2.0);
            ui.add(
                egui::Label::new(
                    RichText::new(m.integrator.label()).size(10.5).color(TEXT_DIM),
                )
                .truncate(),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(format!("O({})", m.integrator.order()))
                        .monospace()
                        .size(10.0)
                        .color(ACCENT),
                );
            });
        });

        // ── Secondary diagnostics ─────────────────────────────────────────────
        ui.add_space(1.0);
        egui::Grid::new("metrics_diag")
            .num_columns(4)
            .spacing([4.0, 2.0])
            .show(ui, |ui| {
                lbl(ui, "vmax");
                val(ui, &format!("{:.3e} {dl}/t", m.max_vel), TEXT_DIM, 9.5);
                lbl(ui, "amax");
                val(ui, &format!("{:.3e}", m.max_acc), TEXT_DIM, 9.5);
                ui.end_row();
            });

        // ── Stability status — always rendered, fixed layout ──────────────────
        ui.add_space(4.0);
        ui.add(egui::Separator::default().spacing(2.0));
        ui.add_space(3.0);

        // Worst severity between energy and angular momentum
        let worst = e_sev.max(lz_sev);
        let sev_col = worst.color();

        // Status line (always present — no layout shift)
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(worst.dot())
                    .size(9.0)
                    .color(sev_col),
            );
            ui.label(
                RichText::new("stability")
                    .size(9.5)
                    .color(TEXT_DIM),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(worst.label())
                        .size(9.5)
                        .color(sev_col)
                        .strong(),
                );
            });
        });

        // Peak row (always present — shows running maximum seen this run)
        if m.steps > 10 {
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new("peak dE/E₀")
                        .size(9.0)
                        .color(TEXT_DIM),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let peak_str = if self.energy_drift_peak == 0.0 {
                        "—".into()
                    } else {
                        format!("{:.2e}", self.energy_drift_peak)
                    };
                    ui.label(
                        RichText::new(peak_str)
                            .monospace()
                            .size(9.0)
                            .color(e_sev.color()),
                    );
                });
            });

            // Lz peak (only if system has non-trivial angular momentum)
            if self.lz_drift_peak > 0.0 {
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new("peak dLz/Lz₀")
                            .size(9.0)
                            .color(TEXT_DIM),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(
                            RichText::new(format!("{:.2e}", self.lz_drift_peak))
                                .monospace()
                                .size(9.0)
                                .color(lz_sev.color()),
                        );
                    });
                });
            }
        } else {
            // Placeholder height while warming up — prevents layout jump
            ui.label(
                RichText::new("warming up…")
                    .size(9.0)
                    .color(TEXT_DIM)
                    .italics(),
            );
        }

        // Actionable hint — only for Warning/Critical, fixed-height placeholder otherwise
        ui.add_space(2.0);
        if let Some(hint) = worst.hint() {
            egui::Frame::NONE
                .fill(Color32::from_rgba_unmultiplied(
                    DANGER.r(), DANGER.g(), DANGER.b(), 18,
                ))
                .stroke(Stroke::new(0.5, sev_col.gamma_multiply(0.4)))
                .corner_radius(3.0)
                .inner_margin(egui::Margin::symmetric(6, 3))
                .show(ui, |ui| {
                    ui.set_width(ui.available_width());
                    ui.add(
                        egui::Label::new(
                            RichText::new(hint).size(9.0).color(sev_col),
                        )
                        .wrap(),
                    );
                });
        } else {
            // Reserve the same vertical space so the panel below doesn't shift
            // when the hint appears or disappears.
            ui.add_space(18.0);
        }
    }
}
