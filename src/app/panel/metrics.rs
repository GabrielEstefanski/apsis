use crate::app::theme::{ACCENT, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, RichText};

/// Column widths for the 2+2 grid: [label, value, label, value]
const LW: f32 = 34.0; // label column
const VW: f32 = 74.0; // value column

fn lbl(ui: &mut egui::Ui, text: &str) {
    ui.add_sized(
        egui::vec2(LW, 0.0),
        egui::Label::new(RichText::new(text).size(10.5).color(TEXT_SEC)),
    );
}

fn val(ui: &mut egui::Ui, text: &str, color: egui::Color32, size: f32) {
    ui.add_sized(
        egui::vec2(VW, 0.0),
        egui::Label::new(RichText::new(text).monospace().size(size).color(color)).truncate(),
    );
}

impl SimulationApp {
    pub(super) fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();
        let tl = &self.physics_cfg.time_label;
        let dl = &self.physics_cfg.dist_label;

        let drift_color = |v: f64| {
            if v.abs() < 1e-8 { SUCCESS }
            else if v.abs() < 1e-5 { ACCENT }
            else { DANGER }
        };

        let de_col  = drift_color(m.rel_energy_error);
        let dlz_col = drift_color(m.rel_angular_momentum_error);

        // ── Conserved quantities ──────────────────────────────────────────
        egui::Grid::new("metrics_main")
            .num_columns(4)
            .spacing([4.0, 3.0])
            .show(ui, |ui| {

            // Row: E | dE/E₀
            lbl(ui, "E");
            val(ui, &fmt_e(m.total_energy), TEXT_PRI, 11.0);
            lbl(ui, "dE/E₀");
            val(ui, &sci(m.rel_energy_error), de_col, 11.0);
            ui.end_row();

            // Row: Lz | dLz/Lz₀
            lbl(ui, "Lz");
            val(ui, &fmt_e(m.angular_momentum_z), TEXT_PRI, 11.0);
            lbl(ui, "dLz");
            val(ui, &sci(m.rel_angular_momentum_error), dlz_col, 11.0);
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
        });

        // ── Integrator — full-width row with badge ────────────────────────
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

        // ── Stability diagnostics (smaller, secondary) ────────────────────
        ui.add_space(1.0);
        egui::Grid::new("metrics_diag")
            .num_columns(4)
            .spacing([4.0, 2.0])
            .show(ui, |ui| {
            lbl(ui, "vmax");
            val(ui, &format!("{:.3e} {dl}/s", m.max_vel), TEXT_DIM, 9.5);
            lbl(ui, "amax");
            val(ui, &format!("{:.3e}", m.max_acc), TEXT_DIM, 9.5);
            ui.end_row();
        });

        // ── Drift alert ───────────────────────────────────────────────────
        let lz_ref_trivial = m.rel_angular_momentum_error.abs() > 1e3
            || m.angular_momentum_z.abs() < 1e-10;
        let energy_bad = m.rel_energy_error.abs() >= 1e-5 && m.total_energy.abs() > 1e-15;
        let lz_bad = !lz_ref_trivial && m.rel_angular_momentum_error.abs() >= 1e-5;

        if energy_bad || lz_bad {
            ui.add_space(2.0);
            let mut parts = Vec::new();
            if energy_bad { parts.push(format!("dE {}", sci(m.rel_energy_error))); }
            if lz_bad     { parts.push(format!("dLz {}", sci(m.rel_angular_momentum_error))); }
            ui.add(
                egui::Label::new(
                    RichText::new(format!("⚠ drift: {}", parts.join("  ")))
                        .size(10.0)
                        .color(DANGER),
                )
                .truncate(),
            );
        }
    }
}

fn sci(v: f64) -> String { format!("{:+.3e}", v) }

/// Format energy/momentum: short fixed-point when small, sci for large/tiny.
fn fmt_e(v: f64) -> String {
    let a = v.abs();
    if a == 0.0           { return "0.0000".into(); }
    if a < 1e-4 || a >= 1e5 { format!("{:+.4e}", v) }
    else                  { format!("{:+.4}", v) }
}
