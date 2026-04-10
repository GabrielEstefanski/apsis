use crate::app::theme::{ACCENT, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::theme::{fix4, sci};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();

        let drift_color = |v: f64| {
            if v.abs() < 1e-8 {
                SUCCESS
            } else if v.abs() < 1e-5 {
                ACCENT
            } else {
                DANGER
            }
        };

        let de_color = drift_color(m.rel_energy_error);
        let dlz_color = drift_color(m.rel_angular_momentum_error);

        egui::Grid::new("metrics_compact")
            .num_columns(4)
            .spacing([4.0, 2.0])
            .show(ui, |ui| {
                // ── Energy ───────────────────────────── //
                ui.label(RichText::new("E").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.total_energy))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_PRI),
                );
                ui.label(RichText::new("dE/E₀").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(sci(m.rel_energy_error))
                        .monospace()
                        .size(10.0)
                        .color(de_color),
                );
                ui.end_row();

                // ── Angular momentum ─────────────────── //
                ui.label(RichText::new("Lz").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.angular_momentum_z))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_PRI),
                );
                ui.label(RichText::new("dLz").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(sci(m.rel_angular_momentum_error))
                        .monospace()
                        .size(10.0)
                        .color(dlz_color),
                );
                ui.end_row();

                // ── Energetics breakdown ─────────────── //
                ui.label(RichText::new("K").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.kinetic))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.label(RichText::new("U").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.potential))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.end_row();

                // ── Simulated time ─────────────────── //
                ui.label(RichText::new("t").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("{:.4e}", m.t))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_PRI),
                );
                ui.label(RichText::new("steps").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("{}", m.steps))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.end_row();

                // ── Integration params ─────────────── //
                ui.label(RichText::new("dt").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("{:.2e}", m.dt))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_SEC),
                );
                ui.label(RichText::new("θ").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("{:.3}", m.theta))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_SEC),
                );
                ui.end_row();

                // ── Integrator ──────────────────────── //
                ui.label(RichText::new("integr.").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("O({})", m.integrator.order()))
                        .monospace()
                        .size(10.0)
                        .color(ACCENT),
                );
                ui.label(RichText::new("").size(10.0)); // spacer
                ui.label(
                    RichText::new(m.integrator.label())
                        .size(9.0)
                        .color(TEXT_DIM),
                );
                ui.end_row();

                // ── Stability diagnostics ───────────── //
                ui.label(RichText::new("vmax").size(10.0).color(TEXT_DIM));
                ui.label(
                    RichText::new(sci(m.max_vel))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.label(RichText::new("amax").size(10.0).color(TEXT_DIM));
                ui.label(
                    RichText::new(sci(m.max_acc))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.end_row();
            });

        // ── Drift alert ───────────────────────────── //
        // Lz relative error is meaningless when initial angular momentum ≈ 0
        // (symmetric configs like figure-8). Detect by huge relative error.
        let lz_ref_trivial = m.rel_angular_momentum_error.abs() > 1e3
            || m.angular_momentum_z.abs() < 1e-10;

        let energy_bad = m.rel_energy_error.abs() >= 1e-5
            && m.total_energy.abs() > 1e-15;
        let lz_bad = !lz_ref_trivial && m.rel_angular_momentum_error.abs() >= 1e-5;

        if energy_bad || lz_bad {
            ui.add_space(2.0);

            let mut parts = Vec::new();
            if energy_bad {
                parts.push(format!("dE {}", sci(m.rel_energy_error)));
            }
            if lz_bad {
                parts.push(format!("dLz {}", sci(m.rel_angular_momentum_error)));
            }

            ui.label(
                RichText::new(format!("⚠ drift: {}", parts.join("  ")))
                    .size(9.5)
                    .color(DANGER),
            );
        }
    }
}
