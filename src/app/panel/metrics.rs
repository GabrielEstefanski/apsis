use crate::app::theme::{ACCENT, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::theme::{fix4, sci};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();
        let dust_clouds = self
            .system
            .bodies()
            .iter()
            .filter(|b| b.is_diffuse_cloud())
            .count();

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

                ui.label(RichText::new("Lz").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.angular_momentum_z))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_PRI),
                );
                ui.label(RichText::new("dLz/L₀").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(sci(m.rel_angular_momentum_error))
                        .monospace()
                        .size(10.0)
                        .color(dlz_color),
                );
                ui.end_row();

                ui.label(RichText::new("K").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.kinetic))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                ui.label(RichText::new("dt").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(format!("{:.2e}", m.dt))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_SEC),
                );
                ui.end_row();

                ui.label(RichText::new("U").size(10.0).color(TEXT_SEC));
                ui.label(
                    RichText::new(fix4(m.potential))
                        .monospace()
                        .size(10.0)
                        .color(TEXT_DIM),
                );
                let max_drift = m.max_rel_energy_error.max(m.max_rel_angular_momentum_error);
                let peak_col = drift_color(max_drift);
                ui.label(RichText::new("peak").size(10.0).color(TEXT_DIM));
                ui.label(
                    RichText::new(sci(max_drift))
                        .monospace()
                        .size(10.0)
                        .color(peak_col),
                );
                ui.end_row();
            });

        // Drift alert banner: shown when either conserved quantity drifts > 1e-5
        let energy_bad = m.rel_energy_error.abs() >= 1e-5;
        let lz_bad = m.rel_angular_momentum_error.abs() >= 1e-5;
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

        if m.fragments_spawned_this_step > 0
            || m.hit_and_runs_this_step > 0
            || dust_clouds > 0
            || m.total_dust_mass > 0.0
        {
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                if m.fragments_spawned_this_step > 0 {
                    ui.label(
                        RichText::new(format!("frags {}", m.fragments_spawned_this_step))
                            .size(9.5)
                            .color(DANGER),
                    );
                }
                if m.hit_and_runs_this_step > 0 {
                    ui.label(
                        RichText::new(format!("H&R {}", m.hit_and_runs_this_step))
                            .size(9.5)
                            .color(ACCENT),
                    );
                }
                if dust_clouds > 0 {
                    ui.label(
                        RichText::new(format!("clouds {}", dust_clouds))
                            .size(9.5)
                            .color(TEXT_SEC),
                    );
                }
                if m.total_dust_mass > 0.0 {
                    ui.label(
                        RichText::new(format!("untracked dust {:.3e}", m.total_dust_mass))
                            .size(9.5)
                            .color(TEXT_DIM),
                    );
                }
            });
        }
    }
}
