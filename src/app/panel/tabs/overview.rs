//! Overview tab — the contextual panel's landing view.
//!
//! Replaces the earlier "Select" tab, which only existed to host selection
//! hints. A Select *tab* is redundant because selection is a canvas gesture:
//! clicking a body opens the right-hand inspector automatically. What was
//! missing was a neutral home surface that's useful when no tool is engaged.
//!
//! The Overview fills that gap. It surfaces:
//!
//! 1. a one-glance system summary (N bodies, total mass, COM);
//! 2. the instantaneous energy breakdown (K, U, E);
//! 3. the five most massive bodies, each clickable to jump to the inspector;
//! 4. a first-use guide when the scene is empty.
//!
//! All of these are derived values — this tab never mutates simulation state
//! except as a navigation shortcut (selection + follow).

use crate::app::icons;
use crate::app::panel::metrics::DriftSeverity;
use crate::app::theme::{
    ACCENT, BORDER, DANGER, SUCCESS, SURFACE_CARD, TEXT_DIM, TEXT_PRI, TEXT_SEC, section,
};
use crate::app::ui::{PanelTab, SelectionForm, SimulationApp};
use eframe::egui::text::{LayoutJob, TextFormat};
use eframe::egui::{self, Align, Color32, FontId, RichText, Stroke};

const TOP_N: usize = 5;
const SEARCH_LIMIT: usize = 50;

impl SimulationApp {
    pub(super) fn panel_tab_overview(&mut self, ui: &mut egui::Ui) {
        let n_bodies = self.system.bodies().len();

        if n_bodies == 0 {
            self.overview_empty_state(ui);
            return;
        }

        // ── Compact status bar ────────────────────────────────────────────── //
        {
            use std::f64::consts::PI;
            let e_sev = DriftSeverity::from_peak(self.energy_drift_peak);
            let lz_sev = DriftSeverity::from_peak(self.lz_drift_peak);
            let worst = e_sev.max(lz_sev);
            let col = worst.color();
            let m = self.system.metrics();
            let t = m.t;
            let yr = t / (2.0 * PI);
            let yr_str = if yr.abs() < 0.01 {
                format!("{:.2e} yr", yr)
            } else if yr.abs() < 10_000.0 {
                format!("{:.2} yr", yr)
            } else {
                format!("{:.2e} yr", yr)
            };
            ui.horizontal(|ui| {
                ui.label(RichText::new(worst.dot()).size(9.0).color(col));
                ui.label(RichText::new(worst.label()).size(10.0).color(col).strong());
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new(yr_str).monospace().size(9.5).color(TEXT_DIM));
                    ui.label(RichText::new("·").size(9.0).color(TEXT_DIM));
                    ui.label(
                        RichText::new(format!("T {:.3e}", t))
                            .monospace()
                            .size(9.5)
                            .color(TEXT_SEC),
                    );
                });
            });
            ui.add_space(2.0);
        }

        // ── System ─────────────────────────────────────────────────────────── //
        section(ui, "SYSTEM");
        let m = self.system.metrics();
        let total_mass = self.system.total_mass();

        kv(ui, "bodies", &format!("{n_bodies}"));
        kv(ui, "total mass", &sci(total_mass));
        kv(ui, "COM", &format!("({:+.3}, {:+.3})", m.com_x, m.com_y));

        // ── Energy ─────────────────────────────────────────────────────────── //
        section(ui, "ENERGY");
        kv(ui, "E (total)", &sci(m.total_energy));
        kv(ui, "K (kinetic)", &sci(m.kinetic));
        kv(ui, "U (potential)", &sci(m.potential));
        kv(ui, "Lz", &sci(m.angular_momentum_z));
        kv_col(
            ui,
            "ΔE / E₀",
            &format!("{:+.2e}", m.rel_energy_error),
            instant_color(m.rel_energy_error),
        );
        let lz_triv = m.angular_momentum_z.abs() < 1e-10;
        if !lz_triv {
            kv_col(
                ui,
                "ΔLz / Lz₀",
                &format!("{:+.2e}", m.rel_angular_momentum_error),
                instant_color(m.rel_angular_momentum_error),
            );
        }

        // ── Solver ─────────────────────────────────────────────────────────── //
        section(ui, "SOLVER");
        kv(ui, "dt", &format!("{:.3e}", m.dt));
        kv(ui, "θ (opening)", &format!("{:.3}", m.theta));
        kv(ui, "steps", &format!("{}", m.steps));
        kv(ui, "vmax", &format!("{:.3e}", m.max_vel));
        kv(ui, "amax", &format!("{:.3e}", m.max_acc));
        if let Some(rec) = m.recommended_dt {
            let ratio = m.dt / rec;
            let col = if ratio <= 2.0 {
                TEXT_DIM
            } else if ratio <= 10.0 {
                ACCENT
            } else {
                DANGER
            };
            kv_col(ui, "suggested dt", &format!("{:.3e}", rec), col);
        }

        // ── Stability (peak drift since reset) ──────────────────────────────── //
        if m.steps > 10 && (self.energy_drift_peak > 0.0 || self.lz_drift_peak > 0.0) {
            section(ui, "STABILITY");
            if self.energy_drift_peak > 0.0 {
                kv_col(
                    ui,
                    "peak ΔE/E₀",
                    &format!("{:.2e}", self.energy_drift_peak),
                    DriftSeverity::from_peak(self.energy_drift_peak).color(),
                );
            }
            if self.lz_drift_peak > 0.0 {
                kv_col(
                    ui,
                    "peak ΔLz/Lz₀",
                    &format!("{:.2e}", self.lz_drift_peak),
                    DriftSeverity::from_peak(self.lz_drift_peak).color(),
                );
            }
        }

        // ── Bodies list with search ────────────────────────────────────────── //
        section(ui, "BODIES");

        // Search bar
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 4.0;
            ui.label(RichText::new(icons::SEARCH).size(11.0).color(TEXT_DIM));
            let search = ui.add(
                egui::TextEdit::singleline(&mut self.overview_search)
                    .desired_width(ui.available_width() - 24.0)
                    .hint_text("Search bodies…")
                    .font(FontId::proportional(10.5))
                    .text_color(TEXT_PRI),
            );
            // Clear button
            if !self.overview_search.is_empty() {
                if ui
                    .add(
                        egui::Button::new(RichText::new("×").size(11.0).color(TEXT_DIM))
                            .fill(Color32::TRANSPARENT)
                            .frame(false),
                    )
                    .on_hover_text("Clear search")
                    .clicked()
                {
                    self.overview_search.clear();
                    search.request_focus();
                }
            }
        });
        ui.add_space(3.0);

        // Snapshot body data — avoids borrow conflicts when we later call
        // self.system.name() and self.selected_body simultaneously.
        let entries: Vec<(usize, String, f64)> = {
            let bodies = self.system.bodies();
            (0..bodies.len())
                .map(|i| {
                    let raw = self.system.name(i);
                    let name =
                        if raw.is_empty() { format!("body #{i}") } else { raw.to_owned() };
                    (i, name, bodies[i].mass)
                })
                .collect()
        };

        let query = self.overview_search.trim().to_lowercase();
        let searching = !query.is_empty();

        // Build display list: filtered search OR top-N by mass
        let display: Vec<(usize, String, f64, Option<usize>)> = if searching {
            entries
                .into_iter()
                .filter(|(_, name, _)| name.to_lowercase().contains(&query))
                .take(SEARCH_LIMIT)
                .map(|(i, name, mass)| (i, name, mass, None))
                .collect()
        } else {
            let mut sorted = entries;
            sorted.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
            sorted
                .into_iter()
                .take(TOP_N)
                .enumerate()
                .map(|(rank, (i, name, mass))| (i, name, mass, Some(rank + 1)))
                .collect()
        };

        let match_count = display.len();
        let mut clicked_idx: Option<usize> = None;

        for (idx, name, mass, rank) in &display {
            let is_selected = self.selected_body == Some(*idx);
            if body_row(ui, *rank, name, *mass, is_selected).clicked() {
                clicked_idx = Some(*idx);
            }
        }

        if let Some(idx) = clicked_idx {
            self.selected_body = Some(idx);
            if let Some(b) = self.system.bodies().get(idx) {
                self.selection_form =
                    Some(SelectionForm::from_body(b, self.system.name(idx)));
            }
        }

        ui.add_space(2.0);
        if searching {
            let label = match match_count {
                0 => "no matches".to_owned(),
                1 => "1 match".to_owned(),
                n => format!("{n} matches"),
            };
            ui.label(RichText::new(label).size(9.5).color(TEXT_DIM).italics());
        } else if n_bodies > TOP_N {
            ui.label(
                RichText::new(format!("… +{} more", n_bodies - TOP_N))
                    .size(9.5)
                    .color(TEXT_DIM)
                    .italics(),
            );
        }
    }

    // ── Empty scene: first-use guidance ──────────────────────────────────── //
    fn overview_empty_state(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);
        egui::Frame::NONE
            .fill(SURFACE_CARD)
            .stroke(Stroke::new(0.5, BORDER))
            .corner_radius(8.0)
            .inner_margin(egui::Margin::symmetric(12, 10))
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                ui.label(
                    RichText::new("No bodies in this scene")
                        .size(13.0)
                        .color(TEXT_PRI)
                        .strong(),
                );
                ui.add_space(4.0);
                ui.label(
                    RichText::new(
                        "Choose a curated preset or place objects by hand — the canvas is ready.",
                    )
                    .size(10.5)
                    .color(TEXT_SEC),
                );
            });

        ui.add_space(6.0);
        section(ui, "QUICK START");

        if quick_start_card(
            ui,
            icons::TOOL_TEMPLATES,
            "Presets",
            "Solar systems, binaries, choreographies…",
            "Tab 3",
        ) {
            self.panel_tab = PanelTab::Templates;
        }
        ui.add_space(6.0);
        if quick_start_card(
            ui,
            icons::TOOL_ADD,
            "Add objects",
            "Single body, rings, or clusters — place on the canvas",
            "Tab 2",
        ) {
            self.panel_tab = PanelTab::Add;
        }

        section(ui, "TIPS");
        tip(ui, "Space", "play / pause");
        tip(ui, "F", "fit view to all bodies");
        tip(ui, "B", "hide / show sidebar");
        tip(ui, "H", "open keyboard shortcuts");
        tip(ui, "Ctrl+Z", "undo last change");

        let _ = ACCENT; // reserved for highlight state
    }
}

/// Wide “simulator tile” row — Universe Sandbox–style entry point.
fn quick_start_card(
    ui: &mut egui::Ui,
    icon: &str,
    title: &str,
    subtitle: &str,
    hint: &str,
) -> bool {
    let w = ui.available_width();
    let mut job = LayoutJob::default();
    job.append(
        icon,
        0.0,
        TextFormat {
            font_id: FontId::proportional(18.0),
            color: TEXT_SEC,
            valign: Align::Center,
            ..Default::default()
        },
    );
    job.append(
        &format!("  {title}\n"),
        0.0,
        TextFormat {
            font_id: FontId::proportional(12.0),
            color: TEXT_PRI,
            valign: Align::Center,
            ..Default::default()
        },
    );
    job.append(
        &format!("  {subtitle}\n"),
        0.0,
        TextFormat {
            font_id: FontId::proportional(10.0),
            color: TEXT_DIM,
            ..Default::default()
        },
    );
    job.append(
        &format!("  {hint}"),
        0.0,
        TextFormat {
            font_id: FontId::proportional(9.0),
            color: TEXT_DIM,
            ..Default::default()
        },
    );

    ui.add(
        egui::Button::new(job)
            .fill(SURFACE_CARD)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(w, 58.0))
            .corner_radius(8.0),
    )
    .clicked()
}

// ── Helpers ──────────────────────────────────────────────────────────────── //

fn kv(ui: &mut egui::Ui, label: &str, value: &str) {
    kv_col(ui, label, value, TEXT_PRI);
}

fn kv_col(ui: &mut egui::Ui, label: &str, value: &str, col: Color32) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(10.0).color(TEXT_SEC));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(value).monospace().size(10.5).color(col));
        });
    });
}

fn instant_color(v: f64) -> Color32 {
    let a = v.abs();
    if a < 1e-8 {
        SUCCESS
    } else if a < 1e-5 {
        TEXT_DIM
    } else if a < 1e-3 {
        ACCENT
    } else {
        DANGER
    }
}

fn body_row(
    ui: &mut egui::Ui,
    rank: Option<usize>,
    name: &str,
    mass: f64,
    selected: bool,
) -> egui::Response {
    let fill = if selected {
        Color32::from_rgba_unmultiplied(ACCENT.r(), ACCENT.g(), ACCENT.b(), 28)
    } else {
        Color32::from_rgb(20, 20, 26)
    };
    let stroke_col = if selected { ACCENT.gamma_multiply(0.6) } else { BORDER };

    // Fixed-width mass column avoids the "button eats all width" problem.
    // Monospace 10.0 pt → ~6 px/char; scientific notation needs ≤ 11 chars.
    const MASS_COL: f32 = 68.0;
    const GAP: f32 = 4.0;

    let mut job = egui::text::LayoutJob::default();
    if let Some(r) = rank {
        job.append(
            &format!("{r}."),
            0.0,
            egui::TextFormat {
                font_id: egui::FontId::monospace(10.0),
                color: TEXT_DIM,
                valign: egui::Align::Center,
                ..Default::default()
            },
        );
        job.append(
            &format!("  {name}"),
            0.0,
            egui::TextFormat {
                font_id: egui::FontId::proportional(10.5),
                color: TEXT_PRI,
                valign: egui::Align::Center,
                ..Default::default()
            },
        );
    } else {
        job.append(
            name,
            0.0,
            egui::TextFormat {
                font_id: egui::FontId::proportional(10.5),
                color: TEXT_PRI,
                valign: egui::Align::Center,
                ..Default::default()
            },
        );
    }

    let btn_w = (ui.available_width() - MASS_COL - GAP).max(60.0);
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        let resp = ui.add(
            egui::Button::new(job)
                .fill(fill)
                .stroke(Stroke::new(0.5, stroke_col))
                .min_size(egui::vec2(btn_w, 20.0))
                .corner_radius(3.0),
        );
        ui.add_space(GAP);
        ui.add_sized(
            egui::vec2(MASS_COL, 20.0),
            egui::Label::new(
                RichText::new(sci(mass)).monospace().size(10.0).color(TEXT_DIM),
            ),
        );
        resp
    })
    .inner
}

fn tip(ui: &mut egui::Ui, key: &str, what: &str) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(key).size(10.0).color(TEXT_PRI).monospace());
        ui.label(RichText::new("·").size(10.0).color(TEXT_DIM));
        ui.label(RichText::new(what).size(10.0).color(TEXT_SEC));
    });
}

fn sci(v: f64) -> String {
    let a = v.abs();
    if a == 0.0 {
        "0".into()
    } else if (1e-3..1e5).contains(&a) {
        format!("{:.3}", v)
    } else {
        format!("{:+.3e}", v)
    }
}
