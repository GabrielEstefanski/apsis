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

use crate::app::theme::{ACCENT, BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC, section};
use crate::app::ui::{PanelTab, SelectionForm, SimulationApp};
use eframe::egui::{self, Color32, RichText, Stroke};

const TOP_N: usize = 5;

impl SimulationApp {
    pub(super) fn panel_tab_overview(&mut self, ui: &mut egui::Ui) {
        let n_bodies = self.system.bodies().len();

        if n_bodies == 0 {
            self.overview_empty_state(ui);
            return;
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

        // ── Most massive bodies ────────────────────────────────────────────── //
        section(ui, "TOP BODIES");

        // Collect (index, mass, name) for the N heaviest. We keep the index so
        // the row stays a navigation target even when names collide or are
        // empty. Small N + one scan avoids any sorting pressure.
        let bodies = self.system.bodies();
        let mut idx_mass: Vec<(usize, f64)> =
            bodies.iter().enumerate().map(|(i, b)| (i, b.mass)).collect();
        idx_mass.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let take = idx_mass.len().min(TOP_N);

        let mut clicked_idx: Option<usize> = None;
        for (rank, (idx, mass)) in idx_mass[..take].iter().enumerate() {
            let name = {
                let n = self.system.name(*idx);
                if n.is_empty() { format!("body #{idx}") } else { n.to_owned() }
            };
            let is_selected = self.selected_body == Some(*idx);

            if body_row(ui, rank + 1, &name, *mass, is_selected).clicked() {
                clicked_idx = Some(*idx);
            }
        }

        if let Some(idx) = clicked_idx {
            // Single click → select + focus the inspector. Double-click UX
            // (follow + zoom) already lives on the canvas; we keep this path
            // minimal so it composes predictably.
            self.selected_body = Some(idx);
            if let Some(b) = self.system.bodies().get(idx) {
                self.selection_form = Some(SelectionForm::from_body(b, self.system.name(idx)));
            }
        }

        if n_bodies > TOP_N {
            ui.add_space(2.0);
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
        ui.label(
            RichText::new("Empty scene")
                .size(13.0)
                .color(TEXT_PRI)
                .strong(),
        );
        ui.label(
            RichText::new("Add a body, or load a preset to get started.")
                .size(10.0)
                .color(TEXT_DIM),
        );

        section(ui, "QUICK START");

        let shortcut_row = |ui: &mut egui::Ui, tab: PanelTab, blurb: &str| -> bool {
            ui.horizontal(|ui| {
                let clicked = ui
                    .add(
                        egui::Button::new(
                            RichText::new(tab.label()).size(10.5).color(TEXT_PRI),
                        )
                        .fill(Color32::from_rgb(20, 20, 26))
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(86.0, 22.0)),
                    )
                    .clicked();
                ui.label(RichText::new(blurb).size(10.0).color(TEXT_SEC));
                clicked
            })
            .inner
        };

        if shortcut_row(ui, PanelTab::Templates, "· load a preset scenario") {
            self.panel_tab = PanelTab::Templates;
        }
        ui.add_space(3.0);
        if shortcut_row(ui, PanelTab::Add, "· place a body manually") {
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

// ── Helpers ──────────────────────────────────────────────────────────────── //

fn kv(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(10.0).color(TEXT_SEC));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(value).monospace().size(10.5).color(TEXT_PRI));
        });
    });
}

fn body_row(
    ui: &mut egui::Ui,
    rank: usize,
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

    let mut job = egui::text::LayoutJob::default();
    job.append(
        &format!("{rank}."),
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

    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        let resp = ui.add(
            egui::Button::new(job)
                .fill(fill)
                .stroke(Stroke::new(0.5, stroke_col))
                .min_size(egui::vec2(ui.available_width() - 72.0, 20.0))
                .corner_radius(3.0),
        );
        ui.add_space(4.0);
        ui.label(RichText::new(sci(mass)).monospace().size(10.0).color(TEXT_DIM));
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
