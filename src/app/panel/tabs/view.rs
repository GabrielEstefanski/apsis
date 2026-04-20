//! Display tool — toggles for visual layers rendered on top of bodies.
//!
//! All controls here are purely presentational (do not affect physics).
//! Migrated from the old top-bar toggle strip to give each layer breathing
//! room and space for per-layer sub-controls (e.g. trail width).

use crate::app::theme::{BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC, section};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn panel_tab_view(&mut self, ui: &mut egui::Ui) {
        ui.add_space(2.0);
        ui.label(
            RichText::new("Display")
                .size(13.0)
                .color(TEXT_PRI)
                .strong(),
        );
        ui.label(
            RichText::new("Visual layers — no effect on physics.")
                .size(10.0)
                .color(TEXT_DIM),
        );

        section(ui, "LAYERS");

        toggle_row(ui, &mut self.show_grid, "Grid", "Reference grid in world units");
        toggle_row(ui, &mut self.show_trails, "Trails", "Body position history");
        if self.show_trails {
            ui.indent("trail_opts", |ui| {
                ui.horizontal(|ui| {
                    ui.label(RichText::new("width").size(10.0).color(TEXT_SEC));
                    ui.add(
                        egui::DragValue::new(&mut self.trail_width)
                            .speed(0.1)
                            .range(0.5_f32..=20.0)
                            .max_decimals(1)
                            .suffix(" px"),
                    );
                });
            });
        }
        toggle_row(
            ui,
            &mut self.show_orbit_ellipses,
            "Orbit ellipses",
            "Keplerian fit of each body's trajectory",
        );

        section(ui, "VECTORS");

        toggle_row(ui, &mut self.show_vectors, "Velocity", "Instantaneous v for each body");
        toggle_row(
            ui,
            &mut self.show_force_vectors,
            "Force",
            "Net gravitational force for each body",
        );

        section(ui, "DIAGNOSTIC");

        toggle_row(
            ui,
            &mut self.show_belts,
            "Tree structure",
            "Barnes-Hut cells & asteroid belt hints",
        );
    }
}

fn toggle_row(ui: &mut egui::Ui, value: &mut bool, label: &str, hover: &str) {
    let col = if *value { TEXT_PRI } else { TEXT_SEC };
    let resp = ui.add(
        egui::Button::new(
            RichText::new(format!("{}  {}", if *value { "●" } else { "○" }, label))
                .size(11.0)
                .color(col),
        )
        .fill(Color32::TRANSPARENT)
        .stroke(Stroke::new(0.5, if *value { BORDER } else { Color32::TRANSPARENT }))
        .min_size(egui::vec2(ui.available_width(), 24.0))
        .corner_radius(4.0),
    );
    if resp.clicked() {
        *value = !*value;
    }
    resp.on_hover_text(hover);
}
