//! Display tool — everything that changes how the simulation is *shown*,
//! never how it is *computed*.
//!
//! Owns three families of controls:
//!
//! 1. **Layers**   — grid, trails, orbit ellipses, belts (toggles + per-layer
//!    sub-controls like trail width / sampling / min-mass-ratio).
//! 2. **Vectors**  — velocity / force overlays.
//! 3. **Colour**   — SPLASH / yt-style data-driven colouring
//!    (`field × normalizer × colormap`). `None` = material colours.
//!
//! Physics parameters (θ, ε, G, seed, integrator) live in the Advanced
//! (Physics) tab. Time/speed/algorithm live in the playbar. Keeping this
//! separation hard prevents the "grab-bag tab" that the old Config panel
//! grew into.

use crate::app::theme::{ACCENT_DIM, BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC, section};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Color32, RichText, Stroke};

// ── Layout constants (shared with Advanced / Physics tab) ────────────────────

const DV_W: f32 = 72.0;
const LBL_W: f32 = 80.0;

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
            RichText::new("Visual layers & colouring — no effect on physics.")
                .size(10.0)
                .color(TEXT_DIM),
        );

        // ── LAYERS ──────────────────────────────────────────────────────────
        section(ui, "LAYERS");

        toggle_row(ui, &mut self.show_grid, "Grid", "Reference grid in world units");

        toggle_row(ui, &mut self.show_trails, "Trails", "Body position history");
        if self.show_trails {
            ui.indent("trail_opts", |ui| {
                // Trail width (presentation).
                kv_drag(
                    ui,
                    "width",
                    "Line width in pixels.",
                    |ui| {
                        ui.add(
                            egui::DragValue::new(&mut self.trail_width)
                                .speed(0.1)
                                .range(0.5_f32..=20.0)
                                .max_decimals(1)
                                .suffix(" px"),
                        );
                    },
                );

                // Sample density — physics-thread sampler config (arc-length
                // threshold multiplier). Presented here because it's a
                // visual-quality knob, not a physics knob.
                let mut trail_every = self.trail_recorder.interval_multiplier();
                let te_tip = "Record one trail point every N arc-length\n\
                    trigger events. 1 = max density; higher = sparser,\n\
                    longer-lived trails.";
                let changed = kv_drag(ui, "sample every", te_tip, |ui| {
                    ui.add(
                        egui::DragValue::new(&mut trail_every)
                            .speed(1)
                            .range(1..=256usize),
                    )
                    .changed()
                });
                if changed {
                    self.trail_recorder.set_interval_multiplier(trail_every);
                    self.system.set_trail_sampler(self.trail_recorder.sampler_kind());
                }

                let mr_tip = "Minimum mass ratio (body / dominant body)\n\
                    for a trail to be shown. Raise to hide asteroid-mass\n\
                    objects; 0 = show all bodies.";
                kv_drag(ui, "min mass ratio", mr_tip, |ui| {
                    ui.add(
                        egui::DragValue::new(&mut self.trail_min_mass_ratio)
                            .speed(1e-8)
                            .range(0.0_f64..=1.0)
                            .custom_formatter(|v, _| format!("{v:.1e}"))
                            .custom_parser(|s| s.parse::<f64>().ok()),
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
        toggle_row(
            ui,
            &mut self.show_belts,
            "Tree structure",
            "Barnes-Hut cells & asteroid belt hints",
        );

        // ── VECTORS ─────────────────────────────────────────────────────────
        section(ui, "VECTORS");

        toggle_row(ui, &mut self.show_vectors, "Velocity", "Instantaneous v for each body");
        toggle_row(
            ui,
            &mut self.show_force_vectors,
            "Force",
            "Net gravitational force for each body",
        );

        // ── COLOUR ──────────────────────────────────────────────────────────
        // Data-driven colouring (SPLASH / yt-style). Off = material colours.
        section(ui, "COLOUR");
        self.colour_section(ui);
    }

    /// Renders the COLOUR subsection: enable toggle + field/colormap/normalizer
    /// dropdowns + resolved-range readout.
    fn colour_section(&mut self, ui: &mut egui::Ui) {
        use crate::render::color::ColorViewSelection;

        let mut enabled = self.color_view.is_some();
        let resp = ui.checkbox(
            &mut enabled,
            RichText::new("Colour by data").size(10.5).color(TEXT_PRI),
        )
        .on_hover_text(
            "Enable data-driven colouring (SPLASH / yt-style).\n\
             Disabled: each body uses its material colour.",
        );
        if resp.changed() {
            if enabled {
                self.color_view = Some(ColorViewSelection::default_velocity());
            } else {
                self.color_view = None;
                self.color_view_range = None;
            }
        }

        let Some(sel) = self.color_view.as_mut() else {
            return;
        };

        ui.add_space(2.0);

        // ── Field dropdown ───────────────────────────────────────────────
        let current_field_name = self
            .field_registry
            .get(&sel.field_id)
            .map(|f| f.name())
            .unwrap_or("(?)")
            .to_string();
        let mut new_field_id: Option<String> = None;
        let mut new_prefers_log: Option<bool> = None;
        kv_combo(
            ui,
            "field",
            "Scalar quantity sampled per body.\n\
             Velocity, mass, acceleration, kinetic energy.",
            "view_field",
            current_field_name,
            |ui| {
                for f in self.field_registry.iter() {
                    let selected = sel.field_id == f.id();
                    if ui
                        .selectable_label(
                            selected,
                            RichText::new(f.name()).size(10.0).color(TEXT_PRI),
                        )
                        .clicked()
                        && !selected
                    {
                        new_field_id = Some(f.id().to_string());
                        new_prefers_log = Some(f.prefers_log());
                    }
                }
            },
        );
        if let Some(id) = new_field_id {
            sel.field_id = id;
            // Auto-pick a sensible normalizer when the field changes.
            sel.normalizer_id = if new_prefers_log.unwrap_or(false) {
                "log".into()
            } else {
                "linear".into()
            };
            sel.range = None;
        }

        // ── Colormap dropdown ─────────────────────────────────────────────
        let current_cm_name = self
            .colormap_registry
            .get(&sel.colormap_id)
            .map(|c| c.name())
            .unwrap_or("(?)")
            .to_string();
        let mut new_cm: Option<String> = None;
        kv_combo(
            ui,
            "colormap",
            "Colour ramp.\n\
             Perceptually-uniform: viridis / inferno / plasma\n\
             Diverging: cool-warm",
            "view_colormap",
            current_cm_name,
            |ui| {
                for c in self.colormap_registry.iter() {
                    let selected = sel.colormap_id == c.id();
                    if ui
                        .selectable_label(
                            selected,
                            RichText::new(c.name()).size(10.0).color(TEXT_PRI),
                        )
                        .clicked()
                        && !selected
                    {
                        new_cm = Some(c.id().to_string());
                    }
                }
            },
        );
        if let Some(id) = new_cm {
            sel.colormap_id = id;
        }

        // ── Normalizer dropdown ───────────────────────────────────────────
        let current_nm_name = self
            .normalizer_registry
            .get(&sel.normalizer_id)
            .map(|n| n.name())
            .unwrap_or("(?)")
            .to_string();
        let mut new_nm: Option<String> = None;
        kv_combo(
            ui,
            "normalize",
            "How values map into [0, 1] before the colour ramp.\n\
             Log is preferred for mass / acceleration (many decades).",
            "view_normalizer",
            current_nm_name,
            |ui| {
                for n in self.normalizer_registry.iter() {
                    let selected = sel.normalizer_id == n.id();
                    if ui
                        .selectable_label(
                            selected,
                            RichText::new(n.name()).size(10.0).color(TEXT_PRI),
                        )
                        .clicked()
                        && !selected
                    {
                        new_nm = Some(n.id().to_string());
                    }
                }
            },
        );
        if let Some(id) = new_nm {
            sel.normalizer_id = id;
        }

        // ── Range readout ─────────────────────────────────────────────────
        if let Some((lo, hi)) = self.color_view_range {
            let unit = self
                .field_registry
                .get(&sel.field_id)
                .map(|f| f.unit_label())
                .unwrap_or("");
            ui.add_space(2.0);
            ui.label(
                RichText::new(format!("range {unit}: {lo:.3e} … {hi:.3e}"))
                    .size(9.0)
                    .color(TEXT_DIM)
                    .monospace(),
            );
        }

        let _ = ACCENT_DIM; // silence unused import if theme changes later
    }
}

// ── Shared helpers ───────────────────────────────────────────────────────────

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

/// Label on the left, a DragValue-style widget on the right. Used inside the
/// indented trail-options block so rows visually align with the rest of the tab.
fn kv_drag<R>(
    ui: &mut egui::Ui,
    label: &str,
    tip: &str,
    add: impl FnOnce(&mut egui::Ui) -> R,
) -> R {
    ui.horizontal(|ui| {
        ui.add_sized(
            egui::vec2(LBL_W, 18.0),
            egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)),
        )
        .on_hover_text(tip);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_sized(egui::vec2(DV_W, 18.0), egui::Label::new("")); // reserve width
            add(ui)
        })
        .inner
    })
    .inner
}

/// Label on the left, ComboBox on the right.
fn kv_combo(
    ui: &mut egui::Ui,
    label: &str,
    tip: &str,
    id_salt: &str,
    current: String,
    contents: impl FnOnce(&mut egui::Ui),
) {
    ui.horizontal(|ui| {
        ui.add_sized(
            egui::vec2(LBL_W, 18.0),
            egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)),
        )
        .on_hover_text(tip);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            egui::ComboBox::from_id_salt(id_salt)
                .selected_text(RichText::new(current).size(10.0).color(TEXT_PRI))
                .show_ui(ui, contents);
        });
    });
}
