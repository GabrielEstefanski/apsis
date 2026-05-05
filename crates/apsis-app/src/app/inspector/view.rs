//! Inspector view — pure consumer of [`super::data::InspectorData`].
//!
//! Renders the structured payload through design primitives. Owns no
//! domain logic; recovers no quantities. NaN values pass through to the
//! formatter and surface as `—`. Auto-shown sections render only when
//! their corresponding payload is present.

use eframe::egui::{self, Align, FontFamily, FontId, Layout, RichText, Sense, Ui};
use egui_phosphor::regular::{CARET_DOWN, CARET_RIGHT};

use super::data::{
    ActionKind, AggregateData, CameraRelativeData, EnergyData, Identity, InspectorData,
    KinematicState, OrbitData, PerturbationData, RelationsData,
};
use super::format::{QuantityType, format_value};
use crate::app::design::primitives::{
    FieldRow, FlashTracker, IconButton, Section, Subgroup, hairline,
};
use crate::app::design::tokens::{border, color, shape, space, typography};

/// State that persists across frames — `More` expander toggles plus the
/// flash tracker for Bloomberg-style value-changed pulses. Lives outside
/// [`InspectorData`] so the data payload itself stays a pure value.
#[derive(Default)]
pub struct InspectorState {
    pub state_more_open: bool,
    pub orbit_more_open: bool,
    pub flash: FlashTracker,
}

/// Render one Inspector frame and return the index of the action that
/// was clicked this frame, if any. The caller dispatches intent against
/// its own action vocabulary; this module knows nothing about what
/// "Focus camera" or "Delete" mean.
pub fn show(ui: &mut Ui, data: &InspectorData, state: &mut InspectorState) -> Option<usize> {
    let mut clicked_action: Option<usize> = None;
    ui.allocate_ui_with_layout(
        egui::vec2(ui.available_width(), ui.available_height()),
        Layout::top_down(Align::LEFT),
        |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.add_space(space::S5);
                show_header(ui, &data.header);
                ui.add_space(space::S4);
                hairline(ui, space::S4);

                Section::new("Identity").show(ui, |ui| identity_rows(ui, &data.identity));
                Section::new("State").show(ui, |ui| {
                    state_rows(ui, &data.state, &mut state.state_more_open, &mut state.flash);
                });

                if let Some(orbit) = &data.orbit {
                    Section::new(&format!("Orbit  ›  {}", orbit.primary_name)).show(ui, |ui| {
                        orbit_rows(ui, orbit, &mut state.orbit_more_open, &mut state.flash);
                    });
                }

                if let Some(rel) = &data.relations {
                    Section::new("Relations").show(ui, |ui| relations_rows(ui, rel));
                }

                if let Some(energy) = &data.energy {
                    Section::new("Energy").show(ui, |ui| {
                        energy_rows(ui, energy, &mut state.flash);
                    });
                }

                if !data.perturbations.is_empty() {
                    Section::new("Perturbations").show(ui, |ui| {
                        for pert in &data.perturbations {
                            perturbation_block(ui, pert, &mut state.flash);
                        }
                    });
                }

                if let Some(cam) = &data.camera_relative {
                    Section::new("Camera-relative").show(ui, |ui| {
                        camera_rows(ui, cam, &mut state.flash);
                    });
                }

                if !data.actions.is_empty() {
                    Section::new("Actions").show(ui, |ui| {
                        clicked_action = action_rows(ui, &data.actions);
                    });
                }

                ui.add_space(space::S5);
            });
        },
    );
    clicked_action
}

/// Render the aggregate (multi-select) inspector frame and return the clicked
/// action index, if any. The caller dispatches actions by index.
pub fn show_aggregate(ui: &mut Ui, data: &AggregateData) -> Option<usize> {
    let mut clicked: Option<usize> = None;
    ui.allocate_ui_with_layout(
        egui::vec2(ui.available_width(), ui.available_height()),
        Layout::top_down(Align::LEFT),
        |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                ui.add_space(space::S5);
                show_aggregate_header(ui, data.count);
                ui.add_space(space::S4);
                hairline(ui, space::S4);

                Section::new("Selection").show(ui, |ui| aggregate_rows(ui, data));

                if !data.actions.is_empty() {
                    Section::new("Actions").show(ui, |ui| {
                        clicked = action_rows(ui, &data.actions);
                    });
                }

                ui.add_space(space::S5);
            });
        },
    );
    clicked
}

fn show_aggregate_header(ui: &mut Ui, count: usize) {
    let medium = FontFamily::Name(typography::font::SANS_MEDIUM.into());
    ui.horizontal(|ui| {
        ui.add_space(space::S4);
        ui.label(
            RichText::new(format!("{count} bodies"))
                .font(FontId::new(typography::text::LG, medium))
                .color(color::fg::PRIMARY),
        );
    });
    ui.horizontal(|ui| {
        ui.add_space(space::S4);
        ui.label(
            RichText::new("selection")
                .font(FontId::new(typography::text::XS, FontFamily::Proportional))
                .color(color::fg::TERTIARY),
        );
    });
}

fn aggregate_rows(ui: &mut Ui, data: &AggregateData) {
    let (s, u) = format_value(data.total_mass_kg, QuantityType::Mass);
    ui.add(FieldRow::new("Total mass", &s, u));

    Subgroup::new("COM position").show(ui, |ui| {
        for (axis, value) in ["x", "y", "z"].iter().zip(data.com_m.iter()) {
            let (val, unit) = format_value(*value, QuantityType::DistanceVector);
            ui.add(FieldRow::new(axis, &val, unit).indented(1));
        }
    });

    let v_speed =
        (data.v_com_m_s[0].powi(2) + data.v_com_m_s[1].powi(2) + data.v_com_m_s[2].powi(2)).sqrt();
    let (val, unit) = format_value(v_speed, QuantityType::VelocityVector);
    ui.add(FieldRow::new("|v_COM|", &val, unit));

    let (val, unit) = format_value(data.bounding_radius_m, QuantityType::DistanceVector);
    ui.add(FieldRow::new("Bounding radius", &val, unit));
}

// ── Header ───────────────────────────────────────────────────────────────────

fn show_header(ui: &mut Ui, header: &super::data::Header) {
    ui.horizontal(|ui| {
        let medium = FontFamily::Name(typography::font::SANS_MEDIUM.into());
        ui.add_space(space::S4);
        ui.label(
            RichText::new(&header.name)
                .font(FontId::new(typography::text::LG, medium))
                .color(color::fg::PRIMARY),
        );
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            ui.add_space(space::S4);
            let (rect, _) =
                ui.allocate_exact_size(egui::vec2(shape::SWATCH, shape::SWATCH), Sense::hover());
            ui.painter().rect_filled(rect, border::radius::NONE, header.swatch);
        });
    });
    ui.horizontal(|ui| {
        ui.add_space(space::S4);
        ui.label(
            RichText::new(&header.breadcrumb)
                .font(FontId::new(typography::text::XS, FontFamily::Proportional))
                .color(color::fg::TERTIARY),
        );
    });
}

// ── Identity ─────────────────────────────────────────────────────────────────

fn identity_rows(ui: &mut Ui, id: &Identity) {
    let (s, u) = format_value(id.mass_kg, QuantityType::Mass);
    ui.add(FieldRow::new("Mass", &s, u));
    if let Some(r) = id.radius_m {
        let (s, u) = format_value(r, QuantityType::DistanceScalar);
        ui.add(FieldRow::new("Radius", &s, u));
    }
}

// ── State ────────────────────────────────────────────────────────────────────

fn state_rows(ui: &mut Ui, s: &KinematicState, more_open: &mut bool, flash: &mut FlashTracker) {
    Subgroup::new("Position").show(ui, |ui| {
        for (axis, value) in ["x", "y", "z"].iter().zip(s.position_m.iter()) {
            let (val, unit) = format_value(*value, QuantityType::DistanceVector);
            let key = format!("state.pos.{axis}");
            let started = flash.observe(&key, &val);
            ui.add(FieldRow::new(axis, &val, unit).indented(1).flash(started));
        }
    });

    Subgroup::new("Velocity").show(ui, |ui| {
        for (axis, value) in ["vx", "vy", "vz"].iter().zip(s.velocity_m_s.iter()) {
            let (val, unit) = format_value(*value, QuantityType::VelocityVector);
            let key = format!("state.vel.{axis}");
            let started = flash.observe(&key, &val);
            ui.add(FieldRow::new(axis, &val, unit).indented(1).flash(started));
        }
    });

    more_disclosure(ui, more_open);
    if *more_open {
        let speed =
            (s.velocity_m_s[0].powi(2) + s.velocity_m_s[1].powi(2) + s.velocity_m_s[2].powi(2))
                .sqrt();
        let (val, unit) = format_value(speed, QuantityType::VelocityVector);
        ui.add(FieldRow::new("|v|", &val, unit));
    }
}

// ── Orbit ────────────────────────────────────────────────────────────────────

fn orbit_rows(ui: &mut Ui, o: &OrbitData, more_open: &mut bool, flash: &mut FlashTracker) {
    let (val, unit) = format_value(o.semi_major_axis_m, QuantityType::DistanceScalar);
    ui.add(FieldRow::new("Semi-major axis (a)", &val, unit));

    let (val, unit) = format_value(o.eccentricity, QuantityType::Eccentricity);
    ui.add(FieldRow::new("Eccentricity (e)", &val, unit));

    let (val, unit) = format_value(o.inclination_rad, QuantityType::AngleStatic);
    ui.add(FieldRow::new("Inclination (i)", &val, unit));

    let (val, unit) = format_value(o.period_s, QuantityType::Time);
    ui.add(FieldRow::new("Period (T)", &val, unit));

    more_disclosure(ui, more_open);
    if *more_open {
        let (val, unit) = format_value(o.lon_ascending_node_rad, QuantityType::AngleStatic);
        ui.add(FieldRow::new("Asc. node (Ω)", &val, unit));

        let (val, unit) = format_value(o.argument_of_pericenter_rad, QuantityType::AngleStatic);
        ui.add(FieldRow::new("Arg. of pericenter (ω)", &val, unit));

        let (val, unit) = format_value(o.true_anomaly_rad, QuantityType::AngleDynamic);
        let started = flash.observe("orbit.nu", &val);
        ui.add(FieldRow::new("True anomaly (ν)", &val, unit).flash(started));

        let (val, unit) = format_value(o.mean_anomaly_rad, QuantityType::AngleDynamic);
        let started = flash.observe("orbit.M", &val);
        ui.add(FieldRow::new("Mean anomaly (M)", &val, unit).flash(started));

        let (val, unit) = format_value(o.eccentric_anomaly_rad, QuantityType::AngleDynamic);
        let started = flash.observe("orbit.E", &val);
        ui.add(FieldRow::new("Eccentric anomaly (E)", &val, unit).flash(started));

        let (val, unit) = format_value(o.pericenter_m, QuantityType::DistanceScalar);
        ui.add(FieldRow::new("Pericenter (q)", &val, unit));

        let (val, unit) = format_value(o.apocenter_m, QuantityType::DistanceScalar);
        ui.add(FieldRow::new("Apocenter (Q)", &val, unit));
    }
}

// ── Relations ────────────────────────────────────────────────────────────────

fn relations_rows(ui: &mut Ui, rel: &RelationsData) {
    ui.add(FieldRow::new("Primary", &rel.primary_name, ""));
    ui.add(FieldRow::new("Relation", rel.kind.label(), ""));
    ui.add(FieldRow::new("Secondary", &rel.secondary_name, ""));
    ui.add(FieldRow::new("Frame", &rel.frame_label, ""));
}

// ── Energy ───────────────────────────────────────────────────────────────────

fn energy_rows(ui: &mut Ui, e: &EnergyData, flash: &mut FlashTracker) {
    let (val, unit) = format_value(e.kinetic_j, QuantityType::Energy);
    let started = flash.observe("energy.K", &val);
    ui.add(FieldRow::new("Kinetic", &val, unit).flash(started));

    let (val, unit) = format_value(e.potential_j, QuantityType::Energy);
    let started = flash.observe("energy.U", &val);
    ui.add(FieldRow::new("Potential", &val, unit).flash(started));

    let (val, unit) = format_value(e.specific_j, QuantityType::Energy);
    let started = flash.observe("energy.eps", &val);
    ui.add(FieldRow::new("Specific", &val, unit).flash(started));
}

// ── Perturbations ────────────────────────────────────────────────────────────

fn perturbation_block(ui: &mut Ui, p: &PerturbationData, flash: &mut FlashTracker) {
    ui.horizontal(|ui| {
        ui.add_space(space::S5);
        let (rect, _) =
            ui.allocate_exact_size(egui::vec2(shape::DOT_LIVE, shape::DOT_LIVE), Sense::hover());
        let dot_color = if p.active { color::signal::LIVE } else { color::fg::TERTIARY };
        if p.active {
            ui.painter().circle_filled(rect.center(), shape::DOT_LIVE * 0.5, dot_color);
        } else {
            ui.painter().circle_stroke(
                rect.center(),
                shape::DOT_LIVE * 0.5,
                eframe::egui::Stroke::new(border::width::HAIRLINE, dot_color),
            );
        }
        ui.add_space(space::S2);
        ui.label(
            RichText::new(&p.name)
                .font(FontId::new(typography::text::BASE, FontFamily::Proportional))
                .color(color::fg::PRIMARY),
        );
    });
    for r in &p.readouts {
        let key = format!("pert.{}.{}", p.name, r.label);
        let started = flash.observe(&key, &r.value_str);
        ui.add(FieldRow::new(&r.label, &r.value_str, &r.unit).indented(1).flash(started));
    }
}

// ── Camera-relative ──────────────────────────────────────────────────────────

fn camera_rows(ui: &mut Ui, c: &CameraRelativeData, flash: &mut FlashTracker) {
    let (val, unit) = format_value(c.distance_m, QuantityType::DistanceScalar);
    let started = flash.observe("cam.dist", &val);
    ui.add(FieldRow::new("Distance", &val, unit).flash(started));

    // Radial velocity in km/s, signed. Negative = approaching.
    let v_kms = c.radial_velocity_m_s / 1.0e3;
    let val = format!("{v_kms:+.1}");
    let started = flash.observe("cam.vrad", &val);
    ui.add(FieldRow::new("Approaching", &val, "km/s").flash(started));

    let (val, unit) = format_value(c.apparent_size_arcsec, QuantityType::Arcsecond);
    ui.add(FieldRow::new("Apparent size", &val, unit));

    let (val, unit) = format_value(c.off_axis_rad, QuantityType::AngleStatic);
    ui.add(FieldRow::new("Off-axis", &val, unit));
}

// ── Actions ──────────────────────────────────────────────────────────────────

fn action_rows(ui: &mut Ui, actions: &[super::data::ActionData]) -> Option<usize> {
    let mut clicked = None;
    for (i, a) in actions.iter().enumerate() {
        let mut btn = IconButton::new(&a.label);
        if let Some(ico) = &a.icon {
            btn = btn.with_icon(ico);
        }
        if let Some(sc) = &a.shortcut {
            btn = btn.with_shortcut(sc);
        }
        if a.kind == ActionKind::Destructive {
            btn = btn.danger();
        }
        if ui.add(btn).clicked() {
            clicked = Some(i);
        }
    }
    clicked
}

// ── Disclosure ───────────────────────────────────────────────────────────────

fn more_disclosure(ui: &mut Ui, open: &mut bool) {
    ui.add_space(space::S2);
    let row_height = space::S5;
    let total_w = ui.available_width();
    let (rect, response) = ui.allocate_exact_size(egui::vec2(total_w, row_height), Sense::click());

    if response.hovered() {
        ui.painter().rect_filled(rect, border::radius::SHARP, color::bg::HOVER);
    }
    if response.clicked() {
        *open = !*open;
    }

    let glyph = if *open { CARET_DOWN } else { CARET_RIGHT };
    let label = if *open { "Less" } else { "More" };
    let text_color = if response.hovered() { color::fg::SECONDARY } else { color::fg::TERTIARY };

    // Painter (non-interactive) instead of ui.label, so clicks land on
    // the outer rect rather than being captured by the label widget.
    ui.painter().text(
        egui::pos2(rect.left() + space::S5, rect.center().y),
        egui::Align2::LEFT_CENTER,
        format!("{glyph}  {label}"),
        FontId::new(typography::text::SM, FontFamily::Proportional),
        text_color,
    );
}
