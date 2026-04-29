use crate::app::theme::{ACCENT, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::theme::{field, metric, primary_btn, secondary_btn, section};
use crate::app::ui::{SelectionForm, SimulationApp, UndoRecord};
use apsis::domain::body::{Body, radius_from_density_mass};
use apsis::physics::orbital::OrbitType;
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn inspector_content(&mut self, ui: &mut egui::Ui, idx: usize) {
        // Inspector stays visible during a Precision Run (read-only
        // view of the selected body is useful) but every interactive
        // widget below is disabled.
        if self.is_editing_locked() {
            ui.disable();
        }
        let body = self.system.bodies()[idx];

        ui.horizontal(|ui| {
            let [cr, cg, cb] = body.color;
            let col = egui::Color32::from_rgb(cr, cg, cb);
            let (dot_rect, _) =
                ui.allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
            ui.painter().circle_filled(dot_rect.center(), 6.0, col);

            let raw_name = self.system.name(idx);
            let display_name =
                if raw_name.is_empty() { format!("body #{idx}") } else { raw_name.to_owned() };
            ui.label(RichText::new(display_name).size(12.5).color(TEXT_PRI).strong());
        });

        // ── Camera shortcuts ─────────────────────────────────────────────── //
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 4.0;
            let follow_col = if self.follow_selected_body { ACCENT } else { TEXT_DIM };
            let follow_label = if self.follow_selected_body { "Following" } else { "Follow" };
            if ui
                .add(
                    egui::Button::new(
                        RichText::new(format!("⊙ {follow_label}")).size(9.5).color(follow_col),
                    )
                    .fill(egui::Color32::TRANSPARENT)
                    .stroke(egui::Stroke::new(
                        0.5,
                        if self.follow_selected_body {
                            ACCENT.gamma_multiply(0.5)
                        } else {
                            crate::app::theme::BORDER
                        },
                    ))
                    .min_size(egui::vec2(72.0, 20.0))
                    .corner_radius(3.0),
                )
                .on_hover_text("Lock camera to this body")
                .clicked()
            {
                self.follow_selected_body = !self.follow_selected_body;
            }

            if ui
                .add(
                    egui::Button::new(RichText::new("⊕ Fit").size(9.5).color(TEXT_DIM))
                        .fill(egui::Color32::TRANSPARENT)
                        .stroke(egui::Stroke::new(0.5, crate::app::theme::BORDER))
                        .min_size(egui::vec2(44.0, 20.0))
                        .corner_radius(3.0),
                )
                .on_hover_text("Zoom to fit this body")
                .clicked()
            {
                self.fit_to_view();
            }

            let pinned = self.pinned_orbits.contains(&idx);
            let pin_col = if pinned { ACCENT } else { TEXT_DIM };
            let pin_label = if pinned { "Pinned" } else { "Pin" };
            if ui
                .add(
                    egui::Button::new(
                        RichText::new(format!("⚲ {pin_label}")).size(9.5).color(pin_col),
                    )
                    .fill(egui::Color32::TRANSPARENT)
                    .stroke(egui::Stroke::new(
                        0.5,
                        if pinned { ACCENT.gamma_multiply(0.5) } else { crate::app::theme::BORDER },
                    ))
                    .min_size(egui::vec2(56.0, 20.0))
                    .corner_radius(3.0),
                )
                .on_hover_text(
                    "Keep this body's orbit visible regardless of the\n\
                     global orbit filters (levels, top-N, degeneracy).",
                )
                .clicked()
            {
                if pinned {
                    self.pinned_orbits.remove(&idx);
                } else {
                    self.pinned_orbits.insert(idx);
                }
            }
        });

        ui.add_space(10.0);

        let ul = &self.physics_cfg.dist_label.clone();
        let um = &self.physics_cfg.mass_label.clone();
        let ut = &self.physics_cfg.time_label.clone();

        section(ui, "LIVE");
        metric(ui, "x", &format!("{:.5e} {ul}", body.x), TEXT_DIM);
        metric(ui, "y", &format!("{:.5e} {ul}", body.y), TEXT_DIM);
        metric(ui, "vx", &format!("{:.5e} {ul}/{ut}", body.vx), TEXT_DIM);
        metric(ui, "vy", &format!("{:.5e} {ul}/{ut}", body.vy), TEXT_DIM);
        metric(ui, "mass", &format!("{:.5e} {um}", body.mass), TEXT_DIM);
        metric(ui, "r", &format!("{:.5e} {ul}", body.physical_radius), TEXT_DIM);
        metric(ui, "soft ε", &format!("{:.5e} {ul}", body.softening), TEXT_DIM);
        metric(ui, "ρ", &format!("{:.4e} {um}/{ul}³", body.density), TEXT_DIM);

        // ── ORBITAL ELEMENTS ─────────────────────────────────────────── //
        section(ui, "ORBITAL");

        let elems = self.system.orbital_elements().get(idx).and_then(|e| *e);

        if let Some(el) = elems {
            let u_dist = &self.physics_cfg.dist_label;
            let u_time = &self.physics_cfg.time_label;

            // Orbit type badge
            let (type_str, type_col) = match el.orbit_type {
                OrbitType::Elliptical => ("elliptical", SUCCESS),
                OrbitType::Parabolic => ("parabolic", ACCENT),
                OrbitType::Hyperbolic => ("hyperbolic", DANGER),
            };
            ui.horizontal(|ui| {
                ui.label(RichText::new("type").size(11.0).color(TEXT_SEC));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(RichText::new(type_str).size(11.0).color(type_col).strong());
                });
            });

            metric(ui, "primary", &format!("body #{}", el.primary_idx), TEXT_DIM);

            // Semi-major axis
            let a_str =
                if el.a.is_finite() { format!("{:.4e} {}", el.a, u_dist) } else { "∞".into() };
            metric(ui, "a  (semi-major)", &a_str, TEXT_PRI);

            // Eccentricity
            metric(ui, "e  (eccentric.)", &format!("{:.6}", el.e), TEXT_PRI);

            // Period
            let t_str = if el.period.is_finite() {
                format!("{:.4e} {}", el.period, u_time)
            } else {
                "∞".into()
            };
            metric(ui, "T  (period)", &t_str, TEXT_PRI);

            // Specific angular momentum
            metric(ui, "h  (ang. mom.)", &format!("{:.4e}", el.h_vec.length()), TEXT_DIM);

            // Specific orbital energy
            metric(ui, "ε  (orb. energy)", &format!("{:.4e}", el.energy), TEXT_DIM);

            // Argument of periapsis (degrees) — only meaningful for eccentric orbits
            if el.e > 1e-4 {
                metric(ui, "ω  (peri. arg.)", &format!("{:.2}°", el.omega.to_degrees()), TEXT_DIM);
            }
        } else {
            ui.label(
                RichText::new("  N/A — system has < 2 bodies").size(10.0).color(TEXT_DIM).italics(),
            );
        }

        section(ui, "COLOR");

        let [r, g, b_] = body.color;
        let mut color_rgb: [f32; 3] = [r as f32 / 255.0, g as f32 / 255.0, b_ as f32 / 255.0];
        let color_changed = ui.color_edit_button_rgb(&mut color_rgb).changed();
        let is_custom = body.color != body.material.props().base_color;
        ui.label(
            RichText::new(if is_custom { "custom" } else { "auto (material)" })
                .size(9.5)
                .color(TEXT_DIM),
        );
        let reset_color = is_custom && secondary_btn(ui, "Reset color");

        if color_changed {
            let old = self.system.bodies()[idx];
            let mut b = old;
            b.color = [
                (color_rgb[0] * 255.0) as u8,
                (color_rgb[1] * 255.0) as u8,
                (color_rgb[2] * 255.0) as u8,
            ];
            self.push_undo(UndoRecord::EditedBody {
                idx,
                old_body: old,
                old_name: self.system.name(idx).to_owned(),
            });
            self.system.update_body(idx, b);
        }
        if reset_color {
            let old = self.system.bodies()[idx];
            let mut b = old;
            b.color = b.material.props().base_color;
            self.push_undo(UndoRecord::EditedBody {
                idx,
                old_body: old,
                old_name: self.system.name(idx).to_owned(),
            });
            self.system.update_body(idx, b);
        }

        section(ui, "EDIT");

        if self.selection_form.is_none() {
            let name = self.system.name(idx).to_owned();
            self.selection_form = Some(SelectionForm::from_body(&body, &name));
        }

        let (apply, delete, error) = {
            let form = self.selection_form.as_mut().unwrap();
            field(ui, "name", &mut form.name);
            field(ui, "x", &mut form.x);
            field(ui, "y", &mut form.y);
            field(ui, "vx", &mut form.vx);
            field(ui, "vy", &mut form.vy);
            field(ui, "mass", &mut form.mass);
            field(ui, "density", &mut form.density);

            let radius_preview = {
                let m: Option<f64> = form.mass.parse().ok().filter(|&v| v > 0.0);
                let d: Option<f64> = form.density.parse().ok().filter(|&v| v > 0.0);
                match (m, d) {
                    (Some(m), Some(d)) => format!("{:.5}", radius_from_density_mass(d, m)),
                    _ => "—".into(),
                }
            };
            ui.add_space(2.0);
            ui.horizontal(|ui| {
                ui.label(RichText::new("r_phys →").size(10.0).color(TEXT_DIM));
                ui.label(RichText::new(&radius_preview).monospace().size(10.0).color(TEXT_SEC));
            });

            let err = form.error.clone();
            ui.add_space(6.0);
            let apply = primary_btn(ui, "Apply");
            ui.add_space(3.0);
            let delete = secondary_btn(ui, "Delete body");
            (apply, delete, err)
        };

        if let Some(err) = error {
            ui.add_space(4.0);
            ui.label(RichText::new(err).size(10.0).color(DANGER));
        }

        if apply {
            let parsed = (|| -> Option<Body> {
                let f = self.selection_form.as_ref().unwrap();
                let mass = f.mass.parse::<f64>().ok().filter(|&v| v > 0.0)?;
                let density = f.density.parse::<f64>().ok().filter(|&v| v > 0.0)?;
                let mut b = Body::rocky(mass)
                    .at(f.x.parse().ok()?, f.y.parse().ok()?)
                    .with_velocity(f.vx.parse().ok()?, f.vy.parse().ok()?);
                b.density = density;
                b.sync_physical_properties();
                b.softening = b.softening.max(b.physical_radius * 2.0);
                Some(b)
            })();
            match parsed {
                Some(b) => {
                    let old_body = self.system.bodies()[idx];
                    let old_name = self.system.name(idx).to_owned();
                    self.push_undo(UndoRecord::EditedBody { idx, old_body, old_name });
                    self.system.update_body(idx, b);
                    let name = self.selection_form.as_ref().unwrap().name.clone();
                    if !name.is_empty() {
                        self.system.set_name(idx, name);
                    }
                    self.selection_form.as_mut().unwrap().error = None;
                },
                None => {
                    self.selection_form.as_mut().unwrap().error = Some("invalid values".into());
                },
            }
        }

        if delete {
            let body = self.system.bodies()[idx];
            let name = self.system.name(idx).to_owned();
            let old_last = self.system.bodies().len().saturating_sub(1);
            self.push_undo(UndoRecord::RemovedBody { body, name });
            self.system.remove_body(idx);
            self.pins_after_swap_remove(idx, old_last);
            self.selected_body = None;
            self.follow_selected_body = false;
            self.selection_form = None;
        }
    }
}
