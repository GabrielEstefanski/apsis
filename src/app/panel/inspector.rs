use crate::app::theme::{DANGER, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::theme::{field, metric, primary_btn, secondary_btn, section};
use crate::app::ui::{SelectionForm, SimulationApp};
use crate::domain::body::{Body, default_moment_inertia, radius_from_density_mass};
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn inspector_content(&mut self, ui: &mut egui::Ui, idx: usize) {
        let body = self.system.bodies()[idx];

        ui.horizontal(|ui| {
            let [cr, cg, cb] = body.color;
            let col = egui::Color32::from_rgb(cr, cg, cb);
            let (dot_rect, _) =
                ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
            ui.painter().circle_filled(dot_rect.center(), 5.0, col);
            ui.label(
                RichText::new(format!("body #{}", idx))
                    .size(12.0)
                    .color(TEXT_PRI)
                    .strong(),
            );
        });

        ui.add_space(10.0);

        section(ui, "LIVE");
        metric(ui, "x", &format!("{:.5}", body.x), TEXT_DIM);
        metric(ui, "y", &format!("{:.5}", body.y), TEXT_DIM);
        metric(ui, "vx", &format!("{:.5}", body.vx), TEXT_DIM);
        metric(ui, "vy", &format!("{:.5}", body.vy), TEXT_DIM);
        metric(ui, "mass", &format!("{:.5}", body.mass), TEXT_DIM);
        metric(ui, "r_phys", &format!("{:.5}", body.physical_radius), TEXT_DIM);
        metric(ui, "r_coll", &format!("{:.5}", body.radius), TEXT_DIM);
        metric(ui, "soft", &format!("{:.5}", body.softening), TEXT_DIM);
        metric(ui, "density", &format!("{:.4e}", body.density), TEXT_DIM);

        section(ui, "COLOR");

        let [r, g, b_] = body.color;
        let mut color_rgb: [f32; 3] =
            [r as f32 / 255.0, g as f32 / 255.0, b_ as f32 / 255.0];
        let color_changed = ui.color_edit_button_rgb(&mut color_rgb).changed();
        let is_custom = body.color != body.material.props().base_color;
        ui.label(
            RichText::new(if is_custom { "custom" } else { "auto (material)" })
                .size(9.5)
                .color(TEXT_DIM),
        );
        let reset_color = is_custom && secondary_btn(ui, "Reset color");

        if color_changed {
            let mut b = self.system.bodies()[idx];
            b.color = [
                (color_rgb[0] * 255.0) as u8,
                (color_rgb[1] * 255.0) as u8,
                (color_rgb[2] * 255.0) as u8,
            ];
            self.system.update_body(idx, b);
        }
        if reset_color {
            let mut b = self.system.bodies()[idx];
            b.color = b.material.props().base_color;
            self.system.update_body(idx, b);
        }

        section(ui, "EDIT");

        if self.selection_form.is_none() {
            self.selection_form = Some(SelectionForm::from_body(&body));
        }

        let (apply, delete, error) = {
            let form = self.selection_form.as_mut().unwrap();
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
                ui.label(
                    RichText::new(&radius_preview)
                        .monospace()
                        .size(10.0)
                        .color(TEXT_SEC),
                );
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
                let mut b = Body::new(
                    f.x.parse().ok()?,
                    f.y.parse().ok()?,
                    f.vx.parse().ok()?,
                    f.vy.parse().ok()?,
                    mass,
                    crate::domain::materials::Material::Rocky,
                );
                b.density = density;
                b.sync_physical_properties();
                b.radius = b.physical_radius;
                b.softening = b.softening.max(b.physical_radius * 2.0);
                b.moment_inertia = default_moment_inertia(mass, b.physical_radius);
                Some(b)
            })();
            match parsed {
                Some(b) => {
                    self.system.update_body(idx, b);
                    self.selection_form.as_mut().unwrap().error = None;
                }
                None => {
                    self.selection_form.as_mut().unwrap().error =
                        Some("invalid values".into());
                }
            }
        }

        if delete {
            self.system.remove_body(idx);
            self.selected_body = None;
            self.selection_form = None;
        }
    }
}
