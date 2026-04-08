use crate::app::theme::{ACCENT_DIM, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_SEC};
use crate::app::theme::{field, primary_btn};
use crate::app::ui::{BodyForm, SimulationApp, SpawnTab};
use crate::domain::body::{Body, radius_from_density_mass};
use crate::domain::materials::Material;
use crate::physics::gravity::G;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn panel_tab_add(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        self.spawn_sub_tab_bar(ui);
        ui.add_space(6.0);
        match self.spawn_tab {
            SpawnTab::Single => self.panel_add_single(ui),
            SpawnTab::Ring => self.panel_add_ring(ui),
            SpawnTab::Cluster => self.panel_add_cluster(ui),
        }
    }

    fn spawn_sub_tab_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let cur = self.spawn_tab;
            let sub_btn = |ui: &mut egui::Ui, tab: SpawnTab, label: &str| -> bool {
                let active = cur == tab;
                let col = if active {
                    crate::app::theme::TEXT_PRI
                } else {
                    TEXT_SEC
                };
                let fill = if active {
                    ACCENT_DIM
                } else {
                    Color32::TRANSPARENT
                };
                ui.add(
                    egui::Button::new(RichText::new(label).size(10.0).color(col))
                        .fill(fill)
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(58.0, 20.0)),
                )
                .clicked()
            };
            if sub_btn(ui, SpawnTab::Single, "single") {
                self.spawn_tab = SpawnTab::Single;
            }
            if sub_btn(ui, SpawnTab::Ring, "ring") {
                self.spawn_tab = SpawnTab::Ring;
            }
            if sub_btn(ui, SpawnTab::Cluster, "cluster") {
                self.spawn_tab = SpawnTab::Cluster;
            }
        });
    }

    fn panel_add_single(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let (col, lbl) = if self.place_mode {
                (SUCCESS, "● place on")
            } else {
                (TEXT_SEC, "○ place off")
            };
            let btn = egui::Button::new(RichText::new(lbl).size(10.5).color(col))
                .fill(Color32::TRANSPARENT)
                .stroke(Stroke::new(
                    0.5,
                    if self.place_mode { SUCCESS } else { BORDER },
                ))
                .min_size(egui::vec2(ui.available_width(), 20.0));
            if ui.add(btn).clicked() {
                self.place_mode = !self.place_mode;
            }
        });

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.label(RichText::new("mass").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.place_mass)
                        .speed(0.1)
                        .clamp_range(1e-6..=1e6),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("density").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.place_density)
                        .speed(0.01)
                        .clamp_range(1e-6..=1e6),
                );
            });
        });
        {
            let r = radius_from_density_mass(self.place_density, self.place_mass);
            ui.label(
                RichText::new(format!("r → {:.4}", r))
                    .size(9.5)
                    .color(TEXT_DIM),
            );
        }

        ui.add_space(6.0);
        ui.label(RichText::new("or enter coords").size(9.5).color(TEXT_DIM));
        ui.add_space(3.0);

        field(ui, "x", &mut self.form.x);
        field(ui, "y", &mut self.form.y);
        field(ui, "vx", &mut self.form.vx);
        field(ui, "vy", &mut self.form.vy);
        field(ui, "mass", &mut self.form.mass);
        field(ui, "dens", &mut self.form.density);

        {
            let r_preview = {
                let m: Option<f64> = self.form.mass.parse().ok().filter(|&v| v > 0.0);
                let d: Option<f64> = self.form.density.parse().ok().filter(|&v| v > 0.0);
                match (m, d) {
                    (Some(m), Some(d)) => format!("{:.5}", radius_from_density_mass(d, m)),
                    _ => "—".into(),
                }
            };
            ui.horizontal(|ui| {
                ui.label(RichText::new("r →").size(9.5).color(TEXT_DIM));
                ui.label(
                    RichText::new(&r_preview)
                        .monospace()
                        .size(9.5)
                        .color(TEXT_SEC),
                );
            });
        }

        ui.add_space(6.0);
        if primary_btn(ui, "+ Add body") {
            match self.form.try_build() {
                Some(body) => {
                    self.system.add_body(body);
                    self.form = BodyForm::default();
                    self.form_error = None;
                }
                None => self.form_error = Some("invalid values".into()),
            }
        }
        if let Some(err) = &self.form_error {
            ui.add_space(4.0);
            ui.label(RichText::new(err).size(10.0).color(DANGER));
        }
    }

    fn panel_add_ring(&mut self, ui: &mut egui::Ui) {
        let total_m = self.system.total_mass();
        let auto_vel = if total_m > 0.0 && self.spawn_ring_radius > 0.0 {
            (G * total_m / self.spawn_ring_radius).sqrt() * self.spawn_ring_vel_scale
        } else {
            0.0
        };

        ui.horizontal(|ui| {
            ui.label(RichText::new("radius").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_ring_radius)
                        .speed(0.1)
                        .clamp_range(0.1..=1000.0),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("count").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_ring_count)
                        .speed(1.0)
                        .clamp_range(2..=2000u32),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("mass each").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_ring_mass)
                        .speed(0.001)
                        .clamp_range(1e-6..=1e6),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("vel scale").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_ring_vel_scale)
                        .speed(0.01)
                        .clamp_range(0.0..=5.0),
                );
            });
        });

        ui.add_space(3.0);
        let m = self.system.metrics();
        ui.label(
            RichText::new(format!(
                "v ≈ {:.3}  around COM ({:.1}, {:.1})",
                auto_vel, m.com_x, m.com_y
            ))
            .size(9.0)
            .color(TEXT_DIM),
        );
        ui.add_space(6.0);

        if primary_btn(ui, "+ Add ring") {
            let n = self.spawn_ring_count as usize;
            if n > 0 && self.spawn_ring_radius > 0.0 {
                let total_m = self.system.total_mass();
                let center = self.system.metrics();

                let v = if total_m > 0.0 {
                    (G * total_m / self.spawn_ring_radius).sqrt() * self.spawn_ring_vel_scale
                } else {
                    0.0
                };

                for i in 0..n {
                    let angle = (i as f64 / n as f64) * std::f64::consts::TAU;

                    let x = center.com_x + self.spawn_ring_radius * angle.cos();
                    let y = center.com_y + self.spawn_ring_radius * angle.sin();

                    // tangential velocity
                    let vx = -v * angle.sin();
                    let vy = v * angle.cos();

                    let mut b = Body::new(x, y, vx, vy, self.spawn_ring_mass, Material::Rocky);

                    b.density = self.place_density;
                    b.sync_physical_properties();

                    self.system.add_body(b);
                }
            }
        }
    }

    fn panel_add_cluster(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label(RichText::new("radius").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_cluster_radius)
                        .speed(0.1)
                        .clamp_range(0.1..=1000.0),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("count").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_cluster_count)
                        .speed(1.0)
                        .clamp_range(1..=5000u32),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("mass each").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_cluster_mass)
                        .speed(0.01)
                        .clamp_range(1e-6..=1e6),
                );
            });
        });
        ui.horizontal(|ui| {
            ui.label(RichText::new("vel disp").size(10.0).color(TEXT_SEC));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.add(
                    egui::DragValue::new(&mut self.spawn_cluster_vel_disp)
                        .speed(0.01)
                        .clamp_range(0.0..=100.0),
                );
            });
        });

        ui.add_space(6.0);

        if primary_btn(ui, "+ Add cluster") {
            let n = self.spawn_cluster_count as usize;
            if n > 0 {
                let center = self.system.metrics();

                for _ in 0..n {
                    let r = self.spawn_cluster_radius * rand::random::<f64>().sqrt();
                    let theta = rand::random::<f64>() * std::f64::consts::TAU;

                    let x = center.com_x + r * theta.cos();
                    let y = center.com_y + r * theta.sin();

                    let vx = (rand::random::<f64>() - 0.5) * self.spawn_cluster_vel_disp;
                    let vy = (rand::random::<f64>() - 0.5) * self.spawn_cluster_vel_disp;

                    let mut b = Body::new(x, y, vx, vy, self.spawn_cluster_mass, Material::Rocky);

                    b.density = self.place_density;
                    b.sync_physical_properties();

                    self.system.add_body(b);
                }
            }
        }
    }
}
