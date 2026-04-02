use crate::app::templates::{spawn_cluster, spawn_ring, template_bodies};
use crate::app::theme::{
    field, fix4, metric, primary_btn, sci, secondary_btn, section, template_btn,
    ACCENT, ACCENT_DIM, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC,
};
use crate::app::ui::{BodyForm, SimulationApp, SpawnTab};
use crate::physics::gravity::G;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    pub(super) fn draw_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("controls")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14.0, 14.0)),
            )
            .min_width(220.0)
            .max_width(220.0)
            .show(ctx, |ui| {
                ui.set_width(192.0);
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.set_width(192.0);
                    self.panel_header(ui);
                    self.panel_simulation(ui);
                    self.panel_templates(ui);
                    self.panel_add(ui);
                    self.panel_metrics(ui);
                    self.panel_system(ui);
                });
            });
    }

    fn panel_header(&self, ui: &mut egui::Ui) {
        ui.label(RichText::new("GRAVITY SIM").size(14.0).color(TEXT_PRI).strong());
        ui.label(
            RichText::new(format!(
                "{} bodies  ·  {}",
                self.system.bodies().len(),
                if self.paused { "paused" } else { "running" }
            ))
            .size(10.0)
            .color(TEXT_DIM),
        );
    }

    fn panel_simulation(&mut self, ui: &mut egui::Ui) {
        section(ui, "SIMULATION");

        ui.horizontal(|ui| {
            let (lbl, col) = if self.paused {
                ("▶  RESUME", SUCCESS)
            } else {
                ("⏸  PAUSE", ACCENT)
            };
            let btn = egui::Button::new(RichText::new(lbl).size(12.0).color(col))
                .fill(Color32::TRANSPARENT)
                .stroke(Stroke::new(1.0, col))
                .min_size(egui::vec2(ui.available_width(), 28.0));
            if ui.add(btn).clicked() {
                self.paused = !self.paused;
            }
        });

        ui.add_space(6.0);
        ui.label(RichText::new("proposed dt").size(10.0).color(TEXT_SEC));
        ui.add(
            egui::Slider::new(&mut self.proposed_dt, 1e-5..=5e-3)
                .logarithmic(true)
                .show_value(true)
                .text(""),
        );

        ui.label(RichText::new("steps / frame").size(10.0).color(TEXT_SEC));
        ui.add(
            egui::Slider::new(&mut self.steps_per_frame, 1..=50)
                .show_value(true)
                .text(""),
        );

        ui.label(RichText::new("zoom").size(10.0).color(TEXT_SEC));
        ui.add(
            egui::Slider::new(&mut self.scale, 1.0..=150.0)
                .logarithmic(true)
                .show_value(true)
                .text(""),
        );

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_grid, RichText::new("grid").size(11.0).color(TEXT_SEC));
            ui.checkbox(&mut self.show_trails, RichText::new("trails").size(11.0).color(TEXT_SEC));
            ui.checkbox(&mut self.show_vectors, RichText::new("vel").size(11.0).color(TEXT_SEC));
        });
    }

    fn panel_templates(&mut self, ui: &mut egui::Ui) {
        section(ui, "TEMPLATES");

        let templates = [
            ("binary", "Binária"),
            ("figure8", "Figure-8"),
            ("solar", "Solar"),
            ("pythagorean", "Pitagórico"),
            ("belt", "Cinturão"),
            ("galaxies", "Galáxias"),
        ];

        egui::Grid::new("tpl_grid")
            .num_columns(2)
            .spacing([4.0, 4.0])
            .show(ui, |ui| {
                for (i, (key, label)) in templates.iter().enumerate() {
                    if template_btn(ui, label) {
                        let bodies = template_bodies(key);
                        self.system.load_bodies(bodies);
                        self.paused = true;
                        self.offset = egui::Vec2::ZERO;
                    }
                    if i % 2 == 1 {
                        ui.end_row();
                    }
                }
            });
    }

    fn panel_add(&mut self, ui: &mut egui::Ui) {
        section(ui, "ADD");

        ui.horizontal(|ui| {
            let cur = self.spawn_tab;
            let tab_btn = |ui: &mut egui::Ui, tab: SpawnTab, label: &str| -> bool {
                let active = cur == tab;
                let col = if active { TEXT_PRI } else { TEXT_SEC };
                let fill = if active { ACCENT_DIM } else { Color32::TRANSPARENT };
                ui.add(
                    egui::Button::new(RichText::new(label).size(10.5).color(col))
                        .fill(fill)
                        .stroke(Stroke::new(0.5, BORDER))
                        .min_size(egui::vec2(58.0, 20.0)),
                )
                .clicked()
            };
            if tab_btn(ui, SpawnTab::Single, "single") {
                self.spawn_tab = SpawnTab::Single;
            }
            if tab_btn(ui, SpawnTab::Ring, "ring") {
                self.spawn_tab = SpawnTab::Ring;
            }
            if tab_btn(ui, SpawnTab::Cluster, "cluster") {
                self.spawn_tab = SpawnTab::Cluster;
            }
        });

        ui.add_space(6.0);

        match self.spawn_tab {
            SpawnTab::Single => self.panel_add_single(ui),
            SpawnTab::Ring => self.panel_add_ring(ui),
            SpawnTab::Cluster => self.panel_add_cluster(ui),
        }
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
        ui.label(RichText::new("mass").size(10.0).color(TEXT_SEC));
        ui.add(
            egui::DragValue::new(&mut self.place_mass)
                .speed(0.1)
                .clamp_range(1e-6..=1e6),
        );

        ui.add_space(6.0);
        ui.label(RichText::new("or enter coords").size(9.5).color(TEXT_DIM));
        ui.add_space(3.0);

        field(ui, "x", &mut self.form.x);
        field(ui, "y", &mut self.form.y);
        field(ui, "vx", &mut self.form.vx);
        field(ui, "vy", &mut self.form.vy);
        field(ui, "mass", &mut self.form.mass);

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
        if primary_btn(ui, &format!("+ Spawn {} ring", self.spawn_ring_count)) {
            let m = self.system.metrics();
            let bodies = spawn_ring(
                m.com_x,
                m.com_y,
                self.spawn_ring_radius,
                self.spawn_ring_count as usize,
                self.spawn_ring_mass,
                auto_vel,
            );
            for b in bodies {
                self.system.add_body(b);
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
        if primary_btn(
            ui,
            &format!("+ Spawn {} cluster", self.spawn_cluster_count),
        ) {
            let m = self.system.metrics();
            let bodies = spawn_cluster(
                m.com_x,
                m.com_y,
                self.spawn_cluster_radius,
                self.spawn_cluster_count as usize,
                self.spawn_cluster_mass,
                self.spawn_cluster_vel_disp,
            );
            for b in bodies {
                self.system.add_body(b);
            }
        }
    }

    fn panel_metrics(&self, ui: &mut egui::Ui) {
        section(ui, "METRICS");

        let m = self.system.metrics();
        let de_color = if m.rel_energy_error.abs() < 1e-8 {
            SUCCESS
        } else if m.rel_energy_error.abs() < 1e-5 {
            ACCENT
        } else {
            DANGER
        };

        metric(ui, "E", &fix4(m.total_energy), TEXT_PRI);
        metric(ui, "dE/E₀", &sci(m.rel_energy_error), de_color);
        metric(ui, "Lz", &fix4(m.angular_momentum_z), TEXT_PRI);
        metric(ui, "θ", &format!("{:.4}", m.theta), TEXT_SEC);
        metric(ui, "dt", &format!("{:.2e}", m.dt), TEXT_SEC);
        metric(ui, "K", &fix4(m.kinetic), TEXT_DIM);
        metric(ui, "U", &fix4(m.potential), TEXT_DIM);
    }

    fn panel_system(&mut self, ui: &mut egui::Ui) {
        section(ui, "SYSTEM");

        if secondary_btn(ui, "Zero COM velocity") {
            self.system.zero_com_velocity();
            self.system.reset_energy_baseline();
        }
        ui.add_space(4.0);
        if secondary_btn(ui, "Reset energy baseline") {
            self.system.reset_energy_baseline();
        }
        ui.add_space(4.0);
        if secondary_btn(ui, "Clear all") {
            self.system.load_bodies(vec![]);
            self.paused = true;
        }
    }
}
