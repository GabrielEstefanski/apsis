use crate::app::config::PhysicsConfig;
use crate::app::templates::{
    spawn_cluster, spawn_ring, template_bodies, TemplateCategory, TEMPLATE_CATALOG,
};
use crate::app::theme::{
    ACCENT, ACCENT_DIM, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC, field,
    fix4, metric, primary_btn, sci, secondary_btn, section, template_btn,
};
use crate::app::ui::{BodyForm, PanelTab, SelectionForm, SimulationApp, SpawnTab};
use crate::domain::body::{Body, radius_from_density_mass};
use crate::physics::gravity::G;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    // ── TOP TOOLBAR ───────────────────────────────────────────────────────── //

    pub(super) fn draw_toolbar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("toolbar")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(12.0, 5.0)),
            )
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new("GRAVITY SIM")
                            .size(12.0)
                            .color(TEXT_PRI)
                            .strong(),
                    );

                    ui.separator();

                    // Play / Pause
                    let (lbl, col) = if self.paused {
                        ("▶  Run", SUCCESS)
                    } else {
                        ("⏸  Pause", ACCENT)
                    };
                    if ui
                        .add(
                            egui::Button::new(RichText::new(lbl).size(11.0).color(col))
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(1.0, col))
                                .min_size(egui::vec2(72.0, 22.0)),
                        )
                        .clicked()
                    {
                        self.paused = !self.paused;
                    }

                    ui.separator();

                    // Simulation parameters
                    let dt_speed = self.proposed_dt * 0.05;

                    ui.label(RichText::new("dt").size(10.0).color(TEXT_SEC));
                    ui.add(
                        egui::DragValue::new(&mut self.proposed_dt)
                            .speed(dt_speed)
                            .clamp_range(1e-5..=0.5)
                            .min_decimals(3)
                            .max_decimals(5),
                    );

                    ui.label(RichText::new("zoom").size(10.0).color(TEXT_SEC));
                    ui.add(
                        egui::DragValue::new(&mut self.scale)
                            .speed(0.5)
                            .clamp_range(1.0..=500.0f32)
                            .max_decimals(1),
                    );

                    ui.separator();

                    // Display toggles
                    ui.checkbox(
                        &mut self.show_grid,
                        RichText::new("grid").size(11.0).color(TEXT_SEC),
                    );
                    ui.checkbox(
                        &mut self.show_trails,
                        RichText::new("trails").size(11.0).color(TEXT_SEC),
                    );
                    ui.checkbox(
                        &mut self.show_vectors,
                        RichText::new("vel").size(11.0).color(TEXT_SEC),
                    );
                    ui.checkbox(
                        &mut self.show_force_vectors,
                        RichText::new("force").size(11.0).color(TEXT_SEC),
                    );
                    ui.checkbox(
                        &mut self.show_impact_normals,
                        RichText::new("nrm").size(11.0).color(TEXT_SEC),
                    );

                    ui.separator();

                    // Right-aligned system buttons
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        let n = self.system.bodies().len();
                        ui.label(
                            RichText::new(format!("{} bodies", n))
                                .size(10.0)
                                .color(TEXT_DIM),
                        );
                        ui.separator();
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("Clear").size(10.0).color(DANGER),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, DANGER))
                                .min_size(egui::vec2(46.0, 20.0)),
                            )
                            .clicked()
                        {
                            self.system.load_bodies(vec![]);
                            self.paused = true;
                        }
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("Reset E₀").size(10.0).color(TEXT_SEC),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, BORDER))
                                .min_size(egui::vec2(60.0, 20.0)),
                            )
                            .clicked()
                        {
                            self.system.reset_energy_baseline();
                        }
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("Zero COM").size(10.0).color(TEXT_SEC),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, BORDER))
                                .min_size(egui::vec2(60.0, 20.0)),
                            )
                            .clicked()
                        {
                            self.system.zero_com_velocity();
                            self.system.reset_energy_baseline();
                        }
                    });
                });
            });
    }

    // ── LEFT PANEL ────────────────────────────────────────────────────────── //

    pub(super) fn draw_panel(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("controls")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(12.0, 10.0)),
            )
            .min_width(220.0)
            .max_width(220.0)
            .show(ctx, |ui| {
                ui.set_width(196.0);

                // ── Metrics (always visible, no scroll) ──────────────── //
                self.panel_metrics_compact(ui);

                // ── Time speed (always visible) ───────────────────────── //
                self.panel_time_speed(ui);

                // ── Tab bar ──────────────────────────────────────────── //
                self.panel_tab_bar(ui);

                ui.add_space(4.0);

                // ── Tab content (scrollable) ─────────────────────────── //
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.set_width(196.0);
                    match self.panel_tab {
                        PanelTab::Add => self.panel_tab_add(ui),
                        PanelTab::Templates => self.panel_tab_templates(ui),
                        PanelTab::Config => self.panel_tab_config(ui),
                    }
                });
            });
    }

    // ── METRICS (pinned) ──────────────────────────────────────────────────── //

    fn panel_metrics_compact(&self, ui: &mut egui::Ui) {
        let m = self.system.metrics();

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
                // Peak drift: max of dE and dLz, shows worst-ever violation
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
                if m.total_dust_mass > 0.0 {
                    ui.label(
                        RichText::new(format!("dust {:.3e}", m.total_dust_mass))
                            .size(9.5)
                            .color(TEXT_DIM),
                    );
                }
            });
        }
    }

    // ── TAB BAR ───────────────────────────────────────────────────────────── //

    // ── TIME SPEED ────────────────────────────────────────────────────────── //

    fn panel_time_speed(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        // Header row: label + current multiplier
        ui.horizontal(|ui| {
            ui.label(RichText::new("TIME SPEED").size(9.5).color(TEXT_DIM).strong());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let col = if self.steps_per_frame > 1 { ACCENT } else { TEXT_DIM };
                ui.label(
                    RichText::new(format!("×{}", self.steps_per_frame))
                        .monospace()
                        .size(10.0)
                        .color(col),
                );
            });
        });

        ui.add_space(2.0);
        ui.add(
            egui::Slider::new(&mut self.steps_per_frame, 1..=1000u32)
                .logarithmic(true)
                .show_value(false),
        );
    }

    // ── TAB BAR ───────────────────────────────────────────────────────────── //

    fn panel_tab_bar(&mut self, ui: &mut egui::Ui) {
        const TABS: &[(PanelTab, &str)] = &[
            (PanelTab::Add, "Add"),
            (PanelTab::Templates, "Library"),
            (PanelTab::Config, "Config"),
        ];
        let w = (196.0 - 4.0 * (TABS.len() as f32 - 1.0)) / TABS.len() as f32;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 4.0;
            for (tab, label) in TABS {
                let active = self.panel_tab == *tab;
                let col = if active { TEXT_PRI } else { TEXT_SEC };
                let fill = if active { ACCENT_DIM } else { Color32::TRANSPARENT };
                if ui
                    .add(
                        egui::Button::new(RichText::new(*label).size(10.5).color(col))
                            .fill(fill)
                            .stroke(Stroke::new(0.5, BORDER))
                            .min_size(egui::vec2(w, 22.0)),
                    )
                    .clicked()
                {
                    self.panel_tab = *tab;
                }
            }
        });
    }

    // ── TAB: ADD ──────────────────────────────────────────────────────────── //

    fn panel_tab_add(&mut self, ui: &mut egui::Ui) {
        ui.add_space(6.0);
        ui.horizontal(|ui| {
            let cur = self.spawn_tab;
            let sub_btn = |ui: &mut egui::Ui, tab: SpawnTab, label: &str| -> bool {
                let active = cur == tab;
                let col = if active { TEXT_PRI } else { TEXT_SEC };
                let fill = if active { ACCENT_DIM } else { Color32::TRANSPARENT };
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
        if primary_btn(ui, &format!("+ Spawn {} cluster", self.spawn_cluster_count)) {
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

    // ── TAB: TEMPLATES ────────────────────────────────────────────────────── //

    fn panel_tab_templates(&mut self, ui: &mut egui::Ui) {
        ui.add_space(2.0);
        ui.label(
            RichText::new("Load a preset — clears current bodies.")
                .size(9.5)
                .color(TEXT_DIM)
                .italics(),
        );

        const CATEGORIES: &[TemplateCategory] = &[
            TemplateCategory::Bodies,
            TemplateCategory::Formations,
            TemplateCategory::Collisions,
        ];

        for &cat in CATEGORIES {
            let entries: Vec<_> = TEMPLATE_CATALOG.iter().filter(|e| e.category == cat).collect();
            if entries.is_empty() {
                continue;
            }

            section(ui, cat.label());

            egui::Grid::new(cat.grid_id())
                .num_columns(2)
                .spacing([4.0, 4.0])
                .show(ui, |ui| {
                    for (i, entry) in entries.iter().enumerate() {
                        if template_btn(ui, entry.label) {
                            let bodies = template_bodies(entry.key);
                            self.system.load_bodies(bodies);
                            self.paused = true;
                            self.offset = egui::Vec2::ZERO;
                        }
                        if i % 2 == 1 {
                            ui.end_row();
                        }
                    }
                    if entries.len() % 2 == 1 {
                        ui.end_row();
                    }
                });
        }
    }

    // ── TAB: CONFIG ───────────────────────────────────────────────────────── //

    fn panel_tab_config(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);

        ui.label(RichText::new("G multiplier").size(10.0).color(TEXT_SEC));
        let g_resp = ui.add(
            egui::Slider::new(&mut self.physics_cfg.g_factor, 0.01..=10.0)
                .logarithmic(true)
                .show_value(true)
                .text(""),
        );
        if g_resp.changed() {
            self.system.set_g_factor(self.physics_cfg.g_factor);
            self.system.reset_energy_baseline();
        }
        ui.label(
            RichText::new(format!("G_eff = {:.4}", self.physics_cfg.g_factor))
                .size(9.0)
                .color(TEXT_DIM),
        );

        ui.add_space(8.0);
        ui.label(
            RichText::new("restitution  (0 = merge)")
                .size(10.0)
                .color(TEXT_SEC),
        );
        if ui
            .add(
                egui::Slider::new(&mut self.collision_cor, 0.0..=1.0)
                    .show_value(true)
                    .text(""),
            )
            .changed()
        {
            self.system.set_cor(self.collision_cor);
        }

        ui.add_space(8.0);
        ui.label(
            RichText::new("UNIT LABELS")
                .size(9.5)
                .color(TEXT_DIM)
                .strong(),
        );
        ui.add_space(3.0);
        egui::Grid::new("unit_labels")
            .num_columns(2)
            .spacing([4.0, 2.0])
            .show(ui, |ui| {
                ui.label(RichText::new("mass").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.mass_label)
                        .desired_width(48.0),
                );
                ui.end_row();
                ui.label(RichText::new("dist").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.dist_label)
                        .desired_width(48.0),
                );
                ui.end_row();
                ui.label(RichText::new("time").size(9.5).color(TEXT_DIM));
                ui.add(
                    egui::TextEdit::singleline(&mut self.physics_cfg.time_label)
                        .desired_width(48.0),
                );
                ui.end_row();
            });

        ui.add_space(8.0);
        if secondary_btn(ui, "Reset to defaults") {
            self.physics_cfg = PhysicsConfig::default();
            self.system.set_g_factor(1.0);
            self.system.set_cor(0.0);
            self.collision_cor = 0.0;
            self.system.reset_energy_baseline();
        }
    }

    // ── INSPECTOR (right panel) ────────────────────────────────────────────── //

    pub(super) fn draw_inspector(&mut self, ctx: &egui::Context) {
        let idx = match self.selected_body {
            Some(i) => i,
            None => return,
        };

        if idx >= self.system.bodies().len() {
            self.selected_body = None;
            self.selection_form = None;
            return;
        }

        egui::SidePanel::right("inspector")
            .frame(
                egui::Frame::none()
                    .fill(crate::app::theme::PANEL_BG)
                    .inner_margin(egui::Margin::symmetric(14.0, 14.0)),
            )
            .min_width(200.0)
            .max_width(200.0)
            .show(ctx, |ui| {
                ui.set_width(172.0);

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
                metric(ui, "radius", &format!("{:.5}", body.radius), TEXT_DIM);
                metric(ui, "density", &format!("{:.4e}", body.density), TEXT_DIM);

                // ── Color picker ─────────────────────────────────────── //
                section(ui, "COLOR");

                let [r, g, b_] = body.color;
                let mut color_rgb: [f32; 3] = [
                    r as f32 / 255.0,
                    g as f32 / 255.0,
                    b_ as f32 / 255.0,
                ];
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
                        ui.label(RichText::new("r →").size(10.0).color(TEXT_DIM));
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
                        let radius = radius_from_density_mass(density, mass);
                        let mut b = Body::new(
                            f.x.parse().ok()?,
                            f.y.parse().ok()?,
                            f.vx.parse().ok()?,
                            f.vy.parse().ok()?,
                            mass,
                            crate::domain::materials::Material::Rocky,
                        );
                        b.radius = radius;
                        b.density = density;
                        b.softening = b.softening.max(radius * 2.0);
                        b.moment_inertia =
                            crate::domain::body::default_moment_inertia(mass, radius);
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
            });
    }
}
