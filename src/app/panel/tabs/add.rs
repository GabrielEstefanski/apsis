use crate::app::theme::primary_btn;
use crate::app::theme::{
    ACCENT, ACCENT_DIM, BORDER, DANGER, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC,
};
use crate::app::ui::{BodyForm, SimulationApp, SpawnTab, UndoRecord};
use crate::domain::body::{Body, radius_from_density_mass};
use crate::domain::materials::{Material, density};
use crate::physics::gravity::G;
use eframe::egui::{self, Color32, RichText, Stroke};

// ── Helpers ───────────────────────────────────────────────────────────────────

/// 3-column material picker grid. Returns `true` if the selection changed.
fn material_grid(ui: &mut egui::Ui, selected: &mut Material) -> bool {
    let mut changed = false;
    let cols = 3;

    egui::Grid::new(ui.id().with("mat_grid")).num_columns(cols).spacing([4.0, 4.0]).show(
        ui,
        |ui| {
            for (i, &mat) in Material::ALL.iter().enumerate() {
                let is_sel = *selected == mat;
                let [r, g, b] = mat.props().base_color;
                let dot_col = Color32::from_rgb(r, g, b);
                let text_col = if is_sel { TEXT_PRI } else { TEXT_SEC };
                let fill = if is_sel { ACCENT_DIM } else { Color32::from_rgb(20, 20, 26) };
                let stroke_col = if is_sel { ACCENT } else { BORDER };

                // Inline colored dot + name via LayoutJob — no post-render painting needed
                let mut job = egui::text::LayoutJob::default();
                job.append(
                    "● ",
                    0.0,
                    egui::TextFormat {
                        font_id: egui::FontId::proportional(9.0),
                        color: dot_col,
                        valign: egui::Align::Center,
                        ..Default::default()
                    },
                );
                job.append(
                    mat.display_name(),
                    0.0,
                    egui::TextFormat {
                        font_id: egui::FontId::proportional(9.5),
                        color: text_col,
                        valign: egui::Align::Center,
                        ..Default::default()
                    },
                );

                let btn = ui.add(
                    egui::Button::new(job)
                        .fill(fill)
                        .stroke(Stroke::new(0.5, stroke_col))
                        .min_size(egui::vec2(0.0, 20.0)),
                );

                if btn.clicked() && !is_sel {
                    *selected = mat;
                    changed = true;
                }

                if (i + 1) % cols == 0 {
                    ui.end_row();
                }
            }
            if Material::ALL.len() % cols != 0 {
                ui.end_row();
            }
        },
    );

    changed
}

/// A label + DragValue row, right-aligned value.
fn drag_row(ui: &mut egui::Ui, label: &str, hint: &str, widget: egui::DragValue) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(10.5).color(TEXT_SEC)).on_hover_text(hint);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add(widget);
        });
    });
}

// ── Section header (inline, without `section()` to avoid extra spacing) ───────

fn sub_section(ui: &mut egui::Ui, label: &str) {
    ui.add_space(8.0);
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(9.0).color(TEXT_DIM).strong());
        ui.add_space(4.0);
        let r = ui.available_rect_before_wrap();
        ui.painter().line_segment(
            [egui::pos2(r.left(), r.center().y), egui::pos2(r.right(), r.center().y)],
            Stroke::new(0.5, BORDER),
        );
    });
    ui.add_space(4.0);
}

// ── Top-level ─────────────────────────────────────────────────────────────────

impl SimulationApp {
    pub(super) fn panel_tab_add(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);
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
            let tab_btn = |ui: &mut egui::Ui, tab: SpawnTab, label: &str| -> bool {
                let active = cur == tab;
                ui.add(
                    egui::Button::new(RichText::new(label).size(10.5).color(if active {
                        TEXT_PRI
                    } else {
                        TEXT_SEC
                    }))
                    .fill(if active { ACCENT_DIM } else { Color32::TRANSPARENT })
                    .stroke(Stroke::new(0.5, if active { ACCENT } else { BORDER }))
                    .min_size(egui::vec2(70.0, 22.0)),
                )
                .clicked()
            };
            if tab_btn(ui, SpawnTab::Single, "Single") {
                self.spawn_tab = SpawnTab::Single;
            }
            if tab_btn(ui, SpawnTab::Ring, "Ring") {
                self.spawn_tab = SpawnTab::Ring;
            }
            if tab_btn(ui, SpawnTab::Cluster, "Cluster") {
                self.spawn_tab = SpawnTab::Cluster;
            }
        });
    }

    // ── Single ────────────────────────────────────────────────────────────────

    fn panel_add_single(&mut self, ui: &mut egui::Ui) {
        // ── MATERIAL ─────────────────────────────────────────────────────────
        sub_section(ui, "MATERIAL");
        let mat_changed = material_grid(ui, &mut self.place_material);
        if mat_changed {
            self.place_mass = self.place_material.default_mass();
            self.place_density = density(self.place_material, self.place_mass);
        }

        // ── PROPERTIES ───────────────────────────────────────────────────────
        sub_section(ui, "PROPERTIES");

        let mass_speed = (self.place_mass * 0.05).max(1e-10);
        let prev_mass = self.place_mass;
        drag_row(
            ui,
            "mass",
            "Body mass in simulation units",
            egui::DragValue::new(&mut self.place_mass)
                .speed(mass_speed)
                .range(1e-12..=1e15_f64)
                .max_decimals(6),
        );
        if self.place_mass != prev_mass {
            self.place_density = density(self.place_material, self.place_mass);
        }

        let dens_speed = (self.place_density * 0.02).max(1e-6);
        drag_row(
            ui,
            "density",
            "Bulk density (controls physical radius)",
            egui::DragValue::new(&mut self.place_density)
                .speed(dens_speed)
                .range(1e-6..=1e15_f64)
                .max_decimals(4),
        );

        let r = radius_from_density_mass(self.place_density, self.place_mass);
        ui.horizontal(|ui| {
            ui.label(RichText::new("r").size(10.5).color(TEXT_DIM));
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(format!("{:.5e}", r)).monospace().size(10.5).color(TEXT_SEC),
                );
            });
        });

        // ── CANVAS PLACEMENT ─────────────────────────────────────────────────
        sub_section(ui, "PLACEMENT");

        let (btn_label, btn_col, btn_fill, btn_stroke) = if self.place_mode {
            (
                "● placing — click canvas",
                SUCCESS,
                Color32::from_rgba_unmultiplied(40, 100, 60, 40),
                Stroke::new(1.0, SUCCESS),
            )
        } else {
            ("○ click canvas to place", TEXT_SEC, Color32::TRANSPARENT, Stroke::new(0.5, BORDER))
        };

        if ui
            .add(
                egui::Button::new(RichText::new(btn_label).size(10.5).color(btn_col))
                    .fill(btn_fill)
                    .stroke(btn_stroke)
                    .min_size(egui::vec2(ui.available_width(), 26.0)),
            )
            .clicked()
        {
            self.place_mode = !self.place_mode;
        }

        if self.place_mode {
            ui.add_space(3.0);
            ui.label(
                RichText::new("drag to set initial velocity").size(9.0).color(TEXT_DIM).italics(),
            );
        }

        // ── MANUAL POSITION (collapsible) ────────────────────────────────────
        ui.add_space(6.0);
        egui::CollapsingHeader::new(RichText::new("Manual coordinates").size(10.0).color(TEXT_SEC))
            .id_salt("single_precise")
            .default_open(false)
            .show(ui, |ui| {
                ui.add_space(4.0);

                // 2-column grid: [label input] [label input]
                egui::Grid::new("precise_form").num_columns(4).spacing([6.0, 4.0]).show(ui, |ui| {
                    let lw = 22.0_f32;
                    let vw = ui.available_width() / 2.0 - lw - 10.0;

                    let mut text_field = |ui: &mut egui::Ui, lbl: &str, val: &mut String| {
                        ui.add_sized(
                            egui::vec2(lw, 18.0),
                            egui::Label::new(RichText::new(lbl).size(10.0).color(TEXT_SEC)),
                        );
                        ui.add_sized(
                            egui::vec2(vw.max(40.0), 18.0),
                            egui::TextEdit::singleline(val)
                                .font(egui::TextStyle::Monospace)
                                .desired_width(vw.max(40.0)),
                        );
                    };

                    text_field(ui, "x", &mut self.form.x);
                    text_field(ui, "vx", &mut self.form.vx);
                    ui.end_row();

                    text_field(ui, "y", &mut self.form.y);
                    text_field(ui, "vy", &mut self.form.vy);
                    ui.end_row();

                    text_field(ui, "m", &mut self.form.mass);
                    text_field(ui, "ρ", &mut self.form.density);
                    ui.end_row();
                });

                // Radius preview
                let r_preview = {
                    let m: Option<f64> = self.form.mass.parse().ok().filter(|&v| v > 0.0);
                    let d: Option<f64> = self.form.density.parse().ok().filter(|&v| v > 0.0);
                    match (m, d) {
                        (Some(m), Some(d)) => format!("{:.4e}", radius_from_density_mass(d, m)),
                        _ => "—".into(),
                    }
                };
                ui.horizontal(|ui| {
                    ui.label(RichText::new("r →").size(9.5).color(TEXT_DIM));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(RichText::new(&r_preview).monospace().size(9.5).color(TEXT_SEC));
                    });
                });

                ui.add_space(4.0);
                if primary_btn(ui, "+ Add body at coords") {
                    match self.form.try_build() {
                        Some(body) => {
                            self.push_undo(UndoRecord::AddedBodies(1));
                            self.system.add_body(body);
                            self.form = BodyForm::default();
                            self.form_error = None;
                        },
                        None => self.form_error = Some("invalid values".into()),
                    }
                }
                if let Some(err) = &self.form_error.clone() {
                    ui.add_space(2.0);
                    ui.label(RichText::new(err).size(9.5).color(DANGER));
                }
            });
    }

    // ── Ring ──────────────────────────────────────────────────────────────────

    fn panel_add_ring(&mut self, ui: &mut egui::Ui) {
        // ── MATERIAL ─────────────────────────────────────────────────────────
        sub_section(ui, "MATERIAL");
        let mat_changed = material_grid(ui, &mut self.spawn_ring_material);
        if mat_changed {
            self.spawn_ring_mass = self.spawn_ring_material.default_mass();
        }

        // ── GEOMETRY ─────────────────────────────────────────────────────────
        sub_section(ui, "GEOMETRY");
        let ring_r_speed = (self.spawn_ring_radius * 0.02).max(0.01);
        let ring_m_speed = (self.spawn_ring_mass * 0.05).max(1e-12);
        drag_row(
            ui,
            "radius",
            "Ring radius around the system's centre of mass",
            egui::DragValue::new(&mut self.spawn_ring_radius)
                .speed(ring_r_speed)
                .range(0.01..=1e6_f64)
                .max_decimals(4),
        );
        drag_row(
            ui,
            "count",
            "Number of bodies evenly spaced around the ring",
            egui::DragValue::new(&mut self.spawn_ring_count).speed(1.0).range(2..=10_000u32),
        );
        drag_row(
            ui,
            "mass each",
            "Mass of each individual body",
            egui::DragValue::new(&mut self.spawn_ring_mass)
                .speed(ring_m_speed)
                .range(1e-12..=1e12_f64)
                .max_decimals(6),
        );

        // ── DYNAMICS ─────────────────────────────────────────────────────────
        sub_section(ui, "DYNAMICS");
        drag_row(
            ui,
            "vel scale",
            "Multiplier on circular orbit speed (1.0 = circular, >1 = eccentric outward, <1 = inward)",
            egui::DragValue::new(&mut self.spawn_ring_vel_scale)
                .speed(0.005)
                .range(0.0..=5.0_f64)
                .max_decimals(3),
        );

        // Computed velocity info
        let total_m = self.system.total_mass();
        let v_circ = if total_m > 0.0 && self.spawn_ring_radius > 0.0 {
            (G * total_m / self.spawn_ring_radius).sqrt()
        } else {
            0.0
        };
        let v_ring = v_circ * self.spawn_ring_vel_scale;
        let m = self.system.metrics();

        ui.add_space(4.0);
        egui::Frame::NONE
            .fill(Color32::from_rgba_unmultiplied(20, 20, 28, 200))
            .stroke(Stroke::new(0.5, BORDER))
            .corner_radius(4.0)
            .inner_margin(egui::Margin::symmetric(8, 5))
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                info_row(ui, "v_circ", &format!("{:.4e}", v_circ));
                info_row(
                    ui,
                    "v_ring",
                    &format!("{:.4e}  (×{:.2})", v_ring, self.spawn_ring_vel_scale),
                );
                info_row(ui, "COM", &format!("({:.3}, {:.3})", m.com_x, m.com_y));
                info_row(
                    ui,
                    "total M",
                    &format!("{:.4e}", self.spawn_ring_mass * self.spawn_ring_count as f64),
                );
            });

        ui.add_space(6.0);
        if primary_btn(ui, &format!("+ Add ring  ({} bodies)", self.spawn_ring_count)) {
            let n = self.spawn_ring_count as usize;
            if n > 0 && self.spawn_ring_radius > 0.0 {
                let total_m = self.system.total_mass();
                let center = self.system.metrics();
                let v = if total_m > 0.0 {
                    (G * total_m / self.spawn_ring_radius).sqrt() * self.spawn_ring_vel_scale
                } else {
                    0.0
                };
                let ring_density = density(self.spawn_ring_material, self.spawn_ring_mass);

                self.push_undo(UndoRecord::AddedBodies(n));
                for i in 0..n {
                    let angle = (i as f64 / n as f64) * std::f64::consts::TAU;
                    let x = center.com_x + self.spawn_ring_radius * angle.cos();
                    let y = center.com_y + self.spawn_ring_radius * angle.sin();
                    let vx = -v * angle.sin();
                    let vy = v * angle.cos();
                    let mut b =
                        Body::new(x, y, vx, vy, self.spawn_ring_mass, self.spawn_ring_material);
                    b.density = ring_density;
                    b.sync_physical_properties();
                    self.system.add_body(b);
                }
            }
        }
    }

    // ── Cluster ───────────────────────────────────────────────────────────────

    fn panel_add_cluster(&mut self, ui: &mut egui::Ui) {
        // ── MATERIAL ─────────────────────────────────────────────────────────
        sub_section(ui, "MATERIAL");
        let mat_changed = material_grid(ui, &mut self.spawn_cluster_material);
        if mat_changed {
            self.spawn_cluster_mass = self.spawn_cluster_material.default_mass();
        }

        // ── GEOMETRY ─────────────────────────────────────────────────────────
        sub_section(ui, "GEOMETRY");
        let clust_r_speed = (self.spawn_cluster_radius * 0.02).max(0.01);
        let clust_m_speed = (self.spawn_cluster_mass * 0.05).max(1e-12);
        drag_row(
            ui,
            "radius",
            "Half-mass radius of the cluster",
            egui::DragValue::new(&mut self.spawn_cluster_radius)
                .speed(clust_r_speed)
                .range(0.01..=1e6_f64)
                .max_decimals(4),
        );
        drag_row(
            ui,
            "count",
            "Number of bodies",
            egui::DragValue::new(&mut self.spawn_cluster_count).speed(1.0).range(1..=50_000u32),
        );
        drag_row(
            ui,
            "mass each",
            "Mass of each body",
            egui::DragValue::new(&mut self.spawn_cluster_mass)
                .speed(clust_m_speed)
                .range(1e-12..=1e12_f64)
                .max_decimals(6),
        );

        // ── DYNAMICS ─────────────────────────────────────────────────────────
        sub_section(ui, "DYNAMICS");
        let clust_v_speed = (self.spawn_cluster_vel_disp * 0.05).max(0.001);
        drag_row(
            ui,
            "σ_v",
            "1D velocity dispersion per component. Set to 0 for a cold collapse.",
            egui::DragValue::new(&mut self.spawn_cluster_vel_disp)
                .speed(clust_v_speed)
                .range(0.0..=1e4_f64)
                .max_decimals(4),
        );

        // Virial ratio estimate: Q = K/|W|; for a uniform sphere W ~ -3GM²/5R
        let total_mass = self.spawn_cluster_mass * self.spawn_cluster_count as f64;
        let sigma2 = self.spawn_cluster_vel_disp * self.spawn_cluster_vel_disp;
        // K = N * (3/2) * m * sigma^2  (3D kinetic, 3 components)
        let kinetic_est = self.spawn_cluster_count as f64 * 1.5 * self.spawn_cluster_mass * sigma2;
        // |W| ~ 3GM²/5R   (uniform sphere)
        let pot_est = if self.spawn_cluster_radius > 0.0 {
            3.0 * G * total_mass * total_mass / (5.0 * self.spawn_cluster_radius)
        } else {
            0.0
        };
        let virial = if pot_est > 0.0 { kinetic_est / pot_est } else { 0.0 };

        let virial_label = if virial < 0.3 {
            "cold collapse"
        } else if virial < 0.7 {
            "sub-virial"
        } else if virial < 1.3 {
            "near virial"
        } else {
            "super-virial"
        };
        let virial_col = if (virial - 0.5).abs() < 0.3 { SUCCESS } else { TEXT_DIM };

        ui.add_space(4.0);
        egui::Frame::NONE
            .fill(Color32::from_rgba_unmultiplied(20, 20, 28, 200))
            .stroke(Stroke::new(0.5, BORDER))
            .corner_radius(4.0)
            .inner_margin(egui::Margin::symmetric(8, 5))
            .show(ui, |ui| {
                ui.set_width(ui.available_width());
                info_row(ui, "total M", &format!("{:.4e}", total_mass));
                info_row(
                    ui,
                    "σ_3D",
                    &format!("{:.4e}", (3.0_f64).sqrt() * self.spawn_cluster_vel_disp),
                );
                ui.horizontal(|ui| {
                    ui.label(RichText::new("Q = K/|W|").size(9.5).color(TEXT_DIM));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(RichText::new(virial_label).size(9.5).color(virial_col));
                        ui.label(
                            RichText::new(format!("{:.2}  ", virial))
                                .monospace()
                                .size(9.5)
                                .color(virial_col),
                        );
                    });
                });
            });

        ui.add_space(6.0);
        if primary_btn(ui, &format!("+ Add cluster  ({} bodies)", self.spawn_cluster_count)) {
            let n = self.spawn_cluster_count as usize;
            if n > 0 {
                let center = self.system.metrics();
                let clust_density = density(self.spawn_cluster_material, self.spawn_cluster_mass);

                self.push_undo(UndoRecord::AddedBodies(n));
                use rand::{SeedableRng, RngExt};
                use rand::rngs::SmallRng;
                let seed = self.system.seed();
                let mut rng: SmallRng = if seed == 0 {
                    rand::make_rng()
                } else {
                    SmallRng::seed_from_u64(seed)
                };
                for _ in 0..n {
                    let r = self.spawn_cluster_radius * rng.random::<f64>().sqrt();
                    let theta = rng.random::<f64>() * std::f64::consts::TAU;
                    let x = center.com_x + r * theta.cos();
                    let y = center.com_y + r * theta.sin();
                    let vx = (rng.random::<f64>() - 0.5) * self.spawn_cluster_vel_disp * 2.0;
                    let vy = (rng.random::<f64>() - 0.5) * self.spawn_cluster_vel_disp * 2.0;
                    let mut b = Body::new(
                        x,
                        y,
                        vx,
                        vy,
                        self.spawn_cluster_mass,
                        self.spawn_cluster_material,
                    );
                    b.density = clust_density;
                    b.sync_physical_properties();
                    self.system.add_body(b);
                }
            }
        }
    }
}

// ── Info row helper ───────────────────────────────────────────────────────────

fn info_row(ui: &mut egui::Ui, label: &str, value: &str) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(9.5).color(TEXT_DIM));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(value).monospace().size(9.5).color(TEXT_SEC));
        });
    });
}
