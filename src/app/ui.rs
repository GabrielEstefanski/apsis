use crate::app::config::PhysicsConfig;
use crate::app::theme::apply_visuals;
use crate::core::system::System;
use crate::domain::body::{Body, default_moment_inertia, radius_from_density_mass};
use eframe::egui;

#[derive(PartialEq, Clone, Copy)]
pub enum SpawnTab {
    Single,
    Ring,
    Cluster,
}

#[derive(PartialEq, Clone, Copy)]
pub enum PanelTab {
    Add,
    Templates,
    Config,
}

pub struct BodyForm {
    pub x: String,
    pub y: String,
    pub vx: String,
    pub vy: String,
    pub mass: String,
    pub density: String,
}

pub struct SelectionForm {
    pub x: String,
    pub y: String,
    pub vx: String,
    pub vy: String,
    pub mass: String,
    pub density: String,
    pub error: Option<String>,
}

impl SelectionForm {
    pub fn from_body(b: &Body) -> Self {
        Self {
            x: format!("{:.6}", b.x),
            y: format!("{:.6}", b.y),
            vx: format!("{:.6}", b.vx),
            vy: format!("{:.6}", b.vy),
            mass: format!("{:.6}", b.mass),
            density: format!("{:.6e}", b.density),
            error: None,
        }
    }
}

impl Default for BodyForm {
    fn default() -> Self {
        Self {
            x: "0.0".into(),
            y: "0.0".into(),
            vx: "0.0".into(),
            vy: "0.0".into(),
            mass: "1.0".into(),
            density: "1.0".into(),
        }
    }
}

impl BodyForm {
    pub fn try_build(&self) -> Option<Body> {
        let mass: f64 = self.mass.parse().ok().filter(|&v| v > 0.0)?;
        let density: f64 = self.density.parse().ok().filter(|&v| v > 0.0)?;
        let radius = radius_from_density_mass(density, mass);
        let mut b = Body::new(
            self.x.parse().ok()?,
            self.y.parse().ok()?,
            self.vx.parse().ok()?,
            self.vy.parse().ok()?,
            mass,
        );
        b.density = density;
        b.radius = radius;
        b.softening = b.softening.max(radius * 2.0);
        b.moment_inertia = default_moment_inertia(mass, radius);
        Some(b)
    }
}

pub struct SimulationApp {
    pub(super) system: System,
    pub(super) proposed_dt: f64,
    pub(super) paused: bool,
    pub(super) scale: f32,
    pub(super) offset: egui::Vec2,
    pub(super) form: BodyForm,
    pub(super) form_error: Option<String>,
    pub(super) show_trails: bool,
    pub(super) show_grid: bool,
    pub(super) show_vectors: bool,
    pub(super) steps_per_frame: u32,
    pub(super) collision_cor: f64,
    pub(super) place_mode: bool,
    pub(super) place_drag_start: Option<egui::Pos2>,
    pub(super) place_mass: f64,
    pub(super) place_density: f64,
    pub(super) spawn_tab: SpawnTab,
    pub(super) spawn_ring_radius: f64,
    pub(super) spawn_ring_count: u32,
    pub(super) spawn_ring_mass: f64,
    pub(super) spawn_ring_vel_scale: f64,
    pub(super) spawn_cluster_radius: f64,
    pub(super) spawn_cluster_count: u32,
    pub(super) spawn_cluster_mass: f64,
    pub(super) spawn_cluster_vel_disp: f64,
    pub(super) selected_body: Option<usize>,
    pub(super) selection_form: Option<SelectionForm>,

    pub(super) physics_cfg: PhysicsConfig,
    pub(super) panel_tab: PanelTab,
}

impl SimulationApp {
    pub fn new(system: System) -> Self {
        Self {
            system,
            proposed_dt: 1e-3,
            paused: true,
            scale: 10.0,
            offset: egui::Vec2::ZERO,
            form: BodyForm::default(),
            form_error: None,
            show_trails: true,
            show_grid: true,
            show_vectors: false,
            steps_per_frame: 1,
            collision_cor: 0.0,
            place_mode: false,
            place_drag_start: None,
            place_mass: 1.0,
            place_density: 1.0,
            spawn_tab: SpawnTab::Single,
            spawn_ring_radius: 10.0,
            spawn_ring_count: 60,
            spawn_ring_mass: 0.01,
            spawn_ring_vel_scale: 1.0,
            spawn_cluster_radius: 5.0,
            spawn_cluster_count: 30,
            spawn_cluster_mass: 1.0,
            spawn_cluster_vel_disp: 0.5,
            selected_body: None,
            selection_form: None,
            physics_cfg: PhysicsConfig::default(),
            panel_tab: PanelTab::Add,
        }
    }
}

impl eframe::App for SimulationApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        apply_visuals(ctx);
        self.draw_toolbar(ctx);
        self.draw_panel(ctx);
        self.draw_inspector(ctx);
        self.draw_canvas(ctx);

        if !self.paused {
            for _ in 0..self.steps_per_frame {
                self.system.step_adaptive(self.proposed_dt);
            }
        }

        ctx.request_repaint();
    }
}
