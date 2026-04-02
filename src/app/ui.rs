use crate::app::theme::apply_visuals;
use crate::core::system::System;
use crate::domain::body::Body;
use eframe::egui;

#[derive(PartialEq, Clone, Copy)]
pub enum SpawnTab {
    Single,
    Ring,
    Cluster,
}

pub struct BodyForm {
    pub x: String,
    pub y: String,
    pub vx: String,
    pub vy: String,
    pub mass: String,
}

impl Default for BodyForm {
    fn default() -> Self {
        Self {
            x: "0.0".into(),
            y: "0.0".into(),
            vx: "0.0".into(),
            vy: "0.0".into(),
            mass: "1.0".into(),
        }
    }
}

impl BodyForm {
    pub fn try_build(&self) -> Option<Body> {
        Some(Body {
            x: self.x.parse().ok()?,
            y: self.y.parse().ok()?,
            vx: self.vx.parse().ok()?,
            vy: self.vy.parse().ok()?,
            mass: self.mass.parse().ok()?,
        })
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
    pub(super) place_mode: bool,
    pub(super) place_drag_start: Option<egui::Pos2>,
    pub(super) place_mass: f64,
    pub(super) spawn_tab: SpawnTab,
    pub(super) spawn_ring_radius: f64,
    pub(super) spawn_ring_count: u32,
    pub(super) spawn_ring_mass: f64,
    pub(super) spawn_ring_vel_scale: f64,
    pub(super) spawn_cluster_radius: f64,
    pub(super) spawn_cluster_count: u32,
    pub(super) spawn_cluster_mass: f64,
    pub(super) spawn_cluster_vel_disp: f64,
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
            place_mode: false,
            place_drag_start: None,
            place_mass: 1.0,
            spawn_tab: SpawnTab::Single,
            spawn_ring_radius: 10.0,
            spawn_ring_count: 60,
            spawn_ring_mass: 0.01,
            spawn_ring_vel_scale: 1.0,
            spawn_cluster_radius: 5.0,
            spawn_cluster_count: 30,
            spawn_cluster_mass: 1.0,
            spawn_cluster_vel_disp: 0.5,
        }
    }
}

impl eframe::App for SimulationApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        apply_visuals(ctx);
        self.draw_panel(ctx);
        self.draw_canvas(ctx);

        if !self.paused {
            for _ in 0..self.steps_per_frame {
                self.system.step_adaptive(self.proposed_dt);
            }
        }

        ctx.request_repaint();
    }
}
