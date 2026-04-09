use crate::app::config::PhysicsConfig;
use crate::app::render_hints::{BodyRenderHints, compute_render_hints};
use crate::app::theme::{BG, apply_visuals};
use crate::core::system::System;
use crate::domain::body::Body;
use crate::domain::materials::Material;
use crate::render::{TrailRenderer, WgpuBackend};
use crate::templates::Template;
use std::sync::{Arc, Mutex};

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

#[derive(PartialEq, Clone, Copy)]
pub enum SemanticScaleMode {
    Physical,
    Comparative,
    Illustrative,
}

impl SemanticScaleMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Physical => "physical",
            Self::Comparative => "compare",
            Self::Illustrative => "illustrate",
        }
    }
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
        let mut b = Body::new(
            self.x.parse().ok()?,
            self.y.parse().ok()?,
            self.vx.parse().ok()?,
            self.vy.parse().ok()?,
            mass,
            crate::domain::materials::Material::Rocky,
        );

        b.density = density;
        b.sync_physical_properties();

        Some(b)
    }
}

pub struct SimulationApp {
    pub(super) system: System,
    pub(super) paused: bool,
    pub(super) scale: f32,
    pub(super) body_size_boost: f32,
    pub(super) semantic_scale_mode: SemanticScaleMode,
    pub(super) offset: egui::Vec2,
    pub(super) form: BodyForm,
    pub(super) form_error: Option<String>,
    pub(super) show_trails: bool,
    pub(super) show_orbit_ellipses: bool,
    pub(super) show_grid: bool,
    pub(super) show_vectors: bool,
    pub(super) steps_per_frame: u32,
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
    pub(super) dragging_body: Option<usize>,
    pub(super) drag_start_world: Option<(f64, f64)>,
    pub(super) selection_form: Option<SelectionForm>,
    pub(super) physics_cfg: PhysicsConfig,
    pub(super) panel_tab: PanelTab,
    pub(super) show_force_vectors: bool,
    pub(super) template_drag: Option<Box<dyn Fn() -> Template>>,
    pub(super) body_angles: Vec<f64>,
    pub(super) render_hints: Vec<BodyRenderHints>,
    pub(super) show_belts: bool,

    pub(super) place_material: Material,

    pub(super) trail: Option<TrailRenderer>,

    // Camera inertia
    pub(super) zoom_vel: f32,
    pub(super) pan_vel: egui::Vec2,

    pub(super) backend: Arc<Mutex<WgpuBackend>>,
    pub(super) device: Option<Arc<wgpu::Device>>,
    pub(super) queue: Option<Arc<wgpu::Queue>>,
    pub(super) format: Option<wgpu::TextureFormat>,
}

impl SimulationApp {
    pub fn new(system: System) -> Self {
        Self {
            system,
            paused: true,
            scale: 10.0,
            body_size_boost: 64.0,
            semantic_scale_mode: SemanticScaleMode::Comparative,
            offset: egui::Vec2::ZERO,
            form: BodyForm::default(),
            form_error: None,
            show_trails: true,
            show_orbit_ellipses: false,
            show_grid: true,
            show_vectors: false,
            steps_per_frame: 1,
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
            dragging_body: None,
            drag_start_world: None,
            selection_form: None,
            physics_cfg: PhysicsConfig::default(),
            panel_tab: PanelTab::Add,
            show_force_vectors: false,
            render_hints: Vec::new(),
            body_angles: Vec::new(),
            template_drag: None,
            show_belts: false,
            place_material: Material::Rocky,
            trail: None,

            zoom_vel: 0.0,
            pan_vel: egui::Vec2::ZERO,

            backend: Arc::new(Mutex::new(WgpuBackend::new())),
            device: None,
            queue: None,
            format: None,
        }
    }

    fn draw_frame(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        apply_visuals(&ctx);

        if !self.paused {
            for _ in 0..self.steps_per_frame {
                self.system.step();
            }

            self.system.push_trail();
            self.render_hints = compute_render_hints(self.system.bodies());
        }

        // 🔴 PAINÉIS NO CONTEXT (ANTES DE QUALQUER CENTRAL)
        self.draw_toolbar(&ctx);
        self.draw_panel(&ctx);
        self.draw_inspector(&ctx);

        // 🔴 CENTRAL POR ÚLTIMO
        egui::CentralPanel::default()
            .frame(egui::Frame::NONE.fill(BG))
            .show(&ctx, |ui| {
                self.draw_canvas(ui);
            });

        if !self.paused {
            ctx.request_repaint();
        }
    }
}

impl eframe::App for SimulationApp {
    fn save(&mut self, _storage: &mut dyn eframe::Storage) {}

    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let render_state = frame.wgpu_render_state().unwrap();

        self.device = Some(render_state.device.clone().into());
        self.queue = Some(render_state.queue.clone().into());
        self.format = Some(render_state.target_format);

        self.draw_frame(ui);
    }
}
