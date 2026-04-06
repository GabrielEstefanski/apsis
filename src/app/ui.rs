use crate::app::config::PhysicsConfig;
use crate::app::theme::apply_visuals;
use crate::core::system::System;
use crate::domain::body::{Body, default_moment_inertia};
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
        b.radius = b.physical_radius;
        b.softening = b.softening.max(b.physical_radius * 2.0);
        b.moment_inertia = default_moment_inertia(mass, b.physical_radius);
        Some(b)
    }
}

/// Short-lived visual effect spawned at a collision site.
pub struct ImpactEffect {
    pub world_x: f64,
    pub world_y: f64,
    /// Impact normal direction (bj → bi).
    pub nx: f32,
    pub ny: f32,
    /// Age in [0, 1]: 0 = just created, 1 = expired.  Lifetime ≈ 0.2 s.
    pub age: f32,
    /// Burst particles: [world_x, world_y, world_vx, world_vy].
    pub particles: Vec<[f64; 4]>,
}

pub struct SimulationApp {
    pub(super) system: System,
    pub(super) proposed_dt: f64,
    pub(super) paused: bool,
    pub(super) scale: f32,
    pub(super) body_size_boost: f32,
    pub(super) semantic_scale_mode: SemanticScaleMode,
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

    pub(super) show_force_vectors: bool,
    pub(super) show_impact_normals: bool,
    /// Per-body accumulated rotation angle (radians), for the spoke indicator.
    pub(super) body_angles: Vec<f64>,
    /// Active visual impact effects.
    pub(super) impact_effects: Vec<ImpactEffect>,
}

impl SimulationApp {
    pub fn new(system: System) -> Self {
        Self {
            system,
            proposed_dt: 1e-3,
            paused: true,
            scale: 10.0,
            body_size_boost: 64.0,
            semantic_scale_mode: SemanticScaleMode::Comparative,
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
            show_force_vectors: false,
            show_impact_normals: false,
            body_angles: Vec::new(),
            impact_effects: Vec::new(),
        }
    }
}

impl eframe::App for SimulationApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        apply_visuals(ctx);

        // ── Advance visual effects (real-wall-clock time) ─────────────────── //
        let dt_real = ctx.input(|i| i.unstable_dt);

        // Spawn burst effects from the previous physics step's collision events.
        for event in self.system.take_impact_events() {
            let kick = (event.v_rel * 0.04).clamp(0.005, 1.5);
            const N: usize = 10;
            let mut particles = Vec::with_capacity(N);
            for k in 0..N {
                let angle = std::f64::consts::TAU * k as f64 / N as f64;
                // Deterministic speed variation: avoids rand dependency in UI layer.
                let r_factor = 0.55 + 0.45 * ((k * 7 + 3) % N) as f64 / N as f64;
                particles.push([
                    event.x,
                    event.y,
                    angle.cos() * kick * r_factor,
                    angle.sin() * kick * r_factor,
                ]);
            }
            self.impact_effects.push(ImpactEffect {
                world_x: event.x,
                world_y: event.y,
                nx: event.nx as f32,
                ny: event.ny as f32,
                age: 0.0,
                particles,
            });
        }

        // Advance and expire effects.
        for effect in &mut self.impact_effects {
            effect.age = (effect.age + dt_real / 0.2).min(1.0);
            for p in &mut effect.particles {
                p[0] += p[2] * dt_real as f64;
                p[1] += p[3] * dt_real as f64;
            }
        }
        self.impact_effects.retain(|e| e.age < 1.0);

        // ── Physics step ──────────────────────────────────────────────────── //
        if !self.paused {
            for _ in 0..self.steps_per_frame {
                self.system.step_adaptive(self.proposed_dt);
            }
        }

        // ── Update rotation angles (must happen after step) ───────────────── //
        {
            let bodies = self.system.bodies();
            if self.body_angles.len() != bodies.len() {
                self.body_angles.resize(bodies.len(), 0.0);
            }
            if !self.paused {
                let phys_dt = self.system.metrics().dt
                    * self.steps_per_frame as f64;
                for (i, b) in bodies.iter().enumerate() {
                    self.body_angles[i] += b.omega_z * phys_dt;
                }
            }
        }

        // ── Render ────────────────────────────────────────────────────────── //
        self.draw_toolbar(ctx);
        self.draw_panel(ctx);
        self.draw_inspector(ctx);
        self.draw_canvas(ctx);

        ctx.request_repaint();
    }
}
