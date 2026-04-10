use crate::app::config::PhysicsConfig;
use crate::app::render_hints::{BodyRenderHints, compute_render_hints};
use crate::app::theme::{BG, apply_visuals};
use crate::core::physics_thread::{PhysicsHandle, spawn as spawn_physics};
use crate::core::recorder::SimRecorder;
use crate::core::snapshot::{SaveEntry, SimSnapshot, list_saves};
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
    pub name: String,
    pub x: String,
    pub y: String,
    pub vx: String,
    pub vy: String,
    pub mass: String,
    pub density: String,
    pub error: Option<String>,
}

impl SelectionForm {
    pub fn from_body(b: &Body, name: &str) -> Self {
        Self {
            name: name.to_owned(),
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
    pub(super) system: PhysicsHandle,
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
    pub(super) trail_width: f32,

    pub(super) place_material: Material,

    pub(super) trail: Option<TrailRenderer>,

    // Camera inertia + animation
    pub(super) zoom_vel: f32,
    pub(super) pan_vel: egui::Vec2,
    /// Smooth-pan target offset; `None` when idle.
    pub(super) camera_anim_target: Option<egui::Vec2>,
    /// When `true`, `draw_frame` will call `fit_to_view` on the next frame
    /// that has a non-empty body list. Used after template/snapshot loads
    /// where bodies arrive asynchronously from the physics thread.
    pub(super) pending_fit: bool,

    // Canvas hover
    pub(super) hovered_body: Option<usize>,

    pub(super) backend: Arc<Mutex<WgpuBackend>>,
    pub(super) device: Option<Arc<wgpu::Device>>,
    pub(super) queue: Option<Arc<wgpu::Queue>>,
    pub(super) format: Option<wgpu::TextureFormat>,

    // ── CSV export ────────────────────────────────────────────────────────────
    /// Active recorder, `None` when not recording.
    pub(super) recorder: Option<SimRecorder>,
    /// Simulated-time gap between successive CSV records.
    pub(super) record_interval: f64,
    /// Base path prefix (no extension); e.g. `"./run01"`.
    pub(super) record_base_path: String,
    /// Last error from starting a recording session, shown in the UI.
    pub(super) record_error: Option<String>,

    // ── Save / Load ───────────────────────────────────────────────────────────
    /// Directory where `.grav` save files are written.
    pub(super) save_dir: String,
    /// Real-time seconds between automatic saves (0 = disabled).
    pub(super) autosave_interval_secs: f64,
    /// Instant of the last successful auto- or manual save.
    pub(super) last_save_instant: std::time::Instant,
    /// `true` while the load-save browser modal is open.
    pub(super) show_save_modal: bool,
    /// Entries shown in the modal; refreshed when the modal opens.
    pub(super) save_modal_entries: Vec<SaveEntry>,
    /// Any error message to display in the modal.
    pub(super) save_modal_error: Option<String>,
    /// Snapshot staged for confirmation before loading (avoids accidental overwrites).
    pub(super) pending_load: Option<SimSnapshot>,
}

impl SimulationApp {
    pub fn new(system: System) -> Self {
        let mut physics_cfg = PhysicsConfig::default();
        physics_cfg.integrator = system.integrator();
        physics_cfg.theta = system.theta();
        physics_cfg.softening_scale = system.softening_scale();
        physics_cfg.trail_every = system.trail_every();

        let physics = spawn_physics(system, true /* start paused */);

        Self {
            system: physics,
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
            physics_cfg,
            panel_tab: PanelTab::Add,
            show_force_vectors: false,
            render_hints: Vec::new(),
            body_angles: Vec::new(),
            template_drag: None,
            show_belts: false,
            trail_width: 1.5,
            place_material: Material::Rocky,
            trail: None,

            zoom_vel: 0.0,
            pan_vel: egui::Vec2::ZERO,
            camera_anim_target: None,
            pending_fit: false,
            hovered_body: None,

            backend: Arc::new(Mutex::new(WgpuBackend::new())),
            device: None,
            queue: None,
            format: None,

            recorder: None,
            record_interval: 0.01,
            record_base_path: "./sim_export".into(),
            record_error: None,

            save_dir: "./saves".into(),
            autosave_interval_secs: 300.0,
            last_save_instant: std::time::Instant::now(),
            show_save_modal: false,
            save_modal_entries: Vec::new(),
            save_modal_error: None,
            pending_load: None,
        }
    }

    fn draw_frame(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        apply_visuals(&ctx);

        // ── Sync latest physics state into local cache ────────────────────────
        self.system.sync();

        // ── Pending fit-to-view (after async template/snapshot load) ──────────
        if self.pending_fit && !self.system.bodies().is_empty() && !self.system.is_loading() {
            self.fit_to_view();
            self.pending_fit = false;
        }

        // Forward UI-controlled parameters to the physics thread every frame.
        // These are cheap fire-and-forget sends; the thread drains them before
        // each batch, so latency is at most one batch period (~100 µs).
        self.system.set_paused(self.paused);
        self.system.set_steps_per_frame(self.steps_per_frame);

        // Recompute render hints from the freshly-synced body list.
        self.render_hints = compute_render_hints(self.system.bodies());

        // ── CSV recording (render-rate sampling) ──────────────────────────────
        if let Some(rec) = self.recorder.as_mut() {
            let t = self.system.t();
            if rec.should_record(t) {
                let metrics = self.system.metrics();
                let _ = rec.record(
                    t,
                    self.system.bodies(),
                    &metrics,
                    self.system.orbital_elements(),
                );
            }
        }

        // ── Auto-save ─────────────────────────────────────────────────────────
        if self.autosave_interval_secs > 0.0
            && !self.system.bodies().is_empty()
            && self.last_save_instant.elapsed().as_secs_f64() >= self.autosave_interval_secs
        {
            let _ = self.do_save();
        }

        self.draw_toolbar(&ctx);
        self.draw_panel(&ctx);
        self.draw_inspector(&ctx);
        self.draw_save_modal(&ctx);

        egui::CentralPanel::default()
            .frame(egui::Frame::NONE.fill(BG))
            .show(&ctx, |ui| {
                self.draw_canvas(ui);
            });

        if !self.paused {
            // Running: repaint every frame to keep the canvas live.
            ctx.request_repaint();
        } else {
            // Paused: still repaint at ~20 Hz so updates from the physics thread
            // (body added, save loaded, etc.) appear without requiring user input.
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
    }

    /// Perform a manual or auto-save. Returns the saved path on success.
    pub(super) fn do_save(&mut self) -> Result<std::path::PathBuf, String> {
        let mut snap = self.system.to_snapshot();
        snap.save_id = SimSnapshot::new_id();
        let dir = std::path::Path::new(&self.save_dir);
        match snap.save_to_dir(dir) {
            Ok(p) => {
                self.last_save_instant = std::time::Instant::now();
                Ok(p)
            }
            Err(e) => Err(e.to_string()),
        }
    }

    /// Open the modal and refresh the file listing.
    pub(super) fn open_save_modal(&mut self) {
        self.save_modal_entries = list_saves(std::path::Path::new(&self.save_dir));
        self.save_modal_error = None;
        self.pending_load = None;
        self.show_save_modal = true;
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
