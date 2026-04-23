use crate::app::config::PhysicsConfig;
use crate::app::render_hints::{BodyRenderHints, compute_render_hints};
use crate::app::theme::{BG, apply_visuals};
use crate::core::physics_thread::{PhysicsHandle, spawn as spawn_physics};
use crate::core::system::System;
use crate::domain::body::Body;
use crate::domain::materials::Material;
use crate::io::recorder::SimRecorder;
use crate::io::snapshot::{SaveEntry, SimSnapshot, list_saves};
use crate::physics::integrator::IntegratorKind;
use crate::render::{TrailRenderer, WgpuBackend};
use crate::templates::{Template, UnitSystem};
use std::collections::HashSet;
use std::sync::{Arc, Mutex};

// ── Undo ──────────────────────────────────────────────────────────────────────

/// A reversible simulation mutation. Stored in a bounded stack; Ctrl+Z pops
/// and reverses the last entry.
pub enum UndoRecord {
    /// N bodies were appended to the end of the body list.
    /// Undo: remove the last N bodies.
    AddedBodies(usize),
    /// A body was removed at `idx`.
    /// Undo: re-append it (index is not preserved, but state is).
    RemovedBody { body: Body, name: String },
    /// A body at `idx` was edited (position, velocity, mass, …).
    /// Undo: restore old values.
    EditedBody { idx: usize, old_body: Body, old_name: String },
}

/// Maximum number of undo records kept in memory.
const UNDO_LIMIT: usize = 20;

#[derive(PartialEq, Clone, Copy)]
pub enum SpawnTab {
    Single,
    Ring,
    Cluster,
}

#[derive(PartialEq, Clone, Copy)]
pub enum PanelTab {
    Overview,
    Add,
    Templates,
    View,
    Camera,
    Config,
}

impl PanelTab {
    pub fn label(self) -> &'static str {
        match self {
            Self::Overview => "Overview",
            Self::Add => "Add Body",
            Self::Templates => "Templates",
            Self::View => "Display",
            Self::Camera => "Camera",
            Self::Config => "Physics",
        }
    }

    pub const ALL: [PanelTab; 6] = [
        PanelTab::Overview,
        PanelTab::Add,
        PanelTab::Templates,
        PanelTab::View,
        PanelTab::Camera,
        PanelTab::Config,
    ];
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
    /// Per-level visibility for the global orbit overlay. Indexed 0..=3
    /// where 3 also covers all deeper sub-satellite levels.
    pub(super) orbit_visible_levels: [bool; 4],
    /// Cap on the number of global orbit overlays drawn per frame.
    /// Candidates are ranked by an influence score combining log-mass and
    /// viewport proximity; only the top-N survive.
    pub(super) orbit_top_n: usize,
    /// When `true`, orbits whose geometry is numerically fragile are
    /// suppressed: near-parabolic (|1 − e| < 0.005) and orbits whose
    /// periapsis lies inside the primary's physical radius.
    pub(super) orbit_hide_degenerate: bool,
    /// Frame-coherent gravitational hierarchy used by the orbit overlay
    /// filter pipeline. Not persisted — rebuilt from live state each tick.
    pub(super) orbit_hierarchy: crate::physics::orbit_hierarchy::OrbitHierarchy,
    /// Bodies whose orbits are drawn unconditionally (bypass level, top-N
    /// and degeneracy filters). Pins are stored by body index; the canvas
    /// prunes out-of-range entries each frame so collision-merges don't
    /// leave dangling pins.
    pub(super) pinned_orbits: HashSet<usize>,
    pub(super) show_grid: bool,
    pub(super) show_vectors: bool,
    /// Target sim-time advance per real second (sim units/s).
    /// Maps directly to `PhysicsCmd::SetSimRateTarget`.
    /// Default: 2π ≈ 1 yr/s in internal units (G=1, AU, solar masses).
    pub(super) sim_rate_target: f64,

    /// IAS15 error tolerance. Forwarded to the physics thread every frame.
    /// Ignored when a different integrator is active.
    pub(super) ias15_epsilon: f64,
    /// Simulation-time duration the user has configured for the next
    /// Precision Run. The run's `t_target` is resolved as
    /// `system.t() + precision_run_duration` at Start time, so changing
    /// the field does not retroactively affect an in-flight run.
    pub(super) precision_run_duration: f64,
    /// When `Some`, the user has selected a Precision-profile
    /// integrator and the confirmation modal is pending resolution.
    /// The variant is the kind to apply if the user accepts. While
    /// the modal is up, `physics_cfg.integrator` and the physics
    /// thread's integrator are NOT changed yet — the intent is held
    /// here until Continue / Cancel resolves the dialog.
    pub(super) precision_confirmation_pending: Option<IntegratorKind>,
    /// When `true`, subsequent selections of a Precision-profile
    /// integrator skip the modal for the rest of the session. Reset
    /// on restart — this is a session-local courtesy, not a durable
    /// preference, so "don't show again" does not silently survive
    /// across app restarts.
    pub(super) precision_confirmation_session_skip: bool,
    pub(super) place_mode: bool,
    pub(super) place_drag_start: Option<egui::Pos2>,
    pub(super) place_mass: f64,
    pub(super) place_density: f64,
    pub(super) spawn_tab: SpawnTab,
    pub(super) spawn_ring_radius: f64,
    pub(super) spawn_ring_count: u32,
    pub(super) spawn_ring_mass: f64,
    pub(super) spawn_ring_vel_scale: f64,
    pub(super) spawn_ring_material: Material,
    pub(super) spawn_cluster_radius: f64,
    pub(super) spawn_cluster_count: u32,
    pub(super) spawn_cluster_mass: f64,
    pub(super) spawn_cluster_vel_disp: f64,
    pub(super) spawn_cluster_material: Material,
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
    /// Visual preset for trails; combined with `trail_width` to produce the
    /// concrete [`crate::render::TrailStyle`] pushed to the backend.
    pub(super) trail_style_preset: crate::render::TrailStylePreset,
    /// Render-side trail recorder. Owns the ring buffer and sampling policy.
    pub(super) trail_recorder: crate::render::TrailRecorder,
    /// Minimum body-mass / dominant-mass ratio required to show a trail.
    /// Bodies below this threshold have their trail alpha zeroed out.
    /// Default 1e-6 suppresses asteroid-mass bodies automatically.
    pub(super) trail_min_mass_ratio: f64,

    pub(super) place_material: Material,

    pub(super) trail: Option<TrailRenderer>,

    // Camera inertia + animation
    pub(super) zoom_vel: f32,
    pub(super) pan_vel: egui::Vec2,
    pub(super) follow_selected_body: bool,
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

    // ── Unit system ───────────────────────────────────────────────────────────
    /// Unit system of the most recently loaded template.
    ///
    /// Updated each time a template is dropped onto the canvas.  Used to
    /// annotate the CSV metadata header with physical unit information.
    /// Defaults to `UnitSystem::dimensionless()` until a template is loaded.
    pub(super) active_units: UnitSystem,

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

    // ── Undo ──────────────────────────────────────────────────────────────────
    pub(super) undo_stack: Vec<UndoRecord>,

    // ── Shortcuts guide ───────────────────────────────────────────────────────
    pub(super) show_shortcuts_modal: bool,

    // ── Settings modal ────────────────────────────────────────────────────────
    pub(super) show_settings_modal: bool,

    // ── Sidebar visibility ───────────────────────────────────────────────────
    /// When `true`, the left contextual panel is hidden and the canvas
    /// expands to fill its space. Toggled by the toolbar button, `B`, or
    /// clicking the active tool tab again.
    pub(super) sidebar_collapsed: bool,

    // ── Drift tracking ────────────────────────────────────────────────────────
    /// Running maximum of |dE/E₀| seen since the last simulation reset.
    /// Never decreases mid-run; reset to 0 on load / clear / snapshot restore.
    pub(super) energy_drift_peak: f64,
    /// Running maximum of |dLz/Lz₀| (only tracked when |Lz₀| is non-trivial).
    pub(super) lz_drift_peak: f64,

    // ── Single-step ───────────────────────────────────────────────────────────
    /// When `true`, the simulation was unpaused for exactly one frame.
    /// `draw_frame` re-pauses on the next tick after physics has run.
    pub(super) step_pending: bool,

    // ── Overview search ───────────────────────────────────────────────────────
    pub(super) overview_search: String,

    // ── Templates modal ───────────────────────────────────────────────────────
    pub(super) show_templates_modal: bool,
    pub(super) templates_search: String,
    /// Pre-computed (description, body_count) for every TEMPLATES entry.
    /// Built once at startup so the modal never rebuilds templates per frame.
    pub(super) templates_meta: Vec<(&'static str, usize)>,

    // ── Simulation identity ───────────────────────────────────────────────────
    /// User-assigned name for the current simulation. Empty until the user names it.
    pub(super) sim_name: String,
    /// When `true`, show the "Name this simulation" prompt before the next save.
    pub(super) show_name_prompt: bool,
    /// Editable buffer for the name prompt text field.
    pub(super) pending_name_input: String,
    /// Reproducibility seed for the current simulation.
    /// Generated fresh when a new simulation starts, preserved when loading a save.
    pub(super) sim_seed: u64,

    // ── Data-driven colour pipeline (SPLASH / yt-style) ──────────────────────
    /// Scalar fields available for body colouring (velocity, mass, …).
    pub(super) field_registry: crate::domain::field::FieldRegistry,
    /// Colour ramps available (viridis, inferno, plasma, cool_warm, grayscale).
    pub(super) colormap_registry: crate::render::color::ColormapRegistry,
    /// Normalizers available (linear, log).
    pub(super) normalizer_registry: crate::render::color::NormalizerRegistry,
    /// Active colour view. `None` falls back to material-based body colours,
    /// which is the default — data-driven colouring is opt-in.
    pub(super) color_view: Option<crate::render::color::ColorViewSelection>,
    /// Last resolved data range from the active colour view. Cached so the
    /// colour bar and numeric readouts can render without re-evaluating.
    pub(super) color_view_range: Option<(f64, f64)>,
}

impl SimulationApp {
    pub fn new(system: System) -> Self {
        let mut physics_cfg = PhysicsConfig::default();
        physics_cfg.integrator = system.integrator_kind();
        physics_cfg.theta = system.theta();
        physics_cfg.softening_scale = system.softening_scale();

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
            orbit_visible_levels: [true, true, true, false],
            orbit_top_n: 24,
            orbit_hide_degenerate: true,
            orbit_hierarchy: crate::physics::orbit_hierarchy::OrbitHierarchy::new(),
            pinned_orbits: HashSet::new(),
            show_grid: true,
            show_vectors: false,
            sim_rate_target: std::f64::consts::TAU,
            ias15_epsilon: 1e-9,
            // Default: roughly one "internal year" at the default
            // unit system (G = 1 gives orbital period = 2π). Users can
            // override before starting a run.
            precision_run_duration: 2.0 * std::f64::consts::PI,
            precision_confirmation_pending: None,
            precision_confirmation_session_skip: false,
            place_mode: false,
            place_drag_start: None,
            place_mass: 1.0,
            place_density: 1.0,
            spawn_tab: SpawnTab::Single,
            spawn_ring_radius: 10.0,
            spawn_ring_count: 60,
            spawn_ring_mass: 0.01,
            spawn_ring_vel_scale: 1.0,
            spawn_ring_material: Material::Rocky,
            spawn_cluster_radius: 5.0,
            spawn_cluster_count: 30,
            spawn_cluster_mass: 1.0,
            spawn_cluster_vel_disp: 0.5,
            spawn_cluster_material: Material::Rocky,
            selected_body: None,
            dragging_body: None,
            drag_start_world: None,
            selection_form: None,
            physics_cfg,
            panel_tab: PanelTab::Overview,
            show_force_vectors: false,
            render_hints: Vec::new(),
            body_angles: Vec::new(),
            template_drag: None,
            show_belts: false,
            trail_width: 1.5,
            trail_style_preset: crate::render::TrailStylePreset::UniverseSandbox,
            trail_recorder: crate::render::TrailRecorder::new(),
            trail_min_mass_ratio: 1e-7,
            place_material: Material::Rocky,
            trail: None,

            zoom_vel: 0.0,
            pan_vel: egui::Vec2::ZERO,
            follow_selected_body: false,
            camera_anim_target: None,
            pending_fit: false,
            hovered_body: None,

            backend: Arc::new(Mutex::new(WgpuBackend::new())),
            device: None,
            queue: None,
            format: None,

            active_units: UnitSystem::dimensionless(),

            recorder: None,
            record_interval: 0.01,
            record_base_path: "./records/sim_export".into(),
            record_error: None,

            save_dir: "./saves".into(),
            autosave_interval_secs: 300.0,
            last_save_instant: std::time::Instant::now(),
            show_save_modal: false,
            save_modal_entries: Vec::new(),
            save_modal_error: None,
            pending_load: None,

            step_pending: false,
            sidebar_collapsed: false,
            undo_stack: Vec::new(),
            show_shortcuts_modal: false,
            show_settings_modal: false,

            energy_drift_peak: 0.0,
            lz_drift_peak: 0.0,

            overview_search: String::new(),

            show_templates_modal: false,
            templates_search: String::new(),
            templates_meta: crate::templates::TEMPLATES
                .iter()
                .map(|e| { let t = (e.build)(0); (t.description, t.body_count()) })
                .collect(),

            sim_name: String::new(),
            show_name_prompt: false,
            pending_name_input: String::new(),
            sim_seed: SimSnapshot::new_seed(),

            field_registry: crate::domain::field::FieldRegistry::standard(),
            colormap_registry: crate::render::color::ColormapRegistry::standard(),
            normalizer_registry: crate::render::color::NormalizerRegistry::standard(),
            color_view: None,
            color_view_range: None,
        }
    }

    fn draw_frame(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        apply_visuals(&ctx);

        // ── Sync latest physics state into local cache ────────────────────────
        self.system.sync();

        // ── Dense-output: interpolate body positions to the render instant ────
        // Advances t_render by sim_rate_target × wall_delta so bodies move
        // smoothly between physics publishes.  Skipped while paused so the
        // display freezes at the last physics position.
        if !self.paused {
            let wall_delta = ctx.input(|i| i.unstable_dt as f64).min(0.2);
            self.system.advance_render_time(wall_delta, self.sim_rate_target);
        }

        // ── Single-step: re-pause after one frame of physics ─────────────────
        if self.step_pending {
            self.step_pending = false;
            self.paused = true;
        }

        // ── Pending fit-to-view (after async template/snapshot load) ──────────
        if self.pending_fit && !self.system.bodies().is_empty() && !self.system.is_loading() {
            self.fit_to_view();
            self.pending_fit = false;
        }

        // Forward UI-controlled parameters to the physics thread every frame.
        // These are cheap fire-and-forget sends; the thread drains them before
        // each batch, so latency is at most one batch period (~100 µs).
        self.system.set_paused(self.paused);
        self.system.set_sim_rate_target(self.sim_rate_target);
        self.system.set_ias15_epsilon(self.ias15_epsilon);

        // Recompute render hints from the freshly-synced body list.
        self.render_hints = compute_render_hints(self.system.bodies());

        // ── Drift peak tracking (runs after sync, before draw) ───────────────
        {
            let m = self.system.metrics();
            // Suppress noise from the very first steps before the system stabilises.
            if m.steps > 10 && m.total_energy.abs() > 1e-30 {
                let ae = m.rel_energy_error.abs();
                if ae > self.energy_drift_peak {
                    self.energy_drift_peak = ae;
                }
            }
            let lz_trivial =
                m.angular_momentum_z.abs() < 1e-10 || m.rel_angular_momentum_error.abs() > 1e3;
            if m.steps > 10 && !lz_trivial {
                let alz = m.rel_angular_momentum_error.abs();
                if alz > self.lz_drift_peak {
                    self.lz_drift_peak = alz;
                }
            }
        }

        // ── CSV recording (render-rate sampling) ──────────────────────────────
        if let Some(rec) = self.recorder.as_mut() {
            let t = self.system.t();
            if rec.should_record(t) {
                let metrics = self.system.metrics();
                let _ =
                    rec.record(t, self.system.bodies(), &metrics, self.system.orbital_elements());
            }
        }

        // ── Auto-save ─────────────────────────────────────────────────────────
        if self.autosave_interval_secs > 0.0
            && !self.system.bodies().is_empty()
            && self.last_save_instant.elapsed().as_secs_f64() >= self.autosave_interval_secs
        {
            let _ = self.do_save_auto();
        }

        // ── Global keyboard shortcuts ─────────────────────────────────────────
        let (ctrl_z, space, key_f, key_h, key_b, tool_keys) = ctx.input_mut(|i| {
            (
                i.consume_key(egui::Modifiers::CTRL, egui::Key::Z),
                i.consume_key(egui::Modifiers::NONE, egui::Key::Space),
                i.consume_key(egui::Modifiers::NONE, egui::Key::F),
                i.consume_key(egui::Modifiers::NONE, egui::Key::H),
                i.consume_key(egui::Modifiers::NONE, egui::Key::B),
                // Tool rail 1..6 — always open the sidebar onto that tool
                // (keyboard = predictable, never toggles closed).
                [
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num1),
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num2),
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num3),
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num4),
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num5),
                    i.consume_key(egui::Modifiers::NONE, egui::Key::Num6),
                ],
            )
        });
        if ctrl_z {
            self.perform_undo();
        }
        if space {
            self.paused = !self.paused;
        }
        if key_f {
            self.fit_to_view();
        }
        if key_h {
            self.show_shortcuts_modal = !self.show_shortcuts_modal;
        }
        if key_b {
            self.sidebar_collapsed = !self.sidebar_collapsed;
        }
        for (i, pressed) in tool_keys.iter().enumerate() {
            if *pressed {
                self.activate_tool(PanelTab::ALL[i]);
            }
        }

        // Registration order carves space: top → bottom → left rail → left
        // panel → right inspector → central canvas. See `panel/mod.rs` for the
        // layout diagram.
        self.draw_toolbar(&ctx);
        // Bottom strip: the Precision Run panel owns the slot whenever
        // either (a) a run is actively in progress or (b) a Precision-
        // profile integrator is selected (panel then shows the Setup
        // view so the user can configure and start a run). The real-
        // time playbar returns only when neither condition is met.
        let integrator_is_precision = self.physics_cfg.integrator.execution_profile()
            == crate::physics::integrator::traits::ExecutionProfile::Precision;
        let run_in_flight = {
            let ctrl = self.system.precision_controller();
            let guard = ctrl.lock().unwrap();
            !matches!(
                guard.state(),
                crate::core::precision_run::RunState::Idle
            )
        };
        if integrator_is_precision || run_in_flight {
            self.draw_precision_panel(&ctx);
        } else {
            self.draw_playbar(&ctx);
        }
        self.draw_tool_rail(&ctx);
        if !self.sidebar_collapsed {
            self.draw_panel(&ctx);
        }
        self.draw_inspector(&ctx);

        self.draw_save_modal(&ctx);
        self.draw_shortcuts_modal(&ctx);
        self.draw_settings_modal(&ctx);
        self.draw_precision_confirmation_modal(&ctx);
        self.draw_templates_modal(&ctx);
        self.draw_name_prompt(&ctx);

        egui::CentralPanel::default().frame(egui::Frame::NONE.fill(BG)).show(&ctx, |ui| {
            self.draw_canvas(ui);
        });

        // Re-apply paused state after UI rendering — button clicks this frame
        // (play/pause, step) may have changed self.paused after the early sync.
        self.system.set_paused(self.paused);
        self.system.set_sim_rate_target(self.sim_rate_target);
        self.system.set_ias15_epsilon(self.ias15_epsilon);

        if !self.paused {
            // Running: repaint every frame to keep the canvas live.
            ctx.request_repaint();
        } else {
            // Paused: still repaint at ~20 Hz so updates from the physics thread
            // (body added, save loaded, etc.) appear without requiring user input.
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
    }

    /// Manual save. If the simulation has no name yet, shows the name prompt
    /// instead of saving immediately (returns an Err to indicate deferred).
    pub(super) fn do_save(&mut self) -> Result<std::path::PathBuf, String> {
        if self.sim_name.is_empty() {
            self.pending_name_input = String::new();
            self.show_name_prompt = true;
            return Err("Prompting for name".into());
        }
        self.do_save_named()
    }

    /// Auto-save: saves with current sim_name (or "Unnamed" if empty) without prompting.
    pub(super) fn do_save_auto(&mut self) -> Result<std::path::PathBuf, String> {
        if self.sim_name.is_empty() {
            self.sim_name = "Unnamed".to_owned();
        }
        self.do_save_named()
    }

    fn do_save_named(&mut self) -> Result<std::path::PathBuf, String> {
        let mut snap = self.system.to_snapshot();
        snap.save_id = SimSnapshot::new_id();
        snap.sim_name = self.sim_name.clone();
        snap.seed = self.sim_seed;
        // Capture trail for visual continuity on reload.
        snap.trail = Some(self.trail_recorder.to_snapshot());
        snap.trail_every = self.trail_recorder.interval_multiplier();
        let dir = std::path::Path::new(&self.save_dir);
        match snap.save_to_dir(dir) {
            Ok(p) => {
                self.last_save_instant = std::time::Instant::now();
                Ok(p)
            },
            Err(e) => Err(e.to_string()),
        }
    }

    /// Draw the "Name this simulation" prompt window.
    pub(super) fn draw_name_prompt(&mut self, ctx: &egui::Context) {
        if !self.show_name_prompt {
            return;
        }

        egui::Window::new("Name this simulation")
            .id(egui::Id::new("name_prompt"))
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .show(ctx, |ui| {
                ui.set_width(280.0);
                ui.label(
                    egui::RichText::new("Give this simulation a name so you can find it later.")
                        .size(10.0)
                        .color(crate::app::theme::TEXT_SEC),
                );
                ui.add_space(6.0);
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.pending_name_input)
                        .desired_width(ui.available_width())
                        .hint_text("e.g. Solar System"),
                );
                resp.request_focus();
                ui.add_space(6.0);
                let enter = resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                ui.horizontal(|ui| {
                    let confirmed = enter
                        || ui
                            .add(
                                egui::Button::new(
                                    egui::RichText::new("Save")
                                        .size(10.5)
                                        .color(crate::app::theme::SUCCESS),
                                )
                                .fill(egui::Color32::TRANSPARENT)
                                .stroke(egui::Stroke::new(1.0, crate::app::theme::SUCCESS))
                                .min_size(egui::vec2(60.0, 22.0)),
                            )
                            .clicked();

                    let skipped = ui
                        .add(
                            egui::Button::new(
                                egui::RichText::new("Skip")
                                    .size(10.5)
                                    .color(crate::app::theme::TEXT_DIM),
                            )
                            .fill(egui::Color32::TRANSPARENT)
                            .stroke(egui::Stroke::new(0.5, crate::app::theme::BORDER))
                            .min_size(egui::vec2(48.0, 22.0)),
                        )
                        .clicked();

                    if confirmed || skipped {
                        let name = self.pending_name_input.trim().to_owned();
                        self.sim_name = if name.is_empty() { "Unnamed".to_owned() } else { name };
                        self.show_name_prompt = false;
                        // Now that we have a name, complete the deferred save.
                        if let Ok(_) = self.do_save_named() {
                            if self.show_save_modal {
                                self.save_modal_entries =
                                    list_saves(std::path::Path::new(&self.save_dir));
                            }
                        }
                    }
                });
            });
    }

    // ── Undo ──────────────────────────────────────────────────────────────────

    /// Remap `pinned_orbits` after a `swap_remove` at `removed_idx`. The body
    /// previously at `old_last` (the pre-removal last index) now sits at
    /// `removed_idx`; update the pin set to follow that move.
    pub(super) fn pins_after_swap_remove(&mut self, removed_idx: usize, old_last: usize) {
        self.pinned_orbits.remove(&removed_idx);
        if old_last != removed_idx && self.pinned_orbits.remove(&old_last) {
            self.pinned_orbits.insert(removed_idx);
        }
    }

    /// Push a record onto the undo stack. Drops the oldest entry if the stack
    /// exceeds `UNDO_LIMIT`.
    pub(super) fn push_undo(&mut self, record: UndoRecord) {
        self.undo_stack.push(record);
        if self.undo_stack.len() > UNDO_LIMIT {
            self.undo_stack.remove(0);
        }
    }

    /// Reverse the last mutation on the undo stack.
    /// Reset drift peak accumulators. Call whenever the simulation is
    /// restarted / a new scenario loaded (peaks are no longer meaningful).
    pub(super) fn reset_drift_peaks(&mut self) {
        self.energy_drift_peak = 0.0;
        self.lz_drift_peak = 0.0;
    }

    pub(super) fn perform_undo(&mut self) {
        if let Some(record) = self.undo_stack.pop() {
            match record {
                UndoRecord::AddedBodies(n) => {
                    let total = self.system.bodies().len();
                    // Remove the last `n` bodies (they were appended, so they sit at the end).
                    for i in (total.saturating_sub(n)..total).rev() {
                        self.system.remove_body(i);
                    }
                    // Appended bodies cannot be pinned (pin requires selection), but
                    // be defensive: drop any pin that lands out of range after undo.
                    let new_len = self.system.bodies().len();
                    self.pinned_orbits.retain(|&i| i < new_len);
                    // If the selected body was one of the removed ones, clear selection.
                    if let Some(sel) = self.selected_body {
                        if sel >= total.saturating_sub(n) {
                            self.selected_body = None;
                            self.selection_form = None;
                        }
                    }
                },
                UndoRecord::RemovedBody { body, name } => {
                    // The body will land at index `bodies().len()` after the add.
                    let future_idx = self.system.bodies().len();
                    self.system.add_body(body);
                    self.system.set_name(future_idx, name);
                },
                UndoRecord::EditedBody { idx, old_body, old_name } => {
                    self.system.update_body(idx, old_body);
                    self.system.set_name(idx, old_name.clone());
                    // Keep the inspector form in sync if this body is selected.
                    if self.selected_body == Some(idx) {
                        self.selection_form = Some(SelectionForm::from_body(&old_body, &old_name));
                    }
                },
            }
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
