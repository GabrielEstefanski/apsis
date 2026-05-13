use crate::app::config::PhysicsConfig;
use crate::app::design::theme as design_theme;
use crate::app::render_hints::{BodyRenderHints, compute_render_hints};
use crate::app::theme::BG;
use crate::render::{TrailRenderer, WgpuBackend};
use apsis::core::physics_thread::{PhysicsHandle, spawn as spawn_physics};
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::domain::body_preset::{self, BodyClass, BodyPreset};
use apsis::io::recorder::SimRecorder;
use apsis::io::snapshot::{SaveEntry, SimSnapshot, list_saves};
use apsis::physics::integrator::IntegratorKind;
use apsis::templates::{Template, UnitSystem};
use std::collections::{BTreeSet, HashSet};
use std::sync::{Arc, Mutex};

// ── Trail class filter ────────────────────────────────────────────────────────

/// Per-class visibility checkboxes for the trail overlay.
///
/// `Star`, `Planet`, and `Moon` start visible; `Asteroid` and `Comet`
/// start hidden so a freshly-loaded `solar_system` template does not
/// drown the canvas in 600+ minor-body trails. Toggle individually
/// from the view tab.
#[derive(Debug, Clone, Copy)]
pub struct TrailClassFilter {
    pub star: bool,
    pub planet: bool,
    pub moon: bool,
    pub asteroid: bool,
    pub comet: bool,
}

impl Default for TrailClassFilter {
    fn default() -> Self {
        Self { star: true, planet: true, moon: true, asteroid: false, comet: false }
    }
}

impl TrailClassFilter {
    /// Whether the given class passes the filter. Bodies whose class
    /// is [`BodyClass::Unknown`] always pass — class filtering only
    /// gates explicitly tagged bodies, leaving authored-set + per-body
    /// override to drive the rest.
    pub fn allows(&self, class: BodyClass) -> bool {
        match class {
            BodyClass::Star => self.star,
            BodyClass::Planet => self.planet,
            BodyClass::Moon => self.moon,
            BodyClass::Asteroid => self.asteroid,
            BodyClass::Comet => self.comet,
            BodyClass::Unknown => true,
        }
    }
}

/// Camera defaults a template suggests for its first frame. Both
/// fields are optional — `None` falls back to bounding-sphere fit and
/// world-Y-up convention respectively.
#[derive(Debug, Clone, Copy)]
pub struct TemplateCameraHints {
    pub up: Option<[f64; 3]>,
    pub distance: Option<f64>,
}

// ── Body selection ────────────────────────────────────────────────────────────

/// Selection state for bodies on the canvas.
///
/// Invariant: `Multi` always contains `len >= 2` indices. Transitions through
/// [`BodySelection::toggle`], [`BodySelection::select_single`], and
/// [`BodySelection::default`] maintain this invariant automatically.
#[derive(Default)]
pub enum BodySelection {
    #[default]
    None,
    /// Exactly one body — the primary/focused body for camera-follow and single
    /// body inspector.
    Single(usize),
    /// Two or more bodies selected simultaneously. `len >= 2` always.
    Multi(BTreeSet<usize>),
}

impl BodySelection {
    pub fn is_some(&self) -> bool {
        !matches!(self, Self::None)
    }

    /// Returns the index when exactly one body is selected; `None` otherwise.
    pub fn single(&self) -> Option<usize> {
        match self {
            Self::Single(i) => Some(*i),
            _ => None,
        }
    }

    /// Returns `true` if `idx` is anywhere in the selection.
    pub fn contains(&self, idx: usize) -> bool {
        match self {
            Self::None => false,
            Self::Single(i) => *i == idx,
            Self::Multi(set) => set.contains(&idx),
        }
    }

    /// Toggle `idx` into/out of the selection, normalising the invariant:
    /// - `None + idx` → `Single(idx)`
    /// - `Single(idx) + idx` → `None`
    /// - `Single(i) + idx` → `Multi({i, idx})`
    /// - `Multi(set) + idx` → removes if present; collapses to `Single` / `None`
    ///   when the set shrinks below two elements.
    pub fn toggle(self, idx: usize) -> Self {
        match self {
            Self::None => Self::Single(idx),
            Self::Single(i) if i == idx => Self::None,
            Self::Single(i) => {
                let mut set = BTreeSet::new();
                set.insert(i);
                set.insert(idx);
                Self::Multi(set)
            },
            Self::Multi(mut set) => {
                if set.contains(&idx) {
                    set.remove(&idx);
                } else {
                    set.insert(idx);
                }
                match set.len() {
                    0 => Self::None,
                    1 => Self::Single(*set.iter().next().unwrap()),
                    _ => Self::Multi(set),
                }
            },
        }
    }

    /// Replace the selection with exactly one body.
    pub fn select_single(idx: usize) -> Self {
        Self::Single(idx)
    }
}

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
    /// Whole body list was replaced (e.g. clicking a template). Stores the
    /// pre-replacement snapshot so undo can restore it.
    /// Undo: reload the snapshot via `load_named_bodies`.
    ReplacedBodies { previous: Vec<apsis::domain::body::NamedBody> },
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
            x: format!("{:.6}", b.pos_x),
            y: format!("{:.6}", b.pos_y),
            vx: format!("{:.6}", b.vel_x),
            vy: format!("{:.6}", b.vel_y),
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
        let mut b = Body::rocky(mass)
            .at(self.x.parse().ok()?, self.y.parse().ok()?)
            .with_velocity(self.vx.parse().ok()?, self.vy.parse().ok()?);

        b.density = density;
        b.sync_physical_properties();

        Some(b)
    }
}

pub struct SimulationApp {
    pub(super) system: PhysicsHandle,
    pub(super) paused: bool,
    pub(super) semantic_scale_mode: SemanticScaleMode,
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
    pub(super) orbit_hierarchy: apsis::physics::orbit_hierarchy::OrbitHierarchy,
    /// EMA filter on osculating elements of overlay-drawn orbits. Removes
    /// the high-frequency jitter induced by indirect N-body perturbations
    /// without falsifying secular dynamics. State-only; no physics impact.
    pub(super) orbit_smoother: crate::render::orbit_smoother::OrbitSmoother,
    /// Simulation time at the previous render frame; used to compute
    /// dt_sim for [`orbit_smoother`]. NaN sentinel = first frame; the
    /// smoother snaps to current invariants in that case.
    pub(super) orbit_smoother_last_t: f64,
    /// Bodies whose orbits are drawn unconditionally (bypass level, top-N
    /// and degeneracy filters). Pins are stored by body index; the canvas
    /// prunes out-of-range entries each frame so collision-merges don't
    /// leave dangling pins.
    pub(super) pinned_orbits: HashSet<usize>,
    pub(super) show_grid: bool,
    pub(super) show_vectors: bool,
    /// User-controlled exposure offset in stops. Multiplies the
    /// reflective HDR plane by `2^exposure_ev` on top of the auto-
    /// exposure scalar. Zero is "no offset"; positive brightens.
    pub(super) exposure_ev: f32,
    /// Target sim-time advance per real second (sim units/s).
    /// Maps directly to `PhysicsCmd::SetSimRateTarget`.
    /// Default: 2π ≈ 1 yr/s in internal units (G=1, AU, solar masses).
    pub(super) sim_rate_target: f64,
    /// Latched state of the playbar's "physics behind target" cue.
    /// Driven by [`crate::app::panel::playbar::shortfall_with_hysteresis`]
    /// — enters at 80 % achieved, exits at 90 %, so the inline
    /// `target → achieved` text doesn't flicker when the achieved
    /// rate hovers around the threshold. Initialised to `false`;
    /// updated each frame by the playbar.
    pub(super) shortfall_active: bool,

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
    /// Consumer-side ring buffer of events published on
    /// [`apsis::core::log`]. Attached to the bus at app construction;
    /// coalesces same-key bursts and feeds the top-bar notification
    /// center.
    pub(super) notifications: Arc<Mutex<crate::app::notifications::NotificationStore>>,
    pub(super) show_notifications_panel: bool,
    pub(super) notifications_filter: crate::app::notifications::NotificationFilter,
    pub(super) place_mode: bool,
    pub(super) place_drag_start: Option<egui::Pos2>,
    pub(super) place_mass: f64,
    pub(super) place_density: f64,
    pub(super) spawn_tab: SpawnTab,
    pub(super) spawn_ring_radius: f64,
    pub(super) spawn_ring_count: u32,
    pub(super) spawn_ring_mass: f64,
    pub(super) spawn_ring_vel_scale: f64,
    pub(super) spawn_ring_preset: &'static BodyPreset,
    pub(super) spawn_cluster_radius: f64,
    pub(super) spawn_cluster_count: u32,
    pub(super) spawn_cluster_mass: f64,
    pub(super) spawn_cluster_vel_disp: f64,
    pub(super) spawn_cluster_preset: &'static BodyPreset,
    pub(super) selection: BodySelection,
    /// Persistent state for the new design-system Inspector — `More`
    /// expander toggles plus the Bloomberg-flash tracker. Lives on
    /// the app rather than inside [`InspectorData`] so the per-frame
    /// payload stays a pure value.
    pub(super) inspector_state: crate::app::inspector::InspectorState,
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
    pub(super) trail_recorder: apsis::core::trail::TrailRecorder,
    /// Per-class visibility for trails (Star, Planet, Moon, Asteroid,
    /// Comet). A body whose class checkbox is unchecked has its trail
    /// hidden regardless of the authored-set rule. Class
    /// [`BodyClass::Unknown`] is not gated — those bodies fall under
    /// the per-body override only.
    pub(super) trail_class_filter: TrailClassFilter,
    /// Per-body explicit overrides keyed by body index. `Some(true)`
    /// forces a trail on, `Some(false)` forces it off. Absent entries
    /// fall back to the authored-set + class default. Index validity
    /// follows the same swap_remove-aware pruning as `pinned_orbits`.
    pub(super) trail_per_body_override: std::collections::HashMap<usize, bool>,

    pub(super) place_preset: &'static BodyPreset,

    pub(super) trail: Option<TrailRenderer>,

    pub(super) follow_selected_body: bool,
    /// Active follow handover. Captured on click-to-focus; decays to
    /// `None` once the camera lands on the body. Layered on top of
    /// `follow_selected_body` so toggle-style call sites elsewhere
    /// (inspector button, Esc) keep their current shape.
    pub(super) follow_transition: Option<crate::app::camera::FollowTransition>,
    /// Template-supplied default-view hints. Consumed alongside
    /// [`pending_fit`](Self::pending_fit) on the next tick where bodies
    /// are loaded; takes precedence over bounding-sphere fit.
    pub(super) pending_camera_hints: Option<TemplateCameraHints>,
    /// Active scenario's orbital plane normal in world coords, unit
    /// length. Drives place-mode's drop plane so bodies created via
    /// click-to-place land on the same plane as the visible orbits
    /// (and the future grid). Defaults to `+Z` — the convention
    /// `state_from_elements` writes heliocentric ecliptic templates
    /// into. Updated when a template loads via its
    /// [`TemplateCameraHints::up`].
    pub(super) orbital_plane_up: glam::DVec3,
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
    pub(super) field_registry: apsis::domain::field::FieldRegistry,
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

    // ── Perturbation catalog ──────────────────────────────────────────────────
    /// User-facing list of available non-gravitational perturbations. Each
    /// entry holds an `enabled` flag; toggling it calls `apply_perturbations`
    /// which rebuilds and sends the active stack to the physics thread.
    pub(super) perturbation_catalog: Vec<crate::app::perturbation::PerturbationCatalogEntry>,

    // ── Render-loop diagnostics ───────────────────────────────────────────────
    pub(super) diagnostics: crate::app::diagnostics::Diagnostics,
    pub(super) show_fps_hud: bool,

    // ── 3D camera ─────────────────────────────────────────────────────────────
    pub(super) camera: crate::app::camera::OrbitCamera,
    pub(super) camera_input_config: crate::app::camera::input::CameraInputConfig,
    pub(super) show_camera_triad: bool,
}

impl SimulationApp {
    /// `true` while a Precision Run owns the simulation (state != Idle).
    /// UI surfaces that mutate physics state — canvas gestures,
    /// inspector fields, config sliders, clear/load actions —
    /// read this to render as disabled. The backend already drops
    /// such commands (see `PhysicsHandle::send`); the UI-level gate
    /// is just communication.
    pub(super) fn is_editing_locked(&self) -> bool {
        let ctrl = self.system.precision_controller();
        let guard = ctrl.lock().unwrap();
        !matches!(guard.state(), apsis::core::precision_run::RunState::Idle)
    }

    /// Short hint shown on hover over disabled edit controls.
    pub(super) fn editing_lock_hint(&self) -> &'static str {
        "Precision run in progress — editing is disabled until the run completes"
    }

    /// Rebuild and push the active perturbation stack to the physics thread.
    /// Call whenever an entry in `perturbation_catalog` changes.
    pub(super) fn apply_perturbations(&mut self) {
        let ps: Vec<Box<dyn apsis::physics::integrator::PerturbationForce>> = self
            .perturbation_catalog
            .iter()
            .filter(|e| e.enabled)
            .map(|e| e.descriptor.build())
            .collect();
        self.system.set_perturbations(ps);
    }

    pub fn new(system: System) -> Self {
        let mut physics_cfg = PhysicsConfig::default();
        physics_cfg.integrator = system.integrator_kind();
        physics_cfg.theta = system.theta();
        physics_cfg.softening_scale = system.softening_scale();

        let physics = spawn_physics(system, true /* start paused */);

        Self {
            system: physics,
            paused: true,
            semantic_scale_mode: SemanticScaleMode::Comparative,
            form: BodyForm::default(),
            form_error: None,
            show_trails: true,
            show_orbit_ellipses: false,
            orbit_visible_levels: [true, true, true, false],
            orbit_top_n: 24,
            orbit_hide_degenerate: true,
            orbit_hierarchy: apsis::physics::orbit_hierarchy::OrbitHierarchy::new(),
            orbit_smoother: crate::render::orbit_smoother::OrbitSmoother::new(),
            orbit_smoother_last_t: f64::NAN,
            pinned_orbits: HashSet::new(),
            show_grid: true,
            show_vectors: false,
            exposure_ev: 0.0,
            sim_rate_target: std::f64::consts::TAU,
            shortfall_active: false,
            ias15_epsilon: 1e-9,
            // Default: roughly one "internal year" at the default
            // unit system (G = 1 gives orbital period = 2π). Users can
            // override before starting a run.
            precision_run_duration: 2.0 * std::f64::consts::PI,
            precision_confirmation_pending: None,
            precision_confirmation_session_skip: false,
            notifications: {
                let store =
                    Arc::new(Mutex::new(crate::app::notifications::NotificationStore::new()));
                let _sub = crate::app::notifications::attach_to_bus(store.clone());
                store
            },
            show_notifications_panel: false,
            notifications_filter: crate::app::notifications::NotificationFilter::default(),
            place_mode: false,
            place_drag_start: None,
            place_mass: 1.0,
            place_density: 1.0,
            spawn_tab: SpawnTab::Single,
            spawn_ring_radius: 10.0,
            spawn_ring_count: 60,
            spawn_ring_mass: 0.01,
            spawn_ring_vel_scale: 1.0,
            spawn_ring_preset: &body_preset::ROCKY,
            spawn_cluster_radius: 5.0,
            spawn_cluster_count: 30,
            spawn_cluster_mass: 1.0,
            spawn_cluster_vel_disp: 0.5,
            spawn_cluster_preset: &body_preset::ROCKY,
            selection: BodySelection::default(),
            inspector_state: crate::app::inspector::InspectorState::default(),
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
            trail_recorder: apsis::core::trail::TrailRecorder::new(),
            trail_class_filter: TrailClassFilter::default(),
            trail_per_body_override: std::collections::HashMap::new(),
            place_preset: &body_preset::ROCKY,
            trail: None,

            follow_selected_body: false,
            follow_transition: None,
            pending_camera_hints: None,
            orbital_plane_up: glam::DVec3::Z,
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
            templates_meta: apsis::templates::TEMPLATES
                .iter()
                .map(|e| {
                    let t = e.build(0);
                    (t.description, t.body_count())
                })
                .collect(),

            sim_name: String::new(),
            show_name_prompt: false,
            pending_name_input: String::new(),
            sim_seed: SimSnapshot::new_seed(),

            field_registry: apsis::domain::field::FieldRegistry::standard(),
            colormap_registry: crate::render::color::ColormapRegistry::standard(),
            normalizer_registry: crate::render::color::NormalizerRegistry::standard(),
            color_view: None,
            color_view_range: None,

            perturbation_catalog: crate::app::perturbation::default_catalog(),

            diagnostics: crate::app::diagnostics::Diagnostics::new(),
            show_fps_hud: true,

            camera: crate::app::camera::OrbitCamera::new(crate::app::camera::CameraPose::new(
                glam::DVec3::ZERO,
                // azimuth = 0: eye in the world-XZ plane (no horizontal
                // rotation around the world-up axis). Combined with the
                // elevation below, the camera lands on the +Z half-space
                // — perpendicular to the XY orbital plane that
                // `state_from_elements` writes solar-system bodies into.
                0.0,
                // elevation ≈ 28°: dead-perpendicular (el = 0) gives a
                // pure top-down ecliptic projection; a small positive
                // tilt reads as "3D" without losing the orbital-plane
                // overview. Matches the default scene-load view of
                // NASA Eyes / Universe Sandbox / Solar System Scope.
                0.5,
                50.0,
            )),
            camera_input_config: crate::app::camera::input::CameraInputConfig::default(),
            show_camera_triad: true,
        }
    }

    fn draw_frame(&mut self, ui: &mut egui::Ui) {
        let ctx = ui.ctx().clone();

        design_theme::install(&ctx);

        // ── Sync latest physics state into local cache ────────────────────────
        self.system.sync();

        // Dense-output interpolation. Animation is gated on both the
        // real-time paused flag AND the Precision Run state so the
        // view freezes cleanly on Paused / Completed (avoiding the
        // teleport artefact when stale snapshots would otherwise be
        // advanced mid-pause). During an active Precision run the
        // sim_rate_target no longer drives physics — the observed
        // throughput does — so the render uses that value instead to
        // keep body motion visually matched to the physics thread.
        let (should_advance, rate) = {
            use apsis::core::precision_run::RunState;
            let ctrl = self.system.precision_controller();
            let guard = ctrl.lock().unwrap();
            match guard.state() {
                RunState::Running { .. } | RunState::Pausing { .. } | RunState::Aborting { .. } => {
                    (true, self.system.sim_rate())
                },
                RunState::Paused { .. } | RunState::Completed { .. } => (false, 0.0),
                RunState::Idle => (!self.paused, self.sim_rate_target),
            }
        };
        if should_advance {
            let wall_delta = ctx.input(|i| i.unstable_dt as f64).min(0.2);
            self.system.advance_render_time(wall_delta, rate);
        }

        let animating = should_advance || self.dragging_body.is_some();
        self.diagnostics.tick(animating);

        // Camera spring integration lives in `draw_canvas`, after the
        // follow loop has set this frame's `target`.

        // ── Single-step: re-pause after one frame of physics ─────────────────
        if self.step_pending {
            self.step_pending = false;
            self.paused = true;
        }

        // ── Pending view (template hints take precedence over fit) ────────────
        if self.pending_fit && !self.system.bodies().is_empty() && !self.system.is_loading() {
            self.apply_pending_view();
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
            == apsis::physics::integrator::traits::ExecutionProfile::Precision;
        let run_in_flight = {
            let ctrl = self.system.precision_controller();
            let guard = ctrl.lock().unwrap();
            !matches!(guard.state(), apsis::core::precision_run::RunState::Idle)
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
        self.draw_notifications_panel(&ctx);
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

        // Same swap_remove dance for the per-body trail override map: drop
        // the entry at the removed slot, then move the entry at `old_last`
        // (which physically migrates into `removed_idx`) so the explicit
        // override survives the body shuffle.
        self.trail_per_body_override.remove(&removed_idx);
        if old_last != removed_idx
            && let Some(v) = self.trail_per_body_override.remove(&old_last)
        {
            self.trail_per_body_override.insert(removed_idx, v);
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
                    // Drop any selected indices that landed in the removed range.
                    let first_removed = total.saturating_sub(n);
                    let new_sel = match &self.selection {
                        BodySelection::Single(sel) if *sel >= first_removed => {
                            Some(BodySelection::default())
                        },
                        BodySelection::Multi(set) => {
                            let valid: BTreeSet<usize> =
                                set.iter().copied().filter(|&i| i < first_removed).collect();
                            (valid.len() < set.len()).then(|| match valid.len() {
                                0 => BodySelection::default(),
                                1 => BodySelection::Single(*valid.iter().next().unwrap()),
                                _ => BodySelection::Multi(valid),
                            })
                        },
                        _ => None,
                    };
                    if let Some(sel) = new_sel {
                        self.selection = sel;
                        self.selection_form = None;
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
                    if self.selection.single() == Some(idx) {
                        self.selection_form = Some(SelectionForm::from_body(&old_body, &old_name));
                    }
                },
                UndoRecord::ReplacedBodies { previous } => {
                    self.system.load_named_bodies(previous);
                    self.pinned_orbits.clear();
                    self.selection = BodySelection::default();
                    self.selection_form = None;
                    self.pending_fit = true;
                    self.reset_drift_peaks();
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
