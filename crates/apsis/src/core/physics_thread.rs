//! Background physics thread.
//!
//! Moves the N-body step loop off the UI thread so heavy simulations (large N,
//! many steps/frame) no longer block rendering.
//!
//! # Communication model
//!
//! ```text
//!  UI thread                        Physics thread
//!  ─────────────────────────────    ──────────────────────────────
//!  PhysicsHandle
//!    .sync()         ←─────────── Arc<Mutex<RenderState>>  (latest frame)
//!    .set_paused()   ──────────→  mpsc::Receiver<PhysicsCmd>
//!    .add_body()     ──────────→
//!    .to_snapshot()  ←─ block ──  mpsc::SyncSender<SimSnapshot>
//! ```
//!
//! The physics thread publishes a new [`RenderState`] at ~60 Hz regardless of
//! how fast the simulation is running. The UI calls [`PhysicsHandle::sync`] once
//! per frame to pull the latest snapshot into its local cache, then reads from
//! that cache through the same method names as the old `System` API — so
//! call-sites throughout the UI are unchanged.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::time::{Duration, Instant};

use crate::core::adaptive::DtMode;
use crate::core::metrics::Metrics;
use crate::core::precision_run::{PrecisionRunController, RunOutcome, RunState, TelemetryBuilder};
use crate::core::system::System;
use crate::core::trail::{TrailSampler, TrailSamplerKind};
use crate::domain::body::{Body, NamedBody};
use crate::io::snapshot::SimSnapshot;
use crate::math::Vec3;
use crate::physics::integrator::{IntegratorKind, PerturbationForce};
use crate::physics::orbital::OrbitalElements;

// ── Render state ──────────────────────────────────────────────────────────────

/// Complete simulation state published by the physics thread for the UI to read.
///
/// Held under `Arc<Mutex<…>>`. The physics thread locks briefly to overwrite the
/// whole struct; the UI locks briefly (or `try_lock`s) to clone what it needs.
#[derive(Clone)]
pub struct RenderState {
    pub bodies: Vec<Body>,
    pub names: Vec<String>,
    pub metrics: Metrics,
    pub orbital_elements: Vec<Option<OrbitalElements>>,
    pub softening_scale: f64,
    pub seed: u64,
    /// Simulation time units advanced per wall-second (rolling 500 ms window).
    /// Zero until the first measurement completes.
    pub sim_rate: f64,
    /// Accumulated world-space COM translation (render units) since the last
    /// frame. The render-side [`TrailRecorder`](crate::core::trail::TrailRecorder)
    /// reads and clears this each tick to keep trail positions aligned with
    /// the shifted body coordinate system.
    pub pending_com_shift: (f32, f32),
    /// Body-position columns sampled this tick by the trail sampler, in the
    /// order they were produced. Each entry has length `bodies.len()`. The
    /// UI thread drains this vector every sync and feeds it to the
    /// [`TrailRecorder`](crate::core::trail::TrailRecorder).
    pub trail_samples: Vec<Vec<[f32; 3]>>,
    /// Gravitational accelerations from the last completed physics step,
    /// published so the render-side
    /// [`AccelerationMagnitudeField`](crate::domain::field::acceleration::AccelerationMagnitudeField)
    /// can paint bodies by |a| without needing to recompute forces.
    pub accelerations: Vec<Vec3>,

    /// Dense-output snapshot from the most recently completed sub-step.
    /// The render thread uses this to interpolate body positions at any
    /// `t ∈ [t₀, t₀ + dt]` without re-running physics.
    pub step_snapshot: Option<crate::physics::integrator::DenseSnapshot>,
}

// ── Commands ──────────────────────────────────────────────────────────────────

/// Operations the UI sends to the physics thread.
///
/// The channel is unbounded so no command is ever dropped. Commands are drained
/// at the start of every physics iteration before stepping.
pub enum PhysicsCmd {
    SetPaused(bool),
    SetDt(f64),
    /// Target sim-time advance per real second (sim units/s).
    /// The physics thread integrates until `system.t() >= t_target` or the
    /// hard CPU cap (`MAX_BATCH_WALL_MS`) is exceeded — whichever comes first.
    SetSimRateTarget(f64),
    /// IAS15 error tolerance (ε). No-op for other integrators.
    SetIas15Epsilon(f64),
    SetExactThreshold(usize),
    SetSeed(u64),
    SetTheta(f64),
    SetSofteningScale(f64),
    SetGFactor(f64),
    SetIntegrator(IntegratorKind),
    AddBody(Body),
    AddNamedBody(NamedBody),
    AddBodies(Vec<Body>),
    AddNamedBodies(Vec<NamedBody>),
    RemoveBody(usize),
    UpdateBody(usize, Body),
    SetName(usize, String),
    LoadBodies(Vec<Body>),
    LoadNamedBodies(Vec<NamedBody>),
    ZeroComVelocity,
    RestoreSnapshot(SimSnapshot),
    SetDtMode(DtMode),
    SetAdaptiveTheta(bool),
    /// Swap the trail-sampler strategy. Anchors (if any) are discarded; the
    /// next step restores them and records an initial sample.
    SetTrailSampler(TrailSamplerKind),
    /// Replace the full perturbation stack atomically. Clears all previously
    /// registered forces, then registers each entry in order.
    SetPerturbations(Vec<Box<dyn PerturbationForce>>),
    Shutdown,
}

// ── PhysicsHandle ─────────────────────────────────────────────────────────────

/// Held by the UI. Mirrors the `System` read API so most call-sites are
/// unchanged; all mutations send commands to the background thread.
pub struct PhysicsHandle {
    cmd_tx: mpsc::Sender<PhysicsCmd>,
    /// Shared render state written by the physics thread.
    render: Arc<Mutex<RenderState>>,
    /// Set to `true` by the physics thread while processing a heavy command
    /// (load/restore). The UI reads this to show the loading overlay.
    loading: Arc<AtomicBool>,
    /// Controls whether orbital elements are updated at 60 Hz (true) or 0.5 Hz
    /// (false). Set to true only while the inspector panel is visible.
    orbital_elements_needed: Arc<AtomicBool>,

    // ── Local cache (updated by sync()) ──────────────────────────────────
    bodies: Vec<Body>,
    names: Vec<String>,
    metrics: Metrics,
    orbital_elements: Vec<Option<OrbitalElements>>,
    softening_scale: f64,
    sim_rate: f64,
    accelerations: Vec<Vec3>,
    /// Pending COM shift published by the physics thread this frame.
    /// Consumed by TrailRecorder on the UI thread.
    pending_com_shift: (f32, f32),
    /// Trail samples published by the physics thread this frame. Drained by
    /// the UI thread on [`sync`] and fed to the [`TrailRecorder`].
    pending_trail_samples: Vec<Vec<[f32; 3]>>,

    /// Latest dense-output snapshot, used by [`advance_render_time`](Self::advance_render_time).
    step_snapshot: Option<crate::physics::integrator::DenseSnapshot>,

    /// Current render time (absolute sim units). Advanced by
    /// [`advance_render_time`](Self::advance_render_time) each frame; always
    /// clamped to the window `[snap.t0, snap.t0 + snap.dt]`.
    t_render: f64,

    /// Shared Precision Run controller. The UI locks this to read the
    /// current run state / telemetry and to issue intent (start,
    /// request_pause, resume, request_abort, acknowledge). The physics
    /// thread locks it at each iteration to observe the state and to
    /// write `mark_paused` / `mark_aborted` / `mark_completed` at
    /// substep boundaries. See [`core::precision_run`](crate::core::precision_run).
    precision: Arc<Mutex<PrecisionRunController>>,

    _thread: thread::JoinHandle<()>,
}

impl PhysicsHandle {
    // ── Frame sync ────────────────────────────────────────────────────────

    /// Pull the latest physics state into the local cache.
    /// Non-blocking: if the physics thread currently holds the lock the cached
    /// values from the previous frame are kept.
    pub fn sync(&mut self) {
        if let Ok(mut rs) = self.render.try_lock() {
            self.bodies = rs.bodies.clone();
            self.names = rs.names.clone();
            self.metrics = rs.metrics;
            self.orbital_elements = rs.orbital_elements.clone();
            self.softening_scale = rs.softening_scale;
            self.sim_rate = rs.sim_rate;
            self.accelerations.clone_from(&rs.accelerations);
            // Drain the COM shift so the physics thread can start fresh.
            let shift = rs.pending_com_shift;
            rs.pending_com_shift = (0.0, 0.0);
            self.pending_com_shift = shift;
            // Drain accumulated trail samples (move; avoids allocation).
            self.pending_trail_samples.append(&mut rs.trail_samples);
            // Pull latest dense-output snapshot.
            if rs.step_snapshot.is_some() {
                self.step_snapshot = rs.step_snapshot.clone();
            }
        }
    }

    /// Advance the render time by `wall_delta` seconds (real wall clock) at the
    /// given sim-rate target, then overwrite the cached body positions with
    /// interpolated values from the latest [`DenseSnapshot`].
    ///
    /// Call this once per render frame, after [`sync`](Self::sync), while the
    /// simulation is running (skip when paused to freeze the display).
    pub fn advance_render_time(&mut self, wall_delta: f64, sim_rate_target: f64) {
        // Drop the frame rather than risk a panic. Three independent guards:
        // `dt > 0` (avoid divide-by-zero in `h`); `n_bodies()` matches body
        // count (avoid index past `x0`); shape-consistent snapshot (avoid
        // index past a shorter internal array; defence in depth against the
        // WH-class producer bug fixed in PR #14).
        let snap = match &self.step_snapshot {
            Some(s)
                if s.dt > 0.0 && s.n_bodies() == self.bodies.len() && s.is_shape_consistent() =>
            {
                s
            },
            _ => return,
        };

        self.t_render += sim_rate_target * wall_delta;
        // Clamp to the valid window — don't extrapolate beyond the accepted step.
        let t0 = snap.t0;
        let t1 = t0 + snap.dt;
        self.t_render = self.t_render.clamp(t0, t1);

        let h = (self.t_render - t0) / snap.dt;
        // Clone snap to release the borrow before mutating bodies.
        let snap = snap.clone();
        for (i, body) in self.bodies.iter_mut().enumerate() {
            let p = snap.interpolate(i, h);
            body.x = p.x;
            body.y = p.y;
            body.z = p.z;
        }
    }

    /// Returns and resets the accumulated COM shift since the last sync.
    ///
    /// Forwarded to [`crate::core::trail::TrailRecorder::apply_com_shift`] each frame.
    pub fn take_pending_com_shift(&mut self) -> (f32, f32) {
        std::mem::replace(&mut self.pending_com_shift, (0.0, 0.0))
    }

    /// Drains the trail samples accumulated since the last call.
    ///
    /// Each returned column is a `Vec<[f32; 3]>` of length `bodies.len()`
    /// produced by the physics-thread's sampler. Consumed by
    /// [`TrailRecorder::ingest`](crate::core::trail::TrailRecorder::ingest).
    pub fn take_trail_samples(&mut self) -> Vec<Vec<[f32; 3]>> {
        std::mem::take(&mut self.pending_trail_samples)
    }

    pub fn sim_rate(&self) -> f64 {
        self.sim_rate
    }

    /// Tell the physics thread whether it should compute orbital elements at
    /// full frame rate (true, inspector open) or at 0.5 Hz (false, inspector closed).
    pub fn set_orbital_elements_needed(&self, needed: bool) {
        self.orbital_elements_needed.store(needed, Ordering::Relaxed);
    }

    // ── Read methods — same signatures as System ──────────────────────────

    /// `true` while the physics thread is processing a heavy command (load / restore).
    /// The UI uses this to display the loading overlay.
    pub fn is_loading(&self) -> bool {
        self.loading.load(Ordering::Relaxed)
    }

    /// Clone of the shared Precision Run controller. The UI uses this
    /// to drive the run lifecycle: call `.lock().unwrap().start(t_target, t0)`
    /// to begin, `.request_pause()` / `.request_abort()` / `.resume()`
    /// for intent, `.state()` and `.telemetry()` to observe.
    ///
    /// The returned `Arc` can be cheaply cloned further if multiple UI
    /// subsystems need their own handle. All consumers see the same
    /// controller — this is the single source of truth for the run
    /// lifecycle, shared between the physics thread and every UI
    /// subscriber.
    pub fn precision_controller(&self) -> Arc<Mutex<PrecisionRunController>> {
        self.precision.clone()
    }

    pub fn bodies(&self) -> &[Body] {
        &self.bodies
    }
    pub fn names(&self) -> &[String] {
        &self.names
    }
    pub fn total_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }
    pub fn name(&self, idx: usize) -> &str {
        self.names.get(idx).map(|s| s.as_str()).unwrap_or("")
    }
    pub fn metrics(&self) -> Metrics {
        self.metrics
    }
    pub fn orbital_elements(&self) -> &[Option<OrbitalElements>] {
        &self.orbital_elements
    }
    /// Accelerations from the last completed physics step (one per body).
    pub fn accelerations(&self) -> &[Vec3] {
        &self.accelerations
    }
    pub fn t(&self) -> f64 {
        self.metrics.t
    }
    pub fn steps(&self) -> u64 {
        self.metrics.steps
    }
    pub fn dt(&self) -> f64 {
        self.metrics.dt
    }
    pub fn theta(&self) -> f64 {
        self.metrics.theta
    }
    pub fn integrator_kind(&self) -> IntegratorKind {
        self.metrics.integrator_kind
    }
    pub fn softening_scale(&self) -> f64 {
        self.softening_scale
    }
    pub fn g_factor(&self) -> f64 {
        self.metrics.g_factor
    }

    // ── Write methods — tunnel through command channel ────────────────────

    fn send(&self, cmd: PhysicsCmd) {
        // During an active Precision Run, physics-state mutations are
        // dropped rather than applied — the UI surfaces should be
        // disabled in this state so this backstop rarely fires, but
        // silently ignoring stray commands keeps the integration
        // deterministic even if a code path slips through. This
        // matches the REBOUND script model: while `reb_integrate` is
        // running, nothing external perturbs the system.
        if self.is_precision_run_active() && cmd_blocks_during_precision(&cmd) {
            return;
        }
        let _ = self.cmd_tx.send(cmd);
    }

    fn is_precision_run_active(&self) -> bool {
        !matches!(
            self.precision.lock().unwrap().state(),
            crate::core::precision_run::RunState::Idle
        )
    }

    pub fn set_paused(&self, paused: bool) {
        self.send(PhysicsCmd::SetPaused(paused));
    }
    /// Set the target sim-time advance per real second (sim units/s).
    /// Must be positive; values ≤ 0 are ignored.
    pub fn set_sim_rate_target(&self, rate: f64) {
        if rate > 0.0 {
            self.send(PhysicsCmd::SetSimRateTarget(rate));
        }
    }

    /// Set the IAS15 error tolerance. No-op for other integrators.
    pub fn set_ias15_epsilon(&self, eps: f64) {
        self.send(PhysicsCmd::SetIas15Epsilon(eps));
    }
    pub fn set_dt(&self, dt: f64) {
        self.send(PhysicsCmd::SetDt(dt));
    }
    pub fn set_exact_threshold(&self, n: usize) {
        self.send(PhysicsCmd::SetExactThreshold(n));
    }
    pub fn set_seed(&self, seed: u64) {
        self.send(PhysicsCmd::SetSeed(seed));
    }
    pub fn seed(&self) -> u64 {
        self.render.lock().unwrap().seed
    }
    pub fn set_theta(&self, theta: f64) {
        self.send(PhysicsCmd::SetTheta(theta));
    }
    pub fn set_softening_scale(&self, s: f64) {
        self.send(PhysicsCmd::SetSofteningScale(s));
    }
    pub fn set_g_factor(&self, g: f64) {
        self.send(PhysicsCmd::SetGFactor(g));
    }
    pub fn set_integrator(&self, kind: IntegratorKind) {
        self.send(PhysicsCmd::SetIntegrator(kind));
    }
    pub fn add_body(&self, body: Body) {
        self.send(PhysicsCmd::AddBody(body));
    }
    pub fn add_named_body(&self, named_body: NamedBody) {
        self.send(PhysicsCmd::AddNamedBody(named_body));
    }
    /// Add a batch of bodies in one operation. Sets the loading flag immediately
    /// so the overlay appears in the same frame as the user action.
    pub fn add_bodies(&self, bodies: Vec<Body>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::AddBodies(bodies));
    }
    /// Add a batch of named bodies while preserving authored template names.
    pub fn add_named_bodies(&self, bodies: Vec<NamedBody>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::AddNamedBodies(bodies));
    }
    pub fn remove_body(&self, idx: usize) {
        self.send(PhysicsCmd::RemoveBody(idx));
    }
    pub fn update_body(&self, idx: usize, body: Body) {
        self.send(PhysicsCmd::UpdateBody(idx, body));
    }
    pub fn set_name(&self, idx: usize, name: String) {
        self.send(PhysicsCmd::SetName(idx, name));
    }
    pub fn load_bodies(&self, bodies: Vec<Body>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::LoadBodies(bodies));
    }
    pub fn load_named_bodies(&self, bodies: Vec<NamedBody>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::LoadNamedBodies(bodies));
    }
    pub fn zero_com_velocity(&self) {
        self.send(PhysicsCmd::ZeroComVelocity);
    }
    pub fn set_dt_mode(&self, mode: DtMode) {
        self.send(PhysicsCmd::SetDtMode(mode));
    }
    pub fn set_adaptive_theta(&self, enabled: bool) {
        self.send(PhysicsCmd::SetAdaptiveTheta(enabled));
    }
    /// Install a new trail sampler. Takes effect on the next physics step.
    pub fn set_trail_sampler(&self, kind: TrailSamplerKind) {
        self.send(PhysicsCmd::SetTrailSampler(kind));
    }
    pub fn set_perturbations(&self, ps: Vec<Box<dyn PerturbationForce>>) {
        self.send(PhysicsCmd::SetPerturbations(ps));
    }

    pub fn restore_from_snapshot(&self, snap: &SimSnapshot) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::RestoreSnapshot(snap.clone()));
    }

    /// Build a snapshot from the locally cached render state.
    ///
    /// This is non-blocking and always returns valid data, using the body
    /// positions and metrics that were last synced from the physics thread.
    /// At worst the data is one frame stale — completely acceptable for saves.
    ///
    /// (The old blocking-RPC path via `RequestSnapshot` was removed because it
    /// could time out when the physics thread was mid-batch, producing an empty
    /// snapshot with 0 bodies / 0 steps.)
    pub fn to_snapshot(&self) -> SimSnapshot {
        use crate::io::snapshot::BodyRecord;
        let m = self.metrics;
        SimSnapshot {
            save_id: 0, // set by caller
            t: m.t,
            steps: m.steps,
            dt: m.dt,
            theta: m.theta,
            softening_scale: self.softening_scale,
            g_factor: m.g_factor,
            integrator_kind: m.integrator_kind,
            trail_every: 1, // trail_every now owned by TrailRecorder; app sets on snapshot
            sim_name: String::new(), // set by caller
            seed: 0,        // set by caller
            trail: None,    // set by caller
            bodies: self.bodies.iter().map(BodyRecord::from_body).collect(),
            names: self.names.clone(),
        }
    }
}

// ── spawn ─────────────────────────────────────────────────────────────────────

/// Spawn the physics thread and return a handle for the UI.
///
/// `paused` sets the initial run state; the UI can change it later via
/// [`PhysicsHandle::set_paused`].
pub fn spawn(mut system: System, paused: bool) -> PhysicsHandle {
    let (cmd_tx, cmd_rx) = mpsc::channel();

    // Snapshot initial state before moving system into the thread.
    system.update_orbital_elements();
    let initial = RenderState {
        bodies: system.bodies().to_vec(),
        names: system.names().to_vec(),
        metrics: system.metrics(),
        orbital_elements: system.orbital_elements().to_vec(),
        softening_scale: system.softening_scale(),
        seed: system.seed(),
        sim_rate: 0.0,
        pending_com_shift: (0.0, 0.0),
        trail_samples: Vec::new(),
        accelerations: system.last_accelerations().to_vec(),
        step_snapshot: None,
    };

    let render = Arc::new(Mutex::new(initial.clone()));
    let render_thr = render.clone();
    let loading = Arc::new(AtomicBool::new(false));
    let loading_thr = loading.clone();
    let orbital_needed = Arc::new(AtomicBool::new(true));
    let orbital_needed_thr = orbital_needed.clone();
    let precision = Arc::new(Mutex::new(PrecisionRunController::new()));
    let precision_thr = precision.clone();

    let thread = thread::spawn(move || {
        physics_loop(
            system,
            cmd_rx,
            render_thr,
            loading_thr,
            paused,
            orbital_needed_thr,
            precision_thr,
        );
    });

    PhysicsHandle {
        cmd_tx,
        render,
        loading,
        orbital_elements_needed: orbital_needed,
        bodies: initial.bodies,
        names: initial.names,
        metrics: initial.metrics,
        orbital_elements: initial.orbital_elements,
        softening_scale: initial.softening_scale,
        sim_rate: 0.0,
        accelerations: initial.accelerations.clone(),
        pending_com_shift: (0.0, 0.0),
        pending_trail_samples: Vec::new(),
        step_snapshot: None,
        t_render: 0.0,
        precision,
        _thread: thread,
    }
}

// ── Physics loop ──────────────────────────────────────────────────────────────

// ── Publish helpers ───────────────────────────────────────────────────────────

/// Cheap publish: only body positions + metrics.
///
/// Uses `copy_from_slice` when N is unchanged (Body is Copy) to avoid
/// heap allocation on the hot path.
fn publish_positions(system: &System, rs: &mut RenderState) {
    let src = system.bodies();
    if rs.bodies.len() == src.len() {
        rs.bodies.copy_from_slice(src);
    } else {
        rs.bodies = src.to_vec();
    }
    rs.metrics = system.metrics();
}

/// Full publish: positions + names + orbital elements + config.
fn publish_full(
    system: &mut System,
    rs: &mut RenderState,
    com_shift_acc: (f32, f32),
    trail_samples: &mut Vec<Vec<[f32; 3]>>,
) {
    publish_positions(system, rs);
    rs.names = system.names().to_vec();
    rs.orbital_elements = system.orbital_elements().to_vec();
    rs.softening_scale = system.softening_scale();
    rs.seed = system.seed();
    let src_acc = system.last_accelerations();
    if rs.accelerations.len() == src_acc.len() {
        rs.accelerations.copy_from_slice(src_acc);
    } else {
        rs.accelerations = src_acc.to_vec();
    }
    // Accumulate COM shift; the UI thread drains it on sync().
    rs.pending_com_shift.0 += com_shift_acc.0;
    rs.pending_com_shift.1 += com_shift_acc.1;
    // Hand off accumulated trail samples — move, don't copy.
    rs.trail_samples.append(trail_samples);
}

// ── Command dispatch ──────────────────────────────────────────────────────────

/// Whether a command must be ignored while a Precision Run is
/// active. Matches the REBOUND mental model: during
/// `reb_integrate`, the simulation is not externally perturbable.
/// Physics-state mutations are blocked; render-only / cosmetic /
/// lifecycle commands still pass through.
fn cmd_blocks_during_precision(cmd: &PhysicsCmd) -> bool {
    match cmd {
        PhysicsCmd::SetDt(_)
        | PhysicsCmd::SetIas15Epsilon(_)
        | PhysicsCmd::SetExactThreshold(_)
        | PhysicsCmd::SetSeed(_)
        | PhysicsCmd::SetTheta(_)
        | PhysicsCmd::SetSofteningScale(_)
        | PhysicsCmd::SetGFactor(_)
        | PhysicsCmd::SetIntegrator(_)
        | PhysicsCmd::SetDtMode(_)
        | PhysicsCmd::SetAdaptiveTheta(_)
        | PhysicsCmd::AddBody(_)
        | PhysicsCmd::AddNamedBody(_)
        | PhysicsCmd::AddBodies(_)
        | PhysicsCmd::AddNamedBodies(_)
        | PhysicsCmd::RemoveBody(_)
        | PhysicsCmd::UpdateBody(_, _)
        | PhysicsCmd::LoadBodies(_)
        | PhysicsCmd::LoadNamedBodies(_)
        | PhysicsCmd::ZeroComVelocity
        | PhysicsCmd::RestoreSnapshot(_)
        | PhysicsCmd::SetPerturbations(_) => true,

        PhysicsCmd::SetPaused(_)
        | PhysicsCmd::SetSimRateTarget(_)
        | PhysicsCmd::SetName(_, _)
        | PhysicsCmd::SetTrailSampler(_)
        | PhysicsCmd::Shutdown => false,
    }
}

/// Effect of applying a single [`PhysicsCmd`]. Returned by [`apply_cmd`] so the
/// caller can exit the thread cleanly on shutdown.
enum CmdEffect {
    Continue,
    Shutdown,
}

/// Applies one command to the physics-thread state. Extracted so the physics
/// loop can drain the command channel from multiple points: once at the top
/// of each batch (for responsiveness while paused), and again mid-batch (so
/// `SetPaused`, `RemoveBody`, `Shutdown` etc. are not delayed by a long batch).
fn apply_cmd(
    cmd: PhysicsCmd,
    system: &mut System,
    paused: &mut bool,
    sim_rate_target: &mut f64,
    needs_full_publish: &mut bool,
    trail_sampler: &mut Box<dyn TrailSampler>,
    loading: &Arc<AtomicBool>,
) -> CmdEffect {
    match cmd {
        PhysicsCmd::Shutdown => return CmdEffect::Shutdown,
        PhysicsCmd::SetPaused(p) => *paused = p,
        PhysicsCmd::SetSimRateTarget(rate) => *sim_rate_target = rate.max(1e-9),
        PhysicsCmd::SetIas15Epsilon(eps) => system.set_ias15_epsilon(eps),
        PhysicsCmd::SetDt(dt) => {
            system.set_dt(dt);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetExactThreshold(n) => {
            system.set_exact_threshold(n);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetSeed(s) => {
            system.set_seed(s);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetTheta(theta) => {
            system.set_theta(theta);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetSofteningScale(s) => {
            system.set_softening_scale(s);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetGFactor(g) => {
            system.set_g_factor(g);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetIntegrator(i) => {
            system.set_integrator(i);
            *needs_full_publish = true;
        },
        PhysicsCmd::AddBody(b) => {
            system.add_body(b);
            *needs_full_publish = true;
        },
        PhysicsCmd::AddNamedBody(named_body) => {
            system.add_named_body(named_body);
            *needs_full_publish = true;
        },
        PhysicsCmd::AddBodies(bodies) => {
            loading.store(true, Ordering::Relaxed);
            system.add_bodies(bodies);
            loading.store(false, Ordering::Relaxed);
            *needs_full_publish = true;
        },
        PhysicsCmd::AddNamedBodies(bodies) => {
            loading.store(true, Ordering::Relaxed);
            system.add_named_bodies(bodies);
            loading.store(false, Ordering::Relaxed);
            *needs_full_publish = true;
        },
        PhysicsCmd::RemoveBody(idx) => {
            system.remove_body(idx);
            *needs_full_publish = true;
        },
        PhysicsCmd::UpdateBody(idx, b) => {
            system.update_body(idx, b);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetName(idx, name) => {
            system.set_name(idx, name);
            *needs_full_publish = true;
        },
        PhysicsCmd::LoadBodies(bodies) => {
            loading.store(true, Ordering::Relaxed);
            system.load_bodies(bodies);
            loading.store(false, Ordering::Relaxed);
            *needs_full_publish = true;
        },
        PhysicsCmd::LoadNamedBodies(bodies) => {
            loading.store(true, Ordering::Relaxed);
            system.load_named_bodies(bodies);
            loading.store(false, Ordering::Relaxed);
            *needs_full_publish = true;
        },
        PhysicsCmd::ZeroComVelocity => {
            system.zero_com_velocity();
            *needs_full_publish = true;
        },
        PhysicsCmd::SetDtMode(mode) => {
            system.set_dt_mode(mode);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetAdaptiveTheta(enabled) => {
            system.set_adaptive_theta(enabled);
            *needs_full_publish = true;
        },
        PhysicsCmd::SetTrailSampler(kind) => {
            *trail_sampler = kind.build();
        },
        PhysicsCmd::SetPerturbations(ps) => {
            system.clear_perturbations();
            for p in ps {
                system.add_perturbation(p);
            }
            *needs_full_publish = true;
        },
        PhysicsCmd::RestoreSnapshot(snap) => {
            loading.store(true, Ordering::Relaxed);
            system.restore_from_snapshot(&snap);
            loading.store(false, Ordering::Relaxed);
            *needs_full_publish = true;
        },
    }
    CmdEffect::Continue
}

// ── Physics loop ──────────────────────────────────────────────────────────────

fn physics_loop(
    mut system: System,
    cmd_rx: mpsc::Receiver<PhysicsCmd>,
    render: Arc<Mutex<RenderState>>,
    loading: Arc<AtomicBool>,
    initial_paused: bool,
    orbital_needed: Arc<AtomicBool>,
    precision: Arc<Mutex<PrecisionRunController>>,
) {
    let mut paused = initial_paused;
    // Target sim-time advance per real second (sim units/s).
    // Default: 2π ≈ 1 yr/s in the internal unit system.
    let mut sim_rate_target: f64 = std::f64::consts::TAU;
    let mut needs_full_publish = false;

    // Accumulated COM shift across the current step batch. Published to
    // RenderState once per frame; consumed by TrailRecorder on the UI thread.
    let mut com_shift_acc: (f32, f32) = (0.0, 0.0);

    // Trail sampler: decides — per physics step — whether to enqueue a
    // trail sample for the UI. Default is arc-length based (see
    // [`TrailSamplerKind::default`]). Swapped via PhysicsCmd::SetTrailSampler.
    let mut trail_sampler: Box<dyn TrailSampler> = TrailSamplerKind::default().build();
    // Samples accumulated since the last publish_full; swapped into the
    // RenderState when that happens. Each column has length `bodies.len()`.
    let mut trail_samples_pending: Vec<Vec<[f32; 3]>> = Vec::new();
    // Safety cap: limit how many samples can accumulate between publishes
    // so a runaway frame cannot produce an unbounded backlog. Sized to
    // cover ~6 full orbits at the default arc-length density (314/orbit),
    // which is well above what the trail ring buffer displays anyway.
    const TRAIL_SAMPLES_PER_BATCH_MAX: usize = 2048;

    let full_interval = Duration::from_millis(16);
    let mut last_full = Instant::now().checked_sub(full_interval).unwrap_or_else(Instant::now);
    // Orbital elements are O(N²) — update at full rate only when inspector is open.
    let orbital_fast = Duration::from_millis(16);
    let orbital_slow = Duration::from_secs(2);
    let mut last_orbital = Instant::now().checked_sub(orbital_slow).unwrap_or_else(Instant::now);

    let pos_interval = Duration::from_millis(8);
    let mut last_pos = Instant::now();

    // Sim-rate measurement (rolling 500 ms window).
    let mut rate_wall = Instant::now();
    let mut rate_sim_acc = 0.0_f64;
    let mut current_sim_rate = 0.0_f64;

    // Precision Run baseline — the `SimSnapshot` captured at the
    // moment a run enters `Running` for the first time. Held entirely
    // within this function; the UI never sees it. Purpose: `Abort`
    // restores the simulation to this exact pre-run state, so an
    // aborted precision run leaves no residue in `System`.
    //
    // Lifecycle:
    //   * `None` when no run is in progress.
    //   * `Some(snapshot)` from the first observation of
    //     `RunState::Running` until the run completes (`mark_completed`,
    //     `mark_aborted`, or `mark_errored`).
    //   * Cleared on `Reached` and `Errored` without use.
    //   * Consumed on `UserAborted` — `system.restore_from_snapshot`
    //     puts the bodies back to the pre-run state.
    //
    // `SimSnapshot` is an owned deep copy of body state (and integrator
    // configuration), not a reference into `System`. Holding it here
    // during a run is safe and does not alias mutable state.
    let mut precision_baseline: Option<SimSnapshot> = None;

    const POS_CHECK_STEPS: u32 = 8;
    // Hard CPU cap per batch. The sim-rate target is the primary control;
    // this cap only fires when physics cannot keep up (large N, tiny dt, high target).
    const MAX_BATCH_WALL_MS: u64 = 500;
    let min_batch_period = Duration::from_micros(100);
    let mut prev_batch_start = Instant::now();

    loop {
        let batch_start = Instant::now();

        while let Ok(cmd) = cmd_rx.try_recv() {
            if matches!(
                apply_cmd(
                    cmd,
                    &mut system,
                    &mut paused,
                    &mut sim_rate_target,
                    &mut needs_full_publish,
                    &mut trail_sampler,
                    &loading,
                ),
                CmdEffect::Shutdown
            ) {
                return;
            }
        }

        // ── Precision Run branch ─────────────────────────────────────────────
        //
        // When a precision run is active, its own run-to-completion logic
        // supersedes the sim-rate-target loop below. The two modes do not
        // interleave within a single iteration: `step_precision_run_tick`
        // returns once the run pauses, aborts, completes, or consumes its
        // per-tick budget, at which point we fall through to the normal
        // publish step.
        let entry_state = precision.lock().unwrap().state();
        let precision_active = matches!(
            entry_state,
            RunState::Running { .. } | RunState::Pausing { .. } | RunState::Aborting { .. }
        );

        // Capture baseline on first observation of `Running` in this run.
        // See `precision_baseline` declaration for lifecycle semantics.
        if precision_baseline.is_none() && matches!(entry_state, RunState::Running { .. }) {
            precision_baseline = Some(system.to_snapshot());
            crate::info_diag!(
                crate::core::log::Source::PhysicsThread,
                "precision run: baseline captured for abort rollback",
                t = system.t(),
                steps = system.steps(),
            );
        }

        if precision_active {
            step_precision_run_tick(&mut system, &precision);

            // Observe post-tick outcome. If abort landed during the tick,
            // restore the baseline *before* publishing so the UI never
            // sees the mid-run state after the user clicked Abort.
            // Commit/Reached and Errored drop the baseline unused.
            let post_state = precision.lock().unwrap().state();
            match post_state {
                RunState::Completed { outcome: RunOutcome::UserAborted } => {
                    if let Some(snap) = precision_baseline.take() {
                        system.restore_from_snapshot(&snap);
                        crate::info_diag!(
                            crate::core::log::Source::PhysicsThread,
                            "precision run aborted: state restored to baseline",
                            t = system.t(),
                            steps = system.steps(),
                        );
                        needs_full_publish = true;
                    }
                },
                RunState::Completed { outcome: RunOutcome::Reached } => {
                    precision_baseline = None;
                    crate::info_diag!(
                        crate::core::log::Source::PhysicsThread,
                        "precision run completed: target simulation time reached",
                        t = system.t(),
                        steps = system.steps(),
                    );
                },
                RunState::Completed { outcome: RunOutcome::Errored } => {
                    // Defensive: drop the baseline rather than leave
                    // stale state around. Today no producer emits
                    // Errored, so this branch is exercised by tests
                    // only; leaving it explicit keeps the state
                    // machine exhaustive.
                    precision_baseline = None;
                },
                _ => {},
            }

            if let Ok(mut rs) = render.try_lock() {
                publish_full(&mut system, &mut rs, com_shift_acc, &mut trail_samples_pending);
                rs.sim_rate = current_sim_rate;
            }
            com_shift_acc = (0.0, 0.0);
            // Short sleep so the UI render thread has a chance to
            // grab the render lock between ticks; precision mode
            // is not bound by sim_rate_target so we do not spin.
            thread::sleep(Duration::from_micros(200));
            continue;
        }

        if !paused {
            // Compute the sim-time we want to reach this batch.
            // wall_delta is clamped to 200 ms to avoid a spiral-of-death on
            // slow frames: if physics fell behind, we don't try to catch up.
            let wall_delta =
                prev_batch_start.elapsed().min(Duration::from_millis(200)).as_secs_f64();
            prev_batch_start = batch_start;
            let t_target = system.t() + sim_rate_target * wall_delta;
            let hard_deadline = batch_start + Duration::from_millis(MAX_BATCH_WALL_MS);
            let mut steps_since_check = 0u32;
            let mut shutdown = false;

            // Hand the batch deadline to adaptive integrators so their
            // retry loop can cooperate with the batch budget rather than
            // spinning past it. Fixed-step integrators ignore this.
            system.set_step_deadline(Some(hard_deadline));

            'batch: loop {
                system.step();
                steps_since_check += 1;

                // Forward the fresh dense-output snapshot to RenderState so the
                // render thread always has the most recent sub-step to interpolate.
                // try_lock is non-blocking; missing one snapshot frame is harmless.
                if system.last_dense_snapshot.is_some() {
                    if let Ok(mut rs) = render.try_lock() {
                        rs.step_snapshot = system.last_dense_snapshot.clone();
                    }
                }

                let dt = system.metrics().dt;
                rate_sim_acc += dt;

                // Drain any COM shift the integrator applied this step.
                let (sx, sy) = system.take_com_shift();
                com_shift_acc.0 += sx;
                com_shift_acc.1 += sy;

                // Arc-length / step-count sampling at physics-step granularity.
                // Capture position columns into a local queue; publish_full
                // hands them off to the UI thread in one move.
                if trail_samples_pending.len() < TRAIL_SAMPLES_PER_BATCH_MAX
                    && trail_sampler.should_sample(system.bodies())
                {
                    let col: Vec<[f32; 3]> = system
                        .bodies()
                        .iter()
                        .map(|b| [b.x as f32, b.y as f32, b.z as f32])
                        .collect();
                    trail_samples_pending.push(col);
                }

                if steps_since_check >= POS_CHECK_STEPS {
                    steps_since_check = 0;

                    let now = Instant::now();

                    if now.duration_since(last_full) >= full_interval {
                        let orb_interval = if orbital_needed.load(Ordering::Relaxed) {
                            orbital_fast
                        } else {
                            orbital_slow
                        };
                        if now.duration_since(last_orbital) >= orb_interval {
                            system.update_orbital_elements();
                            last_orbital = now;
                        }

                        if let Ok(mut rs) = render.try_lock() {
                            publish_full(
                                &mut system,
                                &mut rs,
                                com_shift_acc,
                                &mut trail_samples_pending,
                            );
                            rs.sim_rate = current_sim_rate;
                        }
                        com_shift_acc = (0.0, 0.0);

                        last_full = now;
                        last_pos = now;
                        needs_full_publish = false;
                    } else if now.duration_since(last_pos) >= pos_interval {
                        last_pos = now;

                        if let Ok(mut rs) = render.try_lock() {
                            publish_positions(&system, &mut rs);
                        }
                    }

                    // Mid-batch command drain.
                    let was_paused = paused;
                    while let Ok(cmd) = cmd_rx.try_recv() {
                        if matches!(
                            apply_cmd(
                                cmd,
                                &mut system,
                                &mut paused,
                                &mut sim_rate_target,
                                &mut needs_full_publish,
                                &mut trail_sampler,
                                &loading,
                            ),
                            CmdEffect::Shutdown
                        ) {
                            shutdown = true;
                            break;
                        }
                    }
                    if shutdown
                        || (paused && !was_paused)
                        || system.t() >= t_target
                        || now >= hard_deadline
                    {
                        break 'batch;
                    }
                }
            }

            // Clear the per-batch deadline so any step taken outside the
            // batch loop (single-step command, template load side-effects)
            // does not inherit a stale one.
            system.set_step_deadline(None);

            if shutdown {
                return;
            }

            // Update sim-rate measurement every 500 ms.
            let wall_elapsed = rate_wall.elapsed();
            if wall_elapsed >= Duration::from_millis(500) {
                current_sim_rate = rate_sim_acc / wall_elapsed.as_secs_f64();
                rate_sim_acc = 0.0;
                rate_wall = Instant::now();
            }
        }

        // ── Post-batch full publish ────────────────────────────────
        let now = Instant::now();

        if needs_full_publish || (!paused && now.duration_since(last_full) >= full_interval) {
            let orb_interval =
                if orbital_needed.load(Ordering::Relaxed) { orbital_fast } else { orbital_slow };
            if now.duration_since(last_orbital) >= orb_interval {
                system.update_orbital_elements();
                last_orbital = now;
            }

            if let Ok(mut rs) = render.try_lock() {
                publish_full(&mut system, &mut rs, com_shift_acc, &mut trail_samples_pending);
                rs.sim_rate = current_sim_rate;
            }
            com_shift_acc = (0.0, 0.0);

            last_full = now;
            last_pos = now;
            needs_full_publish = false;
        }

        // ── Throttle ──────────────────────────────────────────────
        if paused {
            thread::sleep(Duration::from_millis(8));
        } else {
            let elapsed = batch_start.elapsed();

            if elapsed < min_batch_period {
                thread::sleep(min_batch_period - elapsed);
            }
        }
    }
}

// ── Precision Run tick ────────────────────────────────────────────────────────
//
// Drives one "tick" of a precision run: observes UI intent at the substep
// boundary (pause/abort), steps the system until `t_target` or the per-tick
// budget is reached, and pushes a fresh telemetry snapshot back into the
// controller. The outer physics loop calls this instead of the real-time
// batch loop whenever the controller is in `Running`, `Pausing`, or
// `Aborting`.
//
// # Budget
//
// `TICK_BUDGET` bounds how long one call can block the physics thread
// without yielding. It is NOT a simulation-rate constraint — precision
// runs ignore `sim_rate_target` on purpose. It exists so the render
// thread can regularly acquire `render` / `precision` locks (progress
// bar at ~20 Hz) and so pause/abort intent is observed promptly.
//
// # Why intent is observed between sub-steps
//
// `PrecisionRunController::request_pause` flips state to `Pausing` from
// the UI thread. The physics thread observes at the top of its
// per-substep loop and calls `mark_paused` *only after* the substep
// currently in flight (none at entry) completes — i.e., at a clean
// boundary. This preserves determinism: no state change ever lands
// mid-Picard-iteration.
fn step_precision_run_tick(system: &mut System, precision: &Arc<Mutex<PrecisionRunController>>) {
    const TICK_BUDGET: Duration = Duration::from_millis(50);

    // Snapshot entry state. If we are in Pausing / Aborting, the "current
    // substep" is already complete (we are between iterations), so the
    // confirmation fires immediately. Otherwise, extract the run bounds.
    let (t_target, t_start, prev_peak_err) = {
        let ctrl = precision.lock().unwrap();
        match ctrl.state() {
            RunState::Pausing { .. } => {
                drop(ctrl);
                precision.lock().unwrap().mark_paused();
                return;
            },
            RunState::Aborting { .. } => {
                drop(ctrl);
                precision.lock().unwrap().mark_aborted();
                return;
            },
            RunState::Running { t_target, t_start, .. } => {
                (t_target, t_start, ctrl.telemetry().peak_energy_err)
            },
            // Defensive: caller (`physics_loop`) only invokes this while the
            // controller is in one of the above states. Any other state is a
            // benign race and we just return.
            _ => return,
        }
    };

    // Clear any step deadline left over from real-time batches — precision
    // runs must not self-terminate on the interactive CPU cap.
    system.set_step_deadline(None);

    let tick_start = Instant::now();
    let substeps_at_tick_start =
        system.adaptive_stats().map(|s| s.substeps).unwrap_or_else(|| system.steps());
    let mut running_peak_err = prev_peak_err;

    while system.t() < t_target && tick_start.elapsed() < TICK_BUDGET {
        // Poll intent between substeps. Lock is held only for the enum read;
        // the UI thread can request pause/abort concurrently.
        match precision.lock().unwrap().state() {
            RunState::Running { .. } => {}, // proceed
            RunState::Pausing { .. } => {
                precision.lock().unwrap().mark_paused();
                break;
            },
            RunState::Aborting { .. } => {
                precision.lock().unwrap().mark_aborted();
                break;
            },
            _ => break,
        }

        system.step();

        let err_abs = system.rel_energy_error().abs();
        if err_abs > running_peak_err {
            running_peak_err = err_abs;
        }
    }

    // Compute telemetry once per tick.
    let wall_delta = tick_start.elapsed().as_secs_f64().max(1e-9);
    let substeps_now =
        system.adaptive_stats().map(|s| s.substeps).unwrap_or_else(|| system.steps());
    let substeps_in_tick = substeps_now.saturating_sub(substeps_at_tick_start);
    let substeps_per_s = (substeps_in_tick as f64) / wall_delta;
    let current_dt = system.dt();
    let sim_time_delta = substeps_in_tick as f64 * current_dt;
    let sim_time_per_s = sim_time_delta / wall_delta;

    let progress = {
        let span = (t_target - t_start).max(f64::MIN_POSITIVE);
        ((system.t() - t_start) / span).clamp(0.0, 1.0) as f32
    };

    let mut builder = TelemetryBuilder::new()
        .with_current_dt(current_dt)
        .with_substeps_per_second(substeps_per_s)
        .with_sim_time_per_second(sim_time_per_s)
        .with_progress_fraction(progress)
        .with_peak_energy_err(running_peak_err)
        .with_current_energy_err(system.rel_energy_error());

    if let Some(stats) = system.adaptive_stats() {
        builder = builder
            .with_substeps(stats.substeps)
            .with_rejections_picard(stats.rejections_picard)
            .with_rejections_truncation(stats.rejections_truncation)
            .with_picard_iters(stats.picard_iters)
            .with_degraded(stats.degraded);
    }

    let telemetry = builder.finish();
    let reached_target = system.t() >= t_target;

    let mut ctrl = precision.lock().unwrap();
    ctrl.update_telemetry(telemetry);
    if reached_target && matches!(ctrl.state(), RunState::Running { .. }) {
        ctrl.mark_completed();
    }
}
