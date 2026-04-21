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
use crate::domain::body::{Body, NamedBody};
use crate::core::metrics::Metrics;
use crate::io::snapshot::SimSnapshot;
use crate::core::system::System;
use crate::physics::integrator::IntegratorKind;
use crate::physics::orbital::OrbitalElements;
use crate::render::trail::{TrailSampler, TrailSamplerKind};

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
    /// frame. The render-side [`TrailRecorder`](crate::render::TrailRecorder)
    /// reads and clears this each tick to keep trail positions aligned with
    /// the shifted body coordinate system.
    pub pending_com_shift: (f32, f32),
    /// Body-position columns sampled this tick by the trail sampler, in the
    /// order they were produced. Each entry has length `bodies.len()`. The
    /// UI thread drains this vector every sync and feeds it to the
    /// [`TrailRecorder`](crate::render::TrailRecorder).
    pub trail_samples: Vec<Vec<[f32; 2]>>,
    /// Gravitational accelerations from the last completed physics step,
    /// published so the render-side
    /// [`AccelerationMagnitudeField`](crate::domain::field::acceleration::AccelerationMagnitudeField)
    /// can paint bodies by |a| without needing to recompute forces.
    pub accelerations: Vec<(f64, f64)>,
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
    ZeroComVelocity,
    RestoreSnapshot(SimSnapshot),
    SetDtMode(DtMode),
    SetAdaptiveTheta(bool),
    /// Swap the trail-sampler strategy. Anchors (if any) are discarded; the
    /// next step restores them and records an initial sample.
    SetTrailSampler(TrailSamplerKind),
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
    accelerations: Vec<(f64, f64)>,
    /// Pending COM shift published by the physics thread this frame.
    /// Consumed by TrailRecorder on the UI thread.
    pending_com_shift: (f32, f32),
    /// Trail samples published by the physics thread this frame. Drained by
    /// the UI thread on [`sync`] and fed to the [`TrailRecorder`].
    pending_trail_samples: Vec<Vec<[f32; 2]>>,

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
            self.pending_trail_samples
                .extend(rs.trail_samples.drain(..));
        }
    }

    /// Returns and resets the accumulated COM shift since the last sync.
    ///
    /// Forwarded to [`crate::render::TrailRecorder::apply_com_shift`] each frame.
    pub fn take_pending_com_shift(&mut self) -> (f32, f32) {
        std::mem::replace(&mut self.pending_com_shift, (0.0, 0.0))
    }

    /// Drains the trail samples accumulated since the last call.
    ///
    /// Each returned column is a `Vec<[f32; 2]>` of length `bodies.len()`
    /// produced by the physics-thread's sampler. Consumed by
    /// [`TrailRecorder::ingest`](crate::render::TrailRecorder::ingest).
    pub fn take_trail_samples(&mut self) -> Vec<Vec<[f32; 2]>> {
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
    pub fn accelerations(&self) -> &[(f64, f64)] {
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
        // Unbounded channel: send never blocks or errors while thread is live.
        let _ = self.cmd_tx.send(cmd);
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
            seed: 0,                 // set by caller
            trail: None,             // set by caller
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
    };

    let render = Arc::new(Mutex::new(initial.clone()));
    let render_thr = render.clone();
    let loading = Arc::new(AtomicBool::new(false));
    let loading_thr = loading.clone();
    let orbital_needed = Arc::new(AtomicBool::new(true));
    let orbital_needed_thr = orbital_needed.clone();

    let thread = thread::spawn(move || {
        physics_loop(system, cmd_rx, render_thr, loading_thr, paused, orbital_needed_thr);
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
    trail_samples: &mut Vec<Vec<[f32; 2]>>,
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
    rs.trail_samples.extend(trail_samples.drain(..));
}

// ── Command dispatch ──────────────────────────────────────────────────────────

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
    let mut trail_samples_pending: Vec<Vec<[f32; 2]>> = Vec::new();
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
    let mut last_orbital = Instant::now()
        .checked_sub(orbital_slow)
        .unwrap_or_else(Instant::now);

    let pos_interval = Duration::from_millis(8);
    let mut last_pos = Instant::now();

    // Sim-rate measurement (rolling 500 ms window).
    let mut rate_wall = Instant::now();
    let mut rate_sim_acc = 0.0_f64;
    let mut current_sim_rate = 0.0_f64;

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

        if !paused {
            // Compute the sim-time we want to reach this batch.
            // wall_delta is clamped to 200 ms to avoid a spiral-of-death on
            // slow frames: if physics fell behind, we don't try to catch up.
            let wall_delta = prev_batch_start
                .elapsed()
                .min(Duration::from_millis(200))
                .as_secs_f64();
            prev_batch_start = batch_start;
            let t_target = system.t() + sim_rate_target * wall_delta;
            let hard_deadline = batch_start + Duration::from_millis(MAX_BATCH_WALL_MS);
            let mut steps_since_check = 0u32;
            let mut shutdown = false;

            'batch: loop {
                system.step();
                steps_since_check += 1;

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
                    let col: Vec<[f32; 2]> = system
                        .bodies()
                        .iter()
                        .map(|b| [b.x as f32, b.y as f32])
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
