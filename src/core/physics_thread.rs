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
use crate::render::trail_buffer::TrailBuffer;
use crate::physics::integrator::IntegratorKind;
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
    pub trail_buf: Arc<TrailBuffer>,
    pub metrics: Metrics,
    pub orbital_elements: Vec<Option<OrbitalElements>>,
    pub softening_scale: f64,
    pub trail_every: usize,
}

// ── Commands ──────────────────────────────────────────────────────────────────

/// Operations the UI sends to the physics thread.
///
/// The channel is unbounded so no command is ever dropped. Commands are drained
/// at the start of every physics iteration before stepping.
pub enum PhysicsCmd {
    SetPaused(bool),
    SetDt(f64),
    SetStepsPerFrame(u32),
    SetTheta(f64),
    SetSofteningScale(f64),
    SetGFactor(f64),
    SetIntegrator(IntegratorKind),
    SetTrailEvery(usize),
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

    // ── Local cache (updated by sync()) ──────────────────────────────────
    bodies: Vec<Body>,
    names: Vec<String>,
    metrics: Metrics,
    orbital_elements: Vec<Option<OrbitalElements>>,
    softening_scale: f64,
    trail_every: usize,
    /// Last successfully cloned trail buffer.
    ///
    /// `clone_trail_buf` tries a fresh `try_lock` first; if the physics thread
    /// currently holds the mutex it falls back here instead of returning an
    /// empty buffer (which caused per-frame trail flicker).
    cached_trail_buf: Arc<TrailBuffer>,

    _thread: thread::JoinHandle<()>,
}

impl PhysicsHandle {
    // ── Frame sync ────────────────────────────────────────────────────────

    /// Pull the latest physics state into the local cache.
    /// Non-blocking: if the physics thread currently holds the lock the cached
    /// values from the previous frame are kept.
    pub fn sync(&mut self) {
        if let Ok(rs) = self.render.try_lock() {
            self.bodies = rs.bodies.clone();
            self.names = rs.names.clone();
            self.metrics = rs.metrics;
            self.orbital_elements = rs.orbital_elements.clone();
            self.softening_scale = rs.softening_scale;
            self.trail_every = rs.trail_every;
            self.cached_trail_buf = rs.trail_buf.clone();
        }
    }

    /// Clone the trail buffer for the current frame.
    ///
    /// Tries a fresh `try_lock` for the lowest-latency snapshot.  If the
    /// physics thread currently holds the mutex, returns the cached copy from
    /// the last successful [`sync`] instead of an empty buffer — this
    /// eliminates the per-frame trail flicker that occurred on lock contention.
    pub fn clone_trail_buf(&self) -> Arc<TrailBuffer> {
        self.render
            .try_lock()
            .map(|rs| rs.trail_buf.clone())
            .unwrap_or_else(|_| self.cached_trail_buf.clone())
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
    pub fn trail_every(&self) -> usize {
        self.trail_every
    }

    // ── Write methods — tunnel through command channel ────────────────────

    fn send(&self, cmd: PhysicsCmd) {
        // Unbounded channel: send never blocks or errors while thread is live.
        let _ = self.cmd_tx.send(cmd);
    }

    pub fn set_paused(&self, paused: bool) {
        self.send(PhysicsCmd::SetPaused(paused));
    }
    pub fn set_steps_per_frame(&self, s: u32) {
        self.send(PhysicsCmd::SetStepsPerFrame(s));
    }
    pub fn set_dt(&self, dt: f64) {
        self.send(PhysicsCmd::SetDt(dt));
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
    pub fn set_trail_every(&self, n: usize) {
        self.send(PhysicsCmd::SetTrailEvery(n));
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
    /// could time out when the physics thread was mid-batch with high
    /// `steps_per_frame`, producing an empty snapshot with 0 bodies / 0 steps.)
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
            trail_every: self.trail_every,
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
        trail_buf: Arc::new(system.trail_buf().clone()),
        metrics: system.metrics(),
        orbital_elements: system.orbital_elements().to_vec(),
        softening_scale: system.softening_scale(),
        trail_every: system.trail_every(),
    };

    let render = Arc::new(Mutex::new(initial.clone()));
    let render_thr = render.clone();
    let loading = Arc::new(AtomicBool::new(false));
    let loading_thr = loading.clone();

    let thread = thread::spawn(move || {
        physics_loop(system, cmd_rx, render_thr, loading_thr, paused);
    });

    PhysicsHandle {
        cmd_tx,
        render,
        loading,
        bodies: initial.bodies,
        names: initial.names,
        metrics: initial.metrics,
        orbital_elements: initial.orbital_elements,
        softening_scale: initial.softening_scale,
        trail_every: initial.trail_every,
        cached_trail_buf: initial.trail_buf,
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

/// Full publish: positions + trail buffer + names + orbital elements + config.
///
/// Trail buffer clone (~2 MB for large N) makes this expensive; call at most
/// at ~60 Hz.
fn publish_full(system: &System, rs: &mut RenderState) {
    publish_positions(system, rs);
    rs.names = system.names().to_vec();
    rs.trail_buf = Arc::new(system.trail_buf().clone());
    rs.orbital_elements = system.orbital_elements().to_vec();
    rs.softening_scale = system.softening_scale();
    rs.trail_every = system.trail_every();
}

// ── Physics loop ──────────────────────────────────────────────────────────────

fn physics_loop(
    mut system: System,
    cmd_rx: mpsc::Receiver<PhysicsCmd>,
    render: Arc<Mutex<RenderState>>,
    loading: Arc<AtomicBool>,
    initial_paused: bool,
) {
    let mut paused = initial_paused;
    let mut steps_per_frame = 1u32;
    let mut needs_full_publish = false;

    let mut trail_time_acc: f64 = 0.0;
    let sample_interval: f64 = 0.01;

    let full_interval = Duration::from_millis(16);
    let mut last_full = Instant::now().checked_sub(full_interval).unwrap_or_else(Instant::now);

    let pos_interval = Duration::from_millis(8);
    let mut last_pos = Instant::now();

    const POS_CHECK_STEPS: u32 = 8;
    let min_batch_period = Duration::from_micros(100);

    loop {
        let batch_start = Instant::now();

        while let Ok(cmd) = cmd_rx.try_recv() {
            match cmd {
                PhysicsCmd::Shutdown => return,
                PhysicsCmd::SetPaused(p) => paused = p,
                PhysicsCmd::SetStepsPerFrame(s) => steps_per_frame = s.max(1),
                PhysicsCmd::SetDt(dt) => {
                    system.set_dt(dt);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetTheta(theta) => {
                    system.set_theta(theta);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetSofteningScale(s) => {
                    system.set_softening_scale(s);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetGFactor(g) => {
                    system.set_g_factor(g);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetIntegrator(i) => {
                    system.set_integrator(i);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetTrailEvery(n) => {
                    system.set_trail_every(n);
                    needs_full_publish = true;
                },
                PhysicsCmd::AddBody(b) => {
                    system.add_body(b);
                    needs_full_publish = true;
                },
                PhysicsCmd::AddNamedBody(named_body) => {
                    system.add_named_body(named_body);
                    needs_full_publish = true;
                },
                PhysicsCmd::AddBodies(bodies) => {
                    loading.store(true, Ordering::Relaxed);
                    system.add_bodies(bodies);
                    loading.store(false, Ordering::Relaxed);
                    needs_full_publish = true;
                },
                PhysicsCmd::AddNamedBodies(bodies) => {
                    loading.store(true, Ordering::Relaxed);
                    system.add_named_bodies(bodies);
                    loading.store(false, Ordering::Relaxed);
                    needs_full_publish = true;
                },
                PhysicsCmd::RemoveBody(idx) => {
                    system.remove_body(idx);
                    needs_full_publish = true;
                },
                PhysicsCmd::UpdateBody(idx, b) => {
                    system.update_body(idx, b);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetName(idx, name) => {
                    system.set_name(idx, name);
                    needs_full_publish = true;
                },
                PhysicsCmd::LoadBodies(bodies) => {
                    loading.store(true, Ordering::Relaxed);
                    system.load_bodies(bodies);
                    loading.store(false, Ordering::Relaxed);
                    needs_full_publish = true;
                },
                PhysicsCmd::ZeroComVelocity => {
                    system.zero_com_velocity();
                    needs_full_publish = true;
                },
                PhysicsCmd::SetDtMode(mode) => {
                    system.set_dt_mode(mode);
                    needs_full_publish = true;
                },
                PhysicsCmd::SetAdaptiveTheta(enabled) => {
                    system.set_adaptive_theta(enabled);
                    needs_full_publish = true;
                },
                PhysicsCmd::RestoreSnapshot(snap) => {
                    loading.store(true, Ordering::Relaxed);
                    system.restore_from_snapshot(&snap);
                    loading.store(false, Ordering::Relaxed);
                    needs_full_publish = true;
                },
            }
        }

        if !paused {
            let mut steps_since_check = 0u32;

            for _ in 0..steps_per_frame {
                system.step();
                steps_since_check += 1;

                trail_time_acc += system.metrics().dt;

                if trail_time_acc >= sample_interval {
                    system.push_trail();

                    trail_time_acc -= sample_interval;
                }

                if steps_since_check >= POS_CHECK_STEPS {
                    steps_since_check = 0;

                    let now = Instant::now();

                    if now.duration_since(last_full) >= full_interval {
                        system.update_orbital_elements();

                        if let Ok(mut rs) = render.try_lock() {
                            publish_full(&system, &mut rs);
                        }

                        last_full = now;
                        last_pos = now;
                        needs_full_publish = false;
                    } else if now.duration_since(last_pos) >= pos_interval {
                        last_pos = now;

                        if let Ok(mut rs) = render.try_lock() {
                            publish_positions(&system, &mut rs);
                        }
                    }
                }
            }
        }

        // ── Post-batch full publish ────────────────────────────────
        let now = Instant::now();

        if needs_full_publish || (!paused && now.duration_since(last_full) >= full_interval) {
            system.update_orbital_elements();

            if let Ok(mut rs) = render.try_lock() {
                publish_full(&system, &mut rs);
            }

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
