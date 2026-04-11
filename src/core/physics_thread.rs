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

use crate::core::metrics::Metrics;
use crate::core::snapshot::SimSnapshot;
use crate::core::system::System;
use crate::core::trail_buffer::TrailBuffer;
use crate::domain::body::Body;
use crate::physics::integrator::Integrator;
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
    pub trail_buf: TrailBuffer,
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
    SetIntegrator(Integrator),
    SetTrailEvery(usize),
    AddBody(Body),
    AddBodies(Vec<Body>),
    RemoveBody(usize),
    UpdateBody(usize, Body),
    SetName(usize, String),
    LoadBodies(Vec<Body>),
    ZeroComVelocity,
    RestoreSnapshot(SimSnapshot),
    /// Request a full snapshot; the physics thread sends it back through `tx`.
    RequestSnapshot(mpsc::SyncSender<SimSnapshot>),
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
    // trail_buf is intentionally excluded — use clone_trail_buf() so we
    // avoid cloning the (potentially large) buffer through an extra level.
    bodies: Vec<Body>,
    names: Vec<String>,
    metrics: Metrics,
    orbital_elements: Vec<Option<OrbitalElements>>,
    softening_scale: f64,
    trail_every: usize,

    _thread: thread::JoinHandle<()>,
}

impl PhysicsHandle {
    // ── Frame sync ────────────────────────────────────────────────────────

    /// Pull the latest physics state into the local cache.
    /// Non-blocking: if the physics thread currently holds the lock the cached
    /// values from the previous frame are kept.
    pub fn sync(&mut self) {
        if let Ok(rs) = self.render.try_lock() {
            self.bodies             = rs.bodies.clone();
            self.names              = rs.names.clone();
            self.metrics            = rs.metrics;
            self.orbital_elements   = rs.orbital_elements.clone();
            self.softening_scale    = rs.softening_scale;
            self.trail_every        = rs.trail_every;
        }
    }

    /// Clone the trail buffer directly from the shared render state.
    ///
    /// Prefer this over `trail_buf()` when passing to the GPU backend: it
    /// bypasses the local cache and saves one clone per frame.
    pub fn clone_trail_buf(&self) -> TrailBuffer {
        self.render
            .try_lock()
            .map(|rs| rs.trail_buf.clone())
            .unwrap_or_else(|_| TrailBuffer::new(0))
    }

    // ── Read methods — same signatures as System ──────────────────────────

    /// `true` while the physics thread is processing a heavy command (load / restore).
    /// The UI uses this to display the loading overlay.
    pub fn is_loading(&self) -> bool { self.loading.load(Ordering::Relaxed) }

    pub fn bodies(&self) -> &[Body]                             { &self.bodies }
    pub fn names(&self)  -> &[String]                           { &self.names }
    pub fn total_mass(&self) -> f64 { self.bodies.iter().map(|b| b.mass).sum() }
    pub fn name(&self, idx: usize) -> &str {
        self.names.get(idx).map(|s| s.as_str()).unwrap_or("")
    }
    pub fn metrics(&self)            -> Metrics                 { self.metrics }
    pub fn orbital_elements(&self)   -> &[Option<OrbitalElements>] { &self.orbital_elements }
    pub fn t(&self)                  -> f64  { self.metrics.t }
    pub fn steps(&self)              -> u64  { self.metrics.steps }
    pub fn dt(&self)                 -> f64  { self.metrics.dt }
    pub fn theta(&self)              -> f64  { self.metrics.theta }
    pub fn integrator(&self)         -> Integrator { self.metrics.integrator }
    pub fn softening_scale(&self)    -> f64  { self.softening_scale }
    pub fn g_factor(&self)           -> f64  { self.metrics.g_factor }
    pub fn trail_every(&self)        -> usize { self.trail_every }

    // ── Write methods — tunnel through command channel ────────────────────

    fn send(&self, cmd: PhysicsCmd) {
        // Unbounded channel: send never blocks or errors while thread is live.
        let _ = self.cmd_tx.send(cmd);
    }

    pub fn set_paused(&self, paused: bool)                      { self.send(PhysicsCmd::SetPaused(paused)); }
    pub fn set_steps_per_frame(&self, s: u32)                   { self.send(PhysicsCmd::SetStepsPerFrame(s)); }
    pub fn set_dt(&self, dt: f64)                               { self.send(PhysicsCmd::SetDt(dt)); }
    pub fn set_theta(&self, theta: f64)                         { self.send(PhysicsCmd::SetTheta(theta)); }
    pub fn set_softening_scale(&self, s: f64)                   { self.send(PhysicsCmd::SetSofteningScale(s)); }
    pub fn set_g_factor(&self, g: f64)                          { self.send(PhysicsCmd::SetGFactor(g)); }
    pub fn set_integrator(&self, i: Integrator)                 { self.send(PhysicsCmd::SetIntegrator(i)); }
    pub fn set_trail_every(&self, n: usize)                     { self.send(PhysicsCmd::SetTrailEvery(n)); }
    pub fn add_body(&self, body: Body)                          { self.send(PhysicsCmd::AddBody(body)); }
    /// Add a batch of bodies in one operation. Sets the loading flag immediately
    /// so the overlay appears in the same frame as the user action.
    pub fn add_bodies(&self, bodies: Vec<Body>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::AddBodies(bodies));
    }
    pub fn remove_body(&self, idx: usize)                       { self.send(PhysicsCmd::RemoveBody(idx)); }
    pub fn update_body(&self, idx: usize, body: Body)           { self.send(PhysicsCmd::UpdateBody(idx, body)); }
    pub fn set_name(&self, idx: usize, name: String)            { self.send(PhysicsCmd::SetName(idx, name)); }
    pub fn load_bodies(&self, bodies: Vec<Body>) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::LoadBodies(bodies));
    }
    pub fn zero_com_velocity(&self)                             { self.send(PhysicsCmd::ZeroComVelocity); }
    pub fn restore_from_snapshot(&self, snap: &SimSnapshot) {
        self.loading.store(true, Ordering::Relaxed);
        self.send(PhysicsCmd::RestoreSnapshot(snap.clone()));
    }

    /// Request a full snapshot from the physics thread.
    /// Blocks up to 2 s; intended only for user-triggered saves (rare).
    pub fn to_snapshot(&self) -> SimSnapshot {
        let (tx, rx) = mpsc::sync_channel(1);
        self.send(PhysicsCmd::RequestSnapshot(tx));
        rx.recv_timeout(Duration::from_secs(2)).unwrap_or_else(|_| SimSnapshot {
            save_id: 0,
            t: 0.0,
            steps: 0,
            dt: 0.01,
            theta: 0.5,
            softening_scale: 1.0,
            g_factor: 1.0,
            integrator: Integrator::VelocityVerlet,
            trail_every: 1,
            sim_name: String::new(),
            bodies: vec![],
            names: vec![],
        })
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
        bodies:           system.bodies().to_vec(),
        names:            system.names().to_vec(),
        trail_buf:        system.trail_buf().clone(),
        metrics:          system.metrics(),
        orbital_elements: system.orbital_elements().to_vec(),
        softening_scale:  system.softening_scale(),
        trail_every:      system.trail_every(),
    };

    let render      = Arc::new(Mutex::new(initial.clone()));
    let render_thr  = render.clone();
    let loading     = Arc::new(AtomicBool::new(false));
    let loading_thr = loading.clone();

    let thread = thread::spawn(move || {
        physics_loop(system, cmd_rx, render_thr, loading_thr, paused);
    });

    PhysicsHandle {
        cmd_tx,
        render,
        loading,
        bodies:           initial.bodies,
        names:            initial.names,
        metrics:          initial.metrics,
        orbital_elements: initial.orbital_elements,
        softening_scale:  initial.softening_scale,
        trail_every:      initial.trail_every,
        _thread:          thread,
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
    rs.names            = system.names().to_vec();
    rs.trail_buf        = system.trail_buf().clone();
    rs.orbital_elements = system.orbital_elements().to_vec();
    rs.softening_scale  = system.softening_scale();
    rs.trail_every      = system.trail_every();
}

// ── Physics loop ──────────────────────────────────────────────────────────────

fn physics_loop(
    mut system: System,
    cmd_rx: mpsc::Receiver<PhysicsCmd>,
    render: Arc<Mutex<RenderState>>,
    loading: Arc<AtomicBool>,
    initial_paused: bool,
) {
    let mut paused          = initial_paused;
    let mut steps_per_frame = 1u32;

    // Full publish (trail + everything) at ~60 Hz.
    let full_interval       = Duration::from_millis(16);
    let mut last_full       = Instant::now().checked_sub(full_interval)
                                .unwrap_or_else(Instant::now);

    // Position-only publish (bodies + metrics, no trail clone) at ~120 Hz.
    // Keeps the canvas smooth even when each step batch takes many milliseconds.
    let pos_interval        = Duration::from_millis(8);
    let mut last_pos        = Instant::now();

    // Check elapsed time every this many steps to avoid calling Instant::now()
    // on every single step in fast simulations.
    const POS_CHECK_STEPS: u32 = 8;

    // Minimum wall time per physics batch — caps at ~10 000 batches/s for
    // trivial simulations so we don't spin at 100 % CPU.
    let min_batch_period    = Duration::from_micros(100);

    loop {
        let batch_start = Instant::now();

        // ── Drain commands ────────────────────────────────────────────────
        while let Ok(cmd) = cmd_rx.try_recv() {
            match cmd {
                PhysicsCmd::Shutdown                  => return,
                PhysicsCmd::SetPaused(p)              => paused = p,
                PhysicsCmd::SetStepsPerFrame(s)       => steps_per_frame = s.max(1),
                PhysicsCmd::SetDt(dt)                 => system.set_dt(dt),
                PhysicsCmd::SetTheta(theta)           => system.set_theta(theta),
                PhysicsCmd::SetSofteningScale(s)      => system.set_softening_scale(s),
                PhysicsCmd::SetGFactor(g)             => system.set_g_factor(g),
                PhysicsCmd::SetIntegrator(i)          => system.set_integrator(i),
                PhysicsCmd::SetTrailEvery(n)          => system.set_trail_every(n),
                PhysicsCmd::AddBody(b)                => system.add_body(b),
                PhysicsCmd::AddBodies(bodies)         => {
                    loading.store(true, Ordering::Relaxed);
                    system.add_bodies(bodies);
                    loading.store(false, Ordering::Relaxed);
                }
                PhysicsCmd::RemoveBody(idx)           => system.remove_body(idx),
                PhysicsCmd::UpdateBody(idx, b)        => system.update_body(idx, b),
                PhysicsCmd::SetName(idx, name)        => system.set_name(idx, name),
                PhysicsCmd::LoadBodies(bodies)        => {
                    loading.store(true, Ordering::Relaxed);
                    system.load_bodies(bodies);
                    loading.store(false, Ordering::Relaxed);
                }
                PhysicsCmd::ZeroComVelocity           => system.zero_com_velocity(),
                PhysicsCmd::RestoreSnapshot(snap)     => {
                    loading.store(true, Ordering::Relaxed);
                    system.restore_from_snapshot(&snap);
                    loading.store(false, Ordering::Relaxed);
                }
                PhysicsCmd::RequestSnapshot(tx)       => {
                    let snap = system.to_snapshot();
                    let _ = tx.try_send(snap);
                }
            }
        }

        // ── Step + mid-batch position publish ────────────────────────────
        if !paused {
            let mut steps_since_check = 0u32;

            for _ in 0..steps_per_frame {
                system.step();
                steps_since_check += 1;

                // Every POS_CHECK_STEPS steps, check if it's time for a fast
                // position-only publish (no trail clone — just bodies + metrics).
                if steps_since_check >= POS_CHECK_STEPS {
                    steps_since_check = 0;
                    let now = Instant::now();
                    if now.duration_since(last_pos) >= pos_interval {
                        last_pos = now;
                        if let Ok(mut rs) = render.try_lock() {
                            publish_positions(&system, &mut rs);
                        }
                    }
                }
            }

            system.push_trail();
        }

        // ── Full publish at ~60 Hz (includes trail + names + orbital elems) ──
        let now = Instant::now();
        if now.duration_since(last_full) >= full_interval {
            system.update_orbital_elements();
            if let Ok(mut rs) = render.try_lock() {
                publish_full(&system, &mut rs);
            }
            last_full = now;
            last_pos  = now; // full publish already updated positions
        }

        // ── Throttle ──────────────────────────────────────────────────────
        if paused {
            // Sleep longer when idle — nothing to compute.
            thread::sleep(Duration::from_millis(8));
        } else {
            let elapsed = batch_start.elapsed();
            if elapsed < min_batch_period {
                thread::sleep(min_batch_period - elapsed);
            }
        }
    }
}
