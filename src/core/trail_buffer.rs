//! CPU-side ring buffer for gravitational trail positions.
//!
//! # Layout
//!
//! Positions are stored **column-major**:
//!
//! ```text
//! positions[col * n_bodies + body_idx]
//! ```
//!
//! where `col` cycles through `0..capacity` as a ring buffer and each column
//! represents one recorded time step.  This layout enables a **single
//! contiguous `write_buffer` call** per trail-append step (one full column of
//! `n_bodies × 8` bytes), while keeping reads for a single body's history
//! approximately sequential in the GPU vertex shader.
//!
//! # Dirty tracking
//!
//! [`TrailBuffer`] tracks three dirty flags so the GPU upload is minimal:
//!
//! | Method | Meaning |
//! |--------|---------|
//! [`take_full_upload`] | Topology changed; upload the entire position matrix |
//! [`take_dirty_col`]   | One column was appended; upload only that column |
//! [`take_colors_dirty`]| Body colours changed; re-upload the colour buffer |
//!
//! Each `take_*` method reads and clears the flag atomically (single-threaded).
//!
//! # Adaptive capacity
//!
//! [`adaptive_capacity`] returns the recommended ring-buffer depth for a given
//! body count, keeping `n_bodies × capacity` (total GPU buffer size) roughly
//! constant across the supported simulation scales.

use crate::core::body::Body;
use crate::core::snapshot::TrailSnapshot;

// ── Tuning ────────────────────────────────────────────────────────────────────

/// Maximum number of individually tracked dirty columns before falling back to
/// a full position upload.  Prevents unbounded `write_buffer` call counts at
/// high `steps_per_frame`.
const INCREMENTAL_LIMIT: usize = 8;

// ── Adaptive capacity ─────────────────────────────────────────────────────────

/// Recommended ring-buffer depth (time steps stored) for `n_bodies`.
///
/// Larger values give finer temporal resolution (more samples per orbit arc)
/// at the cost of GPU memory.  The physics thread caps pushes per batch at
/// `capacity / 128`, so the *simulated-time window* visible in the trail is
/// always `128 × steps_per_frame × dt` regardless of capacity — capacity
/// only affects how many distinct position samples exist within that window.
///
/// Memory footprint: `n_bodies × capacity × 8 bytes`.
pub fn adaptive_capacity(n_bodies: usize) -> usize {
    match n_bodies {
        0..=5 => 16_384,   // ≤ 5 bodies  → ≤ 0.7 MB  — 2/3-body problems, binary stars
        6..=50 => 8_192,   // ≤ 50 bodies → ≤ 3.3 MB  — small solar systems
        51..=200 => 2_048, // ≤ 200       → ≤ 3.3 MB  — medium systems
        201..=1_000 => 512,
        1_001..=5_000 => 128,
        5_001..=20_000 => 48,
        _ => 24,
    }
}

// ── TrailBuffer ───────────────────────────────────────────────────────────────

/// Dirty-state for the position buffer.
#[derive(Default, Clone)]
enum PositionsDirty {
    #[default]
    Clean,
    /// One or more specific columns were written.  The vec holds column indices.
    Columns(Vec<u32>),
    /// Topology changed or column count exceeded the incremental limit;
    /// the entire position matrix must be re-uploaded.
    Full,
}

/// Column-major ring buffer of trail positions for all simulation bodies.
#[derive(Clone)]
///
/// Call [`push`](Self::push) after each physics step that should record a
/// position sample.  Call the `take_*` drain methods once per frame to
/// collect what needs to be uploaded to the GPU.
pub struct TrailBuffer {
    /// Flat position storage, column-major: `positions[col * n_bodies + body_idx]`.
    /// Unwritten slots contain `[f32::NAN; 2]`.
    positions: Vec<[f32; 2]>,

    /// RGBA colour per body (linear f32).  Updated when body colours change.
    colors: Vec<[f32; 4]>,

    /// Next write column in the ring buffer (`0..capacity`).
    head: u32,

    /// Number of valid columns stored so far (`0..=capacity`).
    len: u32,

    /// Ring-buffer depth (number of time steps stored simultaneously).
    capacity: u32,

    /// Number of tracked bodies (matrix width).
    n_bodies: u32,

    /// Total number of samples ever pushed into this buffer instance.
    /// Used by the renderer to detect how many new columns appeared since the
    /// previous frame without cloning or diffing the full matrix on the CPU.
    sample_count: u64,

    /// Pending position dirty state.
    pos_dirty: PositionsDirty,

    /// Whether `colors` needs re-uploading.
    colors_dirty: bool,
}

impl TrailBuffer {
    /// Creates an empty buffer sized for `n_bodies` bodies with capacity chosen
    /// by [`adaptive_capacity`].
    pub fn new(n_bodies: usize) -> Self {
        let capacity = adaptive_capacity(n_bodies);
        Self::new_with_capacity(n_bodies, capacity)
    }

    /// Creates an empty buffer for `n_bodies` bodies with an explicit `capacity`.
    ///
    /// Use this when the optimal capacity cannot be derived from `n_bodies` alone
    /// (e.g. when many bodies are belt members that do not render individual trails).
    pub fn new_with_capacity(n_bodies: usize, capacity: usize) -> Self {
        let mut buf = Self {
            positions: Vec::new(),
            colors: Vec::new(),
            head: 0,
            len: 0,
            capacity: 0,
            n_bodies: 0,
            sample_count: 0,
            pos_dirty: PositionsDirty::Clean,
            colors_dirty: false,
        };
        buf.reset(n_bodies, capacity);
        buf
    }

    /// Reinitialises the buffer for `n_bodies` bodies with `capacity` time steps.
    ///
    /// All stored positions and dirty state are discarded.  A full-upload flag
    /// is raised so the GPU buffer is synchronised on the next render frame.
    pub fn reset(&mut self, n_bodies: usize, capacity: usize) {
        self.n_bodies = n_bodies as u32;
        self.capacity = capacity as u32;
        self.head = 0;
        self.len = 0;
        self.sample_count = 0;

        let total = n_bodies * capacity;
        self.positions = vec![[f32::NAN; 2]; total];
        self.colors = vec![[0.0; 4]; n_bodies];

        self.pos_dirty = PositionsDirty::Full;
        self.colors_dirty = true;
    }

    // ── Mutation ──────────────────────────────────────────────────────────── //

    /// Records the current positions of all bodies as one time step.
    ///
    /// `bodies` must have the same length as [`n_bodies`](Self::n_bodies).
    ///
    /// Advances the ring-buffer head and marks the written column dirty.  If
    /// more than [`INCREMENTAL_LIMIT`] columns accumulate without a GPU sync,
    /// the dirty state escalates to a full upload.
    pub fn push(&mut self, bodies: &[Body]) {
        debug_assert_eq!(bodies.len(), self.n_bodies as usize);

        let col = self.head as usize;
        let n = self.n_bodies as usize;
        let base = col * n;

        for (i, b) in bodies.iter().enumerate() {
            self.positions[base + i] = [b.x as f32, b.y as f32];
        }

        // Dirty tracking — escalate to Full after too many incremental columns.
        self.pos_dirty = match &mut self.pos_dirty {
            PositionsDirty::Full => PositionsDirty::Full,
            PositionsDirty::Columns(v) if v.len() >= INCREMENTAL_LIMIT => PositionsDirty::Full,
            PositionsDirty::Columns(v) => {
                v.push(self.head);
                return self.advance_head();
            },
            PositionsDirty::Clean => {
                let mut v = Vec::with_capacity(INCREMENTAL_LIMIT);
                v.push(self.head);
                PositionsDirty::Columns(v)
            },
        };

        self.advance_head();
        self.sample_count += 1;
    }

    /// Applies a rigid translation to **all** stored positions.
    ///
    /// Called during COM recentering.  NaN slots (unwritten) are left as NaN.
    /// After translation the position matrix is marked for a full GPU upload
    /// because every element changes.
    pub fn translate(&mut self, dx: f32, dy: f32) {
        if dx == 0.0 && dy == 0.0 {
            return;
        }
        for p in self.positions.iter_mut() {
            if p[0].is_finite() {
                p[0] += dx;
                p[1] += dy;
            }
        }
        self.pos_dirty = PositionsDirty::Full;
    }

    /// Updates the colour for every body.
    ///
    /// `bodies` must have the same length as [`n_bodies`](Self::n_bodies).
    pub fn update_colors(&mut self, bodies: &[Body]) {
        debug_assert_eq!(bodies.len(), self.n_bodies as usize);
        for (i, b) in bodies.iter().enumerate() {
            self.colors[i] = [
                b.color[0] as f32 / 255.0,
                b.color[1] as f32 / 255.0,
                b.color[2] as f32 / 255.0,
                1.0,
            ];
        }
        self.colors_dirty = true;
    }

    /// Override the alpha channel for specific bodies.
    ///
    /// `show[i] = false` sets alpha to 0 so the trail shader's discard path
    /// skips all segments for that body — effectively hiding its trail without
    /// removing it from the ring buffer.
    ///
    /// Called by the render layer after cloning the buffer, so it never touches
    /// the physics-thread's own copy.
    pub fn apply_visibility(&mut self, show: &[bool]) {
        let mut changed = false;
        for (i, &visible) in show.iter().enumerate() {
            if let Some(c) = self.colors.get_mut(i) {
                let new_alpha: f32 = if visible { 1.0 } else { 0.0 };
                if (c[3] - new_alpha).abs() > 1e-6 {
                    c[3] = new_alpha;
                    changed = true;
                }
            }
        }
        if changed {
            self.colors_dirty = true;
        }
    }

    // ── Dirty-state drain (call once per render frame) ────────────────────── //

    /// Returns the pending position upload and clears the dirty flag.
    ///
    /// The caller is responsible for uploading the returned data to the GPU
    /// before rendering.  Returns [`PendingPositions::Clean`] when no upload
    /// is needed.
    pub fn take_positions_upload(&mut self) -> PendingPositions {
        match std::mem::replace(&mut self.pos_dirty, PositionsDirty::Clean) {
            PositionsDirty::Clean => PendingPositions::Clean,
            PositionsDirty::Full => PendingPositions::Full,
            PositionsDirty::Columns(cols) => PendingPositions::Columns(cols),
        }
    }

    /// Returns whether colours need re-uploading and clears the flag.
    pub fn take_colors_dirty(&mut self) -> bool {
        std::mem::replace(&mut self.colors_dirty, false)
    }

    // ── Accessors ─────────────────────────────────────────────────────────── //

    /// Full flat position matrix (column-major).
    pub fn positions(&self) -> &[[f32; 2]] {
        &self.positions
    }

    /// The slice of positions for one column (`col * n_bodies` elements).
    pub fn column_slice(&self, col: u32) -> &[[f32; 2]] {
        let n = self.n_bodies as usize;
        let base = col as usize * n;
        &self.positions[base..base + n]
    }

    /// Current body colours (RGBA f32, linear).
    pub fn colors(&self) -> &[[f32; 4]] {
        &self.colors
    }

    /// Next write column (ring-buffer head).
    pub fn head(&self) -> u32 {
        self.head
    }

    /// Number of valid columns stored so far.
    pub fn len(&self) -> u32 {
        self.len
    }

    /// Whether fewer than 2 valid columns are stored (no segments to render).
    pub fn is_renderable(&self) -> bool {
        self.len >= 2 && self.n_bodies > 0
    }

    /// Ring-buffer capacity (maximum stored time steps).
    pub fn capacity(&self) -> u32 {
        self.capacity
    }

    /// Number of tracked bodies.
    pub fn n_bodies(&self) -> u32 {
        self.n_bodies
    }

    /// Total number of pushed samples recorded by this buffer.
    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }

    // ── Snapshot persistence ──────────────────────────────────────────────── //

    /// Capture the full trail state into a [`TrailSnapshot`].
    pub fn to_snapshot(&self) -> TrailSnapshot {
        TrailSnapshot {
            n_bodies: self.n_bodies,
            capacity: self.capacity,
            head: self.head,
            len: self.len,
            positions: self.positions.clone(),
        }
    }

    /// Restore trail state from a [`TrailSnapshot`].
    ///
    /// After this call the entire position matrix is flagged for a full GPU
    /// upload so the renderer picks up the restored data on the next frame.
    pub fn restore_from_snapshot(&mut self, snap: &TrailSnapshot) {
        self.n_bodies = snap.n_bodies;
        self.capacity = snap.capacity;
        self.head = snap.head;
        self.len = snap.len;
        self.sample_count = snap.len as u64;
        self.positions = snap.positions.clone();
        self.pos_dirty = PositionsDirty::Full;
        // colours stay dirty from the reset() that preceded this call
    }

    // ── Private ───────────────────────────────────────────────────────────── //

    fn advance_head(&mut self) {
        self.head = (self.head + 1) % self.capacity;
        self.len = (self.len + 1).min(self.capacity);
    }
}

// ── PendingPositions ──────────────────────────────────────────────────────────

/// Describes what position data the GPU buffer needs this frame.
///
/// Returned by [`TrailBuffer::take_positions_upload`].
pub enum PendingPositions {
    /// No position upload needed.
    Clean,
    /// Upload specific columns (incremental path — N × 8 bytes per column).
    Columns(Vec<u32>),
    /// Upload the entire position matrix (topology change or many appends).
    Full,
}
