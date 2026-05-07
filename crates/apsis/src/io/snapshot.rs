//! Binary save/load for deterministic simulation reproduction.
//!
//! # Format
//!
//! Each save is a single `.grav` file in a compact little-endian binary layout:
//!
//! ```text
//! [4]  magic         = b"GRAV"
//! [2]  schema_ver    u16 LE   — 1, 2, 3, 4, 5, or 6
//! [8]  save_id       u64 LE   — unix-millis at save time (unique, sortable)
//! [8]  t             f64 LE   — simulated time
//! [8]  steps         u64 LE
//! [8]  dt            f64 LE
//! [8]  theta         f64 LE
//! [8]  softening     f64 LE   — softening_scale
//! [8]  g_factor      f64 LE
//! [1]  integrator    u8       — 0=VV, 1=Yoshida4, 2=WisdomHolman
//! [4]  trail_every   u32 LE
//! --- v3+ only ---
//! [4]  sim_name_len  u32 LE
//! [N]  sim_name      UTF-8 bytes
//! --- v4+ only ---
//! [8]  seed          u64 LE   — reproducibility seed
//! ----------------
//! [4]  n_bodies      u32 LE
//! per body (68 bytes, v6+; 84 bytes in v1–5):
//!   [8] x  [8] y  [8] vx  [8] vy
//!   [8] mass  [8] density  [8] softening  [8] physical_radius
//!   [1] material_id  [3] color_rgb
//! (v1–5 stored two extra f64s here: omega_z + moment_inertia — read and discarded on load)
//! v2+ names section: n_bodies × (u32 len + UTF-8 bytes)
//! --- v4+ trail section ---
//! [1]  trail_has     u8       — 0=no trail, 1=trail present
//! if trail_has == 1:
//!   [4]  n_bodies    u32 LE   — must match header n_bodies
//!   [4]  capacity    u32 LE
//!   [4]  head        u32 LE
//!   [4]  len         u32 LE
//!   [n_bodies * capacity * 12] positions  — column-major [f32; 3] triples
//!                                            (v4–v8 stored 8 bytes per sample;
//!                                            v9 widened to 12)
//! ```
//!
//! The `save_id` field doubles as the filename: `{save_id}.grav`.
//! The save browser reads only the header fields of each file via
//! [`SimSnapshot::read_entry`], avoiding a full deserialisation pass.
//!
//! # Schema versions
//!
//! | Ver | Changes |
//! |-----|---------|
//! | 1   | Initial release |
//! | 2   | Added per-body name strings |
//! | 3   | Added `sim_name` to header |
//! | 4   | Added `seed` field and trail section |
//! | 5   | `integrator` byte extended: `2 = WisdomHolman` |
//! | 6   | Removed `omega_z` and `moment_inertia` from per-body record |
//! | 7   | Added `z` and `vz` to per-body record (3D port) |
//! | 8   | Replaced `material` byte with `q_pr` (8 bytes f64); `Body` no longer carries a material taxonomy field |
//! | 9   | Trail positions widened from `[f32; 2]` to `[f32; 3]` for the 3D camera |
//! | 10  | Per-body record carries a `BodyClass` byte appended after `q_pr` |
//! | 11  | Per-body record carries an `albedo` f64 appended after the class byte |
//!
//! Older files (ver < 5) round-trip cleanly; the `WisdomHolman` variant
//! simply cannot be expressed in them and defaults to `VelocityVerlet` on
//! load. v6 and earlier files default `z = 0`, `vz = 0` on load — the
//! 3D port introduces those components but a planar-only file remains
//! mathematically equivalent under the v7 reader. v7 and earlier files
//! reconstruct `q_pr` from the legacy material byte via a fixed lookup
//! that mirrors the pre-refactor `Material::q_pr()` table. v8 trail
//! sections store two floats per sample; the v9 reader materialises the
//! third component as `z = 0`. v9 and earlier per-body records have no
//! class byte; the v10 reader assigns `BodyClass::Unknown` so the body
//! falls outside the class-based filters until the user re-tags it.
//! v10 and earlier records have no albedo field; the v11 reader falls
//! back to a class-typical placeholder (0.30 for Planet, 0.10 for
//! Asteroid, 0.04 for Comet, 0.50 for Moon, 0.0 for Star/Unknown).

use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::domain::body::Body;
use crate::domain::body_preset::BodyClass;
use crate::physics::integrator::IntegratorKind;

// ── Constants ─────────────────────────────────────────────────────────────────

pub const MAGIC: [u8; 4] = *b"GRAV";
pub const SCHEMA_VERSION: u16 = 11;

// ── Snapshot ──────────────────────────────────────────────────────────────────

/// Complete minimal state required to reproduce a simulation deterministically.
///
/// Obtained via [`System::to_snapshot`] and restored via
/// [`System::restore_from_snapshot`]. The [`save_id`](Self::save_id) field
/// doubles as the on-disk filename (`{save_id}.grav`).
#[derive(Clone)]
pub struct SimSnapshot {
    /// Unix milliseconds at save time — unique identifier and sort key.
    pub save_id: u64,
    /// Simulated time elapsed at the moment of the snapshot.
    pub t: f64,
    /// Discrete integration step counter.
    pub steps: u64,
    /// Fixed time step `Δt`.
    pub dt: f64,
    /// Barnes–Hut opening angle `θ`.
    pub theta: f64,
    /// Global Plummer softening scale factor.
    pub softening_scale: f64,
    /// Gravitational constant multiplier `G_eff = G₀ · g_factor`.
    pub g_factor: f64,
    /// Active symplectic integrator.
    pub integrator_kind: IntegratorKind,
    /// Trail ring-buffer sampling interval (frames between recorded points).
    pub trail_every: usize,
    /// User-assigned simulation label (v3+). Empty string for older saves.
    pub sim_name: String,
    /// Reproducibility seed (v4+). Zero for saves predating v4.
    pub seed: u64,
    /// Per-body state — the only fields that evolve during integration.
    pub bodies: Vec<BodyRecord>,
    /// Display names parallel to [`bodies`](Self::bodies).
    /// May be empty for v1 saves; auto-generated on load in that case.
    pub names: Vec<String>,
    /// Serialised trail ring-buffer (v4+). `None` when absent or oversized.
    pub trail: Option<TrailSnapshot>,
}

/// Per-body fields stored in a [`SimSnapshot`].
///
/// Mirrors [`Body`] exactly, but uses `Copy` semantics so snapshots can be
/// cloned without heap allocation per body. v8 dropped the legacy material
/// byte in favour of `q_pr`; v10 added a one-byte [`BodyClass`] tag; v11
/// added an `albedo` f64.
#[derive(Clone, Copy)]
pub struct BodyRecord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub mass: f64,
    pub density: f64,
    pub softening: f64,
    pub physical_radius: f64,
    pub color: [u8; 3],
    /// Radiation-pressure receiver coefficient. v8+ persists it
    /// directly; v1–7 reconstructed it from the legacy material byte.
    pub q_pr: f64,
    /// UX taxonomy bucket. v10+ persists it directly; older versions
    /// load as [`BodyClass::Unknown`].
    pub class: BodyClass,
    /// Bond albedo. v11+ persists it directly; older versions fall
    /// back to a class-typical placeholder on load.
    pub albedo: f64,
}

impl BodyRecord {
    /// Capture the current state of `b` into a record.
    pub fn from_body(b: &Body) -> Self {
        Self {
            x: b.x,
            y: b.y,
            z: b.z,
            vx: b.vx,
            vy: b.vy,
            vz: b.vz,
            mass: b.mass,
            density: b.density,
            softening: b.softening,
            physical_radius: b.physical_radius,
            color: b.color,
            q_pr: b.q_pr,
            class: b.class,
            albedo: b.albedo,
        }
    }

    /// Reconstruct a [`Body`] from this record.
    ///
    /// Uses the low-level [`Body::new`] constructor and overlays each
    /// stored field explicitly — the snapshot is the canonical source
    /// for every quantity, so no preset is consulted on load.
    pub fn into_body(self) -> Body {
        let mut b = Body::new(self.mass, self.density)
            .at_3d(self.x, self.y, self.z)
            .with_velocity_3d(self.vx, self.vy, self.vz);
        b.softening = self.softening;
        b.physical_radius = self.physical_radius;
        b.color = self.color;
        b.q_pr = self.q_pr;
        b.class = self.class;
        b.albedo = self.albedo;
        b
    }
}

/// Serialised state of a [`TrailBuffer`](crate::core::trail::TrailBuffer).
///
/// The `positions` array is stored column-major:
/// `positions[col * n_bodies + body_idx]`. Each sample is a 3D world-space
/// point; unwritten slots are encoded as `[NaN, NaN, NaN]`. The in-memory
/// ring buffer pads each entry to four floats so its bytes match the GPU
/// std430 stride for `vec3<f32>`; the on-disk format drops the pad.
#[derive(Clone)]
pub struct TrailSnapshot {
    /// Number of bodies whose trails are stored.
    pub n_bodies: u32,
    /// Ring-buffer capacity (columns).
    pub capacity: u32,
    /// Write-head position at save time.
    pub head: u32,
    /// Number of valid entries at save time.
    pub len: u32,
    /// Flat position array, column-major. `NaN` entries represent unwritten slots.
    pub positions: Vec<[f32; 3]>,
}

// ── Save-browser metadata ─────────────────────────────────────────────────────

/// Lightweight record populated by reading only the header of a `.grav` file.
///
/// Used by the save browser to list saves without deserialising body data.
#[derive(Clone)]
pub struct SaveEntry {
    pub path: PathBuf,
    /// Unix milliseconds at save time — matches the filename stem.
    pub save_id: u64,
    /// Simulated time at save.
    pub t: f64,
    /// Step counter at save.
    pub steps: u64,
    /// Number of bodies in the snapshot.
    pub n_bodies: u32,
    /// User-assigned simulation name (v3+). Empty for older saves.
    pub sim_name: String,
    /// Reproducibility seed (v4+). Zero for older saves.
    pub seed: u64,
}

impl SaveEntry {
    /// Returns the simulation name, or `"Unnamed"` if none was set.
    pub fn display_name(&self) -> &str {
        if self.sim_name.is_empty() { "Unnamed" } else { &self.sim_name }
    }

    /// Returns a human-readable UTC timestamp derived from [`save_id`](Self::save_id).
    pub fn display_date(&self) -> String {
        unix_millis_to_display(self.save_id)
    }
}

/// Converts Unix milliseconds to a `"YYYY-MM-DD  HH:MM"` UTC string
/// without pulling in a date-time dependency.
fn unix_millis_to_display(millis: u64) -> String {
    let total_secs = millis / 1000;
    let time_of_day = total_secs % 86400;
    let h = time_of_day / 3600;
    let m = (time_of_day % 3600) / 60;

    let days = total_secs / 86400;
    let z = days as i64 + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mo = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mo <= 2 { y + 1 } else { y };

    format!("{y:04}-{mo:02}-{d:02}  {h:02}:{m:02}")
}

// ── Low-level I/O primitives ──────────────────────────────────────────────────

fn wu8(w: &mut impl Write, v: u8) -> io::Result<()> {
    w.write_all(&[v])
}
fn wu16(w: &mut impl Write, v: u16) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn wu32(w: &mut impl Write, v: u32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn wu64(w: &mut impl Write, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn wf32(w: &mut impl Write, v: f32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn wf64(w: &mut impl Write, v: f64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn ru16(r: &mut impl Read) -> io::Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}
fn ru32(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn ru64(r: &mut impl Read) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
fn rf32(r: &mut impl Read) -> io::Result<f32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(f32::from_le_bytes(b))
}
fn rf64(r: &mut impl Read) -> io::Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

// ── IntegratorKind codec ──────────────────────────────────────────────────────────

/// Encodes an [`IntegratorKind`] as a single byte for on-disk storage.
///
/// | Byte | Variant |
/// |------|---------|
/// | 0    | `VelocityVerlet` |
/// | 1    | `Yoshida4` |
/// | 2    | `WisdomHolman` (v5+) |
/// | 3    | `Ias15` (v6+)     |
fn integrator_to_u8(i: IntegratorKind) -> u8 {
    match i {
        IntegratorKind::VelocityVerlet => 0,
        IntegratorKind::Yoshida4 => 1,
        IntegratorKind::WisdomHolman => 2,
        IntegratorKind::Ias15 => 3,
    }
}

/// Decodes an [`IntegratorKind`] from a single byte.
///
/// Unknown values fall back to `VelocityVerlet` for forward compatibility.
fn u8_to_integrator(v: u8) -> IntegratorKind {
    match v {
        1 => IntegratorKind::Yoshida4,
        2 => IntegratorKind::WisdomHolman,
        3 => IntegratorKind::Ias15,
        _ => IntegratorKind::VelocityVerlet,
    }
}

// ── Legacy material codec (v1–v7 only) ───────────────────────────────────────
//
// Schema versions ≤ 7 stored a single byte naming a fixed `Material`
// taxonomy variant. v8 replaced that with a direct `q_pr` field. The
// table below mirrors the pre-refactor `Material::q_pr()` values and
// is consulted only when reading a v1–v7 file.

/// Reconstruct `q_pr` from the legacy material byte. Unknown values
/// fall back to `0.0` (non-receiver) so a forward-incompatible byte
/// never injects radiation pressure on a body that wasn't a receiver.
fn legacy_q_pr_from_material_byte(v: u8) -> f64 {
    match v {
        // 0=Rocky, 6=Star, 7=BrownDwarf, 8=WhiteDwarf, 2=Gas, 3=IceGiant — none receive
        1 => 0.7, // Icy
        4 => 1.0, // Asteroid
        5 => 0.9, // Comet
        _ => 0.0,
    }
}

/// Class-typical Bond albedo placeholder for snapshots written
/// before v11. Values mirror the `default_albedo` on the
/// corresponding `BodyPreset`; bodies that need their published
/// value re-tagged once the user opens the inspector.
fn albedo_fallback_for(class: BodyClass) -> f64 {
    match class {
        BodyClass::Star => 0.0,
        BodyClass::Planet => 0.30,
        BodyClass::Moon => 0.50,
        BodyClass::Asteroid => 0.10,
        BodyClass::Comet => 0.04,
        BodyClass::Unknown => 0.30,
    }
}

// ── Snapshot I/O ──────────────────────────────────────────────────────────────

fn unix_millis() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_millis() as u64).unwrap_or(0)
}

impl SimSnapshot {
    /// Generates a fresh unique save ID (Unix milliseconds).
    pub fn new_id() -> u64 {
        unix_millis()
    }

    /// Generates a reproducibility seed.
    ///
    /// Semantically distinct from [`new_id`](Self::new_id): the seed
    /// identifies the *initial* configuration, not the save instant.
    pub fn new_seed() -> u64 {
        unix_millis()
    }

    /// Serialises this snapshot to `dir/{save_id}.grav`, creating the
    /// directory if it does not exist.
    ///
    /// Assigns a fresh [`save_id`](Self::save_id) if the current value is
    /// zero.  Returns the path of the written file.
    pub fn save_to_dir(&mut self, dir: &Path) -> io::Result<PathBuf> {
        if self.save_id == 0 {
            self.save_id = unix_millis();
        }
        std::fs::create_dir_all(dir)?;
        let path = dir.join(format!("{}.grav", self.save_id));
        self.write_to(&path)?;
        Ok(path)
    }

    /// Serialises this snapshot to an explicit file path.
    pub fn write_to(&self, path: &Path) -> io::Result<()> {
        use std::io::BufWriter;
        let mut w = BufWriter::new(std::fs::File::create(path)?);

        // Header
        w.write_all(&MAGIC)?;
        wu16(&mut w, SCHEMA_VERSION)?;
        wu64(&mut w, self.save_id)?;
        wf64(&mut w, self.t)?;
        wu64(&mut w, self.steps)?;
        wf64(&mut w, self.dt)?;
        wf64(&mut w, self.theta)?;
        wf64(&mut w, self.softening_scale)?;
        wf64(&mut w, self.g_factor)?;
        wu8(&mut w, integrator_to_u8(self.integrator_kind))?;
        wu32(&mut w, self.trail_every as u32)?;

        // v3: simulation name
        let name_bytes = self.sim_name.as_bytes();
        wu32(&mut w, name_bytes.len() as u32)?;
        w.write_all(name_bytes)?;

        // v4: reproducibility seed
        wu64(&mut w, self.seed)?;

        // Body records (v11: x, y, z, vx, vy, vz, mass, density, softening,
        //                     physical_radius, color[3], q_pr, class[1], albedo)
        wu32(&mut w, self.bodies.len() as u32)?;
        for b in &self.bodies {
            wf64(&mut w, b.x)?;
            wf64(&mut w, b.y)?;
            wf64(&mut w, b.z)?;
            wf64(&mut w, b.vx)?;
            wf64(&mut w, b.vy)?;
            wf64(&mut w, b.vz)?;
            wf64(&mut w, b.mass)?;
            wf64(&mut w, b.density)?;
            wf64(&mut w, b.softening)?;
            wf64(&mut w, b.physical_radius)?;
            w.write_all(&b.color)?;
            wf64(&mut w, b.q_pr)?;
            wu8(&mut w, b.class.to_u8())?;
            wf64(&mut w, b.albedo)?;
        }

        // v2: per-body names
        for name in &self.names {
            let bytes = name.as_bytes();
            wu32(&mut w, bytes.len() as u32)?;
            w.write_all(bytes)?;
        }

        // v4: trail ring-buffer (v9 widened sample stride to 3 floats).
        match &self.trail {
            Some(trail) => {
                wu8(&mut w, 1)?;
                wu32(&mut w, trail.n_bodies)?;
                wu32(&mut w, trail.capacity)?;
                wu32(&mut w, trail.head)?;
                wu32(&mut w, trail.len)?;
                for pos in &trail.positions {
                    wf32(&mut w, pos[0])?;
                    wf32(&mut w, pos[1])?;
                    wf32(&mut w, pos[2])?;
                }
            },
            None => wu8(&mut w, 0)?,
        }

        w.flush()
    }

    /// Deserialises a snapshot from a `.grav` file.
    ///
    /// Supports all schema versions 1–[`SCHEMA_VERSION`]. Fields absent or
    /// removed in older/newer versions are defaulted or discarded.
    pub fn load_from(path: &Path) -> io::Result<Self> {
        use std::io::BufReader;
        let mut r = BufReader::new(std::fs::File::open(path)?);

        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a .grav file"));
        }

        let ver = ru16(&mut r)?;
        if !(1..=SCHEMA_VERSION).contains(&ver) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported schema version {ver} (reader supports ≤{SCHEMA_VERSION})"),
            ));
        }

        let save_id = ru64(&mut r)?;
        let t = rf64(&mut r)?;
        let steps = ru64(&mut r)?;
        let dt = rf64(&mut r)?;
        let theta = rf64(&mut r)?;
        let softening_scale = rf64(&mut r)?;
        let g_factor = rf64(&mut r)?;
        let mut integ_byte = [0u8; 1];
        r.read_exact(&mut integ_byte)?;
        let integrator_kind = u8_to_integrator(integ_byte[0]);
        let trail_every = ru32(&mut r)? as usize;

        let sim_name = if ver >= 3 {
            let len = ru32(&mut r)? as usize;
            let mut buf = vec![0u8; len];
            r.read_exact(&mut buf)?;
            String::from_utf8(buf).unwrap_or_default()
        } else {
            String::new()
        };

        let seed = if ver >= 4 { ru64(&mut r)? } else { 0 };

        let n_bodies = ru32(&mut r)?;
        let mut bodies = Vec::with_capacity(n_bodies as usize);
        for _ in 0..n_bodies {
            let x = rf64(&mut r)?;
            let y = rf64(&mut r)?;
            // v7 introduces `z` between `y` and `vx`. Older versions
            // were planar; default `z = 0`.
            let z = if ver >= 7 { rf64(&mut r)? } else { 0.0 };
            let vx = rf64(&mut r)?;
            let vy = rf64(&mut r)?;
            let vz = if ver >= 7 { rf64(&mut r)? } else { 0.0 };
            let mass = rf64(&mut r)?;
            let density = rf64(&mut r)?;
            let softening = rf64(&mut r)?;
            let physical_radius = rf64(&mut r)?;
            if ver < 6 {
                // v1–5 stored omega_z + moment_inertia here — read and discard
                let _ = rf64(&mut r)?;
                let _ = rf64(&mut r)?;
            }
            // v8 dropped the material byte and added a q_pr f64 after color.
            // Older versions encoded both as `[material_id u8][color rgb]`;
            // q_pr is reconstructed via the legacy lookup.
            let q_pr = if ver < 8 {
                let mut mat_byte = [0u8; 1];
                r.read_exact(&mut mat_byte)?;
                legacy_q_pr_from_material_byte(mat_byte[0])
            } else {
                0.0 // placeholder; replaced after the color read below
            };
            let mut color = [0u8; 3];
            r.read_exact(&mut color)?;
            let q_pr = if ver >= 8 { rf64(&mut r)? } else { q_pr };
            // v10 appends a one-byte BodyClass after q_pr. Older
            // versions had no class field — load as Unknown so the
            // body sits outside the class-based filters until the
            // user re-tags it.
            let class = if ver >= 10 {
                let mut class_byte = [0u8; 1];
                r.read_exact(&mut class_byte)?;
                BodyClass::from_u8(class_byte[0])
            } else {
                BodyClass::Unknown
            };
            // v11 appends an albedo f64 after the class byte. Older
            // versions fall back to a class-typical placeholder so the
            // photometry pipeline has a sane number to work with until
            // the user retags.
            let albedo = if ver >= 11 { rf64(&mut r)? } else { albedo_fallback_for(class) };

            bodies.push(BodyRecord {
                x,
                y,
                z,
                vx,
                vy,
                vz,
                mass,
                density,
                softening,
                physical_radius,
                color,
                q_pr,
                class,
                albedo,
            });
        }

        let names = if ver >= 2 {
            let mut ns = Vec::with_capacity(n_bodies as usize);
            for _ in 0..n_bodies {
                let len = ru32(&mut r)? as usize;
                let mut buf = vec![0u8; len];
                r.read_exact(&mut buf)?;
                ns.push(String::from_utf8(buf).unwrap_or_default());
            }
            ns
        } else {
            Vec::new()
        };

        let trail = if ver >= 4 {
            let mut has_byte = [0u8; 1];
            r.read_exact(&mut has_byte)?;
            if has_byte[0] != 0 {
                let tn = ru32(&mut r)?;
                let cap = ru32(&mut r)?;
                let hd = ru32(&mut r)?;
                let ln = ru32(&mut r)?;
                let total = (tn as usize) * (cap as usize);
                let mut positions = Vec::with_capacity(total);
                if ver >= 9 {
                    for _ in 0..total {
                        positions.push([rf32(&mut r)?, rf32(&mut r)?, rf32(&mut r)?]);
                    }
                } else {
                    // v4–v8 stored two floats per sample (planar trail);
                    // pad the third component to 0 so the 3D ring buffer
                    // restores planar saves equivalently.
                    for _ in 0..total {
                        positions.push([rf32(&mut r)?, rf32(&mut r)?, 0.0]);
                    }
                }
                Some(TrailSnapshot { n_bodies: tn, capacity: cap, head: hd, len: ln, positions })
            } else {
                None
            }
        } else {
            None
        };

        Ok(SimSnapshot {
            save_id,
            t,
            steps,
            dt,
            theta,
            softening_scale,
            g_factor,
            integrator_kind,
            trail_every,
            sim_name,
            seed,
            bodies,
            names,
            trail,
        })
    }

    /// Reads only the header fields required by the save browser.
    ///
    /// Significantly faster than [`load_from`](Self::load_from) for directory
    /// listings because it skips the body records, names, and trail data.
    pub fn read_entry(path: &Path) -> io::Result<SaveEntry> {
        use std::io::BufReader;
        let mut r = BufReader::new(std::fs::File::open(path)?);

        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a .grav file"));
        }

        let ver = ru16(&mut r)?;
        let save_id = ru64(&mut r)?;
        let t = rf64(&mut r)?;
        let steps = ru64(&mut r)?;

        // Skip: dt(8) + theta(8) + softening(8) + g_factor(8) + integrator(1) + trail_every(4) = 37 bytes.
        let mut skip = [0u8; 37];
        r.read_exact(&mut skip)?;

        let sim_name = if ver >= 3 {
            let len = ru32(&mut r)? as usize;
            let mut buf = vec![0u8; len];
            r.read_exact(&mut buf)?;
            String::from_utf8(buf).unwrap_or_default()
        } else {
            String::new()
        };

        let seed = if ver >= 4 { ru64(&mut r)? } else { 0 };

        let n_bodies = ru32(&mut r)?;

        Ok(SaveEntry { path: path.to_owned(), save_id, t, steps, n_bodies, sim_name, seed })
    }
}

// ── Directory listing ─────────────────────────────────────────────────────────

/// Scans `dir` for `.grav` files and returns their headers sorted newest-first.
///
/// Files that fail to parse are silently skipped so a single corrupt save
/// does not prevent the browser from loading.
pub fn list_saves(dir: &Path) -> Vec<SaveEntry> {
    let Ok(read_dir) = std::fs::read_dir(dir) else {
        return Vec::new();
    };

    let mut entries: Vec<SaveEntry> = read_dir
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|x| x == "grav").unwrap_or(false))
        .filter_map(|e| SimSnapshot::read_entry(&e.path()).ok())
        .collect();

    entries.sort_by_key(|b| std::cmp::Reverse(b.save_id));
    entries
}
