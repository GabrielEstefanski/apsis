//! Binary save/load for deterministic simulation reproduction.
//!
//! # Format
//!
//! Each save is a single `.grav` file in a compact little-endian binary layout:
//!
//! ```text
//! [4]  magic        = b"GRAV"
//! [2]  schema_ver   = 1  (u16 LE)
//! [8]  save_id      u64 LE  — unix-millis at save time (unique, sortable)
//! [8]  t            f64 LE  — simulated time
//! [8]  steps        u64 LE
//! [8]  dt           f64 LE
//! [8]  theta        f64 LE
//! [8]  softening    f64 LE  — softening_scale
//! [8]  g_factor     f64 LE
//! [1]  integrator   u8      — 0=VV, 1=Yoshida4
//! [4]  trail_every  u32 LE
//! [4]  n_bodies     u32 LE
//! per body (84 bytes):
//!   [8] x  [8] y  [8] vx  [8] vy
//!   [8] mass  [8] density  [8] softening  [8] physical_radius
//!   [8] omega_z  [8] moment_inertia
//!   [1] material_id  [3] color_rgb
//! ```
//!
//! The save_id is used as the filename: `{save_id}.grav`.
//! Listing saves just reads the header (first 63 bytes) of each file.

use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::domain::body::Body;
use crate::domain::materials::Material;
use crate::physics::integrator::Integrator;

// ── Constants ─────────────────────────────────────────────────────────────────

pub const MAGIC: [u8; 4] = *b"GRAV";
pub const SCHEMA_VERSION: u16 = 2;

/// Byte offset where body data starts (after header).
const HEADER_BYTES: usize = 4 + 2 + 8 + 8 + 8 + 8 + 8 + 8 + 8 + 1 + 4 + 4; // 71

/// Bytes per body record.
const BODY_BYTES: usize = 10 * 8 + 1 + 3; // 84

// ── Snapshot ──────────────────────────────────────────────────────────────────

/// Complete minimal state required to reproduce a simulation deterministically.
#[derive(Clone)]
pub struct SimSnapshot {
    /// Unix milliseconds at save time — doubles as unique ID and filename.
    pub save_id: u64,
    /// Simulated time when saved.
    pub t: f64,
    /// Discrete step counter.
    pub steps: u64,
    /// Fixed time step.
    pub dt: f64,
    /// Barnes–Hut opening angle.
    pub theta: f64,
    /// Global Plummer softening scale.
    pub softening_scale: f64,
    /// G multiplier.
    pub g_factor: f64,
    /// Active integrator.
    pub integrator: Integrator,
    /// Trail sampling interval.
    pub trail_every: usize,
    /// Body states — the only things that evolve.
    pub bodies: Vec<BodyRecord>,
    /// Display names, parallel to `bodies`. May be empty for v1 saves (auto-generated on load).
    pub names: Vec<String>,
}

/// Per-body fields stored in a snapshot.
#[derive(Clone, Copy)]
pub struct BodyRecord {
    pub x: f64,
    pub y: f64,
    pub vx: f64,
    pub vy: f64,
    pub mass: f64,
    pub density: f64,
    pub softening: f64,
    pub physical_radius: f64,
    pub omega_z: f64,
    pub moment_inertia: f64,
    pub material: Material,
    pub color: [u8; 3],
}

impl BodyRecord {
    pub fn from_body(b: &Body) -> Self {
        Self {
            x: b.x,
            y: b.y,
            vx: b.vx,
            vy: b.vy,
            mass: b.mass,
            density: b.density,
            softening: b.softening,
            physical_radius: b.physical_radius,
            omega_z: b.omega_z,
            moment_inertia: b.moment_inertia,
            material: b.material,
            color: b.color,
        }
    }

    pub fn into_body(self) -> Body {
        let mut b = Body::new(self.x, self.y, self.vx, self.vy, self.mass, self.material);
        // Override derived fields with exactly-stored values to preserve the
        // simulation state rather than recomputing from the density model.
        b.density = self.density;
        b.softening = self.softening;
        b.physical_radius = self.physical_radius;
        b.omega_z = self.omega_z;
        b.moment_inertia = self.moment_inertia;
        b.color = self.color;
        b
    }
}

// ── Metadata for the browser ──────────────────────────────────────────────────

/// Lightweight header read for the load-save browser (no body data needed).
#[derive(Clone)]
pub struct SaveEntry {
    pub path: PathBuf,
    pub save_id: u64,
    pub t: f64,
    pub steps: u64,
    pub n_bodies: u32,
}

impl SaveEntry {
    pub fn display_name(&self) -> String {
        let secs = self.save_id / 1000;
        // Format: unix+{secs} like the recorder timestamp
        format!("unix+{secs}")
    }
}

// ── Serialisation helpers ─────────────────────────────────────────────────────

fn wu16(w: &mut impl Write, v: u16) -> io::Result<()> { w.write_all(&v.to_le_bytes()) }
fn wu32(w: &mut impl Write, v: u32) -> io::Result<()> { w.write_all(&v.to_le_bytes()) }
fn wu64(w: &mut impl Write, v: u64) -> io::Result<()> { w.write_all(&v.to_le_bytes()) }
fn wf64(w: &mut impl Write, v: f64) -> io::Result<()> { w.write_all(&v.to_le_bytes()) }

fn ru16(r: &mut impl Read) -> io::Result<u16> {
    let mut b = [0u8; 2]; r.read_exact(&mut b)?; Ok(u16::from_le_bytes(b))
}
fn ru32(r: &mut impl Read) -> io::Result<u32> {
    let mut b = [0u8; 4]; r.read_exact(&mut b)?; Ok(u32::from_le_bytes(b))
}
fn ru64(r: &mut impl Read) -> io::Result<u64> {
    let mut b = [0u8; 8]; r.read_exact(&mut b)?; Ok(u64::from_le_bytes(b))
}
fn rf64(r: &mut impl Read) -> io::Result<f64> {
    let mut b = [0u8; 8]; r.read_exact(&mut b)?; Ok(f64::from_le_bytes(b))
}

fn integrator_to_u8(i: Integrator) -> u8 {
    match i { Integrator::VelocityVerlet => 0, Integrator::Yoshida4 => 1 }
}
fn u8_to_integrator(v: u8) -> Integrator {
    match v { 1 => Integrator::Yoshida4, _ => Integrator::VelocityVerlet }
}

fn material_to_u8(m: Material) -> u8 {
    match m {
        Material::Rocky     => 0,
        Material::Icy       => 1,
        Material::Gas       => 2,
        Material::IceGiant  => 3,
        Material::Asteroid  => 4,
        Material::Comet     => 5,
        Material::Star      => 6,
        Material::BrownDwarf => 7,
        Material::WhiteDwarf => 8,
    }
}
fn u8_to_material(v: u8) -> Material {
    match v {
        1 => Material::Icy,
        2 => Material::Gas,
        3 => Material::IceGiant,
        4 => Material::Asteroid,
        5 => Material::Comet,
        6 => Material::Star,
        7 => Material::BrownDwarf,
        8 => Material::WhiteDwarf,
        _ => Material::Rocky,
    }
}

// ── Save / Load ───────────────────────────────────────────────────────────────

fn unix_millis() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis() as u64)
        .unwrap_or(0)
}

impl SimSnapshot {
    /// Assign a fresh unique ID (unix millis) and return it.
    pub fn new_id() -> u64 {
        unix_millis()
    }

    /// Write this snapshot to `dir/{save_id}.grav`.
    /// Creates the directory if needed.
    pub fn save_to_dir(&mut self, dir: &Path) -> io::Result<PathBuf> {
        if self.save_id == 0 {
            self.save_id = unix_millis();
        }
        std::fs::create_dir_all(dir)?;
        let path = dir.join(format!("{}.grav", self.save_id));
        self.write_to(&path)?;
        Ok(path)
    }

    /// Write this snapshot to an explicit path.
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
        w.write_all(&[integrator_to_u8(self.integrator)])?;
        wu32(&mut w, self.trail_every as u32)?;
        wu32(&mut w, self.bodies.len() as u32)?;

        // Bodies
        for b in &self.bodies {
            wf64(&mut w, b.x)?;
            wf64(&mut w, b.y)?;
            wf64(&mut w, b.vx)?;
            wf64(&mut w, b.vy)?;
            wf64(&mut w, b.mass)?;
            wf64(&mut w, b.density)?;
            wf64(&mut w, b.softening)?;
            wf64(&mut w, b.physical_radius)?;
            wf64(&mut w, b.omega_z)?;
            wf64(&mut w, b.moment_inertia)?;
            w.write_all(&[material_to_u8(b.material)])?;
            w.write_all(&b.color)?;
        }

        // Names section (v2+): each name as u32 length + UTF-8 bytes
        for name in &self.names {
            let bytes = name.as_bytes();
            wu32(&mut w, bytes.len() as u32)?;
            w.write_all(bytes)?;
        }

        w.flush()
    }

    /// Load a snapshot from a `.grav` file.
    pub fn load_from(path: &Path) -> io::Result<Self> {
        use std::io::BufReader;
        let mut r = BufReader::new(std::fs::File::open(path)?);

        // Magic
        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a .grav file"));
        }

        // Schema version — accept v1 (no names) and v2 (with names)
        let ver = ru16(&mut r)?;
        if ver != 1 && ver != SCHEMA_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("unsupported schema version {ver} (expected {SCHEMA_VERSION})"),
            ));
        }

        let save_id     = ru64(&mut r)?;
        let t           = rf64(&mut r)?;
        let steps       = ru64(&mut r)?;
        let dt          = rf64(&mut r)?;
        let theta       = rf64(&mut r)?;
        let softening_scale = rf64(&mut r)?;
        let g_factor    = rf64(&mut r)?;
        let mut integ_byte = [0u8; 1];
        r.read_exact(&mut integ_byte)?;
        let integrator  = u8_to_integrator(integ_byte[0]);
        let trail_every = ru32(&mut r)? as usize;
        let n_bodies    = ru32(&mut r)?;

        let mut bodies = Vec::with_capacity(n_bodies as usize);
        for _ in 0..n_bodies {
            let x              = rf64(&mut r)?;
            let y              = rf64(&mut r)?;
            let vx             = rf64(&mut r)?;
            let vy             = rf64(&mut r)?;
            let mass           = rf64(&mut r)?;
            let density        = rf64(&mut r)?;
            let softening      = rf64(&mut r)?;
            let physical_radius = rf64(&mut r)?;
            let omega_z        = rf64(&mut r)?;
            let moment_inertia = rf64(&mut r)?;
            let mut mat_byte = [0u8; 1];
            r.read_exact(&mut mat_byte)?;
            let material = u8_to_material(mat_byte[0]);
            let mut color = [0u8; 3];
            r.read_exact(&mut color)?;

            bodies.push(BodyRecord {
                x, y, vx, vy, mass, density, softening, physical_radius,
                omega_z, moment_inertia, material, color,
            });
        }

        // Names section (v2+): read n_bodies name entries
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
            Vec::new() // caller will auto-generate names
        };

        Ok(SimSnapshot {
            save_id, t, steps, dt, theta, softening_scale,
            g_factor, integrator, trail_every, bodies, names,
        })
    }

    /// Read only the header fields needed for the browser listing.
    pub fn read_entry(path: &Path) -> io::Result<SaveEntry> {
        use std::io::BufReader;
        let mut r = BufReader::new(std::fs::File::open(path)?);

        let mut magic = [0u8; 4];
        r.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "not a .grav file"));
        }
        let _ver     = ru16(&mut r)?;
        let save_id  = ru64(&mut r)?;
        let t        = rf64(&mut r)?;
        let steps    = ru64(&mut r)?;
        // skip dt, theta, softening, g_factor, integrator, trail_every
        let mut skip = [0u8; 8 + 8 + 8 + 8 + 1 + 4];
        r.read_exact(&mut skip)?;
        let n_bodies = ru32(&mut r)?;

        Ok(SaveEntry {
            path: path.to_owned(),
            save_id,
            t,
            steps,
            n_bodies,
        })
    }
}

// ── Directory listing ─────────────────────────────────────────────────────────

/// Scan `dir` for `.grav` files and return entries sorted newest-first.
pub fn list_saves(dir: &Path) -> Vec<SaveEntry> {
    let Ok(read_dir) = std::fs::read_dir(dir) else { return Vec::new() };

    let mut entries: Vec<SaveEntry> = read_dir
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|x| x == "grav").unwrap_or(false))
        .filter_map(|e| SimSnapshot::read_entry(&e.path()).ok())
        .collect();

    // Newest first (largest save_id = latest unix millis)
    entries.sort_by(|a, b| b.save_id.cmp(&a.save_id));
    entries
}
