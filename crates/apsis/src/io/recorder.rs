//! Scientific CSV recorder for N-body simulation state.
//!
//! # Output files
//!
//! Each recording session produces **two** companion CSV files:
//!
//! | File | Rows | Description |
//! |------|------|-------------|
//! | `{prefix}_bodies.csv`  | N rows per record | Per-body kinematics, energetics, osculating orbital elements |
//! | `{prefix}_system.csv`  | 1 row per record  | System-level conservation diagnostics |
//!
//! Both files share the same metadata header block (lines beginning with `#`).
//!
//! # Format design decisions
//!
//! - **Long/tidy format** — one body per row; avoids wide columns that scale with N.
//! - **No redundant columns** — `E_total` is not repeated where `ke + pe` suffices;
//!   `dE_rel` is not in the body file because it is a system-level quantity.
//! - **NaN for undefined values** — unbound orbits get `NaN` for `a` and `T`;
//!   parsers (pandas, R, Julia) handle `NaN` natively.
//! - **Sample rate** — the caller controls `record_interval` (simulated time between
//!   records). The recorder writes a row only when `t ≥ next_record_t`, preventing
//!   a full-simulation dump.
//!
//! # Example metadata header
//!
//! ```text
//! # GRAVITY SIMULATOR — Scientific Dataset
//! # generated:        2026-04-09T12:34:56
//! # N:                3
//! # integrator:       Yoshida 4th-order
//! # order:            4
//! # dt:               1.00e-04
//! # theta:            0.500
//! # softening_scale:  1.000
//! # g_factor:         1.000
//! # record_interval:  1.00e-02
//! ```

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::core::metrics::Metrics;
use crate::domain::body::Body;
use crate::physics::energy::per_body_potential_energy;
use crate::physics::orbital::OrbitalElements;
use crate::templates::UnitSystem;

// ── Metadata snapshot ─────────────────────────────────────────────────────────

/// Simulation parameters captured at the moment recording starts.
/// Written once as the comment header in both CSV files.
pub struct RecordMetadata {
    pub n_bodies: usize,
    pub integrator_label: &'static str,
    pub integrator_order: u32,
    pub dt: f64,
    pub theta: f64,
    pub softening_scale: f64,
    pub g_factor: f64,
    pub record_interval: f64,
    /// Physical unit system used by the template that produced this recording.
    /// Written into the metadata block so downstream tools can interpret the
    /// numerical values without re-deriving the unit mapping.
    pub units: UnitSystem,
}

// ── SimRecorder ───────────────────────────────────────────────────────────────

/// Manages two open CSV files and drives timed sampling.
pub struct SimRecorder {
    interval: f64,
    next_record_t: f64,

    bodies_writer: BufWriter<File>,
    system_writer: BufWriter<File>,

    /// Base path (without `_bodies.csv` / `_system.csv` suffix).
    pub base_path: PathBuf,

    /// Total records written (= rows in `_system.csv`).
    pub records_written: u64,

    /// Cached g_factor for per-body PE computation.
    g_factor: f64,
}

impl SimRecorder {
    // ── Construction ──────────────────────────────────────────────────────────

    /// Open both output files, write metadata + column headers, and return a
    /// ready-to-use recorder.
    ///
    /// `base_path` should be a full path **without** extension, e.g.
    /// `/home/user/results/run01`.  The recorder appends `_bodies.csv` and
    /// `_system.csv` automatically.
    pub fn create(base_path: &Path, interval: f64, meta: &RecordMetadata) -> std::io::Result<Self> {
        let bodies_path = Self::suffixed(base_path, "_bodies.csv");
        let system_path = Self::suffixed(base_path, "_system.csv");

        // Create parent directory if needed
        if let Some(parent) = base_path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }

        let mut bw = BufWriter::new(File::create(&bodies_path)?);
        let mut sw = BufWriter::new(File::create(&system_path)?);

        let ts = chrono_or_fallback();
        let header = build_metadata_block(&ts, meta);

        // Bodies file
        writeln!(bw, "{header}")?;
        writeln!(
            bw,
            "t,body_id,x,y,vx,vy,mass,ke,pe,\
             orb_a,orb_e,orb_T,orb_h,orb_eps,orb_omega_deg,orb_type,orb_primary"
        )?;

        // System file
        writeln!(sw, "{header}")?;
        writeln!(
            sw,
            "t,ke_total,pe_total,E_total,dE_rel,\
             Lz_total,dLz_abs,dLz_rel,\
             com_x,com_y,com_vx,com_vy,\
             theta,steps"
        )?;

        Ok(Self {
            interval: interval.max(1e-30),
            next_record_t: 0.0,
            bodies_writer: bw,
            system_writer: sw,
            base_path: base_path.to_owned(),
            records_written: 0,
            g_factor: meta.g_factor,
        })
    }

    // ── Sampling control ──────────────────────────────────────────────────────

    /// Returns `true` when `t` has reached or passed the next scheduled record.
    #[inline]
    pub fn should_record(&self, t: f64) -> bool {
        t >= self.next_record_t
    }

    // ── Recording ────────────────────────────────────────────────────────────

    /// Write one record (N body rows + 1 system row) and advance `next_record_t`.
    ///
    /// Per-body PE is computed here via an O(N²) pass — acceptable because
    /// recording is infrequent relative to the step rate.
    pub fn record(
        &mut self,
        t: f64,
        bodies: &[Body],
        metrics: &Metrics,
        orbital: &[Option<OrbitalElements>],
    ) -> std::io::Result<()> {
        let pe_per_body = per_body_potential_energy(bodies, self.g_factor);

        // ── Body rows ──────────────────────────────────────────────────────
        for (i, body) in bodies.iter().enumerate() {
            let ke = 0.5 * body.mass * (body.vx * body.vx + body.vy * body.vy);
            let pe = pe_per_body.get(i).copied().unwrap_or(f64::NAN);

            let orb = orbital.get(i).and_then(|e| *e);
            let (a, e, period, h, eps, omega_deg, orb_type, primary) = match orb {
                Some(o) => (
                    nan_if_inf(o.a),
                    o.e,
                    nan_if_inf(o.period),
                    o.h,
                    o.energy,
                    o.omega.to_degrees(),
                    o.orbit_type.label(),
                    o.primary_idx as i64,
                ),
                None => {
                    (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN, "none", -1_i64)
                },
            };

            writeln!(
                self.bodies_writer,
                "{t:.6e},{i},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},\
                 {},{},{},{},{},{:.4},{},{}",
                body.x,
                body.y,
                body.vx,
                body.vy,
                body.mass,
                ke,
                pe,
                fmt_f64(a),
                fmt_f64(e),
                fmt_f64(period),
                fmt_f64(h),
                fmt_f64(eps),
                omega_deg,
                orb_type,
                primary,
            )?;
        }

        // ── System row ─────────────────────────────────────────────────────
        writeln!(
            self.system_writer,
            "{t:.6e},{:.6e},{:.6e},{:.6e},{:.4e},{:.6e},{:.4e},{:.4e},\
             {:.6e},{:.6e},{:.6e},{:.6e},{:.4},{}",
            metrics.kinetic,
            metrics.potential,
            metrics.total_energy,
            metrics.rel_energy_error,
            metrics.angular_momentum_z,
            metrics.abs_angular_momentum_error,
            metrics.rel_angular_momentum_error,
            metrics.com_x,
            metrics.com_y,
            metrics.com_vx,
            metrics.com_vy,
            metrics.theta,
            metrics.steps,
        )?;

        self.next_record_t += self.interval;
        self.records_written += 1;

        // Flush every 500 records to protect against crashes
        if self.records_written % 500 == 0 {
            self.bodies_writer.flush()?;
            self.system_writer.flush()?;
        }

        Ok(())
    }

    /// Flush all buffered data to disk.  Call before dropping.
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.bodies_writer.flush()?;
        self.system_writer.flush()?;
        Ok(())
    }

    /// Paths of the two output files.
    pub fn bodies_path(&self) -> PathBuf {
        Self::suffixed(&self.base_path, "_bodies.csv")
    }
    pub fn system_path(&self) -> PathBuf {
        Self::suffixed(&self.base_path, "_system.csv")
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    fn suffixed(base: &Path, suffix: &str) -> PathBuf {
        let mut s = base.as_os_str().to_owned();
        s.push(suffix);
        PathBuf::from(s)
    }
}

impl Drop for SimRecorder {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}

// ── Formatting helpers ────────────────────────────────────────────────────────

/// Formats a float as `NaN` when it is NaN, otherwise in scientific notation.
fn fmt_f64(v: f64) -> String {
    if v.is_nan() { "NaN".into() } else { format!("{v:.6e}") }
}

/// Converts `f64::INFINITY` and `f64::NEG_INFINITY` to `NaN` for CSV output.
fn nan_if_inf(v: f64) -> f64 {
    if v.is_infinite() { f64::NAN } else { v }
}

// ── Metadata block builder ────────────────────────────────────────────────────

fn build_metadata_block(timestamp: &str, m: &RecordMetadata) -> String {
    // Unit conversion lines: emit SI factors when available, "–" otherwise.
    let mass_to_kg =
        m.units.mass_to_kg.map(|v| format!("{v:.4e} kg")).unwrap_or_else(|| "–".into());
    let length_to_m =
        m.units.length_to_m.map(|v| format!("{v:.4e} m")).unwrap_or_else(|| "–".into());
    let time_to_s = m.units.time_to_s.map(|v| format!("{v:.4e} s")).unwrap_or_else(|| "–".into());

    format!(
        "# GRAVITY SIMULATOR — Scientific Dataset\n\
         # ─────────────────────────────────────────────────────\n\
         # generated:        {timestamp}\n\
         # N:                {}\n\
         # integrator:       {}\n\
         # order:            {}\n\
         # dt:               {:.2e}\n\
         # theta:            {:.3}\n\
         # softening_scale:  {:.3}\n\
         # g_factor:         {:.4}\n\
         # record_interval:  {:.2e}\n\
         # ─────────────────────────────────────────────────────\n\
         # unit_system:      {}\n\
         # unit_mass:        {}  (1 sim mass = {mass_to_kg})\n\
         # unit_length:      {}  (1 sim length = {length_to_m})\n\
         # unit_time:        {}  (1 sim time = {time_to_s})\n\
         # ─────────────────────────────────────────────────────",
        m.n_bodies,
        m.integrator_label,
        m.integrator_order,
        m.dt,
        m.theta,
        m.softening_scale,
        m.g_factor,
        m.record_interval,
        m.units.label,
        m.units.mass_unit,
        m.units.length_unit,
        m.units.time_unit,
    )
}

/// Returns an ISO-8601-ish timestamp without pulling in `chrono`.
fn chrono_or_fallback() -> String {
    // std::time doesn't give calendar time easily; use a fixed-format string
    // from SystemTime epoch offset if available, else a placeholder.
    use std::time::{SystemTime, UNIX_EPOCH};
    match SystemTime::now().duration_since(UNIX_EPOCH) {
        Ok(d) => {
            // Very simple: seconds since epoch as a decimal
            // A proper implementation would use chrono, but we keep deps minimal.
            format!("unix+{}", d.as_secs())
        },
        Err(_) => "unknown".into(),
    }
}
