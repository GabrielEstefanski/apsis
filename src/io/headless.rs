//! Headless batch runner — drives a simulation from a [`RunConfig`] without
//! initialising any GPU or window context.
//!
//! # Entry point
//!
//! ```text
//! gravity-sim --config run.toml
//! ```
//!
//! Parses the TOML, builds the system from the named preset, then advances the
//! integrator loop until `t >= duration`, writing `.grav` snapshots and CSV rows
//! at the configured intervals.
//!
//! # Output
//!
//! All files land in `out_dir`:
//!
//! | File | Trigger |
//! |------|---------|
//! | `{sim_name}_{save_id}.grav`   | every `snapshot_interval` sim-time units |
//! | `{sim_name}_bodies.csv`       | every `csv_interval` sim-time units |
//! | `{sim_name}_system.csv`       | every `csv_interval` sim-time units |
//!
//! # Determinism
//!
//! Given the same `run.toml` and `seed > 0`, two runs on the same platform
//! produce bit-identical trajectories.  `seed = 0` uses OS entropy and is
//! intentionally non-deterministic.

use std::path::{Path, PathBuf};

use crate::core::system::System;
use crate::io::recorder::{RecordMetadata, SimRecorder};
use crate::io::run_config::RunConfig;
use crate::physics::integrator::IntegratorKind;
use crate::templates::{catalog::TEMPLATES, instantiate::instantiate};

/// Run a headless batch simulation described by `config`.
///
/// Progress is reported to stderr so it can be redirected independently of
/// stdout.  Errors are returned as boxed `std::error::Error`.
pub fn run(config: &RunConfig) -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new(&config.out_dir);
    std::fs::create_dir_all(out_dir)?;

    // ── Locate template ────────────────────────────────────────────────────────
    let entry = TEMPLATES
        .iter()
        .find(|e| e.name == config.preset)
        .ok_or_else(|| {
            format!(
                "preset {:?} not found; available: {}",
                config.preset,
                TEMPLATES.iter().map(|e| e.name).collect::<Vec<_>>().join(", ")
            )
        })?;

    // ── Build system ───────────────────────────────────────────────────────────
    let template = (entry.build)(config.seed);
    let named_bodies = instantiate(&template);

    let mut system = System::new(
        vec![],
        0.6,   // θ — standard accuracy
        config.dt,
        32,    // max BH tree depth
        1,
    );
    system.set_seed(config.seed);
    system.set_integrator(parse_integrator(&config.integrator)?);
    system.add_named_bodies(named_bodies);

    // ── Optional CSV recorder ──────────────────────────────────────────────────
    let base_path = out_dir.join(&config.sim_name);
    let mut csv: Option<SimRecorder> = if config.csv_interval > 0.0 {
        let meta = RecordMetadata {
            n_bodies: system.bodies().len(),
            integrator_label: system.integrator_kind().label(),
            integrator_order: system.integrator_kind().order(),
            dt: config.dt,
            theta: system.theta(),
            softening_scale: system.softening_scale(),
            g_factor: system.g_factor(),
            record_interval: config.csv_interval,
            units: template.units.clone(),
        };
        Some(SimRecorder::create(&base_path, config.csv_interval, &meta)?)
    } else {
        None
    };

    // Scheduled targets
    let mut next_snapshot = if config.snapshot_interval > 0.0 { 0.0 } else { f64::INFINITY };
    let snap_interval = if config.snapshot_interval > 0.0 {
        config.snapshot_interval
    } else {
        f64::INFINITY
    };

    let n_bodies = system.bodies().len();
    eprintln!(
        "[headless] preset={:?}  N={}  dt={}  duration={}  seed={}",
        config.preset, n_bodies, config.dt, config.duration, config.seed
    );

    // ── Main loop ──────────────────────────────────────────────────────────────
    let mut last_reported = 0.0;
    let report_every = config.duration / 20.0; // ~5% progress ticks

    while system.t() < config.duration {
        system.step();
        let t = system.t();

        // CSV recording
        if let Some(rec) = csv.as_mut() {
            if rec.should_record(t) {
                system.update_orbital_elements();
                let metrics = system.metrics();
                rec.record(t, system.bodies(), &metrics, system.orbital_elements())?;
            }
        }

        // Snapshot saving
        if t >= next_snapshot {
            let path = save_snapshot(&mut system, out_dir, &config.sim_name)?;
            eprintln!("[headless] snapshot → {}", path.display());
            next_snapshot += snap_interval;
        }

        // Progress reporting
        if t - last_reported >= report_every {
            let pct = (t / config.duration * 100.0).min(100.0);
            eprintln!("[headless] {pct:.1}%  t={t:.4e}  steps={}", system.steps());
            last_reported = t;
        }
    }

    // Final snapshot
    let path = save_snapshot(&mut system, out_dir, &config.sim_name)?;
    eprintln!("[headless] final snapshot → {}", path.display());

    if let Some(mut rec) = csv {
        rec.flush()?;
        eprintln!(
            "[headless] CSV written → {}  ({} records)",
            rec.system_path().display(),
            rec.records_written,
        );
    }

    eprintln!("[headless] done  t={:.4e}  steps={}", system.t(), system.steps());
    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn save_snapshot(
    system: &mut System,
    out_dir: &Path,
    sim_name: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let mut snap = system.to_snapshot();
    snap.sim_name = sim_name.to_owned();
    let path = snap.save_to_dir(out_dir)?;
    Ok(path)
}

fn parse_integrator(s: &str) -> Result<IntegratorKind, Box<dyn std::error::Error>> {
    match s {
        "velocity_verlet" => Ok(IntegratorKind::VelocityVerlet),
        "yoshida4" => Ok(IntegratorKind::Yoshida4),
        "wisdom_holman" => Ok(IntegratorKind::WisdomHolman),
        other => Err(format!(
            "unknown integrator {:?}; use velocity_verlet | yoshida4 | wisdom_holman",
            other
        )
        .into()),
    }
}
