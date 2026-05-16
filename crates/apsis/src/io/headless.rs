//! Headless batch runner — drives a simulation from a [`RunConfig`] without
//! initialising any GPU or window context.
//!
//! # Entry point
//!
//! ```text
//! apsis --config run.toml
//! ```
//!
//! Parses the TOML, builds the system from the named preset, then advances the
//! integrator loop until `t >= duration`, writing one Apsis Record + CSV rows
//! at the configured intervals.
//!
//! # Output
//!
//! All files land in `out_dir`:
//!
//! | File | Trigger |
//! |------|---------|
//! | `{sim_name}.apsis`            | written incrementally per `snapshot_interval`; trailer on completion |
//! | `{sim_name}_bodies.csv`       | every `csv_interval` sim-time units |
//! | `{sim_name}_system.csv`       | every `csv_interval` sim-time units |
//!
//! # Determinism
//!
//! Given the same `run.toml` and `seed > 0`, two runs on the same platform
//! produce a bit-identical record frame stream. `seed = 0` uses OS entropy
//! and is intentionally non-deterministic.

use std::path::Path;

use crate::core::system::System;
use crate::io::recorder::{RecordMetadata, SimRecorder};
use crate::io::run_config::RunConfig;
use crate::records::{RecordHook, RecordPolicy, provenance::header_from_system};
use crate::templates::{catalog::TEMPLATES, instantiate::instantiate};
use crate::units::UnitSystem;

/// Run a headless batch simulation described by `config`.
///
/// Progress is reported to stderr so it can be redirected independently of
/// stdout. Errors are returned as boxed `std::error::Error`.
pub fn run(config: &RunConfig) -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new(&config.out_dir);
    std::fs::create_dir_all(out_dir)?;

    // ── Locate template ────────────────────────────────────────────────────────
    let entry = TEMPLATES.iter().find(|e| e.name == config.preset).ok_or_else(|| {
        format!(
            "preset {:?} not found; available: {}",
            config.preset,
            TEMPLATES.iter().map(|e| e.name).collect::<Vec<_>>().join(", ")
        )
    })?;

    // ── Build system ───────────────────────────────────────────────────────────
    let template = entry.build(config.seed);
    let named_bodies = instantiate(&template);

    // Templates carry their own UI/CSV unit metadata in `template.units`;
    // the runtime contract is the dimensionless Hénon system, since every
    // preset's body velocities are calibrated for `G = 1`.
    let mut system = System::new(vec![], UnitSystem::canonical())
        .with_theta(0.6)
        .with_dt(config.dt)
        .with_max_depth(32);
    system.set_seed(config.seed);
    system.set_integrator(config.integrator);
    system.add_named_bodies(named_bodies);

    // ── Apsis Record (replaces the legacy .grav snapshot stream) ───────────────
    let record_path = out_dir.join(format!("{}.apsis", config.sim_name));
    let policy = if config.snapshot_interval > 0.0 {
        RecordPolicy::EveryTime(config.snapshot_interval)
    } else {
        RecordPolicy::BookendsAndEvents
    };
    let header = header_from_system(&system, config.seed, None)?;
    let hook = RecordHook::with_header(&record_path, header, policy)?;
    system.hooks_mut().register(0, Box::new(hook));
    eprintln!("[headless] record → {}", record_path.display());

    // ── Optional CSV recorder ──────────────────────────────────────────────────
    let base_path = out_dir.join(&config.sim_name);
    let mut csv: Option<SimRecorder> = if config.csv_interval > 0.0 {
        let meta = RecordMetadata {
            n_bodies: system.bodies().len(),
            integrator_label: system.integrator_kind().label(),
            integrator_order: system.integrator_kind().order(),
            dt: config.dt,
            theta: system.theta(),
            g_factor: system.g_factor(),
            record_interval: config.csv_interval,
            units: template.units,
        };
        Some(SimRecorder::create(&base_path, config.csv_interval, &meta)?)
    } else {
        None
    };

    let n_bodies = system.bodies().len();
    eprintln!(
        "[headless] preset={:?}  N={}  dt={}  duration={}  seed={}",
        config.preset, n_bodies, config.dt, config.duration, config.seed
    );

    // ── Main loop ──────────────────────────────────────────────────────────────
    let mut last_reported = 0.0;
    let report_every = config.duration / 20.0;

    while system.t() < config.duration {
        system.step();
        let t = system.t();

        if let Some(rec) = csv.as_mut()
            && rec.should_record(t)
        {
            system.update_orbital_elements();
            let metrics = system.metrics();
            rec.record(t, system.bodies(), &metrics, system.orbital_elements())?;
        }

        if t - last_reported >= report_every {
            let pct = (t / config.duration * 100.0).min(100.0);
            eprintln!("[headless] {pct:.1}%  t={t:.4e}  steps={}", system.steps());
            last_reported = t;
        }
    }

    // Dropping the system drops its hook registry, which fires
    // `RecordHook::drop` and writes the trailer + final bookend.
    drop(system);

    if let Some(mut rec) = csv {
        rec.flush()?;
        eprintln!(
            "[headless] CSV written → {}  ({} records)",
            rec.system_path().display(),
            rec.records_written,
        );
    }

    eprintln!("[headless] done");
    Ok(())
}
