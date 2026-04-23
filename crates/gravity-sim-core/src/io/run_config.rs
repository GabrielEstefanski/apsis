//! Declarative run configuration for headless batch simulations.
//!
//! # File format (`run.toml`)
//!
//! ```toml
//! preset    = "Solar System"       # name from TEMPLATES catalog
//! integrator = "velocity_verlet"   # "velocity_verlet" | "yoshida4" | "wisdom_holman"
//! dt        = 0.001
//! duration  = 10.0
//!
//! snapshot_interval = 1.0          # simulated time between .grav saves; 0 = disable
//! csv_interval      = 0.1          # simulated time between CSV rows; 0 = disable
//!
//! seed      = 42                   # 0 = OS entropy (non-deterministic)
//! out_dir   = "out/solar"
//! sim_name  = "solar-run-01"
//! ```
//!
//! All fields are required.  Unknown keys are silently ignored by the TOML parser.
//!
//! # Field notes
//!
//! - `preset` must exactly match the `name` field of a [`TemplateEntry`] in the
//!   catalog.  A mismatch is a hard error at run time.
//! - `seed = 0` → OS entropy; preset will behave non-deterministically.
//! - `snapshot_interval = 0` or `csv_interval = 0` disables that output stream.
//! - `out_dir` is created automatically (including parents) if it does not exist.

use std::path::Path;

use crate::physics::integrator::IntegratorKind;

/// Fully specified parameters for a headless batch run.
#[derive(serde::Deserialize, Debug, Clone)]
pub struct RunConfig {
    /// Template name — must match [`TemplateEntry::name`] exactly.
    pub preset: String,

    /// Integration algorithm slug — see [`IntegratorKind::slug`] for valid values.
    #[serde(deserialize_with = "deserialize_integrator_kind")]
    pub integrator: IntegratorKind,

    /// Fixed timestep (simulation units matching the preset's unit system).
    pub dt: f64,

    /// Total simulated time to advance.
    pub duration: f64,

    /// Simulated time between snapshot saves. `0` disables snapshots.
    pub snapshot_interval: f64,

    /// Simulated time between CSV records. `0` disables CSV output.
    pub csv_interval: f64,

    /// Reproducibility seed.  `0` draws from OS entropy.
    pub seed: u64,

    /// Directory where all output files are written.
    pub out_dir: String,

    /// Human-readable label embedded in snapshot `sim_name` and CSV metadata.
    pub sim_name: String,
}

impl RunConfig {
    /// Parse a `run.toml` file from disk.
    pub fn from_file(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let text = std::fs::read_to_string(path)?;
        let cfg: Self = toml::from_str(&text)?;
        Ok(cfg)
    }
}

fn deserialize_integrator_kind<'de, D>(d: D) -> Result<IntegratorKind, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = <String as serde::Deserialize>::deserialize(d)?;
    s.parse::<IntegratorKind>().map_err(serde::de::Error::custom)
}
