//! TOML header schema for Apsis Records. See ADR-011 §"Header TOML schema".

use crate::physics::gravity::kernel::{Continuity, Exactness};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Header {
    pub apsis: Apsis,
    pub reproducibility: Reproducibility,
    pub unit_system: UnitSystemMeta,
    pub integrator: IntegratorMeta,
    pub kernel: KernelMeta,
    #[serde(default)]
    pub operators: Vec<OperatorMeta>,
    pub bodies: BodiesMeta,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Apsis {
    pub version: String,
    pub git_sha: String,
    pub created_utc: String,
    /// `rustc --version` output captured at build time. Empty when the
    /// build script couldn't invoke rustc (vendored, sandboxed build).
    /// f64 codegen varies between rustc releases; the field is part of
    /// the reproducibility envelope, not just informational.
    #[serde(default)]
    pub rustc_version: String,
    /// Tool that emitted this record. Defaults to `"apsis <version>"`
    /// for records produced by the in-tree writer; downstream wrappers
    /// (alternative bindings, custom binaries) override.
    #[serde(default = "default_generated_by")]
    pub generated_by: String,
}

fn default_generated_by() -> String {
    format!("apsis {}", env!("CARGO_PKG_VERSION"))
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reproducibility {
    pub cargo_lock_blake3: String,
    pub seed: u64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnitSystemMeta {
    pub g: f64,
    pub length: String,
    pub mass: String,
    pub time: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IntegratorMeta {
    pub kind: String,
    pub dt_mode: String,
    pub initial_dt: f64,
    #[serde(default, skip_serializing_if = "serde_json::Map::is_empty")]
    pub params: serde_json::Map<String, serde_json::Value>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KernelMeta {
    pub variant: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub softening: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OperatorMeta {
    pub name: String,
    pub version: String,
    pub crate_hash: String,
    #[serde(default)]
    pub requirements: KernelRequirementsMeta,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct KernelRequirementsMeta {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kernel_exactness: Option<Exactness>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kernel_continuity: Option<Continuity>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodiesMeta {
    pub count: usize,
    pub list: Vec<BodyMeta>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyMeta {
    pub name: String,
    pub mass: f64,
    pub density: f64,
    pub physical_radius: f64,
    pub color: [u8; 3],
    #[serde(default)]
    pub q_pr: f64,
    #[serde(default = "default_albedo")]
    pub albedo: f64,
    pub class: String,
}

fn default_albedo() -> f64 {
    0.5
}

impl Header {
    pub fn to_toml(&self) -> Result<String, toml::ser::Error> {
        toml::to_string_pretty(self)
    }

    pub fn from_toml(s: &str) -> Result<Self, toml::de::Error> {
        toml::from_str(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Header {
        Header {
            apsis: Apsis {
                version: "0.1.0".into(),
                git_sha: "abc123".into(),
                created_utc: "2026-05-16T11:23:45Z".into(),
                rustc_version: "rustc 1.95.0".into(),
                generated_by: "apsis 0.1.0".into(),
            },
            reproducibility: Reproducibility { cargo_lock_blake3: "deadbeef".into(), seed: 42 },
            unit_system: UnitSystemMeta {
                g: 1.0,
                length: "AU".into(),
                mass: "M_sun".into(),
                time: "yr/2pi".into(),
            },
            integrator: IntegratorMeta {
                kind: "IAS15".into(),
                dt_mode: "Adaptive".into(),
                initial_dt: 0.01,
                params: serde_json::json!({"epsilon": 1.0e-9}).as_object().unwrap().clone(),
            },
            kernel: KernelMeta { variant: "Newton".into(), softening: None },
            operators: vec![OperatorMeta {
                name: "apsis-1pn".into(),
                version: "0.1.0".into(),
                crate_hash: "f".repeat(64),
                requirements: KernelRequirementsMeta {
                    kernel_exactness: Some(Exactness::Exact),
                    kernel_continuity: Some(Continuity::Smooth),
                },
            }],
            bodies: BodiesMeta {
                count: 1,
                list: vec![BodyMeta {
                    name: "sun".into(),
                    mass: 1.0,
                    density: 1.408,
                    physical_radius: 4.65e-3,
                    color: [255, 233, 100],
                    q_pr: 0.0,
                    albedo: 0.5,
                    class: "Star".into(),
                }],
            },
        }
    }

    #[test]
    fn header_round_trip_via_toml() {
        let h = sample();
        let s = h.to_toml().unwrap();
        let back: Header = Header::from_toml(&s).unwrap();
        assert_eq!(h, back);
    }

    #[test]
    fn kernel_softening_round_trip() {
        let mut h = sample();
        h.kernel = KernelMeta { variant: "Plummer".into(), softening: Some(1.0e-4) };
        let s = h.to_toml().unwrap();
        let back: Header = Header::from_toml(&s).unwrap();
        assert_eq!(h.kernel, back.kernel);
    }

    #[test]
    fn requirements_optional_per_operator() {
        let mut h = sample();
        h.operators[0].requirements = KernelRequirementsMeta::default();
        let s = h.to_toml().unwrap();
        let back: Header = Header::from_toml(&s).unwrap();
        assert!(back.operators[0].requirements.kernel_exactness.is_none());
        assert!(back.operators[0].requirements.kernel_continuity.is_none());
    }
}
