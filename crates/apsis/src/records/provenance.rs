//! Builds a record [`Header`] from a [`System`] + the workspace `Cargo.lock`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use crate::core::system::System;
use crate::records::header::{
    Apsis, BodiesMeta, BodyMeta, Header, IntegratorMeta, KernelMeta, KernelRequirementsMeta,
    OperatorMeta, Reproducibility, UnitSystemMeta,
};

#[derive(Debug)]
pub enum ProvenanceError {
    LockNotFound(PathBuf),
    LockRead(PathBuf, std::io::Error),
    LockParse(String),
}

impl std::fmt::Display for ProvenanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LockNotFound(p) => {
                write!(f, "Cargo.lock not found searching upward from {}", p.display())
            },
            Self::LockRead(p, e) => write!(f, "failed to read {}: {e}", p.display()),
            Self::LockParse(msg) => write!(f, "Cargo.lock parse error: {msg}"),
        }
    }
}

impl std::error::Error for ProvenanceError {}

/// Build a complete record [`Header`] from a `System`. Reads `Cargo.lock`
/// from `lock_path` (or, if `None`, walks up from the current working
/// directory until a `Cargo.lock` is found).
pub fn header_from_system(
    sys: &System,
    seed: u64,
    lock_path: Option<&Path>,
) -> Result<Header, ProvenanceError> {
    let lock_path = resolve_lock_path(lock_path)?;
    let lock_bytes =
        std::fs::read(&lock_path).map_err(|e| ProvenanceError::LockRead(lock_path.clone(), e))?;
    let lock_hash = blake3::hash(&lock_bytes).to_hex().to_string();
    let lock_index = parse_lock(&lock_bytes)?;

    let units = sys.units();
    let kernel = sys.kernel();
    let kernel_variant = kernel_variant_name(kernel.as_ref());

    let mut operators = Vec::new();
    operators.extend(sys.hamiltonian_perturbations().iter().filter_map(|op| {
        let cit = op.citation()?;
        let req = op.kernel_requirements();
        let crate_hash = lock_index.get(cit.crate_name).map(|(_, h)| h.clone()).unwrap_or_default();
        Some(OperatorMeta {
            name: cit.crate_name.to_string(),
            version: cit.crate_version.to_string(),
            crate_hash,
            requirements: KernelRequirementsMeta {
                kernel_exactness: req.required_exactness,
                kernel_continuity: req.min_continuity,
            },
        })
    }));
    let seen: std::collections::HashSet<String> =
        operators.iter().map(|o| o.name.clone()).collect();
    operators.extend(sys.non_conservative_perturbations().iter().filter_map(|op| {
        let cit = op.citation()?;
        // Avoid duplicate entries when the same crate publishes both
        // Hamiltonian and non-conservative operators.
        if seen.contains(cit.crate_name) {
            return None;
        }
        let req = op.kernel_requirements();
        let crate_hash = lock_index.get(cit.crate_name).map(|(_, h)| h.clone()).unwrap_or_default();
        Some(OperatorMeta {
            name: cit.crate_name.to_string(),
            version: cit.crate_version.to_string(),
            crate_hash,
            requirements: KernelRequirementsMeta {
                kernel_exactness: req.required_exactness,
                kernel_continuity: req.min_continuity,
            },
        })
    }));

    let body_meta: Vec<BodyMeta> = sys
        .bodies()
        .iter()
        .zip(sys.names().iter())
        .map(|(b, name)| BodyMeta {
            name: name.clone(),
            mass: b.mass,
            density: b.density,
            physical_radius: b.physical_radius,
            color: b.color,
            q_pr: b.q_pr,
            albedo: b.albedo,
            class: format!("{:?}", b.class),
        })
        .collect();

    Ok(Header {
        apsis: Apsis {
            version: env!("CARGO_PKG_VERSION").to_string(),
            git_sha: option_env!("APSIS_GIT_SHA").unwrap_or("unknown").to_string(),
            created_utc: rfc3339_utc_now(),
        },
        reproducibility: Reproducibility { cargo_lock_blake3: lock_hash, seed },
        unit_system: UnitSystemMeta {
            g: units.g(),
            length: units.length_label().to_string(),
            mass: units.mass_label().to_string(),
            time: units.time_label().to_string(),
        },
        integrator: IntegratorMeta {
            kind: format!("{:?}", sys.integrator_kind()),
            dt_mode: format!("{:?}", sys.dt_mode()),
            initial_dt: sys.dt(),
            params: serde_json::Map::new(),
        },
        kernel: KernelMeta { variant: kernel_variant, softening: None },
        operators,
        bodies: BodiesMeta { count: body_meta.len(), list: body_meta },
    })
}

fn resolve_lock_path(explicit: Option<&Path>) -> Result<PathBuf, ProvenanceError> {
    if let Some(p) = explicit {
        return Ok(p.to_path_buf());
    }
    let start =
        std::env::current_dir().map_err(|e| ProvenanceError::LockRead(PathBuf::from("."), e))?;
    let mut cursor: Option<&Path> = Some(&start);
    while let Some(dir) = cursor {
        let candidate = dir.join("Cargo.lock");
        if candidate.exists() {
            return Ok(candidate);
        }
        cursor = dir.parent();
    }
    Err(ProvenanceError::LockNotFound(start))
}

/// Index `Cargo.lock` by crate name → (version, checksum). v0.1 uses a
/// minimal manual TOML walk; pulling in `cargo_metadata` would drag
/// considerable transitive deps for a one-shot read.
fn parse_lock(bytes: &[u8]) -> Result<HashMap<String, (String, String)>, ProvenanceError> {
    let s = std::str::from_utf8(bytes)
        .map_err(|e| ProvenanceError::LockParse(format!("non-UTF8: {e}")))?;
    let parsed: toml::Value =
        toml::from_str(s).map_err(|e| ProvenanceError::LockParse(e.to_string()))?;
    let packages = parsed
        .get("package")
        .and_then(|v| v.as_array())
        .ok_or_else(|| ProvenanceError::LockParse("missing [[package]]".into()))?;
    let mut out = HashMap::new();
    for pkg in packages {
        let Some(name) = pkg.get("name").and_then(|v| v.as_str()) else {
            continue;
        };
        let Some(version) = pkg.get("version").and_then(|v| v.as_str()) else {
            continue;
        };
        let checksum = pkg.get("checksum").and_then(|v| v.as_str()).unwrap_or("");
        out.insert(name.to_string(), (version.to_string(), checksum.to_string()));
    }
    Ok(out)
}

/// Extract a short kernel variant name from the concrete type. v0.1
/// uses `std::any::type_name`; a `Kernel::variant_name(&self)` trait
/// method is the natural follow-up when more than two kernels ship.
fn kernel_variant_name(kernel: &dyn crate::physics::gravity::kernel::Kernel) -> String {
    let name = std::any::type_name_of_val(kernel);
    let leaf = name.rsplit("::").next().unwrap_or(name);
    leaf.trim_end_matches("Kernel").to_string()
}

/// Minimal RFC 3339 timestamp without pulling chrono. Format:
/// `YYYY-MM-DDTHH:MM:SSZ`.
fn rfc3339_utc_now() -> String {
    let now =
        SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap_or_default().as_secs();
    let (year, month, day, hour, min, sec) = unix_to_civil(now as i64);
    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{min:02}:{sec:02}Z")
}

// Civil-from-days algorithm from H. S. Hinnant, "Howard Hinnant's date
// algorithms" (public domain). Avoids pulling chrono into the dep tree
// for a single timestamp.
fn unix_to_civil(secs: i64) -> (i32, u32, u32, u32, u32, u32) {
    let days = secs.div_euclid(86_400);
    let time_of_day = secs.rem_euclid(86_400) as u32;
    let hour = time_of_day / 3600;
    let min = (time_of_day / 60) % 60;
    let sec = time_of_day % 60;

    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i32 + era as i32 * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d, hour, min, sec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::physics::integrator::IntegratorKind;
    use crate::units::UnitSystem;

    #[test]
    fn parse_lock_extracts_name_version_checksum() {
        let lock = br#"
[[package]]
name = "alpha"
version = "1.2.3"
checksum = "abcdef"

[[package]]
name = "beta"
version = "0.1.0"
"#;
        let idx = parse_lock(lock).unwrap();
        assert_eq!(idx.len(), 2);
        assert_eq!(idx["alpha"], ("1.2.3".to_string(), "abcdef".to_string()));
        assert_eq!(idx["beta"], ("0.1.0".to_string(), String::new()));
    }

    #[test]
    fn resolve_lock_path_walks_up() {
        let resolved = resolve_lock_path(None).unwrap();
        assert!(resolved.file_name() == Some(std::ffi::OsStr::new("Cargo.lock")));
        assert!(resolved.exists());
    }

    #[test]
    fn rfc3339_format_shape() {
        let s = rfc3339_utc_now();
        // YYYY-MM-DDTHH:MM:SSZ → 20 chars
        assert_eq!(s.len(), 20);
        assert!(s.ends_with('Z'));
        assert_eq!(&s[4..5], "-");
        assert_eq!(&s[10..11], "T");
    }

    #[test]
    fn header_from_system_captures_bodies_and_units() {
        let sun = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
        let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
        let sys = System::new(vec![sun, earth], UnitSystem::canonical())
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(0.01);
        let h = header_from_system(&sys, 42, None).unwrap();
        assert_eq!(h.bodies.count, 2);
        assert_eq!(h.bodies.list.len(), 2);
        assert_eq!(h.reproducibility.seed, 42);
        assert_eq!(h.reproducibility.cargo_lock_blake3.len(), 64);
        assert!(h.operators.is_empty(), "no perturbations registered");
    }
}
