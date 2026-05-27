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
    let kernel_softening = {
        let eps_sq = kernel.epsilon_squared();
        if eps_sq > 0.0 { Some(eps_sq.sqrt()) } else { None }
    };
    let kernel_props = kernel.properties();

    let apsis_sha: &str = option_env!("APSIS_GIT_COMMIT").unwrap_or("");

    let mut operators = Vec::new();
    operators.extend(sys.hamiltonian_perturbations().iter().filter_map(|op| {
        let cit = op.citation()?;
        let req = op.kernel_requirements();
        Some(OperatorMeta {
            name: cit.crate_name.to_string(),
            version: cit.crate_version.to_string(),
            crate_hash: operator_crate_hash(cit.crate_name, &lock_index, apsis_sha),
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
        Some(OperatorMeta {
            name: cit.crate_name.to_string(),
            version: cit.crate_version.to_string(),
            crate_hash: operator_crate_hash(cit.crate_name, &lock_index, apsis_sha),
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
            git_sha: if apsis_sha.is_empty() { "unknown" } else { apsis_sha }.to_string(),
            created_utc: rfc3339_utc_now(),
            rustc_version: option_env!("APSIS_RUSTC_VERSION").unwrap_or("").to_string(),
            generated_by: format!("apsis {}", env!("CARGO_PKG_VERSION")),
        },
        reproducibility: Reproducibility { cargo_lock_blake3: lock_hash, seed },
        unit_system: UnitSystemMeta {
            g: units.g(),
            length: units.length_label().to_string(),
            mass: units.mass_label().to_string(),
            time: units.time_label().to_string(),
        },
        integrator: IntegratorMeta {
            kind: sys.integrator_kind().label().to_string(),
            dt_mode: sys.dt_mode().label().to_string(),
            initial_dt: sys.dt(),
            params: serde_json::Map::new(),
        },
        kernel: KernelMeta {
            variant: kernel_variant,
            softening: kernel_softening,
            exactness: Some(kernel_props.exactness),
            continuity: Some(kernel_props.continuity),
        },
        operators,
        bodies: BodiesMeta { count: body_meta.len(), list: body_meta },
    })
}

/// Resolve an operator's `crate_hash` entry for the record header.
///
/// Three sources, in priority order:
///
/// 1. **Registry checksum.** Operators pulled from crates.io / a git
///    registry have a `checksum` field in `Cargo.lock`. This is the
///    canonical content hash; emit as-is.
/// 2. **Workspace path dep.** Operators that live in the same Cargo
///    workspace as apsis core have no `checksum` in the lockfile.
///    Fall back to the workspace's git SHA prefixed with `workspace:`
///    so the field still identifies the source state — reproducing the
///    run requires checking out apsis at that SHA, which is exactly
///    what the lockfile + workspace path dep imply.
/// 3. **No source state known.** Path dep + no git (tarball, vendored,
///    sandboxed CI). Emit empty string; the runtime treats it as
///    "source unknown" and the reproducibility claim is weakened to
///    "the lockfile + the file tree at build time".
///
/// The `workspace:<sha>` form is machine-parseable: a verifier sees the
/// prefix and knows the hash is a workspace-scope SHA, not a per-crate
/// content hash.
fn operator_crate_hash(
    crate_name: &str,
    lock_index: &HashMap<String, (String, String)>,
    apsis_sha: &str,
) -> String {
    match lock_index.get(crate_name) {
        Some((_, checksum)) if !checksum.is_empty() => checksum.clone(),
        Some(_) if !apsis_sha.is_empty() => format!("workspace:{apsis_sha}"),
        _ => String::new(),
    }
}

/// Read the workspace `Cargo.lock` (walking up from CWD when
/// `lock_path` is `None`) and return its BLAKE3 hash as a 64-char
/// hex string. Same lookup contract as `header_from_system`; callers
/// that only need the hash (e.g. [`crate::core::system::System::cite`])
/// avoid building a whole [`Header`] just to read it.
pub fn lock_blake3(lock_path: Option<&Path>) -> Result<String, ProvenanceError> {
    let lock_path = resolve_lock_path(lock_path)?;
    let lock_bytes =
        std::fs::read(&lock_path).map_err(|e| ProvenanceError::LockRead(lock_path.clone(), e))?;
    Ok(blake3::hash(&lock_bytes).to_hex().to_string())
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

/// Short kernel-implementation label for the record header. Delegates
/// to [`Kernel::variant_name`], which each impl defines as a stable
/// static string. The label identifies the implementation; the regime
/// (softened vs exact) is read from `kernel.softening` populated from
/// [`Kernel::epsilon_squared`].
fn kernel_variant_name(kernel: &dyn crate::physics::gravity::kernel::Kernel) -> String {
    kernel.variant_name().to_string()
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

    fn lock(entries: &[(&str, &str, &str)]) -> HashMap<String, (String, String)> {
        entries
            .iter()
            .map(|(n, v, c)| ((*n).to_string(), ((*v).to_string(), (*c).to_string())))
            .collect()
    }

    #[test]
    fn operator_crate_hash_prefers_registry_checksum() {
        let idx = lock(&[("apsis-1pn", "0.1.0", "deadbeef".repeat(8).as_str())]);
        let h = operator_crate_hash("apsis-1pn", &idx, "wsha");
        assert_eq!(h, "deadbeef".repeat(8));
    }

    #[test]
    fn operator_crate_hash_falls_back_to_workspace_sha_when_no_checksum() {
        let idx = lock(&[("apsis-1pn", "0.1.0", "")]);
        let h = operator_crate_hash("apsis-1pn", &idx, "abc123");
        assert_eq!(h, "workspace:abc123");
    }

    #[test]
    fn operator_crate_hash_empty_when_no_checksum_and_no_workspace_sha() {
        let idx = lock(&[("apsis-1pn", "0.1.0", "")]);
        let h = operator_crate_hash("apsis-1pn", &idx, "");
        assert_eq!(h, "");
    }

    #[test]
    fn operator_crate_hash_empty_when_crate_not_in_lockfile() {
        let idx = lock(&[]);
        let h = operator_crate_hash("unlisted-crate", &idx, "abc123");
        assert_eq!(h, "");
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
