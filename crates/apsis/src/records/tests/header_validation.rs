//! Invalid headers / bad magic / format mismatch / truncation → `Record::open`
//! returns the appropriate `RecordError`.

use crate::records::reader::{Record, RecordError};

fn write_bytes(name: &str, bytes: &[u8]) -> std::path::PathBuf {
    let p = std::env::temp_dir().join(format!("apsis-validation-{name}.apsis"));
    std::fs::write(&p, bytes).unwrap();
    p
}

#[test]
fn bad_magic_errors() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"NOTAPSR0");
    bytes.extend_from_slice(&[0u8; 16]);
    let p = write_bytes("magic", &bytes);
    let err = Record::open(&p).unwrap_err();
    assert!(matches!(err, RecordError::BadMagic), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn unsupported_format_version_errors() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"APSR");
    bytes.extend_from_slice(&999u16.to_le_bytes());
    bytes.extend_from_slice(&0u16.to_le_bytes());
    bytes.extend_from_slice(&0u64.to_le_bytes());
    let p = write_bytes("ver", &bytes);
    let err = Record::open(&p).unwrap_err();
    assert!(matches!(err, RecordError::UnsupportedFormatVersion(999)), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn truncated_header_errors() {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"APSR");
    bytes.extend_from_slice(&1u16.to_le_bytes());
    bytes.extend_from_slice(&0u16.to_le_bytes());
    bytes.extend_from_slice(&100u64.to_le_bytes()); // claims 100 bytes; file ends
    let p = write_bytes("trunc", &bytes);
    let err = Record::open(&p).unwrap_err();
    assert!(matches!(err, RecordError::Io(_)), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn missing_trailer_errors() {
    use crate::records::frame::{Frame, Snapshot};
    use crate::records::header::{
        Apsis, BodiesMeta, Header, IntegratorMeta, KernelMeta, Reproducibility, UnitSystemMeta,
    };

    let header = Header {
        apsis: Apsis {
            version: "0.1.0".into(),
            git_sha: "x".into(),
            created_utc: "2026-05-16T00:00:00Z".into(),
        },
        reproducibility: Reproducibility { cargo_lock_blake3: "00".repeat(32), seed: 0 },
        unit_system: UnitSystemMeta {
            g: 1.0,
            length: "AU".into(),
            mass: "M_sun".into(),
            time: "yr/2pi".into(),
        },
        integrator: IntegratorMeta {
            kind: "Ias15".into(),
            dt_mode: "Fixed".into(),
            initial_dt: 1e-3,
            params: Default::default(),
        },
        kernel: KernelMeta { variant: "Newton".into(), softening: None },
        operators: vec![],
        bodies: BodiesMeta { count: 0, list: vec![] },
    };
    let toml = header.to_toml().unwrap();
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"APSR");
    bytes.extend_from_slice(&1u16.to_le_bytes());
    bytes.extend_from_slice(&0u16.to_le_bytes());
    bytes.extend_from_slice(&(toml.len() as u64).to_le_bytes());
    bytes.extend_from_slice(toml.as_bytes());
    let snap = Frame::Snapshot(Snapshot { t: 0.0, bodies: vec![] });
    snap.write(&mut bytes).unwrap();
    // No trailer.
    let p = write_bytes("no_trailer", &bytes);
    let err = Record::open(&p).unwrap_err();
    assert!(matches!(err, RecordError::MissingTrailer), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}
