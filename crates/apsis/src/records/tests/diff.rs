//! `Record::diff` semantic categorisation.

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::records::diff::HeaderChange;
use crate::records::{Record, RecordHook, RecordPolicy, provenance::header_from_system};
use crate::units::UnitSystem;

fn write_run(path: &std::path::Path, seed: u64, integrator: IntegratorKind, dt: f64) -> Record {
    let sun = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
    let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
    let mut sys = System::new(vec![sun, earth], UnitSystem::canonical())
        .with_integrator(integrator)
        .with_dt(dt);
    let header = header_from_system(&sys, seed, None).unwrap();
    let hook = RecordHook::with_header(path, header, RecordPolicy::BookendsAndEvents).unwrap();
    sys.hooks_mut().register(0, Box::new(hook));
    for _ in 0..50 {
        sys.step();
    }
    drop(sys);
    Record::open(path).unwrap()
}

#[test]
fn identical_runs_diff_empty() {
    let a = std::env::temp_dir().join("apsis-diff-id-a.apsis");
    let b = std::env::temp_dir().join("apsis-diff-id-b.apsis");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
    let ra = write_run(&a, 7, IntegratorKind::Ias15, 1e-3);
    let rb = write_run(&b, 7, IntegratorKind::Ias15, 1e-3);
    let d = ra.diff(&rb).unwrap();
    assert!(d.header.is_empty(), "unexpected header changes: {:?}", d.header);
    assert!(d.frames.trailer_blake3_match);
    assert_eq!(d.frames.trajectory_rms_at_final, Some(0.0));
    assert!(d.is_empty());
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
}

#[test]
fn seed_change_is_categorised() {
    let a = std::env::temp_dir().join("apsis-diff-seed-a.apsis");
    let b = std::env::temp_dir().join("apsis-diff-seed-b.apsis");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
    let ra = write_run(&a, 7, IntegratorKind::Ias15, 1e-3);
    let rb = write_run(&b, 11, IntegratorKind::Ias15, 1e-3);
    let d = ra.diff(&rb).unwrap();
    assert!(
        d.header.iter().any(|c| matches!(c, HeaderChange::SeedChanged { before: 7, after: 11 })),
        "expected SeedChanged, got {:?}",
        d.header
    );
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
}

#[test]
fn integrator_change_categorised_and_trajectory_diverges() {
    let a = std::env::temp_dir().join("apsis-diff-int-a.apsis");
    let b = std::env::temp_dir().join("apsis-diff-int-b.apsis");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
    let ra = write_run(&a, 7, IntegratorKind::Ias15, 1e-3);
    let rb = write_run(&b, 7, IntegratorKind::VelocityVerlet, 1e-3);
    let d = ra.diff(&rb).unwrap();
    assert!(
        d.header.iter().any(|c| matches!(c, HeaderChange::IntegratorKindChanged { .. })),
        "expected IntegratorKindChanged, got {:?}",
        d.header
    );
    assert!(!d.frames.trailer_blake3_match);
    let rms = d.frames.trajectory_rms_at_final.expect("equal body counts → rms defined");
    assert!(rms > 0.0, "trajectory rms must be positive when integrators differ");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
}
