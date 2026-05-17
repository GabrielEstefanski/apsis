//! Mid-run snapshot resume produces a bit-equal continuation.

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::records::{
    DiagnosticCadence, Record, RecordHook, RecordPolicy, provenance::header_from_system,
    restore_into, resume::RestoreError,
};
use crate::units::UnitSystem;

fn sun_earth_jupiter() -> Vec<Body> {
    vec![
        Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
        Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0),
        Body::gas_giant(1e-3).at(5.2, 0.0).with_velocity(0.0, 0.439),
    ]
}

fn assert_bodies_bit_equal(a: &System, b: &System) {
    assert_eq!(a.bodies().len(), b.bodies().len());
    for (i, (ba, bb)) in a.bodies().iter().zip(b.bodies().iter()).enumerate() {
        for (label, x, y) in [
            ("pos_x", ba.pos_x, bb.pos_x),
            ("pos_y", ba.pos_y, bb.pos_y),
            ("pos_z", ba.pos_z, bb.pos_z),
            ("vel_x", ba.vel_x, bb.vel_x),
            ("vel_y", ba.vel_y, bb.vel_y),
            ("vel_z", ba.vel_z, bb.vel_z),
        ] {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "body {i} {label} mismatch: {x} != {y} ({} ulp)",
                (x.to_bits() as i64 - y.to_bits() as i64).abs()
            );
        }
    }
}

/// Run a recorded session for `total_steps` and return the opened
/// record. Snapshot policy emits every `snap_every` steps with resume
/// state captured alongside. `total_steps` MUST stay below 97 so the
/// System's periodic COM recentering never fires and confuses the
/// per-snapshot bit-equality contract.
fn record_session(
    kind: IntegratorKind,
    path: &std::path::Path,
    total_steps: usize,
    snap_every: u32,
) -> Record {
    assert!(total_steps < 97, "guard against COM-recentering at step 97");
    let _ = std::fs::remove_file(path);
    let mut sys = System::new(sun_earth_jupiter(), UnitSystem::canonical())
        .with_integrator(kind)
        .with_dt(1e-3);
    let header = header_from_system(&sys, 0, None).unwrap();
    let hook = RecordHook::with_header(path, header, RecordPolicy::EveryNSteps(snap_every))
        .unwrap()
        .with_diagnostics(DiagnosticCadence::Off)
        .with_resume_capture(true);
    sys.hooks_mut().register(0, Box::new(hook));
    for _ in 0..total_steps {
        sys.step();
    }
    drop(sys);
    Record::open(path).unwrap()
}

fn fresh_system(kind: IntegratorKind) -> System {
    System::new(sun_earth_jupiter(), UnitSystem::canonical()).with_integrator(kind).with_dt(1e-3)
}

#[test]
fn whfast_resume_yields_bit_equal_continuation() {
    let p = std::env::temp_dir().join("apsis-resume-whfast.apsis");
    let record = record_session(IntegratorKind::WHFast, &p, 80, 10);
    // Snapshots at steps 0/10/20/30/40/50/60/70/80 — idx 4 = step 40.
    let remaining = 80 - 40;

    let mut resumed = fresh_system(IntegratorKind::WHFast);
    restore_into(&mut resumed, &record, 4).unwrap();
    for _ in 0..remaining {
        resumed.step();
    }

    let mut straight = fresh_system(IntegratorKind::WHFast);
    for _ in 0..80 {
        straight.step();
    }

    assert_bodies_bit_equal(&resumed, &straight);
    let _ = std::fs::remove_file(&p);
}

#[test]
fn ias15_resume_yields_bit_equal_continuation() {
    let p = std::env::temp_dir().join("apsis-resume-ias15.apsis");
    let record = record_session(IntegratorKind::Ias15, &p, 50, 10);
    // Snapshots at 0/10/20/30/40/50 — idx 2 = step 20.
    let remaining = 50 - 20;

    let mut resumed = fresh_system(IntegratorKind::Ias15);
    restore_into(&mut resumed, &record, 2).unwrap();
    for _ in 0..remaining {
        resumed.step();
    }

    let mut straight = fresh_system(IntegratorKind::Ias15);
    for _ in 0..50 {
        straight.step();
    }

    assert_bodies_bit_equal(&resumed, &straight);
    let _ = std::fs::remove_file(&p);
}

#[test]
fn mercurius_resume_yields_bit_equal_continuation() {
    let p = std::env::temp_dir().join("apsis-resume-mercurius.apsis");
    let record = record_session(IntegratorKind::Mercurius, &p, 60, 10);
    // Snapshots at 0/10/20/30/40/50/60 — idx 3 = step 30.
    let remaining = 60 - 30;

    let mut resumed = fresh_system(IntegratorKind::Mercurius);
    restore_into(&mut resumed, &record, 3).unwrap();
    for _ in 0..remaining {
        resumed.step();
    }

    let mut straight = fresh_system(IntegratorKind::Mercurius);
    for _ in 0..60 {
        straight.step();
    }

    assert_bodies_bit_equal(&resumed, &straight);
    let _ = std::fs::remove_file(&p);
}

#[test]
fn body_count_mismatch_errors() {
    let p = std::env::temp_dir().join("apsis-resume-mismatch.apsis");
    let record = record_session(IntegratorKind::WHFast, &p, 30, 10);

    let bodies =
        vec![Body::star(1.0).at(0.0, 0.0), Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0)];
    let mut sys = System::new(bodies, UnitSystem::canonical())
        .with_integrator(IntegratorKind::WHFast)
        .with_dt(1e-3);
    let err = restore_into(&mut sys, &record, 0).unwrap_err();
    assert!(
        matches!(err, RestoreError::BodyCountMismatch { expected: 2, found: 3 }),
        "got {err:?}"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn integrator_kind_mismatch_errors() {
    let p = std::env::temp_dir().join("apsis-resume-kindmismatch.apsis");
    let record = record_session(IntegratorKind::WHFast, &p, 30, 10);

    let mut sys = fresh_system(IntegratorKind::Ias15);
    let err = restore_into(&mut sys, &record, 0).unwrap_err();
    assert!(matches!(err, RestoreError::IntegratorMismatch { .. }), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn missing_resume_capture_errors() {
    let p = std::env::temp_dir().join("apsis-resume-nocapture.apsis");
    let _ = std::fs::remove_file(&p);
    let mut sys = fresh_system(IntegratorKind::WHFast);
    let header = header_from_system(&sys, 0, None).unwrap();
    let hook = RecordHook::with_header(&p, header, RecordPolicy::EveryNSteps(10)).unwrap();
    sys.hooks_mut().register(0, Box::new(hook));
    for _ in 0..30 {
        sys.step();
    }
    drop(sys);
    let record = Record::open(&p).unwrap();

    let mut fresh = fresh_system(IntegratorKind::WHFast);
    let err = restore_into(&mut fresh, &record, 0).unwrap_err();
    assert!(matches!(err, RestoreError::MissingResumeState { .. }), "got {err:?}");
    let _ = std::fs::remove_file(&p);
}
