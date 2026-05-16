//! Each `RecordPolicy` variant emits the expected number of Snapshot frames.

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::records::{Record, RecordHook, RecordPolicy, provenance::header_from_system};
use crate::units::UnitSystem;

fn run_policy(policy: RecordPolicy, steps: usize, path: &std::path::Path) -> Record {
    let sun = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
    let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
    let mut sys = System::new(vec![sun, earth], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3);
    let header = header_from_system(&sys, 0, None).unwrap();
    let hook = RecordHook::with_header(path, header, policy).unwrap();
    sys.hooks_mut().register(0, Box::new(hook));
    for _ in 0..steps {
        sys.step();
    }
    drop(sys);
    Record::open(path).unwrap()
}

#[test]
fn bookends_only_emits_two_snapshots() {
    let p = std::env::temp_dir().join("apsis-policy-bookends.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run_policy(RecordPolicy::BookendsAndEvents, 50, &p);
    let count = rec.dense().unwrap().count();
    assert_eq!(count, 2, "BookendsAndEvents emitted {count} snapshots; expected 2");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn every_n_steps_emits_expected_count() {
    let p = std::env::temp_dir().join("apsis-policy-everyn.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run_policy(RecordPolicy::EveryNSteps(10), 50, &p);
    // Initial bookend + 5 periodic (steps 10, 20, 30, 40, 50) + final bookend.
    // The final bookend may coincide with the step=50 periodic snapshot;
    // RecordHook dedups by t to avoid emitting twice at the same instant.
    let count = rec.dense().unwrap().count();
    assert!(
        (6..=7).contains(&count),
        "EveryNSteps(10)/50steps emitted {count} snapshots; expected 6 or 7"
    );
    let _ = std::fs::remove_file(&p);
}

#[test]
fn dense_emits_one_per_step() {
    let p = std::env::temp_dir().join("apsis-policy-dense.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run_policy(RecordPolicy::Dense, 20, &p);
    let count = rec.dense().unwrap().count();
    // Initial bookend at t=0 + per-step snapshots (20).
    assert!(count >= 20, "Dense/20 steps emitted {count}; expected >= 20");
    let _ = std::fs::remove_file(&p);
}
