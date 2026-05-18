//! `DiagnosticCadence` end-to-end: writer emits Diagnostic frames at the
//! configured cadence; reader observes them with the expected drift
//! values.

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::records::{
    DiagnosticCadence, Record, RecordHook, RecordPolicy, provenance::header_from_system,
};
use crate::units::UnitSystem;

fn run(cadence: DiagnosticCadence, steps: usize, path: &std::path::Path) -> Record {
    let sun = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
    let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
    let mut sys = System::new(vec![sun, earth], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3);
    let header = header_from_system(&sys, 0, None).unwrap();
    let hook = RecordHook::with_header(path, header, RecordPolicy::BookendsAndEvents)
        .unwrap()
        .with_diagnostics(cadence);
    sys.hooks_mut().register(0, Box::new(hook));
    for _ in 0..steps {
        sys.step();
    }
    drop(sys);
    Record::open(path).unwrap()
}

#[test]
fn off_emits_no_diagnostics() {
    let p = std::env::temp_dir().join("apsis-diag-off.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run(DiagnosticCadence::Off, 50, &p);
    assert_eq!(rec.diagnostics().unwrap().count(), 0);
    let _ = std::fs::remove_file(&p);
}

#[test]
fn every_n_steps_emits_anchor_plus_periodic() {
    let p = std::env::temp_dir().join("apsis-diag-everyn.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run(DiagnosticCadence::EveryNSteps(10), 50, &p);
    // t=0 anchor + steps 10/20/30/40/50.
    let count = rec.diagnostics().unwrap().count();
    assert_eq!(count, 6, "EveryNSteps(10)/50steps emitted {count} diagnostics; expected 6");
    let _ = std::fs::remove_file(&p);
}

#[test]
fn anchor_has_zero_drift() {
    let p = std::env::temp_dir().join("apsis-diag-anchor.apsis");
    let _ = std::fs::remove_file(&p);
    let rec = run(DiagnosticCadence::EveryNSteps(10), 20, &p);
    let first =
        rec.diagnostics().unwrap().next().expect("anchor diagnostic missing").expect("read err");
    assert_eq!(first.t, 0.0);
    assert_eq!(first.d_energy_rel, 0.0);
    assert_eq!(first.d_lz_rel, 0.0);
    let _ = std::fs::remove_file(&p);
}
