//! Paper claim gate: two runs with the same seed + config produce
//! byte-equal record frame streams (everything past the TOML header,
//! which carries a per-run `created_utc` timestamp).

use crate::core::system::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::records::{Record, RecordHook, RecordPolicy, provenance::header_from_system};
use crate::units::UnitSystem;

fn run_with_seed(seed: u64, path: &std::path::Path, steps: usize) {
    let sun = Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0);
    let earth = Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0);
    let mut sys = System::new(vec![sun, earth], UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3)
        .with_seed(seed);

    let header = header_from_system(&sys, seed, None).unwrap();
    let hook = RecordHook::with_header(path, header, RecordPolicy::BookendsAndEvents).unwrap();
    sys.hooks_mut().register(0, Box::new(hook));

    for _ in 0..steps {
        sys.step();
    }
    drop(sys);
}

#[test]
fn two_runs_same_seed_have_equal_frame_streams() {
    let a = std::env::temp_dir().join("apsis-repro-a.apsis");
    let b = std::env::temp_dir().join("apsis-repro-b.apsis");
    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);

    run_with_seed(42, &a, 100);
    run_with_seed(42, &b, 100);

    let rec_a = Record::open(&a).unwrap();
    let rec_b = Record::open(&b).unwrap();

    // Header equality except for the per-run `created_utc` timestamp.
    let mut ha = rec_a.header().clone();
    let mut hb = rec_b.header().clone();
    ha.apsis.created_utc.clear();
    hb.apsis.created_utc.clear();
    assert_eq!(ha, hb, "headers differ in non-timestamp fields");

    // Frame-stream equality from frames_start onward.
    let bytes_a = std::fs::read(&a).unwrap();
    let bytes_b = std::fs::read(&b).unwrap();
    let hl_a = u64::from_le_bytes(bytes_a[8..16].try_into().unwrap()) as usize;
    let hl_b = u64::from_le_bytes(bytes_b[8..16].try_into().unwrap()) as usize;
    assert_eq!(
        &bytes_a[16 + hl_a..],
        &bytes_b[16 + hl_b..],
        "frame stream differs between two same-seed runs"
    );

    let _ = std::fs::remove_file(&a);
    let _ = std::fs::remove_file(&b);
}
