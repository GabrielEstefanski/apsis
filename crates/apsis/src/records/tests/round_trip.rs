//! Write a small record, read it back, assert structural integrity.

use crate::core::hooks::{HookContext, HookPhase, HookPhaseKind, SimHook};
use crate::domain::body::Body;
use crate::records::header::{
    Apsis, BodiesMeta, IntegratorMeta, KernelMeta, Reproducibility, UnitSystemMeta,
};
use crate::records::{Header, Record, RecordHook, RecordPolicy};

fn tmp(name: &str) -> std::path::PathBuf {
    let p = std::env::temp_dir().join(format!("apsis-rt-{name}.apsis"));
    let _ = std::fs::remove_file(&p);
    p
}

fn minimal_header() -> Header {
    Header {
        apsis: Apsis {
            version: "0.1.0".into(),
            git_sha: "test".into(),
            created_utc: "2026-05-16T00:00:00Z".into(),
            rustc_version: "".into(),
            generated_by: "apsis-test".into(),
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
    }
}

fn make_ctx<'a>(bodies: &'a [Body], names: &'a [String], t: f64, steps: u64) -> HookContext<'a> {
    HookContext {
        bodies,
        names,
        t,
        dt: 1e-3,
        steps,
        rel_energy_error: 0.0,
        rel_angular_momentum_error: 0.0,
        phase: HookPhase(HookPhaseKind::PreStep),
    }
}

#[test]
fn write_then_open_round_trip() {
    let path = tmp("basic");
    {
        let mut hook =
            RecordHook::with_header(&path, minimal_header(), RecordPolicy::BookendsAndEvents)
                .unwrap();
        let bodies = vec![Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0)];
        let names = vec!["sun".to_string()];
        let ctx = make_ctx(&bodies, &names, 0.0, 0);
        hook.pre_step(&ctx);
        hook.on_finish(&ctx);
    }
    let rec = Record::open(&path).unwrap();
    assert_eq!(rec.header().apsis.version, "0.1.0");
    assert_eq!(rec.header().reproducibility.cargo_lock_blake3, "00".repeat(32));
    let _ = std::fs::remove_file(&path);
}

#[test]
fn round_trip_with_post_step_advances() {
    let path = tmp("advance");
    {
        let mut hook =
            RecordHook::with_header(&path, minimal_header(), RecordPolicy::BookendsAndEvents)
                .unwrap();
        let bodies = vec![Body::star(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0)];
        let names = vec!["sun".to_string()];
        hook.pre_step(&make_ctx(&bodies, &names, 0.0, 0));
        for s in 1..=10 {
            hook.post_step(&make_ctx(&bodies, &names, s as f64 * 1e-3, s));
        }
        hook.on_finish(&make_ctx(&bodies, &names, 10e-3, 10));
    }
    let rec = Record::open(&path).unwrap();
    let (init, fin) = rec.bookends().unwrap();
    assert_eq!(init.t, 0.0);
    assert!((fin.t - 10e-3).abs() < 1e-12);
    let _ = std::fs::remove_file(&path);
}
