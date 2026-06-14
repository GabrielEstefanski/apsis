//! Exact-finish-time semantics of `System::integrate_for` /
//! `integrate_until`: the loop clips the final step so the returned
//! state is at exactly the requested time. Before this contract,
//! fixed-time measurements sampled the state up to one step past
//! `t_end` — on the Mercury 1PN gate that endpoint-sampling error was
//! the dominant residual (ADR-015).

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

/// Two-body Kepler scenario (e ≈ 0.3, a = 1) valid for every
/// integrator, including the hierarchy-gated Wisdom-Holman family.
fn kepler_system(kind: IntegratorKind, dt: f64) -> System {
    let sun = Body::star(1.0);
    let r0 = 0.7;
    let v0 = (2.0 / r0 - 1.0_f64).sqrt();
    let planet = Body::rocky(1.0e-6).at(r0, 0.0).with_velocity(0.0, v0);
    System::new(vec![sun, planet], UnitSystem::solar_canonical()).with_integrator(kind).with_dt(dt)
}

/// Duration deliberately incommensurate with every dt used below, so a
/// non-clipped final step always overshoots.
const DURATION: f64 = 2.37;

#[test]
fn ias15_lands_exactly_on_t_end() {
    let mut sys = kepler_system(IntegratorKind::Ias15, 1e-4);
    sys.integrate_for(DURATION);
    assert_eq!(
        sys.t(),
        DURATION,
        "IAS15 endpoint t = {} differs from requested {}",
        sys.t(),
        DURATION
    );
}

#[test]
fn fixed_step_integrators_land_exactly_on_t_end() {
    for kind in [
        IntegratorKind::VelocityVerlet,
        IntegratorKind::Yoshida4,
        IntegratorKind::ImplicitMidpoint,
        IntegratorKind::WisdomHolman,
        IntegratorKind::WHFast,
        IntegratorKind::Mercurius,
    ] {
        let mut sys = kepler_system(kind, 1e-3);
        sys.integrate_for(DURATION);
        assert_eq!(
            sys.t(),
            DURATION,
            "{:?} endpoint t = {} differs from requested {}",
            kind,
            sys.t(),
            DURATION
        );
    }
}

#[test]
fn clipped_final_step_keeps_energy_at_roundoff_floor() {
    // IAS15 on Kepler two-body: |dE/E| ~ 1e-13..1e-15 (physics
    // toolkit floor); 10x headroom. The clipped final step must not
    // disturb it.
    let mut sys = kepler_system(IntegratorKind::Ias15, 1e-4);
    sys.integrate_for(DURATION);
    let de = sys.rel_energy_error().expect("energy tracked after stepping");
    assert!(de.abs() < 1e-12, "post-clip |dE/E| = {:.3e} above 10x IAS15 floor", de.abs());
}

#[test]
fn opt_out_restores_overshoot_semantics() {
    let mut sys = kepler_system(IntegratorKind::Ias15, 1e-4);
    sys.set_exact_finish_time(false);
    sys.integrate_for(DURATION);
    assert!(sys.t() > DURATION, "opt-out should run whole steps past t_end (t = {})", sys.t());
}

#[test]
fn segmented_sampling_preserves_controller_rhythm() {
    // A clipped step must not poison the IAS15 controller: `dt_next`
    // is restored after the clip, so each segment boundary costs at
    // most the clipped step itself plus one cold-Picard step. Without
    // the restore, `DT_GROWTH_LIMIT` (7x per accept) forces a
    // multi-step crawl back from the tiny clipped dt at every
    // boundary, and this bound breaks.
    let segments = 50;
    let mut one_shot = kepler_system(IntegratorKind::Ias15, 1e-4);
    one_shot.integrate_for(DURATION);

    // Absolute targets: relative `integrate_for` durations would
    // accumulate the fp round-off of the duration sum itself.
    let mut segmented = kepler_system(IntegratorKind::Ias15, 1e-4);
    for k in 1..=segments {
        segmented.integrate_until(DURATION * k as f64 / segments as f64);
    }

    assert_eq!(one_shot.t(), segmented.t(), "both paths must land on the same t");
    assert!(
        segmented.steps() <= one_shot.steps() + 2 * segments,
        "segmented sampling took {} steps vs {} one-shot (+{} boundaries allowed)",
        segmented.steps(),
        one_shot.steps(),
        2 * segments
    );
}

#[test]
fn integrate_until_past_target_is_a_no_op() {
    let mut sys = kepler_system(IntegratorKind::Ias15, 1e-4);
    sys.integrate_for(DURATION);
    let steps_before = sys.steps();
    sys.integrate_until(DURATION / 2.0);
    assert_eq!(sys.t(), DURATION, "t moved on a past-target call");
    assert_eq!(sys.steps(), steps_before, "steps taken on a past-target call");
}
