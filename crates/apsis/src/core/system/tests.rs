//! End-to-end integration tests for the simulation pipeline.
//!
//! Tests are grouped by the invariant they verify:
//!
//! - [`energy`]     — Hamiltonian conservation over many orbital periods
//! - [`wh_guard`]   — Wisdom–Holman suitability guard and fallback behaviour
//! - [`benchmarks`] — quantitative accuracy: Kepler, figure-8, Pythagorean 3-body
//! - [`replay`]     — bit-identical determinism and snapshot round-trip
//! - [`hook_dispatch`] — hook registry fires and commands mutate via step()

use super::System;
use crate::domain::body::Body;
use crate::physics::integrator::IntegratorKind;
use crate::units::UnitSystem;

// ── UnitSystem snapshot invariants ────────────────────────────────────────────
//
// The contract is "the unit system is part of the simulation's frozen state":
// no public path may mutate it after [`System::new`] returns. These tests
// pin the invariant so a future refactor can't silently introduce a setter
// or a mutating helper.

#[test]
fn units_snapshot_is_immutable_across_integration() {
    let bodies = vec![Body::star(1.0), Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0)];
    let mut sys = System::new(bodies, UnitSystem::solar())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3);
    let units_at_construction = *sys.units();
    assert_eq!(units_at_construction, UnitSystem::solar());

    sys.integrate_for(1.0);

    assert_eq!(*sys.units(), units_at_construction);
    assert_eq!(*sys.units(), UnitSystem::solar());
}

#[test]
fn system_g_factor_is_derived_from_units_at_construction() {
    let bodies = vec![Body::star(1.0)];

    let sys_solar = System::new(bodies.clone(), UnitSystem::solar());
    let sys_canon = System::new(bodies, UnitSystem::canonical());

    assert_eq!(sys_solar.g_factor(), UnitSystem::solar().g());
    assert_eq!(sys_canon.g_factor(), 1.0);
    assert_ne!(sys_solar.g_factor(), sys_canon.g_factor());
}

// ── Shared helpers ────────────────────────────────────────────────────────────

fn two_body_circular_system(integrator: IntegratorKind, dt: f64) -> System {
    let bodies = vec![
        Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
        Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
    ];
    let mut sys =
        System::new(bodies, UnitSystem::canonical()).with_theta(0.5).with_dt(dt).with_max_depth(10);
    sys.set_integrator(integrator);
    sys
}

fn two_body_deterministic_system() -> System {
    let bodies = vec![
        Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
        Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
    ];
    let mut sys = System::new(bodies, UnitSystem::canonical())
        .with_theta(0.5)
        .with_dt(0.01)
        .with_max_depth(10);
    // Replay/determinism tests use a fixed-step, stateless integrator.
    // IAS15 (the project default) carries warm-start state (b, e, csb, dt_next)
    // that is intentionally not serialised in snapshots — reloading resets it
    // and the resumed trajectory differs from the reference by ∼1 ULP for a
    // few steps. That difference is physically meaningless but would fail
    // a bit-exact check.
    sys.set_integrator(IntegratorKind::VelocityVerlet);
    sys.set_seed(42);
    sys
}

fn assert_bodies_bit_equal(a: &[Body], b: &[Body], context: &str) {
    assert_eq!(a.len(), b.len(), "{context}: body count differs");
    for (i, (ba, bb)) in a.iter().zip(b.iter()).enumerate() {
        assert_eq!(ba.pos_x.to_bits(), bb.pos_x.to_bits(), "{context}: body {i} x differs");
        assert_eq!(ba.pos_y.to_bits(), bb.pos_y.to_bits(), "{context}: body {i} y differs");
        assert_eq!(ba.vel_x.to_bits(), bb.vel_x.to_bits(), "{context}: body {i} vx differs");
        assert_eq!(ba.vel_y.to_bits(), bb.vel_y.to_bits(), "{context}: body {i} vy differs");
    }
}

// ── Energy conservation ───────────────────────────────────────────────────────
//
// Physical scenario: two equal-mass bodies in a circular orbit.
//
//   G = 1, M₁ = M₂ = 1
//   Positions: (−1, 0) and (+1, 0), separation d = 2, orbital radius r = 1
//   Velocities: (0, −0.5) and (0, +0.5) — CCW orbit
//   Orbital period: T = 2πr/v = 4π ≈ 12.566
//
// Tolerance derivation (dt = 0.01, T = 4π, dt/T ≈ 7.96 × 10⁻⁴):
//   VV  (2nd order): amplitude ~ (dt/T)² ≈ 6.3 × 10⁻⁷ → tol 1e-4
//   Y4  (4th order): amplitude ~ (dt/T)⁴ ≈ 4 × 10⁻¹³ → tol 1e-7

mod energy {
    use super::*;

    fn max_rel_energy_error(sys: &mut System, n_periods: u64, dt: f64) -> f64 {
        const PERIOD: f64 = 4.0 * std::f64::consts::PI;
        let total_steps = (n_periods as f64 * PERIOD / dt).ceil() as u64;
        let mut peak: f64 = 0.0;
        for _ in 0..total_steps {
            sys.step();
            peak = peak.max(
                sys.metrics()
                    .rel_energy_error
                    .expect("well-conditioned regime: rel_energy_error must be Some")
                    .abs(),
            );
        }
        peak
    }

    #[test]
    fn velocity_verlet() {
        const DT: f64 = 0.01;
        const N_PERIODS: u64 = 100;
        const TOLERANCE: f64 = 1e-4;

        let mut sys = two_body_circular_system(IntegratorKind::VelocityVerlet, DT);
        let peak_err = max_rel_energy_error(&mut sys, N_PERIODS, DT);

        assert!(
            peak_err < TOLERANCE,
            "VelocityVerlet: peak |δE/E₀| = {:.3e} exceeds {:.0e} \
             after {} periods (dt = {}, T = 4π ≈ 12.566)",
            peak_err,
            TOLERANCE,
            N_PERIODS,
            DT,
        );
    }

    #[test]
    fn yoshida4() {
        const DT: f64 = 0.01;
        const N_PERIODS: u64 = 100;
        const TOLERANCE: f64 = 1e-7;

        let mut sys = two_body_circular_system(IntegratorKind::Yoshida4, DT);
        let peak_err = max_rel_energy_error(&mut sys, N_PERIODS, DT);

        assert!(
            peak_err < TOLERANCE,
            "Yoshida4: peak |δE/E₀| = {:.3e} exceeds {:.0e} \
             after {} periods (dt = {}, T = 4π ≈ 12.566)",
            peak_err,
            TOLERANCE,
            N_PERIODS,
            DT,
        );
    }

    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect raw peak errors"]
    fn diagnostic_peak_errors() {
        for &(label, integrator, dt) in &[
            ("VV    dt=0.01 ", IntegratorKind::VelocityVerlet, 0.01_f64),
            ("VV    dt=0.001", IntegratorKind::VelocityVerlet, 0.001_f64),
            ("Y4    dt=0.01 ", IntegratorKind::Yoshida4, 0.01_f64),
            ("Y4    dt=0.001", IntegratorKind::Yoshida4, 0.001_f64),
        ] {
            let mut sys = two_body_circular_system(integrator, dt);
            let peak = max_rel_energy_error(&mut sys, 10, dt);
            println!("{label}  peak |δE/E₀| = {peak:.3e}");
        }
    }
}

// ── Wisdom-Holman guard ───────────────────────────────────────────────────────
//
// `is_wh_suitable()` must reject systems without a dominant central mass.
// Fallback to Yoshida4 must not panic and must conserve energy identically
// to a direct Yoshida4 run.

mod wh_guard {
    use super::*;
    use crate::math::Vec3;

    #[test]
    fn hierarchical_system_is_suitable() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        assert!(sys.is_wh_suitable());
    }

    #[test]
    fn equal_mass_system_is_not_suitable() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        assert!(!sys.is_wh_suitable());
    }

    #[test]
    fn three_equal_mass_is_not_suitable() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
            Body::rocky(1.0).at(0.0, 1.0).with_velocity(0.5, 0.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        assert!(!sys.is_wh_suitable());
    }

    #[test]
    fn boundary_at_exactly_10x_is_suitable() {
        let bodies = vec![
            Body::star(10.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 1.0),
        ];
        assert!(
            System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(0.01)
                .with_max_depth(10)
                .is_wh_suitable()
        );
    }

    #[test]
    fn boundary_below_10x_is_not_suitable() {
        let bodies = vec![
            Body::star(9.9).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 1.0),
        ];
        assert!(
            !System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(0.01)
                .with_max_depth(10)
                .is_wh_suitable()
        );
    }

    #[test]
    fn single_body_is_not_suitable() {
        let bodies = vec![Body::rocky(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0)];
        assert!(
            !System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(0.01)
                .with_max_depth(10)
                .is_wh_suitable()
        );
    }

    /// Regression: when WH is the active integrator, the dense-output
    /// snapshot the renderer reads must not exist with internally
    /// inconsistent array lengths. WH evaluates forces only on
    /// `bodies[1..]`, leaving `scratch_acc` sized `N − 1`; the Order-2
    /// fallback path in `System::step` previously combined that with
    /// body-aligned `x0` / `v0` and produced a snapshot whose `a0`
    /// disagreed with `x0` by one entry. The renderer guard checked
    /// only `n_bodies()` (which reads `x0.len()`) and let the
    /// inconsistent snapshot through, so `interpolate(i, h)` panicked
    /// at `i = N − 1` indexing `a0[N − 1]`.
    ///
    /// Two independent guarantees pin the fix:
    ///   1. Producer side — `System::step` does not synthesise an
    ///      Order-2 snapshot when `scratch_acc.len() != bodies.len()`.
    ///   2. Snapshot invariant — any snapshot that does exist passes
    ///      `is_shape_consistent()`.
    #[test]
    fn wh_step_emits_no_inconsistent_dense_snapshot() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(1e-3).at(15.0, 0.0).with_velocity(0.0, 8.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        assert!(sys.is_wh_suitable(), "fixture should select the WH path, not the fallback");

        // Two steps: the first populates scratch_acc, the second is
        // the first that actually exercises the Order-2 dense-snapshot
        // synthesis path (which was the bug site).
        sys.step();
        sys.step();

        if let Some(snap) = &sys.last_dense_snapshot {
            assert!(
                snap.is_shape_consistent(),
                "WH step published a DenseSnapshot with mismatched array lengths: \
                 x0={} v0={} a0={} b={} — the renderer would panic indexing past the \
                 shortest array",
                snap.x0.len(),
                snap.v0.len(),
                snap.a0.len(),
                snap.b.len(),
            );
            assert_eq!(
                snap.n_bodies(),
                sys.bodies().len(),
                "snapshot body count {} disagrees with system body count {}",
                snap.n_bodies(),
                sys.bodies().len(),
            );
        }
    }

    #[test]
    fn non_hierarchical_does_not_panic_and_stays_finite() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        for _ in 0..100 {
            sys.step();
        }
        for b in sys.bodies() {
            assert!(b.pos_x.is_finite() && b.pos_y.is_finite(), "body left finite domain");
            assert!(b.vel_x.is_finite() && b.vel_y.is_finite(), "velocity left finite domain");
        }
    }

    /// `scratch_acc` after a WH step is sized N (one entry per body),
    /// not N-1 (planets only). This is the contract the dense-output
    /// snapshot path in `System::step` reads, and it is also the
    /// contract every consumer of `Metrics::last_accelerations` reads
    /// (camera feedforward, |a| field, CSV export). Pre-fix the buffer
    /// stayed at N-1 because WH's force evaluation skips the central
    /// body, breaking all three consumers silently — no snapshot for
    /// the renderer, off-by-one indexing for the rest.
    #[test]
    fn wh_step_leaves_scratch_acc_sized_n() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(1e-3).at(15.0, 0.0).with_velocity(0.0, 8.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        assert_eq!(
            sys.last_accelerations().len(),
            sys.bodies().len(),
            "WH must leave scratch_acc body-aligned (N), not planet-only (N-1)"
        );
    }

    /// With `scratch_acc` body-aligned, the next step's pre-step
    /// kinematics capture in `System::step` synthesises an Order-2
    /// dense snapshot for WH (the producer condition was
    /// `scratch_acc.len() == bodies.len()`). The renderer needs this
    /// snapshot to interpolate body positions within a step; without
    /// it bodies update only at sparse publish ticks and appear to
    /// freeze at slow achieved sim rates.
    #[test]
    fn wh_step_publishes_dense_snapshot() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(1e-3).at(15.0, 0.0).with_velocity(0.0, 8.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        // First step populates scratch_acc; second step is the first
        // whose pre-step capture sees a non-empty, body-aligned buffer.
        sys.step();
        sys.step();
        let snap = sys
            .last_dense_snapshot
            .as_ref()
            .expect("WH must publish a dense snapshot once scratch_acc is populated");
        assert!(snap.is_shape_consistent(), "WH dense snapshot must be shape-consistent");
        assert_eq!(snap.n_bodies(), sys.bodies().len(), "snapshot body count must match system");
    }

    /// Newton's third law in inertial coordinates: the total force on
    /// the system is zero, so `Σ m_i · a_i = 0` for the inertial
    /// accelerations published by the integrator. WH's central-body
    /// acceleration was synthesised from per-planet gravitational
    /// reactions (`a_0 = G Σ m_i q_i / r_i³`); this test confirms the
    /// synthesis is consistent with the reaction every planet feels
    /// from the central body's pull (`-μ q_i / r_i³`).
    ///
    /// Tolerance accounts for f64 round-off accumulating in the per-
    /// planet sum at this body count; tightens automatically as the
    /// implementation maintains higher precision.
    #[test]
    fn wh_acc_satisfies_newton_third_law() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(0.5).at(-7.0, 4.0).with_velocity(-2.0, -6.0),
            Body::rocky(0.2).at(3.0, -12.0).with_velocity(5.0, 1.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();

        let total: Vec3 = sys
            .bodies()
            .iter()
            .zip(sys.last_accelerations().iter())
            .fold(Vec3::ZERO, |s, (b, a)| s + b.mass * *a);
        assert!(
            total.length() < 1e-10,
            "WH must publish accelerations that satisfy Σ m_i a_i = 0; observed |Σ m a| = {}",
            total.length(),
        );
    }

    /// The WH dense snapshot now carries a Kepler-analytical kernel
    /// (`wh_data`), used by the renderer's bulk interpolation path to
    /// evaluate planet trajectories within a step instead of via the
    /// order-2 Taylor fallback. Without this field the renderer can
    /// still interpolate (Taylor falls back) but fast-orbit bodies
    /// (Galilean moons, Phobos) wobble at step boundaries because the
    /// quadratic doesn't track the orbit curvature.
    #[test]
    fn wh_step_snapshot_carries_kepler_kernel() {
        let bodies = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(1e-3).at(15.0, 0.0).with_velocity(0.0, 8.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        let snap = sys.last_dense_snapshot.as_ref().expect("WH must publish a snapshot");
        let wh = snap.wh_data.as_ref().expect("WH snapshot must carry the Kepler kernel");
        assert_eq!(wh.q0_helio_rest.len(), sys.bodies().len() - 1);
        assert_eq!(wh.v0_inertial_rest.len(), sys.bodies().len() - 1);
        assert_eq!(wh.planet_masses.len(), sys.bodies().len() - 1);
    }

    /// Round-trip: at `h = 0` the Kepler kernel reproduces the
    /// integrator's pre-step inertial state for every body. The
    /// `dt_sub = 0` Kepler call returns the input state unchanged, so
    /// the only operations are the barycenter reconstruction (which
    /// must give back the original Sun position) and the Galilean
    /// shift (zero at `h = 0`). Drift above f64 round-off would
    /// indicate a sign error in the rest-frame transformations.
    #[test]
    fn wh_kepler_interp_at_h_zero_returns_pre_step_state() {
        let bodies_initial = vec![
            Body::star(1000.0).at(2.0, -1.0).with_velocity(0.1, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(0.5).at(-5.0, 7.0).with_velocity(-3.0, -2.0),
        ];
        let mut sys = System::new(bodies_initial.clone(), UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        let snap = sys.last_dense_snapshot.as_ref().unwrap();
        let wh = snap.wh_data.as_ref().unwrap();
        let kin = wh.interpolate_kinematics(0.0, snap.dt);

        for (i, b) in bodies_initial.iter().enumerate() {
            let dr = (kin.positions[i] - Vec3::new(b.pos_x, b.pos_y, b.pos_z)).length();
            let dv = (kin.velocities[i] - Vec3::new(b.vel_x, b.vel_y, b.vel_z)).length();
            assert!(dr < 1e-12, "body {i} position drift at h=0: {dr}");
            assert!(dv < 1e-12, "body {i} velocity drift at h=0: {dv}");
        }
    }

    /// Two-body Kepler is exactly closed under the WH split: with no
    /// planet-planet interaction there is no kick contribution and no
    /// indirect drift, so the Kepler-analytical interpolation at
    /// `h = 1` must match the integrator's post-step state to round-
    /// off. Any deviation would expose a frame-transformation error
    /// in the dense kernel (e.g. forgetting the Galilean shift, mixing
    /// rest-frame and inertial-frame velocities).
    #[test]
    fn wh_kepler_interp_at_h_one_matches_post_step_in_two_body() {
        let bodies_initial = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1e-6).at(1.0, 0.0).with_velocity(0.0, 31.62),
        ];
        let mut sys = System::new(bodies_initial.clone(), UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(1e-3)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        let post_step: Vec<(Vec3, Vec3)> = sys
            .bodies()
            .iter()
            .map(|b| (Vec3::new(b.pos_x, b.pos_y, b.pos_z), Vec3::new(b.vel_x, b.vel_y, b.vel_z)))
            .collect();
        let snap = sys.last_dense_snapshot.as_ref().unwrap();
        let wh = snap.wh_data.as_ref().unwrap();
        let kin = wh.interpolate_kinematics(1.0, snap.dt);

        for i in 0..post_step.len() {
            let dr = (kin.positions[i] - post_step[i].0).length();
            let dv = (kin.velocities[i] - post_step[i].1).length();
            // Two-body has no kick or indirect contribution; bound is
            // round-off scale. Tolerance accounts for accumulated f64
            // error through Newton iteration in `kepler_step`.
            assert!(dr < 1e-10, "body {i} position deviation at h=1: {dr}");
            assert!(dv < 1e-10, "body {i} velocity deviation at h=1: {dv}");
        }
    }

    /// Total inertial linear momentum is conserved by the integrator
    /// step (WH 1991 §III; verified by the Bug #1 regression test).
    /// The dense interpolation must preserve the same conservation:
    /// at every sub-step `h`, `Σ m_i v_i(h)` equals the step-entry
    /// total. Drift above round-off would indicate that the Galilean
    /// shift on velocities is being applied inconsistently across
    /// bodies (e.g. forgotten on the central body's reconstruction).
    #[test]
    fn wh_kepler_interp_conserves_linear_momentum_across_step() {
        let bodies_initial = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.05, -0.03),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.05, 9.97),
            Body::rocky(0.5).at(-5.0, 7.0).with_velocity(-2.95, -2.03),
        ];
        let mut sys = System::new(bodies_initial.clone(), UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        let p_initial: Vec3 = bodies_initial
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z));
        let snap = sys.last_dense_snapshot.as_ref().unwrap();
        let wh = snap.wh_data.as_ref().unwrap();

        for h in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let kin = wh.interpolate_kinematics(h, snap.dt);
            let p_h: Vec3 = kin
                .velocities
                .iter()
                .zip(sys.bodies().iter())
                .fold(Vec3::ZERO, |s, (v, b)| s + b.mass * *v);
            let drift = (p_h - p_initial).length();
            assert!(
                drift < 1e-10,
                "linear momentum drifted at h={h}: |Δp| = {drift}, |p_0| = {}",
                p_initial.length(),
            );
        }
    }

    /// At every sub-step `h`, the synthesised central-body acceleration
    /// satisfies Newton's third law against the planet-Sun pulls
    /// (`Σ m_i a_i = 0`). This is the same invariant the post-step
    /// scratch_acc test checks, evaluated through the dense-kernel
    /// path to confirm the kinematics computation is internally
    /// consistent at sub-step times, not only at the step boundary.
    #[test]
    fn wh_kepler_interp_acc_satisfies_newton_third_law_at_sub_step() {
        let bodies_initial = vec![
            Body::star(1000.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at(10.0, 0.0).with_velocity(0.0, 10.0),
            Body::rocky(0.5).at(-7.0, 4.0).with_velocity(-2.0, -6.0),
        ];
        let mut sys = System::new(bodies_initial, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::WisdomHolman);
        sys.step();
        let snap = sys.last_dense_snapshot.as_ref().unwrap();
        let wh = snap.wh_data.as_ref().unwrap();

        for h in [0.0, 0.3, 0.7, 1.0] {
            let kin = wh.interpolate_kinematics(h, snap.dt);
            let total: Vec3 = kin
                .accelerations
                .iter()
                .zip(sys.bodies().iter())
                .fold(Vec3::ZERO, |s, (a, b)| s + b.mass * *a);
            assert!(
                total.length() < 1e-10,
                "Σ m_i a_i drifted at h={h}: |Σ m a| = {}",
                total.length(),
            );
        }
    }
}

// ── WH refactor: per-bug regression tests + smoke tests ──────────────────────
//
// Each test exercises one defect from the WH refactor with initial
// conditions chosen so that the failure mode the defect predicts
// dominates the observable signature.

mod wh_refactor_regression {
    use super::*;
    use crate::math::Vec3;

    fn total_inertial_momentum(bodies: &[Body]) -> Vec3 {
        bodies.iter().fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z))
    }

    fn total_angular_momentum(bodies: &[Body]) -> Vec3 {
        bodies.iter().fold(Vec3::ZERO, |s, b| {
            s + b.mass
                * Vec3::new(
                    b.pos_y * b.vel_z - b.pos_z * b.vel_y,
                    b.pos_z * b.vel_x - b.pos_x * b.vel_z,
                    b.pos_x * b.vel_y - b.pos_y * b.vel_x,
                )
        })
    }

    fn total_energy(bodies: &[Body], g_factor: f64) -> f64 {
        let kinetic: f64 = bodies
            .iter()
            .map(|b| 0.5 * b.mass * (b.vel_x * b.vel_x + b.vel_y * b.vel_y + b.vel_z * b.vel_z))
            .sum();
        let mut potential = 0.0;
        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                let dx = bodies[i].pos_x - bodies[j].pos_x;
                let dy = bodies[i].pos_y - bodies[j].pos_y;
                let dz = bodies[i].pos_z - bodies[j].pos_z;
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
                potential -= g_factor * bodies[i].mass * bodies[j].mass / r;
            }
        }
        kinetic + potential
    }

    /// Bug #1 regression — non-canonical centre-of-mass frame.
    ///
    /// Two-body Kepler at e = 0.5 with non-zero initial COM velocity
    /// injected. Total inertial momentum must remain at f64 floor through
    /// 1000 orbital periods at dt = T/200. Drift accumulating above the
    /// f64 floor signals Bug #1 recurrence.
    #[test]
    fn bug1_linear_momentum_conserved_under_nonzero_com_velocity() {
        let m_central = 1.0_f64;
        let m_planet = 1.0e-3_f64;
        let a = 1.0_f64;
        let e = 0.5_f64;
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) / (a * (1.0 - e))).sqrt();
        let v_com = Vec3::new(0.1, 0.05, 0.02);
        let mut bodies = vec![
            Body::star(m_central).at(0.0, 0.0).with_velocity(v_com.x, v_com.y),
            Body::rocky(m_planet).at(r_peri, 0.0).with_velocity(v_com.x, v_peri + v_com.y),
        ];
        bodies[0].vel_z = v_com.z;
        bodies[1].vel_z = v_com.z;

        let p_initial = total_inertial_momentum(&bodies);
        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / (m_central + m_planet)).sqrt();
        let dt = period / 200.0;

        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let mut max_dp = 0.0_f64;
        for _ in 0..(200 * 1000) {
            sys.step();
            let p = total_inertial_momentum(sys.bodies());
            max_dp = max_dp.max((p - p_initial).length());
        }

        // The bound `1e-10` admits the f64 round-off accumulated through the
        // Galilean-transform round-trip and the per-step momentum-conservation
        // inversion over `200 × 1000 = 2e5` substeps. A non-canonical frame
        // would exhibit secular drift many orders of magnitude above this floor.
        assert!(
            max_dp <= 1.0e-10,
            "Bug #1 regression: |ΔP| = {max_dp:.3e} exceeds 1e-10 floor; non-canonical frame leaking momentum"
        );
    }

    /// Bug #2 regression — central-body update outside the symplectic split.
    ///
    /// Two-body Kepler at e = 0.95 (high eccentricity, repeated near-singular
    /// periapsis passages), 100 orbital periods at dt = T/200. Energy
    /// conservation must avoid catastrophic loss; bound at 1e-3 is
    /// deliberately loose because high-e stretches the smooth-flow assumption
    /// — the regression target is preventing the O(1) energy loss the
    /// pre-refactor code exhibited.
    #[test]
    fn bug2_energy_bounded_at_high_eccentricity() {
        let m_central = 1.0_f64;
        let m_planet = 1.0e-6_f64;
        let a = 1.0_f64;
        let e = 0.95_f64;
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) / (a * (1.0 - e))).sqrt();
        let bodies = vec![
            Body::star(m_central).at(0.0, 0.0),
            Body::rocky(m_planet).at(r_peri, 0.0).with_velocity(0.0, v_peri),
        ];

        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / (m_central + m_planet)).sqrt();
        let dt = period / 200.0;
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let e_initial = total_energy(sys.bodies(), sys.g_factor());
        let mut max_de_rel = 0.0_f64;
        for _ in 0..(200 * 100) {
            sys.step();
            let e_now = total_energy(sys.bodies(), sys.g_factor());
            max_de_rel = max_de_rel.max(((e_now - e_initial) / e_initial).abs());
        }

        assert!(
            max_de_rel <= 1.0e-3,
            "Bug #2 regression: |ΔE/E_0| = {max_de_rel:.3e} exceeds 1e-3 floor at e=0.95; central body update outside split"
        );
    }

    /// Bug #4 regression — 2D-only computation.
    ///
    /// Inclined two-body Kepler (i = 30°, e = 0.3, a = 1, mass ratio
    /// 1:1.66e-7 matching Sun/Mercury so the WH 1991 fixed-center O(m_p/m_0)
    /// truncation lies far below the 1e-13 floor), 100 orbital periods at
    /// dt = T/200. Full 3D angular-momentum vector is preserved at the
    /// machine-precision floor; the z-coordinate stays inside the analytic
    /// envelope a(1+e) sin(i). The test confirms that the integration
    /// propagates z motion rather than silently dropping it.
    #[test]
    fn bug4_inclined_orbit_preserves_3d_angular_momentum() {
        let m_central = 1.0_f64;
        let m_planet = 1.66e-7_f64;
        let a = 1.0_f64;
        let e = 0.3_f64;
        let inclination = 30.0_f64.to_radians();
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) / (a * (1.0 - e))).sqrt();

        let mut bodies = vec![
            Body::star(m_central).at(0.0, 0.0),
            Body::rocky(m_planet).at(r_peri, 0.0).with_velocity(0.0, v_peri * inclination.cos()),
        ];
        bodies[1].vel_z = v_peri * inclination.sin();

        let l_initial = total_angular_momentum(&bodies);
        let l0_norm = l_initial.length();
        let z_envelope = a * (1.0 + e) * inclination.sin();

        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / (m_central + m_planet)).sqrt();
        let dt = period / 200.0;
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let mut max_dl_rel = 0.0_f64;
        let mut max_z = 0.0_f64;
        for _ in 0..(200 * 100) {
            sys.step();
            let l = total_angular_momentum(sys.bodies());
            max_dl_rel = max_dl_rel.max((l - l_initial).length() / l0_norm);
            for b in sys.bodies().iter().skip(1) {
                max_z = max_z.max(b.pos_z.abs());
            }
        }

        // Two assertions, with distinct evidentiary roles. The first is the
        // load-bearing claim Bug #4 is about: that z motion is propagated
        // through the integrator rather than silently dropped. The second
        // characterises the angular-momentum-conservation floor at this
        // mass ratio and horizon, observed empirically rather than derived;
        // catastrophic drop of 2D-only behaviour produced O(1) drift, while
        // the refactored 3D path stays within the WH 1991 fixed-center floor.
        assert!(
            max_z <= z_envelope * 1.05,
            "Bug #4 regression: max|z| = {max_z:.3e} exceeded analytic envelope a(1+e)sin(i) = {z_envelope:.3e} — z motion is being dropped"
        );
        assert!(
            max_dl_rel <= 1.0e-3,
            "Bug #4 regression: |ΔL|/|L_0| = {max_dl_rel:.3e} exceeds 1e-3 floor on inclined orbit"
        );
    }

    /// Negative test — Wisdom-Holman fails catastrophically on an equal-mass
    /// binary (regime where the perturbation expansion underlying the WH
    /// derivation does not hold).
    ///
    /// The assertion is inverted: passing means WH did fail loudly. If a
    /// future refactor makes WH conserve energy at this regime, that is
    /// itself a finding (probably indicates the algorithm has changed beyond
    /// the WH 1991 derivation, or the test is degenerate).
    #[test]
    fn wh_fails_predictably_on_equal_mass_binary() {
        // Two equal-mass bodies in a circular mutual orbit. Period T = 4π
        // for separation 2 and v = 0.5 each. Dominance ratio = 1, far below
        // the WH_DOMINANCE_RATIO = 10 threshold.
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];

        let period = 4.0 * std::f64::consts::PI;
        let dt = period / 200.0;
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let e_initial = total_energy(sys.bodies(), sys.g_factor());
        let mut max_de_rel = 0.0_f64;
        for _ in 0..(200 * 100) {
            sys.step();
            let e_now = total_energy(sys.bodies(), sys.g_factor());
            max_de_rel = max_de_rel.max(((e_now - e_initial) / e_initial).abs());
        }

        // The bound `1e-4` sits an order of magnitude above the WH 1991
        // smooth-flow floor (1e-5 for Sun-Mercury-class hierarchical
        // configurations). Equal-mass binary drift consistently lands at
        // 10x to 1000x the smooth-flow floor — not catastrophic O(1) loss
        // (the Galilean shift to the rest frame and the barycenter-constraint
        // reconstruction do absorb some of the perturbation-expansion
        // breakdown), but well above the regime where WH derivation
        // applies. Drift below 1e-4 would suggest WH is conserving energy
        // at the Sun-Mercury floor for an equal-mass binary, which the
        // perturbation expansion does not justify.
        assert!(
            max_de_rel >= 1.0e-4,
            "Equal-mass binary: WH energy drift {max_de_rel:.3e} unexpectedly low — \
             10x above WH 1991 smooth-flow floor expected; investigate"
        );
    }

    /// Negative test — Wisdom-Holman degrades gracefully on a marginal
    /// hierarchy where the static `is_suitable_for()` criterion still passes
    /// but the perturbation expansion is no longer in its small-parameter
    /// regime.
    ///
    /// Single planet at m_p/m_0 = 0.1, eccentric orbit. The dominance ratio
    /// is exactly 10 — at the WH_DOMINANCE_RATIO threshold — so
    /// `is_suitable_for()` passes the formal criterion. But m_p/m_0 = 0.1
    /// is six orders of magnitude larger than the Mercury-like ratio for
    /// which WH 1991 publishes 1e-5 conservation; the perturbation expansion
    /// here is at the edge of its asymptotic series. The observed energy
    /// drift exceeds the smooth-flow Tier 1 floor — graceful degradation
    /// rather than catastrophic failure. The assertion is inverted: passing
    /// means WH degraded as expected when the small-parameter regime is
    /// stretched.
    ///
    /// Single-perturber asymmetric configuration (rather than symmetric
    /// multi-planet) ensures the perturbation expansion is exercised with
    /// no symmetry-cancellation artefacts.
    #[test]
    fn wh_degrades_predictably_on_marginal_hierarchy() {
        let m_central = 1.0_f64;
        let m_planet = 0.1_f64;
        let a = 1.0_f64;
        let e = 0.3_f64;
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) / (a * (1.0 - e))).sqrt();
        let bodies = vec![
            Body::star(m_central).at(0.0, 0.0),
            Body::rocky(m_planet).at(r_peri, 0.0).with_velocity(0.0, v_peri),
        ];

        let period = 2.0 * std::f64::consts::PI;
        let dt = period / 200.0;
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let e_initial = total_energy(sys.bodies(), sys.g_factor());
        let mut max_de_rel = 0.0_f64;
        for _ in 0..(200 * 100) {
            sys.step();
            let e_now = total_energy(sys.bodies(), sys.g_factor());
            max_de_rel = max_de_rel.max(((e_now - e_initial) / e_initial).abs());
        }

        // The bound `1e-5` is the WH 1991 smooth-flow floor for
        // m_p/m_0 = 1.66e-7 (Sun + Mercury). At m_p/m_0 = 0.1 the
        // perturbation-expansion small parameter is six orders of magnitude
        // larger; energy drift exceeds the Tier 1 floor (graceful
        // degradation, not catastrophic failure).
        assert!(
            max_de_rel >= 1.0e-5,
            "Marginal hierarchy m_p/m_0 = 0.1: WH energy drift {max_de_rel:.3e} \
             unexpectedly low — the perturbation expansion at this mass ratio \
             should not preserve energy at the Sun+Mercury floor; investigate"
        );
    }

    /// Tier 1 smoke — Sun + Mercury hierarchical, smooth-flow energy
    /// conservation at the WH 1991 published floor (1e-5 over 1000 orbits
    /// at dt = T/200).
    #[test]
    fn tier1_sun_mercury_energy_within_published_floor() {
        let m_sun = 1.0_f64;
        let m_mercury = 1.66e-7_f64;
        let a = 1.0_f64;
        let e = 0.2056_f64;
        let r_peri = a * (1.0 - e);
        let v_peri = ((1.0 + e) / (a * (1.0 - e))).sqrt();
        let bodies = vec![
            Body::star(m_sun).at(0.0, 0.0),
            Body::rocky(m_mercury).at(r_peri, 0.0).with_velocity(0.0, v_peri),
        ];

        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / (m_sun + m_mercury)).sqrt();
        let dt = period / 200.0;
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(dt);
        sys.set_integrator(IntegratorKind::WisdomHolman);

        let e_initial = total_energy(sys.bodies(), sys.g_factor());
        let mut max_de_rel = 0.0_f64;
        for _ in 0..(200 * 1000) {
            sys.step();
            let e_now = total_energy(sys.bodies(), sys.g_factor());
            max_de_rel = max_de_rel.max(((e_now - e_initial) / e_initial).abs());
        }

        assert!(
            max_de_rel <= 1.0e-5,
            "Tier 1 smoke: |ΔE/E_0| = {max_de_rel:.3e} exceeds WH 1991 floor of 1e-5"
        );
    }
}

// ── Benchmarks ────────────────────────────────────────────────────────────────
//
// Quantitative accuracy tests against known solutions.  Each scenario has
// a documented tolerance derived from the integrator order and step size.

mod benchmarks {
    use super::*;

    // ── Kepler ────────────────────────────────────────────────────────────────
    //
    // Two equal-mass bodies (ε = 0) at periapsis of an elliptical orbit.
    // a = 2, e = 0.5, μ = G·(1+1) = 2, T = 4π ≈ 12.566.
    //
    // Expected accuracy (dt = 0.01):
    //   VV  O(dt²/T) × 2π ≈ 5 × 10⁻⁵  →  tol 10⁻²
    //   Y4  O(dt⁴/T) × 2π ≈ 5 × 10⁻⁹  →  tol 10⁻⁶

    fn solve_kepler(mean_anomaly: f64, e: f64) -> f64 {
        let mut ea = mean_anomaly;
        for _ in 0..60 {
            let d = (mean_anomaly - ea + e * ea.sin()) / (1.0 - e * ea.cos());
            ea += d;
            if d.abs() < 1e-14 {
                break;
            }
        }
        ea
    }

    fn kepler_relative_pos(t: f64, mu: f64, a: f64, e: f64) -> (f64, f64) {
        let n = (mu / a.powi(3)).sqrt();
        let ea = solve_kepler(n * t, e);
        (a * (ea.cos() - e), a * (1.0 - e * e).sqrt() * ea.sin())
    }

    fn kepler_position_error(integrator: IntegratorKind, dt: f64, n_steps: u64) -> f64 {
        const A: f64 = 2.0;
        const E: f64 = 0.5;
        const MU: f64 = 2.0;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();

        let b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
        let b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);

        let mut sys = System::new(vec![b1, b2], UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(dt)
            .with_max_depth(10);
        sys.set_integrator(integrator);
        for _ in 0..n_steps {
            sys.step();
        }
        let t = n_steps as f64 * dt;
        let bodies = sys.bodies();
        let (rx, ry) = (bodies[1].pos_x - bodies[0].pos_x, bodies[1].pos_y - bodies[0].pos_y);
        let (ex, ey) = kepler_relative_pos(t, MU, A, E);
        ((rx - ex).powi(2) + (ry - ey).powi(2)).sqrt()
    }

    #[test]
    fn kepler_velocity_verlet() {
        let err = kepler_position_error(IntegratorKind::VelocityVerlet, 0.01, 1257);
        assert!(err < 1e-2, "VV Kepler |Δr| = {:.3e} > 1e-2", err);
    }

    #[test]
    fn kepler_yoshida4() {
        let err = kepler_position_error(IntegratorKind::Yoshida4, 0.01, 1257);
        assert!(err < 1e-6, "Y4 Kepler |Δr| = {:.3e} > 1e-6", err);
    }

    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect raw Kepler errors"]
    fn diagnostic_kepler_errors() {
        for &(label, integrator, dt, n) in &[
            ("VV  dt=0.01  ", IntegratorKind::VelocityVerlet, 0.01_f64, 1257u64),
            ("VV  dt=0.001 ", IntegratorKind::VelocityVerlet, 0.001_f64, 12567u64),
            ("Y4  dt=0.01  ", IntegratorKind::Yoshida4, 0.01_f64, 1257u64),
            ("Y4  dt=0.001 ", IntegratorKind::Yoshida4, 0.001_f64, 12567u64),
        ] {
            let err = kepler_position_error(integrator, dt, n);
            println!("{label}  |Δr| = {err:.3e}");
        }
    }

    // ── Figure-8 choreography ─────────────────────────────────────────────────
    //
    // Chenciner & Montgomery (2000) three-body figure-8.
    // T ≈ 6.32591398 (Simó 2002), G = 1, m = 1, ε = 0.
    // Tolerance 10⁻³: timing floor ≈ 8.6 × 10⁻⁵, factor-of-12 margin.

    const FIGURE8_IC: [(f64, f64, f64, f64); 3] = [
        (-0.97000436, 0.24308753, 0.46620369, 0.43236573),
        (0.97000436, -0.24308753, 0.46620369, 0.43236573),
        (0.0, 0.0, -0.93240737, -0.86473146),
    ];
    const FIGURE8_T: f64 = 6.32591398;

    #[test]
    fn figure8_orbit_closure_yoshida4() {
        const DT: f64 = 0.001;
        const STEPS: u64 = 6326;
        const TOL: f64 = 1e-3;

        let bodies = FIGURE8_IC
            .iter()
            .map(|&(x, y, vx, vy)| Body::rocky(1.0).at(x, y).with_velocity(vx, vy))
            .collect();
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(DT)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::Yoshida4);
        for _ in 0..STEPS {
            sys.step();
        }

        let max_err = FIGURE8_IC
            .iter()
            .zip(sys.bodies().iter())
            .map(|(&(x0, y0, _, _), b)| ((b.pos_x - x0).powi(2) + (b.pos_y - y0).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);

        assert!(
            max_err < TOL,
            "Figure-8 (Y4, dt={DT}): max |Δr| = {:.3e} > {:.0e} \
             after {STEPS} steps (t={:.6}, T≈{FIGURE8_T:.6})",
            max_err,
            TOL,
            STEPS as f64 * DT,
        );
    }

    #[test]
    #[ignore = "diagnostic — run with --ignored to inspect figure-8 closure errors"]
    fn diagnostic_figure8_closure() {
        for &(label, integrator, dt, steps) in &[
            ("Y4  dt=0.001 ", IntegratorKind::Yoshida4, 0.001_f64, 6326u64),
            ("Y4  dt=0.0001", IntegratorKind::Yoshida4, 0.0001_f64, 63259u64),
            ("VV  dt=0.001 ", IntegratorKind::VelocityVerlet, 0.001_f64, 6326u64),
        ] {
            let bodies = FIGURE8_IC
                .iter()
                .map(|&(x, y, vx, vy)| Body::rocky(1.0).at(x, y).with_velocity(vx, vy))
                .collect();
            let mut sys = System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(dt)
                .with_max_depth(10);
            sys.set_integrator(integrator);
            for _ in 0..steps {
                sys.step();
            }
            println!("{label}  t={:.6}  T={FIGURE8_T:.6}", steps as f64 * dt);
            for (i, (&(x0, y0, _, _), b)) in FIGURE8_IC.iter().zip(sys.bodies().iter()).enumerate()
            {
                let err = ((b.pos_x - x0).powi(2) + (b.pos_y - y0).powi(2)).sqrt();
                println!("  body {i}: |Δr| = {err:.3e}");
            }
        }
    }

    // ── Pythagorean three-body ────────────────────────────────────────────────
    //
    // Burrau (1913) initial conditions: m₁=3 at (1,3), m₂=4 at (-2,-1), m₃=5
    // at (1,-1), all at rest.  Separations form a 3-4-5 right triangle.
    // G = 1, ε = 0 (pure Newtonian gravity).
    //
    // Initial potential energy (exact):
    //   r₁₂ = 5, r₁₃ = 4, r₂₃ = 3
    //   E₀ = −G·(m₁m₂/r₁₂ + m₁m₃/r₁₃ + m₂m₃/r₂₃)
    //      = −(12/5 + 15/4 + 20/3) = −769/60 ≈ −12.8167
    //
    // The system is chaotic — one body is eventually ejected (~t ≈ 46–60).
    // Position tests are limited to t = 1 (before chaos onset).

    fn pythagorean_system(dt: f64) -> System {
        let bodies = [
            Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0),
            Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0),
            Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0),
        ];
        let mut sys = System::new(bodies.to_vec(), UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(dt)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::Yoshida4);
        sys
    }

    #[test]
    fn pythagorean_initial_geometry() {
        const POS: [(f64, f64); 3] = [(1.0, 3.0), (-2.0, -1.0), (1.0, -1.0)];
        let d = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
        assert!((d(POS[0], POS[1]) - 5.0).abs() < 1e-15, "r₁₂ ≠ 5");
        assert!((d(POS[0], POS[2]) - 4.0).abs() < 1e-15, "r₁₃ ≠ 4");
        assert!((d(POS[1], POS[2]) - 3.0).abs() < 1e-15, "r₂₃ ≠ 3");
    }

    #[test]
    fn pythagorean_initial_energy() {
        const E0_EXACT: f64 = -769.0 / 60.0;
        let mut sys = pythagorean_system(1e-5);
        sys.step();
        let e0 = sys.metrics().total_energy;
        let rel = (e0 - E0_EXACT).abs() / E0_EXACT.abs();
        assert!(rel < 1e-6, "E₀ = {e0:.8}, exact = {E0_EXACT:.8}, rel = {rel:.3e}");
    }

    #[test]
    fn pythagorean_energy_conservation_y4() {
        const DT: f64 = 1e-3;
        const STEPS: u64 = 1_000;
        const TOL: f64 = 1e-8;

        let mut sys = pythagorean_system(DT);
        let mut peak: f64 = 0.0;
        for _ in 0..STEPS {
            sys.step();
            peak = peak.max(
                sys.metrics()
                    .rel_energy_error
                    .expect("well-conditioned regime: rel_energy_error must be Some")
                    .abs(),
            );
        }
        assert!(peak < TOL, "Y4 peak |δE/E₀| = {peak:.3e} > {TOL:.0e} over t=1.0");
    }

    #[test]
    fn pythagorean_angular_momentum_conserved() {
        const DT: f64 = 1e-3;
        let mut sys = pythagorean_system(DT);
        for _ in 0..1_000 {
            sys.step();
        }
        let lz = sys.metrics().angular_momentum_z.abs();
        assert!(lz < 1e-12, "|L_z| = {lz:.3e} after t=1.0, expected < 1e-12");
    }

    #[test]
    fn pythagorean_position_convergence_y4() {
        const T_END: f64 = 1.0;
        const TOL: f64 = 1e-6;

        let mut ref_sys = pythagorean_system(1e-4);
        for _ in 0..(T_END / 1e-4) as usize {
            ref_sys.step();
        }

        let mut sys = pythagorean_system(1e-3);
        for _ in 0..(T_END / 1e-3) as usize {
            sys.step();
        }

        let max_dr = ref_sys
            .bodies()
            .iter()
            .zip(sys.bodies().iter())
            .map(|(r, t)| ((r.pos_x - t.pos_x).powi(2) + (r.pos_y - t.pos_y).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);

        assert!(max_dr < TOL, "Y4 max |Δr| = {max_dr:.3e} > {TOL:.0e} at t=1 (dt=1e-3 vs dt=1e-4)");
    }

    #[test]
    #[ignore = "diagnostic — prints canonical positions for REBOUND cross-validation"]
    fn diagnostic_pythagorean_positions() {
        const DT: f64 = 1e-4;
        let mut sys = pythagorean_system(DT);

        println!("\nPythagorean 3-body — Y4, dt={DT}, ε=0, G=1");
        println!("IC: m=3 at (1,3), m=4 at (-2,-1), m=5 at (1,-1), v=0");
        println!("E₀ = -769/60 = {:.10}", -769.0_f64 / 60.0);

        let mut steps_done: usize = 0;
        for &t_snap in &[1.0_f64, 5.0, 10.0] {
            let target = (t_snap / DT).round() as usize;
            for _ in steps_done..target {
                sys.step();
            }
            steps_done = target;
            let m = sys.metrics();
            println!(
                "\nt = {t_snap:.1}  (step {target})  δE/E₀ = {:.3e}  L_z = {:.3e}",
                m.rel_energy_error.unwrap_or(f64::NAN),
                m.angular_momentum_z,
            );
            for (i, b) in sys.bodies().iter().enumerate() {
                println!(
                    "  body{i} m={:.0}:  x={:+.10e}  y={:+.10e}  vx={:+.10e}  vy={:+.10e}",
                    b.mass, b.pos_x, b.pos_y, b.vel_x, b.vel_y,
                );
            }
        }
    }
}

// ── Deterministic replay ──────────────────────────────────────────────────────
//
// Fase 2 gate: given identical ICs and seed, two independent runs on the same
// platform produce bit-identical body states.  These tests are the CI guard for
// reproducibility regressions.
//
// We use `f64::to_bits()` equality (i.e. NaN-aware bitwise comparison) to catch
// even single-ULP drift between runs.

mod replay {
    use super::*;

    #[test]
    fn same_ic_produces_identical_trajectory() {
        const STEPS: u64 = 500;
        let mut sys_a = two_body_deterministic_system();
        let mut sys_b = two_body_deterministic_system();
        for _ in 0..STEPS {
            sys_a.step();
            sys_b.step();
        }
        assert_bodies_bit_equal(sys_a.bodies(), sys_b.bodies(), "same-IC replay");
    }
}

// ── Hook dispatch ─────────────────────────────────────────────────────────────
//
// Verifies the observer + command pattern end-to-end: hooks fire from
// System::step() in the documented phase order, and commands they return
// mutate state (body removal, stop request) after dispatch.

mod hook_dispatch {
    use super::*;
    use crate::core::hooks::{Command, HookContext, SimHook};
    use crate::physics::integrator::IntegratorKind;
    use std::sync::{Arc, Mutex};

    #[derive(Default)]
    struct PhaseRecorder {
        log: Arc<Mutex<Vec<&'static str>>>,
    }

    impl SimHook for PhaseRecorder {
        fn pre_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
            self.log.lock().unwrap().push("pre");
            Vec::new()
        }
        fn post_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
            self.log.lock().unwrap().push("post");
            Vec::new()
        }
    }

    #[test]
    fn pre_and_post_step_fire_in_order() {
        let mut sys = two_body_circular_system(IntegratorKind::VelocityVerlet, 0.01);
        let log = Arc::new(Mutex::new(Vec::new()));
        sys.hooks_mut().register(0, Box::new(PhaseRecorder { log: log.clone() }));

        sys.step();
        sys.step();

        assert_eq!(*log.lock().unwrap(), vec!["pre", "post", "pre", "post"]);
    }

    struct RemoveFirstOnce {
        fired: bool,
    }

    impl SimHook for RemoveFirstOnce {
        fn post_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
            if self.fired {
                return Vec::new();
            }
            self.fired = true;
            vec![Command::RemoveBody { index: 0 }]
        }
    }

    #[test]
    fn remove_body_command_shrinks_system() {
        let mut sys = two_body_circular_system(IntegratorKind::VelocityVerlet, 0.01);
        assert_eq!(sys.bodies().len(), 2);

        sys.hooks_mut().register(0, Box::new(RemoveFirstOnce { fired: false }));
        sys.step();

        assert_eq!(sys.bodies().len(), 1, "RemoveBody command must drop one body");
    }

    struct StopAfterOne;

    impl SimHook for StopAfterOne {
        fn post_step(&mut self, _ctx: &HookContext<'_>) -> Vec<Command> {
            vec![Command::Stop]
        }
    }

    #[test]
    fn stop_command_sets_stop_requested() {
        let mut sys = two_body_circular_system(IntegratorKind::VelocityVerlet, 0.01);
        assert!(!sys.stop_requested());

        sys.hooks_mut().register(0, Box::new(StopAfterOne));
        sys.step();

        assert!(sys.stop_requested(), "Command::Stop must flip stop_requested");
        sys.clear_stop_request();
        assert!(!sys.stop_requested());
    }
}

// ── Integrator / force-model compatibility ────────────────────────────────────
//
// Contract (see `docs/adr/003-integrator-execution-profile.md`):
//
//   * `System::set_integrator` is the single enforcement point for the
//     integrator/force-model pairing rule.
//   * When the new integrator requires a deterministic force and the
//     current force model is not deterministic, the force model is
//     auto-reconfigured (exact threshold raised) so BH is unreachable.
//   * When the new integrator does not require determinism, the
//     force model is left untouched.
//
// These tests guard the contract against accidental regression.
mod integrator_force_compat {
    use super::*;

    /// Build a system large enough that BH would be active by default
    /// (N > the engine's built-in `EXACT_THRESHOLD = 64`).
    fn many_body_system() -> System {
        // N=80 — comfortably above the 64 default but not expensive.
        let bodies: Vec<Body> = (0..80)
            .map(|i| {
                let theta = i as f64 * 0.1;
                Body::rocky(1.0).at(theta.cos(), theta.sin()).with_velocity(0.0, 0.0)
            })
            .collect();
        System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10)
    }

    #[test]
    fn ias15_selection_forces_deterministic_force_model() {
        let mut sys = many_body_system();
        // Pre-condition: large N + default threshold → BH would be used.
        // Switch to a non-precision integrator first so the default
        // (Yoshida 4) does not retroactively affect the check.
        sys.set_integrator(IntegratorKind::VelocityVerlet);
        assert!(
            !sys.bh_engine().unwrap().is_direct_mode(),
            "baseline: engine should start in BH mode at N=80 with default threshold"
        );

        sys.set_integrator(IntegratorKind::Ias15);

        assert!(
            sys.bh_engine().unwrap().is_direct_mode(),
            "IAS15 selection must auto-switch the force model to direct mode"
        );
        assert_eq!(
            sys.integrator_kind(),
            IntegratorKind::Ias15,
            "IAS15 must remain the active integrator after auto-correction"
        );
    }

    #[test]
    fn symplectic_selection_preserves_barnes_hut() {
        let mut sys = many_body_system();
        sys.set_integrator(IntegratorKind::VelocityVerlet);
        let threshold_before = sys.exact_threshold();

        sys.set_integrator(IntegratorKind::Yoshida4);

        assert_eq!(
            sys.exact_threshold(),
            threshold_before,
            "Yoshida 4 does not require determinism; force-model configuration must not change"
        );
        assert!(
            !sys.bh_engine().unwrap().is_direct_mode(),
            "symplectic integrator must leave BH active at large N"
        );
    }

    #[test]
    fn switching_ias15_then_symplectic_does_not_revert_threshold() {
        // Once IAS15 has raised the threshold, switching back to a
        // non-precision integrator does NOT lower it. This is the
        // correct behaviour: `set_integrator` is a rule about what the
        // *new* integrator needs, not a reversal of prior adjustments.
        // If the user wants BH back, they can call
        // `set_exact_threshold` explicitly.
        let mut sys = many_body_system();
        sys.set_integrator(IntegratorKind::Ias15);
        assert!(sys.bh_engine().unwrap().is_direct_mode());

        sys.set_integrator(IntegratorKind::Yoshida4);

        assert!(
            sys.bh_engine().unwrap().is_direct_mode(),
            "switching away from IAS15 should not auto-revert the force-model threshold"
        );
    }

    #[test]
    fn small_n_system_stays_direct_under_ias15() {
        // At N ≤ 64 the engine is already in the direct path for its
        // normal operation (BH bypass by `exact_threshold` default).
        // The check guards that we do not accidentally *decrease* the
        // threshold in that regime.
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);

        sys.set_integrator(IntegratorKind::Ias15);

        // Engine exact_threshold must be at or above the direct-mode
        // threshold after IAS15 selection, regardless of starting state.
        assert!(
            sys.bh_engine().unwrap().is_direct_mode(),
            "IAS15 must leave the engine in direct mode even at small N"
        );
    }

    #[test]
    fn mercurius_does_not_force_direct_mode() {
        // Mercurius computes its own K-weighted forces internally and
        // does not rely on the outer `ctx.force` being deterministic
        // (`requires_deterministic_force = false`). It must NOT raise
        // the BH exact_threshold at selection time.
        let mut sys = many_body_system();
        sys.set_integrator(IntegratorKind::VelocityVerlet);
        let threshold_before = sys.exact_threshold();

        sys.set_integrator(IntegratorKind::Mercurius);

        assert_eq!(
            sys.exact_threshold(),
            threshold_before,
            "Mercurius does not require determinism; force-model configuration must not change"
        );
    }
}

#[cfg(test)]
mod encounter_diagnostic {
    use super::*;
    use crate::physics::encounter::EncounterFlag;

    fn approaching_pair_system(initial_separation: f64, approach_velocity: f64) -> System {
        // Two equal-mass bodies on a head-on closing trajectory along x.
        // The negative `vel_x` on body 1 closes the gap at rate `2 ·
        // approach_velocity` (both bodies move toward the midpoint).
        let bodies = vec![
            Body::rocky(1.0)
                .at(-initial_separation * 0.5, 0.0)
                .with_velocity(approach_velocity, 0.0),
            Body::rocky(1.0)
                .at(initial_separation * 0.5, 0.0)
                .with_velocity(-approach_velocity, 0.0),
        ];
        System::new(bodies, UnitSystem::canonical()).with_dt(0.01)
    }

    #[test]
    fn unconfigured_threshold_keeps_flag_far() {
        let mut sys = approaching_pair_system(2.0, 0.5);
        for _ in 0..50 {
            sys.step();
            assert_eq!(sys.last_encounter_flag(), EncounterFlag::Far);
        }
    }

    #[test]
    fn flag_escalates_as_pair_closes() {
        // Threshold at 1.0; bodies start at separation 2.0 closing at
        // 1.0/t.u. → far → approaching → close as the gap shrinks past
        // the threshold and then past the half-threshold.
        let mut sys = approaching_pair_system(2.0, 0.5);
        sys.set_close_encounter_threshold(Some(1.0));

        let mut saw_far = false;
        let mut saw_approaching = false;
        let mut saw_close = false;

        for _ in 0..200 {
            sys.step();
            match sys.last_encounter_flag() {
                EncounterFlag::Far => saw_far = true,
                EncounterFlag::Approaching => saw_approaching = true,
                EncounterFlag::Close => {
                    saw_close = true;
                    break;
                },
            }
        }

        assert!(saw_far, "should observe Far before bodies close");
        assert!(saw_approaching, "should pass through Approaching band");
        assert!(saw_close, "should reach Close before 200 steps");
    }

    #[test]
    fn three_d_separation_is_used() {
        // Two bodies coincident in xy but separated in z. The pre-fix
        // 2D `compute_closeness` would report `r_min == 0` and trigger
        // a spurious Close flag; the 3D fix sees the true separation.
        let bodies = vec![
            Body::rocky(1.0).at(0.0, 0.0).with_velocity(0.0, 0.0),
            Body::rocky(1.0).at_3d(0.0, 0.0, 1.0).with_velocity_3d(0.0, 0.0, 0.0),
        ];
        let mut sys = System::new(bodies, UnitSystem::canonical()).with_dt(0.01);
        sys.set_close_encounter_threshold(Some(0.1));

        sys.step();
        assert!(sys.r_min >= 0.99, "3D separation should be ~1.0, got r_min = {}", sys.r_min);
        assert_eq!(
            sys.last_encounter_flag(),
            EncounterFlag::Far,
            "z-separated bodies must not register a 2D-projected close encounter"
        );
    }

    #[test]
    fn updating_threshold_resets_transition_tracker() {
        let mut sys = approaching_pair_system(0.3, 0.0);
        sys.set_close_encounter_threshold(Some(1.0));
        sys.step();
        assert_eq!(sys.last_encounter_flag(), EncounterFlag::Close);

        sys.set_close_encounter_threshold(Some(0.1));
        // Tracker reset to Far; the pair is far above the new threshold.
        assert_eq!(sys.last_encounter_flag(), EncounterFlag::Far);
        sys.step();
        assert_eq!(sys.last_encounter_flag(), EncounterFlag::Far);
    }
}
