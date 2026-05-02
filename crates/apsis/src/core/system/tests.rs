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
    let bodies = vec![
        Body::star(1.0).unsoftened(),
        Body::rocky(3e-6).at(1.0, 0.0).with_velocity(0.0, 1.0).unsoftened(),
    ];
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
        assert_eq!(ba.x.to_bits(), bb.x.to_bits(), "{context}: body {i} x differs");
        assert_eq!(ba.y.to_bits(), bb.y.to_bits(), "{context}: body {i} y differs");
        assert_eq!(ba.vx.to_bits(), bb.vx.to_bits(), "{context}: body {i} vx differs");
        assert_eq!(ba.vy.to_bits(), bb.vy.to_bits(), "{context}: body {i} vy differs");
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
            peak = peak.max(sys.metrics().rel_energy_error.abs());
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
            assert!(b.x.is_finite() && b.y.is_finite(), "body left finite domain");
            assert!(b.vx.is_finite() && b.vy.is_finite(), "velocity left finite domain");
        }
    }

    #[test]
    fn fallback_energy_matches_yoshida4_directly() {
        let bodies = vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ];
        let mut sys_wh = System::new(bodies.clone(), UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys_wh.set_integrator(IntegratorKind::WisdomHolman);
        let mut sys_y4 = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(0.01)
            .with_max_depth(10);
        sys_y4.set_integrator(IntegratorKind::Yoshida4);

        for _ in 0..100 {
            sys_wh.step();
            sys_y4.step();
        }

        let err_wh = sys_wh.metrics().rel_energy_error.abs();
        let err_y4 = sys_y4.metrics().rel_energy_error.abs();
        assert!(
            (err_wh - err_y4).abs() < 1e-15,
            "WH fallback energy error {err_wh:.3e} ≠ direct Y4 {err_y4:.3e}"
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

        let mut b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
        b1.softening = 0.0;
        let mut b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);
        b2.softening = 0.0;

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
        let (rx, ry) = (bodies[1].x - bodies[0].x, bodies[1].y - bodies[0].y);
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
            .map(|&(x, y, vx, vy)| {
                let mut b = Body::rocky(1.0).at(x, y).with_velocity(vx, vy);
                b.softening = 0.0;
                b
            })
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
            .map(|(&(x0, y0, _, _), b)| ((b.x - x0).powi(2) + (b.y - y0).powi(2)).sqrt())
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
                .map(|&(x, y, vx, vy)| {
                    let mut b = Body::rocky(1.0).at(x, y).with_velocity(vx, vy);
                    b.softening = 0.0;
                    b
                })
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
                let err = ((b.x - x0).powi(2) + (b.y - y0).powi(2)).sqrt();
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
        let mut bodies = [
            Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0),
            Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0),
            Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0),
        ];
        for b in &mut bodies {
            b.softening = 0.0;
        }
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
            peak = peak.max(sys.metrics().rel_energy_error.abs());
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
            .map(|(r, t)| ((r.x - t.x).powi(2) + (r.y - t.y).powi(2)).sqrt())
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
                m.rel_energy_error, m.angular_momentum_z,
            );
            for (i, b) in sys.bodies().iter().enumerate() {
                println!(
                    "  body{i} m={:.0}:  x={:+.10e}  y={:+.10e}  vx={:+.10e}  vy={:+.10e}",
                    b.mass, b.x, b.y, b.vx, b.vy,
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

    #[test]
    fn snapshot_midpoint_produces_identical_trajectory() {
        const STEPS_BEFORE: u64 = 200;
        const STEPS_AFTER: u64 = 300;

        let mut reference = two_body_deterministic_system();
        for _ in 0..(STEPS_BEFORE + STEPS_AFTER) {
            reference.step();
        }

        let mut replayed = two_body_deterministic_system();
        for _ in 0..STEPS_BEFORE {
            replayed.step();
        }
        let snap = replayed.to_snapshot();
        replayed.restore_from_snapshot(&snap);
        for _ in 0..STEPS_AFTER {
            replayed.step();
        }

        assert_bodies_bit_equal(reference.bodies(), replayed.bodies(), "snapshot replay");
    }

    #[test]
    fn snapshot_file_roundtrip() {
        use crate::io::snapshot::SimSnapshot;

        let mut sys = two_body_deterministic_system();
        for _ in 0..100 {
            sys.step();
        }

        let mut snap = sys.to_snapshot();
        snap.sim_name = "roundtrip-test".to_owned();

        let dir = std::env::temp_dir();
        let path = snap.save_to_dir(&dir).expect("snapshot write failed");

        let loaded = SimSnapshot::load_from(&path).expect("snapshot load failed");
        std::fs::remove_file(&path).ok();

        assert_eq!(loaded.t.to_bits(), snap.t.to_bits(), "t");
        assert_eq!(loaded.steps, snap.steps, "steps");
        assert_eq!(loaded.dt.to_bits(), snap.dt.to_bits(), "dt");
        assert_eq!(loaded.seed, snap.seed, "seed");
        assert_eq!(loaded.sim_name, snap.sim_name, "sim_name");
        assert_eq!(loaded.bodies.len(), snap.bodies.len(), "body count");
        for (i, (l, s)) in loaded.bodies.iter().zip(snap.bodies.iter()).enumerate() {
            assert_eq!(l.x.to_bits(), s.x.to_bits(), "body {i} x roundtrip");
            assert_eq!(l.y.to_bits(), s.y.to_bits(), "body {i} y roundtrip");
            assert_eq!(l.vx.to_bits(), s.vx.to_bits(), "body {i} vx roundtrip");
            assert_eq!(l.vy.to_bits(), s.vy.to_bits(), "body {i} vy roundtrip");
            assert_eq!(l.mass.to_bits(), s.mass.to_bits(), "body {i} mass roundtrip");
        }
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
}
