//! End-to-end integration tests for the simulation pipeline.
//!
//! These tests verify that the full pipeline — force evaluation, integrator,
//! energy tracking — correctly conserves the Hamiltonian over many orbital
//! periods.  They test the *integrated system*, not individual primitives.
//!
//! Physical scenario: two equal-mass bodies in a circular orbit.
//!
//!   G = 1, M₁ = M₂ = 1
//!   Positions: (−1, 0) and (+1, 0), separation d = 2, orbital radius r = 1
//!   Velocities: (0, −0.5) and (0, +0.5) — CCW orbit
//!   Orbital period: T = 2πr/v = 4π ≈ 12.566
//!
//! Tolerance derivation (dt = 0.01, T = 4π, dt/T ≈ 7.96 × 10⁻⁴):
//!   VV  (2nd order): amplitude ~ (dt/T)² ≈ 6.3 × 10⁻⁷ → tol 1e-4
//!   Y4  (4th order): amplitude ~ (dt/T)⁴ ≈ 4 × 10⁻¹³ → tol 1e-7

use super::System;
use crate::domain::body::Body;
use crate::domain::materials::Material;
use crate::physics::integrator::IntegratorKind;

// ── Shared helpers ────────────────────────────────────────────────────────────

fn two_body_circular_system(integrator: IntegratorKind, dt: f64) -> System {
    let bodies = vec![
        Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
        Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
    ];
    let mut sys = System::new(bodies, 0.5, dt, 10, 1);
    sys.set_integrator(integrator);
    sys
}

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

// ── Energy conservation ───────────────────────────────────────────────────────

#[test]
fn energy_conservation_velocity_verlet() {
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
fn energy_conservation_yoshida4() {
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
fn print_peak_errors_diagnostic() {
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

// ── Wisdom-Holman guard ───────────────────────────────────────────────────────

#[test]
fn hierarchical_system_is_suitable() {
    let bodies = vec![
        Body::new(0.0, 0.0, 0.0, 0.0, 1000.0, Material::Star),
        Body::new(10.0, 0.0, 0.0, 10.0, 1.0, Material::Rocky),
    ];
    let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
    sys.set_integrator(IntegratorKind::WisdomHolman);
    assert!(sys.is_wh_suitable());
}

#[test]
fn equal_mass_system_is_not_suitable() {
    let bodies = vec![
        Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
        Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
    ];
    let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
    sys.set_integrator(IntegratorKind::WisdomHolman);
    assert!(!sys.is_wh_suitable());
}

#[test]
fn three_equal_mass_is_not_suitable() {
    let bodies = vec![
        Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
        Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
        Body::new(0.0, 1.0, 0.5, 0.0, 1.0, Material::Rocky),
    ];
    let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
    sys.set_integrator(IntegratorKind::WisdomHolman);
    assert!(!sys.is_wh_suitable());
}

#[test]
fn boundary_at_exactly_10x_is_suitable() {
    let bodies = vec![
        Body::new(0.0, 0.0, 0.0, 0.0, 10.0, Material::Star),
        Body::new(10.0, 0.0, 0.0, 1.0, 1.0, Material::Rocky),
    ];
    assert!(System::new(bodies, 0.5, 0.01, 10, 1).is_wh_suitable());
}

#[test]
fn boundary_below_10x_is_not_suitable() {
    let bodies = vec![
        Body::new(0.0, 0.0, 0.0, 0.0, 9.9, Material::Star),
        Body::new(10.0, 0.0, 0.0, 1.0, 1.0, Material::Rocky),
    ];
    assert!(!System::new(bodies, 0.5, 0.01, 10, 1).is_wh_suitable());
}

#[test]
fn single_body_is_not_suitable() {
    let bodies = vec![Body::new(0.0, 0.0, 0.0, 0.0, 1.0, Material::Rocky)];
    assert!(!System::new(bodies, 0.5, 0.01, 10, 1).is_wh_suitable());
}

#[test]
fn wh_on_non_hierarchical_does_not_panic_and_stays_finite() {
    let bodies = vec![
        Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
        Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
    ];
    let mut sys = System::new(bodies, 0.5, 0.01, 10, 1);
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
fn wh_fallback_energy_matches_yoshida4_directly() {
    let bodies = vec![
        Body::new(-1.0, 0.0, 0.0, -0.5, 1.0, Material::Rocky),
        Body::new(1.0, 0.0, 0.0, 0.5, 1.0, Material::Rocky),
    ];
    let mut sys_wh = System::new(bodies.clone(), 0.5, 0.01, 10, 1);
    sys_wh.set_integrator(IntegratorKind::WisdomHolman);
    let mut sys_y4 = System::new(bodies, 0.5, 0.01, 10, 1);
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

// ── Kepler benchmark ──────────────────────────────────────────────────────────
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

    let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
    b1.softening = 0.0;
    let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
    b2.softening = 0.0;

    let mut sys = System::new(vec![b1, b2], 0.5, dt, 10, 1);
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
fn kepler_position_accuracy_velocity_verlet() {
    let err = kepler_position_error(IntegratorKind::VelocityVerlet, 0.01, 1257);
    assert!(err < 1e-2, "VV Kepler |Δr| = {:.3e} > 1e-2", err);
}

#[test]
fn kepler_position_accuracy_yoshida4() {
    let err = kepler_position_error(IntegratorKind::Yoshida4, 0.01, 1257);
    assert!(err < 1e-6, "Y4 Kepler |Δr| = {:.3e} > 1e-6", err);
}

#[test]
#[ignore = "diagnostic — run with --ignored to inspect raw Kepler errors"]
fn print_kepler_errors_diagnostic() {
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

// ── Figure-8 choreography ─────────────────────────────────────────────────────
//
// Chenciner & Montgomery (2000) three-body figure-8.
// T ≈ 6.32591398 (Simó 2002), G = 1, m = 1, ε = 0.
// Tolerance 10⁻³: timing floor ≈ 8.6 × 10⁻⁵, factor-of-12 margin.

const IC: [(f64, f64, f64, f64); 3] = [
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

    let bodies = IC
        .iter()
        .map(|&(x, y, vx, vy)| {
            let mut b = Body::new(x, y, vx, vy, 1.0, Material::Rocky);
            b.softening = 0.0;
            b
        })
        .collect();
    let mut sys = System::new(bodies, 0.5, DT, 10, 1);
    sys.set_integrator(IntegratorKind::Yoshida4);

    for _ in 0..STEPS {
        sys.step();
    }

    let max_err = IC
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
fn print_figure8_closure_diagnostic() {
    for &(label, integrator, dt, steps) in &[
        ("Y4  dt=0.001 ", IntegratorKind::Yoshida4, 0.001_f64, 6326u64),
        ("Y4  dt=0.0001", IntegratorKind::Yoshida4, 0.0001_f64, 63259u64),
        ("VV  dt=0.001 ", IntegratorKind::VelocityVerlet, 0.001_f64, 6326u64),
    ] {
        let bodies = IC
            .iter()
            .map(|&(x, y, vx, vy)| {
                let mut b = Body::new(x, y, vx, vy, 1.0, Material::Rocky);
                b.softening = 0.0;
                b
            })
            .collect();
        let mut sys = System::new(bodies, 0.5, dt, 10, 1);
        sys.set_integrator(integrator);
        for _ in 0..steps {
            sys.step();
        }
        println!("{label}  t={:.6}  T={FIGURE8_T:.6}", steps as f64 * dt);
        for (i, (&(x0, y0, _, _), b)) in IC.iter().zip(sys.bodies().iter()).enumerate() {
            let err = ((b.x - x0).powi(2) + (b.y - y0).powi(2)).sqrt();
            println!("  body {i}: |Δr| = {err:.3e}");
        }
    }
}
