//! REBOUND parity — Mercurius, apsis side.
//!
//! Sun + 4 outer planets + 1 Jupiter-crossing test particle integrated
//! under apsis Mercurius for 10⁴ years (~840 Jupiter orbits) at outer
//! `dt = 0.01 yr`. Output is sampled at 1-year cadence (10001 samples
//! per body) for cross-implementation comparison against REBOUND
//! MERCURIUS via the matching Python harness.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example rebound_parity_mercurius -p apsis
//! cargo run --release --example rebound_parity_mercurius -p apsis -- --output path/to/apsis.csv
//! ```
//!
//! Default output path: `validation/rebound-parity/mercurius/out/apsis.csv`
//! (relative to the workspace root — `cargo run` should be invoked from
//! the repo root or from the orchestrator in that directory).
//!
//! ## Protocol
//!
//! Initial conditions, integrator settings, metrics, and tolerances declared
//! *a priori* in
//! [`docs/experiments/2026-05-13-rebound-parity-mercurius.md`](../../../../docs/experiments/2026-05-13-rebound-parity-mercurius.md).
//! Constants in this file (`A_*`, `M_*`, `N_YEARS`, `DT_YEARS`,
//! `ALPHA_HILL`) are the protocol's IC values — changes here are protocol
//! changes; update the notebook in lockstep.

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

// ── Protocol constants (mirrored in the notebook) ───────────────────────── //

/// Solar mass (canonical 1 in solar AU-year units).
const M_SUN: f64 = 1.0;

/// Jupiter / Sun mass ratio (Murray & Dermott Table A.2).
const M_JUPITER: f64 = 9.55e-4;
const M_SATURN: f64 = 2.86e-4;
const M_URANUS: f64 = 4.37e-5;
const M_NEPTUNE: f64 = 5.15e-5;
/// Test-particle mass: massless-class probe (1 nano-solar mass).
const M_TEST: f64 = 1.0e-9;

/// Heliocentric semi-major axes (AU).
const A_JUPITER: f64 = 5.20;
const A_SATURN: f64 = 9.58;
const A_URANUS: f64 = 19.18;
const A_NEPTUNE: f64 = 30.07;

/// Test particle: eccentric Jupiter-crossing orbit with non-zero inclination.
const A_TEST: f64 = 4.20;
const E_TEST: f64 = 0.40;
const I_TEST: f64 = 0.05;

/// Total integration horizon (years).
const N_YEARS: u64 = 10_000;
/// Outer integrator step size (years).
const DT_YEARS: f64 = 0.01;
/// Hill-radius multiplier for Mercurius changeover (REBOUND default).
const ALPHA_HILL: f64 = 3.0;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Build the body list in the rest frame ───────────────────────────── //
    //
    // All planets on circular coplanar orbits at named heliocentric
    // distances, true anomalies spread out so the system is not
    // axis-aligned. The test particle starts at periapsis with non-zero
    // inclination — apoapsis sits above Jupiter's orbit, so the orbit
    // crosses Jupiter's twice per test-particle period.
    //
    // The Sun starts at the origin with zero velocity; the per-body
    // velocities below are computed in the heliocentric frame (Sun-rest)
    // and the COM shift below moves the whole system into the COM frame
    // before integration starts.

    let units = UnitSystem::solar();
    let g = units.g();

    let mut bodies = vec![
        Body::star(M_SUN).at_3d(0.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 0.0).unsoftened(),
        circular_planet(M_JUPITER, A_JUPITER, 0.0, g),
        circular_planet(M_SATURN, A_SATURN, std::f64::consts::FRAC_PI_2, g),
        circular_planet(M_URANUS, A_URANUS, std::f64::consts::PI, g),
        circular_planet(M_NEPTUNE, A_NEPTUNE, 1.5 * std::f64::consts::PI, g),
        eccentric_inclined_planet(M_TEST, A_TEST, E_TEST, I_TEST, g),
    ];

    // COM-shift to the rest frame so the integrator starts on the
    // canonical IC (REBOUND will do `sim.move_to_com()` to land at the
    // same configuration).
    com_shift(&mut bodies);

    // ── Integrator setup ────────────────────────────────────────────────── //
    let mut sys =
        System::new(bodies, units).with_integrator(IntegratorKind::Mercurius).with_dt(DT_YEARS);
    sys.set_mercurius_alpha(ALPHA_HILL);

    // ── CSV output ──────────────────────────────────────────────────────── //
    //
    // Wide format: one row per yearly sample, including all body state
    // and conservation diagnostics. Initial state at year=0; subsequent
    // samples at yearly intervals. Total = N_YEARS + 1 = 10001 samples.
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# REBOUND parity — Mercurius — apsis side").unwrap();
    writeln!(w, "# protocol: docs/experiments/2026-05-13-rebound-parity-mercurius.md").unwrap();
    writeln!(w, "# integrator: Mercurius (apsis), alpha={ALPHA_HILL}").unwrap();
    writeln!(w, "# units: solar AU-year (G ≈ {g:.18e})").unwrap();
    writeln!(w, "# n_years={N_YEARS}, dt_years={DT_YEARS}, n_bodies={}", sys.bodies().len())
        .unwrap();
    write_header(&mut w, sys.bodies().len());

    write_sample(&mut w, 0, &sys);
    for year in 1..=N_YEARS {
        let t_target = year as f64;
        sys.integrate_until(t_target);
        write_sample(&mut w, year, &sys);
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {}", N_YEARS + 1, output_path.display());
}

// ── IC helpers ──────────────────────────────────────────────────────────── //

/// Circular planet at heliocentric distance `a`, true anomaly `nu`, in
/// the (x, y) plane around a Sun at the origin. Heliocentric circular
/// speed `v_c = sqrt(G·M_sun / a)`.
fn circular_planet(mass: f64, a: f64, nu: f64, g: f64) -> Body {
    let v_c = (g * M_SUN / a).sqrt();
    let x = a * nu.cos();
    let y = a * nu.sin();
    let vx = -v_c * nu.sin();
    let vy = v_c * nu.cos();
    Body::rocky(mass).at_3d(x, y, 0.0).with_velocity_3d(vx, vy, 0.0).unsoftened()
}

/// Eccentric inclined planet starting at periapsis. Position in the
/// orbit plane is `(r_peri, 0)`; tangent velocity is `v_peri`. The orbit
/// plane is then rotated by inclination `i` around the x-axis: position
/// stays `(r_peri, 0, 0)`; velocity becomes `(0, v_peri·cos(i),
/// v_peri·sin(i))`.
fn eccentric_inclined_planet(mass: f64, a: f64, e: f64, i: f64, g: f64) -> Body {
    let r_peri = a * (1.0 - e);
    // Standard vis-viva at periapsis around a 1 M_sun primary:
    //   v_peri² = G M (2/r_peri - 1/a) = (G M / a) · (1 + e) / (1 - e)
    let v_peri = (g * M_SUN / a * (1.0 + e) / (1.0 - e)).sqrt();
    Body::rocky(mass)
        .at_3d(r_peri, 0.0, 0.0)
        .with_velocity_3d(0.0, v_peri * i.cos(), v_peri * i.sin())
        .unsoftened()
}

/// Subtract the centre-of-mass position and velocity from every body so
/// the system starts at rest in the COM frame.
fn com_shift(bodies: &mut [Body]) {
    let mut com_pos = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut com_vel = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut mtot = 0.0_f64;
    for b in bodies.iter() {
        com_pos.0 += b.mass * b.pos_x;
        com_pos.1 += b.mass * b.pos_y;
        com_pos.2 += b.mass * b.pos_z;
        com_vel.0 += b.mass * b.vel_x;
        com_vel.1 += b.mass * b.vel_y;
        com_vel.2 += b.mass * b.vel_z;
        mtot += b.mass;
    }
    com_pos.0 /= mtot;
    com_pos.1 /= mtot;
    com_pos.2 /= mtot;
    com_vel.0 /= mtot;
    com_vel.1 /= mtot;
    com_vel.2 /= mtot;
    for b in bodies.iter_mut() {
        b.pos_x -= com_pos.0;
        b.pos_y -= com_pos.1;
        b.pos_z -= com_pos.2;
        b.vel_x -= com_vel.0;
        b.vel_y -= com_vel.1;
        b.vel_z -= com_vel.2;
    }
}

// ── Output helpers ──────────────────────────────────────────────────────── //

fn write_header(w: &mut BufWriter<File>, n_bodies: usize) {
    let mut header = String::from("year,t");
    for i in 0..n_bodies {
        header.push_str(&format!(",x{i},y{i},z{i},vx{i},vy{i},vz{i}"));
    }
    header.push_str(",e_total,lz_total");
    writeln!(w, "{header}").unwrap();
}

fn write_sample(w: &mut BufWriter<File>, year: u64, sys: &System) {
    let bodies = sys.bodies();
    let g = sys.units().g();
    let mut row = format!("{year},{:.18e}", sys.t());
    for b in bodies {
        row.push_str(&format!(
            ",{:.18e},{:.18e},{:.18e},{:.18e},{:.18e},{:.18e}",
            b.pos_x, b.pos_y, b.pos_z, b.vel_x, b.vel_y, b.vel_z
        ));
    }
    let e_total = total_energy(bodies, g);
    let lz_total = total_lz(bodies);
    row.push_str(&format!(",{e_total:.18e},{lz_total:.18e}"));
    writeln!(w, "{row}").unwrap();
}

/// Total mechanical energy `KE + PE`, with `KE = ½ Σ m_i v_i²` and
/// `PE = −Σ_{i<j} G m_i m_j / r_ij` (no softening, all bodies are
/// `Body::unsoftened()`). 3D inner products. Matches REBOUND's
/// `sim.energy()` for the Mercurius/IAS15 N-body Hamiltonian.
fn total_energy(bodies: &[Body], g: f64) -> f64 {
    let ke: f64 = bodies
        .iter()
        .map(|b| 0.5 * b.mass * (b.vel_x * b.vel_x + b.vel_y * b.vel_y + b.vel_z * b.vel_z))
        .sum();
    let mut pe = 0.0_f64;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let dz = bodies[i].pos_z - bodies[j].pos_z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            pe -= g * bodies[i].mass * bodies[j].mass / r;
        }
    }
    ke + pe
}

/// Total z-component of angular momentum: `Lz = Σ m_i (x_i v_y - y_i v_x)`.
/// Used as the conservation diagnostic that should agree with REBOUND's
/// `angular_momentum()` z-component to machine precision via the
/// analytical Kepler drift.
fn total_lz(bodies: &[Body]) -> f64 {
    bodies.iter().map(|b| b.mass * (b.pos_x * b.vel_y - b.pos_y * b.vel_x)).sum()
}

// ── CLI ─────────────────────────────────────────────────────────────────── //

fn parse_output_path() -> PathBuf {
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--output" || arg == "-o" {
            return PathBuf::from(args.next().expect("--output requires a path argument"));
        }
    }
    PathBuf::from("validation/rebound-parity/mercurius/out/apsis.csv")
}
