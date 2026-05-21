//! Cross-platform parity — WHFast outer Solar System.
//!
//! Mirrors the [`rebound_parity_mercurius`] scenario (Sun + 4 outer
//! planets + Jupiter-crossing test particle, 10⁴ years at `dt = 0.01 yr`,
//! yearly sampling) but integrates under WHFast instead of Mercurius.
//!
//! ## Purpose
//!
//! Exercise the [`kepler`] universal-variable solver — the only WHFast
//! sub-step path that calls libc-bound transcendentals (`sin`, `cos`,
//! `cosh`, `sinh`, `tanh` via Stumpff functions). PR #165 routed those
//! calls through the `libm` crate; this scenario asserts the bit-equal
//! cross-platform property holds for WHFast end-to-end, not just the
//! IAS15 + 1PN configuration validated by PR #162.
//!
//! ## Physics correctness is NOT asserted here
//!
//! The test particle starts on an eccentric, inclined, Jupiter-crossing
//! orbit. WHFast assumes hierarchical Keplerian motion and does not
//! switch integration scheme during close encounters — the particle's
//! trajectory through Jupiter's neighbourhood is therefore numerically
//! incorrect. This is **deliberate**: the scenario stress-tests the
//! Kepler solver under a regime that hits its boundary, exercising
//! libm calls at non-trivial Stumpff arguments. Cross-platform bit
//! equality is the contract; physics fidelity is the Mercurius parity
//! scenario's job (where the IAS15 sub-step handles the encounter).
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example whfast_outer_solar -p apsis -- --output path/to/whfast.csv
//! ```
//!
//! [`kepler`]: apsis::physics::integrator::kepler

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

// ── Protocol constants (mirror rebound_parity_mercurius) ────────────────── //

const M_SUN: f64 = 1.0;
const M_JUPITER: f64 = 9.55e-4;
const M_SATURN: f64 = 2.86e-4;
const M_URANUS: f64 = 4.37e-5;
const M_NEPTUNE: f64 = 5.15e-5;
const M_TEST: f64 = 1.0e-9;

const A_JUPITER: f64 = 5.20;
const A_SATURN: f64 = 9.58;
const A_URANUS: f64 = 19.18;
const A_NEPTUNE: f64 = 30.07;

const A_TEST: f64 = 4.20;
const E_TEST: f64 = 0.40;
const I_TEST: f64 = 0.05;

const N_YEARS: u64 = 10_000;
const DT_YEARS: f64 = 0.01;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    let units = UnitSystem::solar();
    let g = units.g();

    let mut bodies = vec![
        Body::star(M_SUN).at_3d(0.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 0.0),
        circular_planet(M_JUPITER, A_JUPITER, 0.0, g),
        circular_planet(M_SATURN, A_SATURN, std::f64::consts::FRAC_PI_2, g),
        circular_planet(M_URANUS, A_URANUS, std::f64::consts::PI, g),
        circular_planet(M_NEPTUNE, A_NEPTUNE, 1.5 * std::f64::consts::PI, g),
        eccentric_inclined_planet(M_TEST, A_TEST, E_TEST, I_TEST, g),
    ];

    com_shift(&mut bodies);

    let mut sys =
        System::new(bodies, units).with_integrator(IntegratorKind::WHFast).with_dt(DT_YEARS);

    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# Cross-platform parity — WHFast outer Solar System — apsis side").unwrap();
    writeln!(w, "# integrator: WHFast (apsis)").unwrap();
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

// ── IC helpers (parallel to rebound_parity_mercurius) ───────────────────── //

fn circular_planet(mass: f64, a: f64, nu: f64, g: f64) -> Body {
    let v_c = (g * M_SUN / a).sqrt();
    let x = a * nu.cos();
    let y = a * nu.sin();
    let vx = -v_c * nu.sin();
    let vy = v_c * nu.cos();
    Body::rocky(mass).at_3d(x, y, 0.0).with_velocity_3d(vx, vy, 0.0)
}

fn eccentric_inclined_planet(mass: f64, a: f64, e: f64, i: f64, g: f64) -> Body {
    let r_peri = a * (1.0 - e);
    let v_peri = (g * M_SUN / a * (1.0 + e) / (1.0 - e)).sqrt();
    Body::rocky(mass).at_3d(r_peri, 0.0, 0.0).with_velocity_3d(
        0.0,
        v_peri * i.cos(),
        v_peri * i.sin(),
    )
}

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

// ── Output helpers (parallel) ───────────────────────────────────────────── //

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
    PathBuf::from("validation/cross-platform/windows/whfast_outer_solar.csv")
}
