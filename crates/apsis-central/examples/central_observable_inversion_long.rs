//! Cross-platform parity — apsis-central observable-inversion long horizon.
//!
//! Sun + Mercury-like body, with [`CentralForce::from_apsidal_rate`]
//! registered as a Hamiltonian perturbation at `γ = -3` (the
//! Schwarzschild-effective regime, per Nobili & Roxburgh 1986)
//! targeting an apsidal-rate Mercury would feel under a strong central
//! correction. Integrated under IAS15 for 500 Mercury orbits with
//! per-orbit sampling.
//!
//! ## Purpose
//!
//! Exercise the apsis-central operator's libc-bound transcendentals
//! (`f64::powf` at four sites in force/potential/inversion, `f64::ln`
//! in the `γ = -1` branch — though `γ = -3` here exercises only the
//! `powf` path). PR #167 routed those calls through the `libm` crate;
//! this scenario asserts the bit-equal cross-platform property holds
//! for `CentralForce` end-to-end, not just the IAS15 + 1PN
//! configuration validated by PR #162.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example central_observable_inversion_long -p apsis-central -- --output path/to/central.csv
//! ```

use std::env;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_central::CentralForce;

// ── Protocol constants ──────────────────────────────────────────────────── //

const M_SUN: f64 = 1.0;
const M_MERCURY: f64 = 1.660_114e-7;
const A_MERCURY: f64 = 0.387_098; // AU
const E_MERCURY: f64 = 0.205_63;

/// Target apsidal precession rate (rad / yr) the observable-inversion constructor
/// inverts into the operator's coupling `A`. Order of magnitude matches
/// the GR 1PN Mercury rate (~43 arcsec/century ≈ 2.1 × 10⁻⁶ rad/yr) but
/// the value is otherwise nominal — the scenario asserts cross-platform
/// determinism, not GR fidelity.
const OMEGA_DOT_TARGET: f64 = 2.0e-6;

/// Schwarzschild-effective regime: `a = A · r⁻³`.
const GAMMA: f64 = -3.0;

const N_ORBITS: u64 = 500;
const DT_YEARS: f64 = 1.0e-4;

// ── Main ────────────────────────────────────────────────────────────────── //

fn main() {
    let output_path = parse_output_path();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    let units = UnitSystem::solar();
    let g = units.g();

    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = (g * M_SUN / A_MERCURY * (1.0 + E_MERCURY) / (1.0 - E_MERCURY)).sqrt();

    let bodies = vec![
        Body::star(M_SUN).at_3d(0.0, 0.0, 0.0).with_velocity_3d(0.0, 0.0, 0.0),
        Body::rocky(M_MERCURY).at_3d(r_peri, 0.0, 0.0).with_velocity_3d(0.0, v_peri, 0.0),
    ];

    let force = CentralForce::from_apsidal_rate(
        0, // source: Sun
        1, // target: Mercury
        OMEGA_DOT_TARGET,
        GAMMA,
        &bodies,
        units,
    )
    .expect("observable-inversion: gamma = -3 is non-degenerate, indices valid, Mercury bound");

    let mut sys =
        System::new(bodies, units).with_integrator(IntegratorKind::Ias15).with_dt(DT_YEARS);
    sys.add_hamiltonian_perturbation(Box::new(force))
        .expect("matching units; CentralForce must register");

    // Mercury orbital period in solar AU-year units: T = 2π · a^(3/2) / √(G M_sun)
    let period = 2.0 * std::f64::consts::PI * A_MERCURY.powf(1.5) / (g * M_SUN).sqrt();

    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(
        w,
        "# Cross-platform parity — apsis-central observable-inversion long horizon — apsis side"
    )
    .unwrap();
    writeln!(w, "# integrator: IAS15 (apsis) + CentralForce (apsis-central)").unwrap();
    writeln!(w, "# units: solar AU-year (G ≈ {g:.18e})").unwrap();
    writeln!(
        w,
        "# gamma={GAMMA}, omega_dot_target={OMEGA_DOT_TARGET:e} rad/yr, n_orbits={N_ORBITS}, period_yr={period:.6}"
    )
    .unwrap();
    writeln!(w, "# n_bodies={}", sys.bodies().len()).unwrap();
    write_header(&mut w, sys.bodies().len());

    write_sample(&mut w, 0, &sys);
    for orbit in 1..=N_ORBITS {
        let t_target = orbit as f64 * period;
        sys.integrate_until(t_target);
        write_sample(&mut w, orbit, &sys);
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {}", N_ORBITS + 1, output_path.display());
}

// ── Output helpers (parallel to other parity examples) ──────────────────── //

fn write_header(w: &mut BufWriter<File>, n_bodies: usize) {
    let mut header = String::from("orbit,t");
    for i in 0..n_bodies {
        header.push_str(&format!(",x{i},y{i},z{i},vx{i},vy{i},vz{i}"));
    }
    header.push_str(",e_total,lz_total");
    writeln!(w, "{header}").unwrap();
}

fn write_sample(w: &mut BufWriter<File>, orbit: u64, sys: &System) {
    let bodies = sys.bodies();
    let g = sys.units().g();
    let mut row = format!("{orbit},{:.18e}", sys.t());
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
    PathBuf::from("validation/cross-platform/windows/central_observable_inversion_long.csv")
}
