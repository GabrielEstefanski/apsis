//! REBOUNDx parity — Sun–Mercury 1PN, apsis-1pn side.
//!
//! Runs Sun–Mercury under apsis IAS15 with the `apsis-1pn` first
//! post-Newtonian operator (test-particle Schwarzschild form, applied
//! pairwise) for `N_ORBITS` orbital periods, sampling state and total
//! (Newtonian) energy at the end of each orbit. Output is a CSV consumable
//! by the matching `reboundx_side.py` harness, which runs the REBOUNDx `gr`
//! effect (single-dominant-mass formulation) on identical initial
//! conditions for cross-implementation parity.
//!
//! The two 1PN formulations differ in gauge, coordinates (inertial vs
//! Jacobi) and solve (explicit vs iterative); for Sun–Mercury (mass ratio
//! ~1.7e-7) both reduce to the test-particle Schwarzschild limit, so the
//! comparison is a measurement of the formulation/gauge difference, NOT a
//! bit-parity check. The analytic 43"/century is the primary anchor.
//!
//! ## Run
//!
//! ```text
//! # 1PN on (the measurement):
//! cargo run --release --example reboundx_parity_gr -p apsis-1pn -- \
//!     --output validation/reboundx-parity/gr-mercury/out/apsis.csv
//! # 1PN off (the harness control — should match REBOUND to the ULP floor):
//! cargo run --release --example reboundx_parity_gr -p apsis-1pn -- \
//!     --no-1pn --output validation/reboundx-parity/gr-mercury/out/apsis_kepler.csv
//! ```

use std::env;
use std::f64::consts::PI;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_1pn::{C_SOLAR_UNITS, PostNewtonian1PN};

// ── Protocol constants (Sun–Mercury; mirrored in reboundx_side.py) ──────── //

/// Semi-major axis (AU, solar-canonical units; G = 1).
const A_MERCURY: f64 = 0.387_098;
/// Eccentricity.
const E_MERCURY: f64 = 0.205_63;
/// Primary (Sun) mass (M_sun).
const M_SUN: f64 = 1.0;
/// Secondary (Mercury) mass (M_sun).
const M_MERCURY: f64 = 1.660_114e-7;
/// Number of orbital periods integrated.
const N_ORBITS: u64 = 500;
/// Initial timestep, as a fraction of the orbital period.
const DT_FRACTION_OF_PERIOD: f64 = 1.0e-3;

fn main() {
    let (output_path, with_1pn) = parse_args();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    // ── Initial conditions (COM frame, zero net momentum) ───────────────── //
    //
    // Place Mercury at periapsis; offset both bodies so the centre of mass is
    // at the origin and net momentum is zero — eliminates COM drift as a |Δr|
    // source. Both implementations must evaluate this exact f64 expression.
    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = ((1.0 + E_MERCURY) / (A_MERCURY * (1.0 - E_MERCURY))).sqrt();

    let m_total = M_SUN + M_MERCURY;
    let sun_x = -(M_MERCURY / m_total) * r_peri;
    let sun_vy = -(M_MERCURY / m_total) * v_peri;
    let mercury_x = (M_SUN / m_total) * r_peri;
    let mercury_vy = (M_SUN / m_total) * v_peri;

    let sun = Body::star(M_SUN).at(sun_x, 0.0).with_velocity(0.0, sun_vy);
    let mercury = Body::rocky(M_MERCURY).at(mercury_x, 0.0).with_velocity(0.0, mercury_vy);

    let units = UnitSystem::solar_canonical();
    let period = 2.0 * PI * A_MERCURY.powf(1.5); // G = 1, M ≈ 1
    let dt0 = period * DT_FRACTION_OF_PERIOD;

    let mut sys =
        System::new(vec![sun, mercury], units).with_integrator(IntegratorKind::Ias15).with_dt(dt0);

    if with_1pn {
        sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(
            C_SOLAR_UNITS,
            units,
        )))
        .expect("apsis-1pn shares UnitSystem::solar_canonical() with the system");
    }

    // ── CSV output (schema identical to the rebound-parity scenarios) ────── //
    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    let mode = if with_1pn {
        "apsis-1pn (test-particle pairwise)"
    } else {
        "Newtonian (1PN off, control)"
    };
    writeln!(w, "# REBOUNDx parity — Sun-Mercury 1PN — apsis side").unwrap();
    writeln!(w, "# protocol: paper/notebooks/2026-05-29-reboundx-parity-gr.md").unwrap();
    writeln!(w, "# integrator: IAS15 (apsis); force: {mode}").unwrap();
    writeln!(w, "# units: solar-canonical (AU, yr/2pi, Msun, G=1)").unwrap();
    writeln!(w, "# a={A_MERCURY}, e={E_MERCURY}, m_sun={M_SUN}, m_mercury={M_MERCURY:e}").unwrap();
    writeln!(w, "# c={C_SOLAR_UNITS:.18e}").unwrap();
    writeln!(w, "# period={period:.18e}, dt0={dt0:.18e}, n_orbits={N_ORBITS}, with_1pn={with_1pn}")
        .unwrap();
    writeln!(w, "orbit,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,e_total").unwrap();

    write_sample(&mut w, 0, &sys);
    for orbit in 1..=N_ORBITS {
        sys.integrate_until(period * (orbit as f64));
        write_sample(&mut w, orbit, &sys);
    }

    w.flush().unwrap();
    eprintln!("wrote {} samples to {} (with_1pn={with_1pn})", N_ORBITS + 1, output_path.display());
}

fn write_sample(w: &mut BufWriter<File>, orbit: u64, sys: &System) {
    let bodies = sys.bodies();
    let (b0, b1) = (&bodies[0], &bodies[1]);
    let e_total = total_energy(bodies);
    writeln!(
        w,
        "{orbit},{t:.18e},{x0:.18e},{y0:.18e},{vx0:.18e},{vy0:.18e},{x1:.18e},{y1:.18e},{vx1:.18e},{vy1:.18e},{e:.18e}",
        t = sys.t(),
        x0 = b0.pos_x, y0 = b0.pos_y, vx0 = b0.vel_x, vy0 = b0.vel_y,
        x1 = b1.pos_x, y1 = b1.pos_y, vx1 = b1.vel_x, vy1 = b1.vel_y,
        e = e_total,
    )
    .unwrap();
}

/// Newtonian total energy (matches REBOUND's `sim.energy()`, which is also
/// Newtonian — REBOUNDx effects do not modify it). Under 1PN this is not
/// conserved; both sides drift the same way if the dynamics match.
fn total_energy(bodies: &[Body]) -> f64 {
    let ke: f64 =
        bodies.iter().map(|b| 0.5 * b.mass * (b.vel_x * b.vel_x + b.vel_y * b.vel_y)).sum();
    let mut pe = 0.0;
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let r = (dx * dx + dy * dy).sqrt();
            pe -= bodies[i].mass * bodies[j].mass / r;
        }
    }
    ke + pe
}

fn parse_args() -> (PathBuf, bool) {
    let mut output = PathBuf::from("validation/reboundx-parity/gr-mercury/out/apsis.csv");
    let mut with_1pn = true;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" | "-o" => {
                output = PathBuf::from(args.next().expect("--output requires a path argument"));
            },
            "--no-1pn" => with_1pn = false,
            _ => {},
        }
    }
    (output, with_1pn)
}
