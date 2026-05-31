//! Softened-Plummer apsidal precession measured by apsis across a softening
//! sweep — the measured-data side of the §3.2 convergence figure.
//!
//! For each softening ε, integrates Sun–Mercury under a softened
//! `NewtonKernel(ε)` (pure Newtonian, NO 1PN) with IAS15 and measures the
//! GEOMETRIC apsidal precession per radial period — the same observable the
//! full-potential apsidal-angle quadrature computes, so apsis is compared to
//! the oracle in the oracle's own convention.
//!
//! Why geometric (periapsis-passage) and not osculating-ω drift: the osculating
//! argument of periapsis oscillates within an orbit, and sampling it at integer
//! Kepler periods (T_radial ≠ T_Kepler) folds that wiggle into a sign-
//! oscillating, N-dependent stroboscopic artifact — not a property of apsis's
//! fidelity. Measuring the angle swept between successive true periapses (ṙ: −→+)
//! removes both: the radial-period boundary is the physical one and no osculating
//! element enters.
//!
//! Method: integrate in fine sub-steps, accumulate the continuously-unwrapped
//! position angle Θ(t) (the total angle swept), and detect periapsis passages.
//! Over K radial periods the body sweeps K·(2π + Δϖ), so
//!     Δϖ_per_radial = Θ(periapsis K) / K − 2π.
//! Θ at the IC periapsis is 0 by construction (Mercury starts at periapsis on
//! +x); the endpoint Θ is linearly interpolated to the ṙ=0 crossing, an error
//! that enters divided by K and is therefore negligible.
//!
//! 1PN is OFF: the quadrature oracle is pure softened-Plummer, and the fixed GR
//! signal would contaminate the ε²-shrinking softening precession at small ε.
//!
//! ## Run
//!
//! ```text
//! cargo run --release --example softened_plummer_sweep -p apsis -- \
//!     --output paper/figures/data/apsis_softened_sweep.csv
//! ```

use std::env;
use std::f64::consts::PI;
use std::fs::{File, create_dir_all};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::NewtonKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;

const A_MERCURY: f64 = 0.387_098;
const E_MERCURY: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;
const DEFAULT_N_RADIAL: u64 = 300;
const DT: f64 = 1.0e-4;
/// Sampling sub-steps per Kepler period for the periapsis search and the
/// continuous angle unwrap. 360 keeps the per-step swept angle well under π even
/// at periapsis (fastest point), so the unwrap never skips a turn.
const SAMPLES_PER_PERIOD: f64 = 360.0;

/// Relative-orbit angle, radial velocity, and separation (Mercury w.r.t. Sun).
fn state(sys: &System) -> (f64, f64) {
    let b = sys.bodies();
    let dx = b[1].pos_x - b[0].pos_x;
    let dy = b[1].pos_y - b[0].pos_y;
    let dvx = b[1].vel_x - b[0].vel_x;
    let dvy = b[1].vel_y - b[0].vel_y;
    let r = (dx * dx + dy * dy).sqrt();
    let rdot = (dx * dvx + dy * dvy) / r;
    (dy.atan2(dx), rdot)
}

fn wrap_to_pi(mut d: f64) -> f64 {
    while d > PI {
        d -= 2.0 * PI;
    }
    while d <= -PI {
        d += 2.0 * PI;
    }
    d
}

/// Geometric apsidal precession per radial period (rad) under a softened
/// `NewtonKernel(ε)`, from periapsis-passage detection. Pure Newtonian, no 1PN.
fn measure_geometric_precession_per_radial_rad(epsilon: f64, n_radial: u64) -> f64 {
    let sun = Body::star(1.0);
    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = (2.0 / r_peri - 1.0 / A_MERCURY).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_kernel(Arc::new(NewtonKernel::new(epsilon)))
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(DT);

    let period = compute_elements(sys.bodies(), 1, 0, 1.0).expect("initial elements").period;
    let h = period / SAMPLES_PER_PERIOD;

    // IC is the 0th periapsis: Mercury on +x with ṙ=0, so Θ=0, θ=0.
    let (mut theta_last, mut rdot_prev) = state(&sys);
    let mut theta_cont = 0.0_f64;
    let mut peri_count = 0_u64;
    let mut t = 0.0_f64;

    loop {
        t += h;
        sys.integrate_until(t);
        let (theta, rdot) = state(&sys);
        let dtheta = wrap_to_pi(theta - theta_last);
        theta_cont += dtheta;
        theta_last = theta;

        if rdot_prev < 0.0 && rdot >= 0.0 {
            peri_count += 1;
            if peri_count == n_radial {
                // Linear interp of Θ to the ṙ=0 crossing inside [prev, now]:
                let frac = rdot_prev / (rdot_prev - rdot); // ∈ [0, 1]
                let theta_cont_prev = theta_cont - dtheta;
                let theta_peri = theta_cont_prev + frac * dtheta;
                return theta_peri / (n_radial as f64) - 2.0 * PI;
            }
        }
        rdot_prev = rdot;
    }
}

/// Log-spaced softening sweep spanning both honest edges: small ε where the
/// geometric precession signal (∝ε²) approaches the measurement resolution, and
/// large ε where the leading closed form departs from the exact oracle. ε in AU;
/// r_peri = a(1−e) ≈ 0.307 AU, so the largest ε is ~1/3 of periapsis.
fn epsilon_sweep() -> Vec<f64> {
    const N: usize = 17;
    const LOG_LO: f64 = -3.0; // 1e-3 AU
    const LOG_HI: f64 = -1.0; // 1e-1 AU
    (0..N)
        .map(|i| 10f64.powf(LOG_LO + (LOG_HI - LOG_LO) * (i as f64) / ((N - 1) as f64)))
        .collect()
}

fn main() {
    let (output_path, n_radial) = parse_args();
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent).expect("failed to create output directory");
    }

    let file = File::create(&output_path).expect("failed to open output file");
    let mut w = BufWriter::new(file);

    writeln!(w, "# Softened-Plummer apsidal precession — apsis measured (geometric, no 1PN)")
        .unwrap();
    writeln!(w, "# observable: geometric apsidal angle per radial period (periapsis-passage)")
        .unwrap();
    writeln!(w, "# protocol: IAS15, dt={DT:e}, n_radial={n_radial}, Sun-Mercury softened kernel")
        .unwrap();
    writeln!(w, "# units: solar-canonical (AU, yr/2pi, Msun, G=1)").unwrap();
    writeln!(w, "# a={A_MERCURY}, e={E_MERCURY}, m_mercury={M_MERCURY:e}").unwrap();
    writeln!(w, "eps,precession_per_radial_rad").unwrap();

    for eps in epsilon_sweep() {
        let prec = measure_geometric_precession_per_radial_rad(eps, n_radial);
        writeln!(w, "{eps:.18e},{prec:.18e}").unwrap();
        eprintln!("  eps={eps:.4e} AU -> {prec:.9e} rad/radial-period");
    }

    w.flush().unwrap();
    eprintln!("wrote softened-Plummer geometric sweep to {}", output_path.display());
}

fn parse_args() -> (PathBuf, u64) {
    let mut output = PathBuf::from("paper/figures/data/apsis_softened_sweep.csv");
    let mut n_radial = DEFAULT_N_RADIAL;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--output" | "-o" => {
                output = PathBuf::from(args.next().expect("--output requires a path argument"));
            },
            "--n-radial" => {
                n_radial = args
                    .next()
                    .expect("--n-radial requires a value")
                    .parse()
                    .expect("--n-radial must be a positive integer");
            },
            _ => {},
        }
    }
    (output, n_radial)
}
