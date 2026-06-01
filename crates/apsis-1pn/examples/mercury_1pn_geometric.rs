//! Geometric cross-check of the Mercury 1PN perihelion advance: the §3.2
//! Plummer convergence figure's observable (geometric apsidal angle from
//! periapsis-passage detection) applied to the 1PN-*satisfied* case, as a
//! counterpart to the osculating-ω end-vs-start gate in
//! `mercury_precession_gate.rs`.
//!
//! Sun–Mercury, exact kernel (ε=0), apsis-1pn 1PN with the gate's
//! `C_SOLAR_UNITS` c, IAS15 dt=1e-4. The closed-form Schwarzschild advance per
//! revolution 6πGM/(c²a(1−e²)) and the geometric apsidal precession per radial
//! period are the same convention, so the ratio is a direct test.
//!
//! Unlike the Plummer case (per-orbit signal ~1e-2 rad, geometric saturates the
//! f64 floor at ~1e-7), Mercury's per-orbit precession is tiny (~5e-7 rad), so
//! the geometric (per-orbit) observable is integration-noise-limited: rel_err
//! falls from ~2.6e-2 at N=500 to ~5.5e-4 at N=2e4, oscillating around the
//! prediction (consistent with a regression-precision floor ∝N^-3/2 below the
//! ~N=2000 crossover and an IAS15 random-walk floor ∝N^-1/2 above), and never
//! reaches the test-particle physical ceiling m_Mercury/M_sun ≈ 1.7e-7 at
//! practical N. The osculating cumulative observable (the gate) is
//! correspondingly tighter (1e-4 at N=500).
//!
//! Diagnostic, not a gate; edit N_ORBITS to reproduce the N-series.
//! Run: cargo run --release --example mercury_1pn_geometric -p apsis-1pn

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::geometric_apsidal_precession_per_radial;
use apsis::units::UnitSystem;
use apsis_1pn::{C_SOLAR_UNITS, PostNewtonian1PN};

const A: f64 = 0.387_098;
const E: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;
const N_ORBITS: u64 = 500;
const RAD_TO_ARCSEC: f64 = 180.0 * 3600.0 / PI;

fn main() {
    let units = UnitSystem::solar_canonical();
    let sun = Body::star(1.0);
    let r_peri = A * (1.0 - E);
    let v_peri = (2.0 / r_peri - 1.0 / A).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys =
        System::new(vec![sun, mercury], units).with_integrator(IntegratorKind::Ias15).with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::from_raw_c(C_SOLAR_UNITS, units)))
        .expect("Sun-Mercury 1PN: System and operator share solar_canonical");

    // Geometric apsidal precession per radial period (rad). Exact Kepler kernel,
    // so the only precession source is 1PN.
    let measured_per_radial =
        geometric_apsidal_precession_per_radial(&mut sys, 1, 0, 1.0, N_ORBITS);

    let c = C_SOLAR_UNITS;
    let predicted_per_radial = 6.0 * PI / (c * c * A * (1.0 - E * E));

    let rel_err = (measured_per_radial / predicted_per_radial - 1.0).abs();
    let n = N_ORBITS as f64;

    println!("Mercury 1PN — geometric apsidal precession (periapsis-passage)");
    println!("  measured  (per radial period) = {measured_per_radial:.9e} rad");
    println!("  predicted (Schwarzschild)     = {predicted_per_radial:.9e} rad");
    println!("  rel_err                       = {rel_err:.6e}");
    println!(
        "  cumulative measured  / {N_ORBITS} orbits = {:.6} arcsec",
        measured_per_radial * n * RAD_TO_ARCSEC
    );
    println!(
        "  cumulative predicted / {N_ORBITS} orbits = {:.6} arcsec",
        predicted_per_radial * n * RAD_TO_ARCSEC
    );
}
