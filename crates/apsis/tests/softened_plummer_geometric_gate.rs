//! Release gate for the §3.2 geometric apsidal-precession claim.
//!
//! paper.md §3.2 argues the osculating-ω gate's 0.04 % residual is an
//! osculating-vs-geometric *definition* difference, not a dynamical error,
//! because apsis's geometric apsidal precession reproduces the independent
//! full-potential quadrature oracle to ~1e-7. That ~1e-7 is a load-bearing
//! premise of the paper's argument, so it carries its own gate here.
//!
//! The oracle is the integrator-independent Gauss–Legendre apsidal-angle
//! quadrature (`paper/notebooks/scripts/plummer_apsidal_quadrature.py`), pinned
//! as a constant exactly as `QUADRATURE_DRIFT_PER_KEPLER_RAD` is in the
//! osculating gate; the apsis side is measured fresh through the same shared
//! `geometric_apsidal_precession_per_radial` the figure sweep uses. Pure
//! `NewtonKernel` (no 1PN), so this lives in the apsis crate, not apsis-1pn.

use std::sync::Arc;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::gravity::NewtonKernel;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::geometric_apsidal_precession_per_radial;
use apsis::units::UnitSystem;

const A_MERCURY: f64 = 0.387_098;
const E_MERCURY: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;
const DT: f64 = 1.0e-4;
const N_RADIAL: u64 = 300;

// Geometric apsidal precession per radial period (rad) from the full-potential
// Gauss–Legendre quadrature, at ε = 10^-1.75 and 10^-1.625 — mid-band points of
// the figure sweep where the measured deviation floors at ~1.2–1.4e-7,
// bit-identical Windows/Linux.
const ORACLE_EPS_NEG_1_75: f64 = -0.021_372_842_436_778_505;
const ORACLE_EPS_NEG_1_625: f64 = -0.037_587_403_022_007_14;

fn measure_geometric_precession(epsilon: f64) -> f64 {
    let sun = Body::star(1.0);
    let r_peri = A_MERCURY * (1.0 - E_MERCURY);
    let v_peri = (2.0 / r_peri - 1.0 / A_MERCURY).sqrt();
    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_kernel(Arc::new(NewtonKernel::new(epsilon)))
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(DT);

    geometric_apsidal_precession_per_radial(&mut sys, 1, 0, 1.0, N_RADIAL)
}

/// apsis's geometric apsidal precession reproduces the independent quadrature
/// oracle to < 1e-6. A *dynamical* error would show at the closed-form scale
/// (~1e-2); the 1e-6 bound cleanly separates "geometric agreement holds — hence
/// the §3.2 osculating residual is definitional, not a dynamics error" from a
/// real discrepancy. Observed dev ~1.2–1.4e-7 (~8× headroom), bit-identical
/// across Windows and Linux despite the harness's non-hardened `atan2`.
#[test]
#[ignore = "release-mode integration test; run with `cargo test --release -- --ignored`"]
fn geometric_apsidal_precession_matches_quadrature() {
    for (log10_eps, oracle) in [(-1.75_f64, ORACLE_EPS_NEG_1_75), (-1.625, ORACLE_EPS_NEG_1_625)] {
        let eps = 10f64.powf(log10_eps);
        let measured = measure_geometric_precession(eps);
        let dev = (measured / oracle - 1.0).abs();

        eprintln!(
            "[geometric-gate] eps={eps:.9e}: measured={measured:.12e} \
             oracle={oracle:.12e} dev={dev:.3e}",
        );
        assert!(
            dev < 1e-6,
            "geometric apsidal precession dev {dev:.3e} exceeds 1e-6 at eps={eps:.6e} \
             (measured {measured:.9e} vs quadrature oracle {oracle:.9e})",
        );
    }
}
