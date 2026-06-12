//! Phase-B error-budget run — one parameterised Mercury 1PN gate measurement.
//!
//! Outputs one CSV line to stdout:
//!
//! ```text
//! orbits,ulp,constructor,measured_rad,predicted_rad,rel_err
//! ```
//!
//! Run (release mode required for gate fidelity):
//!
//! ```text
//! cargo run --release --example error_budget_run -p apsis-1pn -- \
//!     --orbits 500 --ulp 0 --constructor raw_c
//! ```
//!
//! # Arguments
//!
//! * `--orbits N` — number of Mercury orbits to integrate (default 500)
//! * `--ulp K` — signed integer; perturbs Mercury's initial x-position
//!   by K ULPs (default 0; K may be negative)
//! * `--constructor` — `for_units` or `raw_c` (default `raw_c`)
//!
//! # ULP perturbation
//!
//! Mercury starts at `r_peri > 0`, so the IEEE-754 bit pattern is an
//! ordinary positive double and bit-level monotonicity holds. For `K ≥ 0`
//! the bits are incremented by K; for `K < 0` the bits are decremented by
//! |K|. Both directions keep `r_peri` positive for any |K| ≤ ~1e15.
//!
//! # Constructor conventions
//!
//! * `raw_c` — `PostNewtonian1PN::from_raw_c(C_SOLAR_UNITS, …)`:
//!   IAU julian-year literal; matches the gate.
//! * `for_units` — `PostNewtonian1PN::for_units(UnitSystem::solar_canonical())`:
//!   derives `c` from Gaussian time (`sqrt(AU³/GM_sun)`);
//!   differs from `C_SOLAR_UNITS` by ~19 ppm (ADR-014).
//!
//! # Predicted advance
//!
//! `predicted = 6π / (c² · A · (1 − E²)) · N`
//! using the constructor's own `c` (`pn.c()`), matching the gate's oracle.

use std::f64::consts::PI;

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;
use apsis_1pn::{C_SOLAR_UNITS, PostNewtonian1PN};

// ── Gate constants (mirrors mercury_precession_gate.rs exactly) ────────────────
const A: f64 = 0.387_098;
const E: f64 = 0.205_63;
const M_MERCURY: f64 = 1.660_114e-7;

fn main() {
    let (n_orbits, ulp_k, use_for_units) = parse_args();

    // ── Initial conditions ────────────────────────────────────────────────────
    let sun = Body::star(1.0);
    let r_peri_base = A * (1.0 - E);
    let v_peri = (2.0 / r_peri_base - 1.0 / A).sqrt();

    // ULP perturbation: r_peri > 0 so bit-monotonicity holds throughout.
    let r_peri = apply_ulp(r_peri_base, ulp_k);

    let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

    // ── Build PostNewtonian1PN ────────────────────────────────────────────────
    let pn = if use_for_units {
        PostNewtonian1PN::for_units(UnitSystem::solar_canonical())
    } else {
        PostNewtonian1PN::from_raw_c(C_SOLAR_UNITS, UnitSystem::solar_canonical())
    };
    let c = pn.c();
    let constructor_label = if use_for_units { "for_units" } else { "raw_c" };

    // ── Build System ──────────────────────────────────────────────────────────
    let mut sys = System::new(vec![sun, mercury], UnitSystem::solar_canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);
    sys.add_hamiltonian_perturbation(Box::new(pn))
        .expect("error_budget_run: shared UnitSystem::solar_canonical()");

    // ── Reference state ───────────────────────────────────────────────────────
    // mu = 1.0, matching the gate's compute_elements call.
    let el0 = compute_elements(sys.bodies(), 1, 0, 1.0)
        .expect("error_budget_run: bound elements at t = 0");

    // ── Integrate ─────────────────────────────────────────────────────────────
    sys.integrate_for(el0.period * (n_orbits as f64));

    // ── Measure perihelion advance ────────────────────────────────────────────
    let el_end = compute_elements(sys.bodies(), 1, 0, 1.0)
        .expect("error_budget_run: bound elements at t = end");

    // Wrap to (−π, π] — identical logic to the gate.
    let measured = {
        let mut d = el_end.omega - el0.omega;
        while d > PI {
            d -= 2.0 * PI;
        }
        while d <= -PI {
            d += 2.0 * PI;
        }
        d
    };

    // ── Predicted advance (first-order GR formula) ────────────────────────────
    // Uses pn.c(), not the literal C_SOLAR_UNITS, so for_units and raw_c
    // produce their own self-consistent predictions.
    let predicted = 6.0 * PI / (c * c * A * (1.0 - E * E)) * (n_orbits as f64);

    let rel_err = (measured - predicted).abs() / predicted.abs();

    // ── Output ────────────────────────────────────────────────────────────────
    println!(
        "{},{},{},{:.17e},{:.17e},{:.17e}",
        n_orbits, ulp_k, constructor_label, measured, predicted, rel_err
    );
}

/// Perturb an f64 value by `k` ULPs. The value must be positive so that
/// bit-level monotonicity holds (IEEE-754 positive doubles have the same
/// ordering as their bit patterns). Negative `k` decrements the bit pattern.
fn apply_ulp(x: f64, k: i64) -> f64 {
    debug_assert!(x > 0.0, "apply_ulp: x must be positive for bit-monotonicity");
    let bits = x.to_bits();
    let perturbed =
        if k >= 0 { bits.wrapping_add(k as u64) } else { bits.wrapping_sub(k.unsigned_abs()) };
    f64::from_bits(perturbed)
}

/// Parse `--orbits N`, `--ulp K`, `--constructor for_units|raw_c` from argv.
/// Unknown / extra arguments are silently ignored so the binary stays terse.
fn parse_args() -> (u64, i64, bool) {
    let args: Vec<String> = std::env::args().collect();
    let mut n_orbits: u64 = 500;
    let mut ulp_k: i64 = 0;
    let mut use_for_units = false; // default: raw_c

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--orbits" => {
                i += 1;
                n_orbits = args[i].parse().expect("--orbits requires a non-negative integer");
            },
            "--ulp" => {
                i += 1;
                ulp_k = args[i].parse().expect("--ulp requires a signed integer");
            },
            "--constructor" => {
                i += 1;
                match args[i].as_str() {
                    "for_units" => use_for_units = true,
                    "raw_c" => use_for_units = false,
                    other => panic!("--constructor must be `for_units` or `raw_c`, got `{other}`"),
                }
            },
            _ => {},
        }
        i += 1;
    }
    (n_orbits, ulp_k, use_for_units)
}
