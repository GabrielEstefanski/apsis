//! Solar-system preset integrated over a long horizon — template + IAS15.
//!
//! Run with:
//!
//! ```text
//! cargo run --example solar_system_long --release
//! ```
//!
//! Demonstrates the high-level `from_template` path: a researcher reproducing
//! a known system in one line, picking an integrator, and printing the two
//! conservation diagnostics that a methods paper would cite.

use apsis::core::system::System;
use apsis::physics::integrator::IntegratorKind;
use apsis::templates::TemplateKind;

fn main() {
    let mut sys = System::from_template(TemplateKind::SolarSystem)
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-3);

    let n0 = sys.bodies().len();
    let e0 = sys.energy();

    // ~16 Jupiter orbits (T_jup ≈ 11.86 yr, simulation year = 2π).
    const T_END: f64 = 1200.0;

    println!("Solar System @ IAS15,  N = {},  target t = {}", n0, T_END);
    println!("  starting E  = {:+.6e}", e0);

    // Dump drift every tenth of the window so regressions show up
    // as a quality trend, not a single final number.
    for k in 1..=10 {
        let t_target = T_END * (k as f64) / 10.0;
        sys.integrate_until(t_target);
        println!(
            "  t = {:>7.1}   dE/E = {:+.3e}   dLz/Lz = {:+.3e}",
            sys.t(),
            sys.energy_delta(),
            sys.lz_delta(),
        );
    }
}
