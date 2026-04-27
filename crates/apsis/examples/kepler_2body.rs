//! Two-body Keplerian orbit — low-level API demonstration.
//!
//! Run with:
//!
//! ```text
//! cargo run --example kepler_2body --release
//! ```
//!
//! Demonstrates direct body construction, integrator choice, and the
//! `integrate_for` + energy-drift query path without going through the
//! template catalog.

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

fn main() {
    // Unit-mass sun at the origin, unit-mass "earth" in a circular orbit
    // of radius 1 with v_circ = 1 (G = 1 in simulation units).
    let sun = Body::star(1.0);
    let earth = Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 1.0);

    let mut sys =
        System::new(vec![sun, earth], UnitSystem::canonical()).with_integrator(IntegratorKind::Ias15).with_dt(1e-3);

    // Integrate for ~16 orbital periods (T = 2π).
    const T_END: f64 = 100.0;
    let steps = sys.integrate_for(T_END);

    println!("Kepler 2-body @ IAS15");
    println!("  t_end       = {:.4e}", sys.t());
    println!("  steps       = {}", steps);
    println!("  energy      = {:+.6e}", sys.energy());
    println!("  dE/E        = {:+.3e}", sys.energy_delta());
    println!("  Lz          = {:+.6e}", sys.lz());
    println!("  dLz/Lz      = {:+.3e}", sys.lz_delta());
}
