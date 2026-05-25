//! Pythagorean three-body problem (Burrau 1913) — stress test for the
//! adaptive step controller during violent close encounters.
//!
//! Run with:
//!
//! ```text
//! cargo run --example pythagorean_close_encounter --release
//! ```
//!
//! Masses 3, 4, 5 at rest on the vertices of a 3-4-5 triangle. The subsequent
//! evolution is chaotic with multiple tight encounters in quick succession —
//! exactly the regime where symplectic integrators fail and IAS15's step-size
//! controller earns its keep.

use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;

fn main() {
    let bodies = vec![
        Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0),
        Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0),
        Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0),
    ];

    let mut sys = System::new(bodies, UnitSystem::canonical())
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(0.1);

    println!("Pythagorean 3-body @ IAS15");
    println!("  starting E  = {:+.6e}", sys.energy());

    // Burrau's window — the strongest close-encounter chain happens in [2, 5].
    const T_END: f64 = 10.0;
    for t_target in [1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, T_END] {
        sys.integrate_until(t_target);
        println!("  t = {:>5.2}   dE/E = {:+.3e}", sys.t(), sys.energy_delta().unwrap_or(f64::NAN),);
    }
}
