//! Three-body figure-eight choreographic orbit (Chenciner & Montgomery 2000).
//!
//! Run with:
//!
//! ```text
//! cargo run --example figure_eight --release
//! ```
//!
//! The famous stable periodic three-body configuration. IAS15 with strict dt
//! preserves the orbit's closed-loop signature over hundreds of periods —
//! a good sanity check on the integrator-plus-core stack.

use gravity_sim_core::core::system::System;
use gravity_sim_core::physics::integrator::IntegratorKind;
use gravity_sim_core::templates::TemplateKind;

fn main() {
    let mut sys = System::from_template(TemplateKind::ThreeBodyFigureEight)
        .with_integrator(IntegratorKind::Ias15)
        .with_dt(1e-4);

    // Figure-eight period ≈ 6.3259; integrate for 100 periods.
    const T_END: f64 = 6.3259 * 100.0;

    println!("Figure-eight @ IAS15,  target t = {:.2}", T_END);
    println!("  starting E = {:+.6e}", sys.energy());

    sys.integrate_for(T_END);

    println!(
        "  ended at t = {:.4}   dE/E = {:+.3e}   dLz/Lz = {:+.3e}",
        sys.t(),
        sys.energy_delta(),
        sys.lz_delta(),
    );
}
