//! Pure Newtonian Kepler 2-body sanity check — no 1PN.
//!
//! If LRL precesses here, the bug is upstream of the perihelion example.

use gravity_sim_core::core::system::System;
use gravity_sim_core::domain::body::Body;
use gravity_sim_core::physics::integrator::IntegratorKind;

fn main() {
    // Circular orbit: e = 0, no reference axis for omega but LRL magnitude ≈ 0.
    // Eccentric orbit: e = 0.5, well-defined LRL.
    //
    // Run both; the LRL vector must be bit-stable across orbits for a
    // pure Newtonian 2-body with IAS15.
    for (label, e, m_secondary) in [
        ("circular (e=0, m_2=1e-7)", 0.0, 1e-7_f64),
        ("eccentric (e=0.5, m_2=1e-7)", 0.5, 1e-7),
        ("eccentric (e=0.5, m_2=1e-3)", 0.5, 1e-3),
    ] {
        println!("\n── {label} ──");
        let a = 1.0;
        let r_peri = a * (1.0 - e);
        let gm = 1.0 + m_secondary;
        let v_peri = (gm * (2.0 / r_peri - 1.0 / a)).sqrt();

        let sun = Body::star(1.0);
        let secondary = Body::rocky(m_secondary)
            .at(r_peri, 0.0)
            .with_velocity(0.0, v_peri);

        let mut sys = System::new(vec![sun, secondary])
            .with_integrator(IntegratorKind::Ias15)
            .with_dt(1e-4);

        let h0: f64 = r_peri * v_peri;
        let period = 2.0 * std::f64::consts::PI * (a * a * a / gm).sqrt();

        println!("  T = {period:.6}, v_peri = {v_peri:.6}, h0 = {h0:.6e}");

        for k in [1u64, 5, 10, 50, 100] {
            sys.integrate_until(period * (k as f64));
            let s = &sys.bodies()[0];
            let m = &sys.bodies()[1];
            let rx = m.x - s.x;
            let ry = m.y - s.y;
            let vrx = m.vx - s.vx;
            let vry = m.vy - s.vy;
            let r = (rx * rx + ry * ry).sqrt();
            let h = rx * vry - ry * vrx;
            let ex = vry * h / gm - rx / r;
            let ey = -vrx * h / gm - ry / r;
            let e_mag = (ex * ex + ey * ey).sqrt();
            let omega = ey.atan2(ex);
            println!(
                "  orbit {:>4}   h = {:+.6e}   |e| = {:.6e}   ex = {:+.4e}   ey = {:+.4e}   ω = {:+.4e}",
                k, h, e_mag, ex, ey, omega
            );
        }
    }
}
