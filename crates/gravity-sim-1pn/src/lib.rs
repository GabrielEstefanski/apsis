//! First post-Newtonian (1PN) gravitational correction — out-of-tree plugin
//! demonstrating the [`PerturbationForce`] extension point of `gravity-sim-core`.
//!
//! # What this crate is for
//!
//! This crate's reason to exist is to *prove*, against the workspace's
//! published API surface, that an external consumer can extend the physics
//! of the simulator **without any `pub(crate)` access, without touching the
//! core sources, and without any dependency other than `gravity-sim-core`
//! itself**. It is the concrete answer to the paper's Phase 3 claim of a
//! "public API with an out-of-tree demonstration".
//!
//! If a future change to the core breaks the API this crate depends on, the
//! CI build of `gravity-sim-1pn` fails — surfacing the contract violation
//! loudly instead of quietly.
//!
//! # Physics
//!
//! The implementation uses the **test-particle 1PN (Schwarzschild) formula**
//! applied pairwise: for every receiver `i` and every source `j ≠ i`,
//!
//! ```text
//!   a_1PN(i←j) = G m_j / (c² r²) · [ (4 G m_j / r − v_i²) · r̂ + 4 (r̂ · v_i) v_i ]
//! ```
//!
//! where `r = r_j − r_i`, `r̂` is the corresponding unit vector, `v_i` is the
//! receiver's velocity, and `c` is the speed of light in simulation units.
//!
//! This is the Schwarzschild limit of the full Einstein–Infeld–Hoffmann (EIH)
//! equations. It is **exact** for a test mass around a dominant source
//! (e.g. Mercury–Sun with `m_Mercury / m_Sun ≈ 2 × 10⁻⁷`) and recovers the
//! GR perihelion precession `Δφ = 6π G M / (c² a (1−e²))` per orbit at
//! leading order.
//!
//! For equal-mass binaries the full EIH cross-terms matter and this crate's
//! result is approximate. That regime is out of scope for the current
//! demonstration; the paper's Mercury test sits squarely in the
//! test-particle regime where the simplified form is canonical.
//!
//! # Usage
//!
//! ```ignore
//! use gravity_sim_core::core::system::System;
//! use gravity_sim_core::domain::body::Body;
//! use gravity_sim_core::physics::integrator::IntegratorKind;
//! use gravity_sim_1pn::PostNewtonian1PN;
//!
//! let sun     = Body::star(1.0);
//! let mercury = Body::rocky(1.66e-7).at(0.307, 0.0).with_velocity(0.0, 2.078);
//! let mut sys = System::new(vec![sun, mercury])
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-4);
//!
//! sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));
//!
//! sys.integrate_for(200.0 * std::f64::consts::PI);
//! println!("dE/E = {:.3e}", sys.energy_delta());
//! ```

#![deny(unsafe_code)]

use gravity_sim_core::domain::body::Body;
use gravity_sim_core::physics::integrator::PerturbationForce;

/// Speed of light expressed in the simulator's canonical solar-system units:
///
/// | Unit      | Value                   |
/// |-----------|-------------------------|
/// | Length    | 1 AU                    |
/// | Mass      | 1 M☉                    |
/// | Time      | 1 year / (2π)           |
/// | G         | 1                       |
///
/// Computed at compile time from SI constants so the derivation, not a
/// hand-transcribed literal, is the source of truth:
///
/// ```text
///   C_SOLAR_UNITS = c_SI · (year_s / 2π) / AU_SI
/// ```
///
/// Current value: approximately `10_065.130` AU per (year / 2π).
pub const C_SOLAR_UNITS: f64 = {
    const C_SI: f64 = 299_792_458.0; // m/s, exact by SI definition
    const AU_SI: f64 = 149_597_870_700.0; // m, IAU 2012
    const YEAR_S: f64 = 365.25 * 86_400.0;
    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
    C_SI * (YEAR_S / TWO_PI) / AU_SI
};

/// First post-Newtonian gravitational correction (Schwarzschild, test-particle
/// form applied pairwise).
///
/// Register via [`System::add_perturbation`](gravity_sim_core::core::system::System::add_perturbation):
///
/// ```ignore
/// sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));
/// ```
///
/// Stateless. Safe to share across threads.
#[derive(Debug, Clone, Copy)]
pub struct PostNewtonian1PN {
    /// Speed of light in the caller's unit system.
    c: f64,
}

impl PostNewtonian1PN {
    /// Construct with an explicit speed of light.
    ///
    /// Use this when the simulation runs in a unit system other than the
    /// canonical solar one — e.g. geometric units (c = 1), or SI.
    pub const fn with_c(c: f64) -> Self {
        Self { c }
    }

    /// Construct for the simulator's canonical solar-system units
    /// ([`C_SOLAR_UNITS`]).
    pub const fn solar_units() -> Self {
        Self::with_c(C_SOLAR_UNITS)
    }

    /// Speed of light this instance was configured with.
    pub const fn c(&self) -> f64 {
        self.c
    }
}

impl PerturbationForce for PostNewtonian1PN {
    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)]) {
        debug_assert_eq!(
            bodies.len(),
            scratch_acc.len(),
            "PerturbationForce contract: scratch_acc must be sized to bodies"
        );

        let c2 = self.c * self.c;

        // Schwarzschild-gauge 1PN acceleration on receiver i due to source j:
        //
        //     a_1PN = (G m_j / c² r²) · [ (4 G m_j / r − v²) n̂ + 4 (n̂·v) v ]
        //
        // where n̂ points FROM SOURCE TO RECEIVER. The vector we compute below
        // (`rhat_*`) points FROM RECEIVER TO SOURCE, i.e. r̂ = −n̂. Substituting
        // r̂ = −n̂ flips the sign of the whole bracket (both n̂ and n̂·v change
        // sign; the outer product n̂·v × v inherits one flip), so in this
        // crate's convention the expression is applied with an overall
        // minus sign.
        for i in 0..bodies.len() {
            let b_i = &bodies[i];
            let vx_i = b_i.vx;
            let vy_i = b_i.vy;
            let v2_i = vx_i * vx_i + vy_i * vy_i;

            let mut ax = 0.0_f64;
            let mut ay = 0.0_f64;

            for j in 0..bodies.len() {
                if i == j {
                    continue;
                }
                let b_j = &bodies[j];
                let dx = b_j.x - b_i.x;
                let dy = b_j.y - b_i.y;
                let r2 = dx * dx + dy * dy;
                if r2 < 1e-30 {
                    continue;
                }
                let r = r2.sqrt();
                let inv_r = r.recip();
                let rhat_x = dx * inv_r; // points receiver → source
                let rhat_y = dy * inv_r;

                let gm_over_r = b_j.mass * inv_r; // G = 1
                let pref = b_j.mass / (c2 * r2); //  G m_j / (c² r²)

                let rhat_dot_v = rhat_x * vx_i + rhat_y * vy_i;
                let scalar_rhat = 4.0 * gm_over_r - v2_i; // (4GM/r − v²)
                let scalar_v = 4.0 * rhat_dot_v; // 4 (r̂·v)

                // Minus sign: see block comment above — our r̂ is −n̂.
                ax -= pref * (scalar_rhat * rhat_x + scalar_v * vx_i);
                ay -= pref * (scalar_rhat * rhat_y + scalar_v * vy_i);
            }

            scratch_acc[i].0 += ax;
            scratch_acc[i].1 += ay;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sanity-check the computed constant against an independently-written
    /// runtime derivation. Guards against someone rewriting the const block
    /// with a bogus formula.
    #[test]
    fn c_solar_units_matches_si_derivation() {
        let c_si: f64 = 299_792_458.0;
        let au_si: f64 = 149_597_870_700.0;
        let year_s: f64 = 365.25 * 86_400.0;
        let derived = c_si * (year_s / (2.0 * std::f64::consts::PI)) / au_si;
        assert!(
            (derived - C_SOLAR_UNITS).abs() < 1e-9,
            "C_SOLAR_UNITS drift: derived {derived}, const {C_SOLAR_UNITS}"
        );
        // Rough order-of-magnitude check — c is ~10⁴ in these units.
        assert!(
            (9000.0..11000.0).contains(&C_SOLAR_UNITS),
            "C_SOLAR_UNITS out of expected range: {C_SOLAR_UNITS}"
        );
    }

    /// Acceleration on a lone body must be zero — no self-interaction.
    #[test]
    fn isolated_body_feels_no_pn_force() {
        let bodies = vec![Body::star(1.0)];
        let mut acc = vec![(0.0, 0.0)];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);
        assert_eq!(acc[0], (0.0, 0.0));
    }

    /// 1PN perturbation is additive: register twice and the effect doubles.
    /// Guards the "must add, not overwrite" clause of the trait contract.
    #[test]
    fn accumulate_is_additive() {
        let bodies = vec![
            Body::star(1.0),
            Body::rocky(1e-6).at(0.5, 0.0).with_velocity(0.0, 1.414),
        ];
        let pn = PostNewtonian1PN::solar_units();

        let mut once = vec![(0.0, 0.0); 2];
        pn.accumulate(&bodies, &mut once);

        let mut twice = vec![(0.0, 0.0); 2];
        pn.accumulate(&bodies, &mut twice);
        pn.accumulate(&bodies, &mut twice);

        for i in 0..2 {
            assert!(
                (twice[i].0 - 2.0 * once[i].0).abs() < 1e-14,
                "body {i} ax not additive"
            );
            assert!(
                (twice[i].1 - 2.0 * once[i].1).abs() < 1e-14,
                "body {i} ay not additive"
            );
        }
    }

    /// At the speed-of-light limit c → ∞, the 1PN correction vanishes.
    /// This is the physical sanity check that the prefactor scales as 1/c².
    #[test]
    fn infinite_c_gives_zero_correction() {
        let bodies = vec![
            Body::star(1.0),
            Body::rocky(1e-6).at(0.387, 0.0).with_velocity(0.0, 2.07),
        ];
        let pn = PostNewtonian1PN::with_c(1e20);
        let mut acc = vec![(0.0, 0.0); 2];
        pn.accumulate(&bodies, &mut acc);
        for (ax, ay) in acc {
            assert!(ax.abs() < 1e-30, "ax={ax} should be ~0 at c→∞");
            assert!(ay.abs() < 1e-30, "ay={ay} should be ~0 at c→∞");
        }
    }

    /// At perihelion with purely tangential motion (n̂·v = 0), the 1PN
    /// correction must be radially **outward** — away from the Sun.
    /// An inward-pointing correction would flip the sign of the perihelion
    /// precession (making it retrograde) and was the bug caught on the
    /// first end-to-end Mercury run. Guards the sign convention across
    /// future refactors.
    #[test]
    fn mercury_perihelion_pn_points_outward() {
        const A_MERCURY: f64 = 0.387_098;
        const E_MERCURY: f64 = 0.205_63;
        const M_MERCURY: f64 = 1.660_114e-7;

        let sun = Body::star(1.0);
        let r_peri = A_MERCURY * (1.0 - E_MERCURY);
        let v_peri = (2.0 / r_peri - 1.0 / A_MERCURY).sqrt();
        let mercury = Body::rocky(M_MERCURY)
            .at(r_peri, 0.0)
            .with_velocity(0.0, v_peri);

        let bodies = vec![sun, mercury];
        let mut acc = vec![(0.0, 0.0); 2];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);

        // Mercury is at (+r_peri, 0); "outward from Sun" = +x direction.
        assert!(
            acc[1].0 > 0.0,
            "1PN correction on Mercury at perihelion must point outward (+x), got ax = {}",
            acc[1].0,
        );
        // Tangential velocity is in +y; with n̂·v = 0 the v-term vanishes,
        // so the correction is purely radial.
        assert!(
            acc[1].1.abs() < 1e-20,
            "ay should be ~0 at purely-tangential perihelion, got {}",
            acc[1].1,
        );
    }

    /// Order-of-magnitude check: Mercury at perihelion sees a 1PN
    /// acceleration roughly (v/c)² smaller than Newtonian. With
    /// v_peri ≈ 2.0 and c ≈ 10 065, β² ≈ 4 × 10⁻⁸, so the ratio of
    /// magnitudes must sit in the 10⁻⁸–10⁻⁶ band.
    #[test]
    fn mercury_pn_has_expected_magnitude() {
        const A_MERCURY: f64 = 0.387_098;
        const E_MERCURY: f64 = 0.205_63;
        const M_MERCURY: f64 = 1.660_114e-7;

        let sun = Body::star(1.0);
        let r_peri = A_MERCURY * (1.0 - E_MERCURY);
        let v_peri = (2.0 / r_peri - 1.0 / A_MERCURY).sqrt(); // G M_sun = 1
        let mercury = Body::rocky(M_MERCURY)
            .at(r_peri, 0.0)
            .with_velocity(0.0, v_peri);

        let bodies = vec![sun, mercury];
        let mut acc = vec![(0.0, 0.0); 2];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);

        let a_pn = (acc[1].0.powi(2) + acc[1].1.powi(2)).sqrt();
        let a_newt = 1.0 / (r_peri * r_peri); // G M / r²
        let ratio = a_pn / a_newt;

        assert!(
            ratio > 1e-8 && ratio < 1e-6,
            "Mercury a_pn/a_newt = {ratio:.3e}, expected band 10⁻⁸–10⁻⁶"
        );
    }
}
