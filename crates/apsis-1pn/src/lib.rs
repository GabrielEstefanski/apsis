//! Out-of-tree perturbation crate for `apsis`.
//!
//! **This crate proves that the perturbation extension contract is
//! buildable, not just documented.** It compiles against the public API
//! alone — no `pub(crate)` access, no patches to core sources, no
//! dependency other than `apsis` itself. A future change to that API
//! that breaks this crate fails CI loudly rather than quietly.
//!
//! Treat this crate as the **template** when writing new perturbation
//! crates (radiation pressure, J2 oblateness, drag, …). See the
//! crate-level README for the full extension-contract specification
//! and the rationale.
//!
//! # ⚠ Critical precondition
//!
//! Attaching 1PN to a softened-gravity system **invalidates the physical
//! model**. For Mercury-like orbits, the numerical apsidal precession
//! induced by Plummer softening alone is ~2 × 10³ larger than the
//! relativistic signal *and inverts its sign*. Energy and angular
//! momentum stay conserved at machine precision while the trajectory
//! is physically wrong.
//!
//! **This is not a numerical error — it is a model violation.**
//!
//! Call [`Body::unsoftened`](apsis::domain::body::Body::unsoftened) on
//! every body or
//! [`System::with_exact_gravity`](apsis::core::system::System::with_exact_gravity)
//! system-wide. The contract is enforced once, in the core: a violation
//! emits a structured warning at registration naming the failed invariant.
//!
//! # Validation signal
//!
//! With the contract enforced, this crate reproduces Mercury's textbook
//! 43 arcsec/century rate to **4.4 ppm**, gated in CI under `mercury-gate`.
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
//! use apsis::core::system::System;
//! use apsis::domain::body::Body;
//! use apsis::physics::integrator::IntegratorKind;
//! use apsis_1pn::PostNewtonian1PN;
//!
//! let sun     = Body::star(1.0);
//! let mercury = Body::rocky(1.66e-7).at(0.307, 0.0).with_velocity(0.0, 2.078);
//! let mut sys = System::new(vec![sun, mercury], UnitSystem::canonical())
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-4);
//!
//! sys.add_perturbation(Box::new(PostNewtonian1PN::solar_units()));
//!
//! sys.integrate_for(200.0 * std::f64::consts::PI);
//! println!("dE/E = {:.3e}", sys.energy_delta());
//! ```

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::PerturbationForce;

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
/// Register via [`System::add_perturbation`](apsis::core::system::System::add_perturbation):
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
    /// The 1PN correction is derived by expanding the geodesic equation
    /// around the Newtonian Hamiltonian `H_N = p²/2m − GMm/r`. The
    /// expansion therefore requires:
    ///
    /// - **Exactness::Exact** — the unperturbed base must be the bit-exact
    ///   1/r potential. Plummer softening would substitute a different
    ///   unperturbed system whose apsidal precession alone is ~2 × 10³
    ///   larger than the 1PN signal for a Mercury-like orbit, swamping
    ///   the physical effect and inverting its sign.
    /// - **Continuity::Smooth** — symplectic integration of the correction
    ///   relies on a smooth Hamiltonian flow; force discontinuities cannot
    ///   be represented within any symplectic splitting scheme regardless
    ///   of integrator order.
    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::exact_and_smooth()
    }

    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [Vec3]) {
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
        // (`rhat`) points FROM RECEIVER TO SOURCE, i.e. r̂ = −n̂. Substituting
        // r̂ = −n̂ flips the sign of the whole bracket (both n̂ and n̂·v change
        // sign; the outer product n̂·v × v inherits one flip), so in this
        // crate's convention the expression is applied with an overall
        // minus sign.
        for i in 0..bodies.len() {
            let b_i = &bodies[i];
            let v_i = Vec3::new(b_i.vx, b_i.vy, b_i.vz);
            let v2_i = v_i.length_squared();

            let mut a = Vec3::ZERO;

            for j in 0..bodies.len() {
                if i == j {
                    continue;
                }
                let b_j = &bodies[j];
                let dx = b_j.x - b_i.x;
                let dy = b_j.y - b_i.y;
                let dz = b_j.z - b_i.z;
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < 1e-30 {
                    continue;
                }
                let r = r2.sqrt();
                let inv_r = r.recip();
                let rhat = Vec3::new(dx * inv_r, dy * inv_r, dz * inv_r); // receiver → source

                let gm_over_r = b_j.mass * inv_r; // G = 1
                let pref = b_j.mass / (c2 * r2); //  G m_j / (c² r²)

                let rhat_dot_v = rhat.dot(v_i);
                let scalar_rhat = 4.0 * gm_over_r - v2_i; // (4GM/r − v²)
                let scalar_v = 4.0 * rhat_dot_v; // 4 (r̂·v)

                // Minus sign: see block comment above — our r̂ is −n̂.
                //
                // Per-axis scalar form `pref · (scalar_rhat · r̂ +
                // scalar_v · v)` is load-bearing: re-associating into
                // Vec3 ops shifts ULPs and moves the Mercury 1PN gate
                // floor by ~10² ppm.
                a.x -= pref * (scalar_rhat * rhat.x + scalar_v * v_i.x);
                a.y -= pref * (scalar_rhat * rhat.y + scalar_v * v_i.y);
                a.z -= pref * (scalar_rhat * rhat.z + scalar_v * v_i.z);
            }

            scratch_acc[i] += a;
        }
    }
}

/// Federation entry point — a [`PerturbationDescriptor`] that consumers
/// register without ever naming [`PostNewtonian1PN`] directly.
///
/// The descriptor delegates `kernel_requirements` to the produced
/// perturbation; metadata (`name`, `description`) is the single
/// authoritative source for any UI surface that lists this plugin.
pub struct Descriptor;

impl apsis::physics::integrator::PerturbationDescriptor for Descriptor {
    fn name(&self) -> &str {
        "General Relativity (1PN)"
    }

    fn description(&self) -> &str {
        "Schwarzschild perihelion advance — Mercury 43 arcsec/century"
    }

    fn kernel_requirements(&self) -> KernelRequirements {
        PostNewtonian1PN::solar_units().kernel_requirements()
    }

    fn build(&self) -> Box<dyn apsis::physics::integrator::PerturbationForce> {
        Box::new(PostNewtonian1PN::solar_units())
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
        let mut acc = vec![Vec3::ZERO];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);
        assert_eq!(acc[0], Vec3::ZERO);
    }

    /// 1PN perturbation is additive: register twice and the effect doubles.
    /// Guards the "must add, not overwrite" clause of the trait contract.
    #[test]
    fn accumulate_is_additive() {
        let bodies =
            vec![Body::star(1.0), Body::rocky(1e-6).at(0.5, 0.0).with_velocity(0.0, 1.414)];
        let pn = PostNewtonian1PN::solar_units();

        let mut once = vec![Vec3::ZERO; 2];
        pn.accumulate(&bodies, &mut once);

        let mut twice = vec![Vec3::ZERO; 2];
        pn.accumulate(&bodies, &mut twice);
        pn.accumulate(&bodies, &mut twice);

        for i in 0..2 {
            assert!((twice[i].x - 2.0 * once[i].x).abs() < 1e-14, "body {i} ax not additive");
            assert!((twice[i].y - 2.0 * once[i].y).abs() < 1e-14, "body {i} ay not additive");
            assert!((twice[i].z - 2.0 * once[i].z).abs() < 1e-14, "body {i} az not additive");
        }
    }

    /// At the speed-of-light limit c → ∞, the 1PN correction vanishes.
    /// This is the physical sanity check that the prefactor scales as 1/c².
    #[test]
    fn infinite_c_gives_zero_correction() {
        let bodies =
            vec![Body::star(1.0), Body::rocky(1e-6).at(0.387, 0.0).with_velocity(0.0, 2.07)];
        let pn = PostNewtonian1PN::with_c(1e20);
        let mut acc = vec![Vec3::ZERO; 2];
        pn.accumulate(&bodies, &mut acc);
        for a in acc {
            assert!(a.x.abs() < 1e-30, "ax={} should be ~0 at c→∞", a.x);
            assert!(a.y.abs() < 1e-30, "ay={} should be ~0 at c→∞", a.y);
            assert!(a.z.abs() < 1e-30, "az={} should be ~0 at c→∞", a.z);
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
        let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

        let bodies = vec![sun, mercury];
        let mut acc = vec![Vec3::ZERO; 2];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);

        // Mercury is at (+r_peri, 0); "outward from Sun" = +x direction.
        assert!(
            acc[1].x > 0.0,
            "1PN correction on Mercury at perihelion must point outward (+x), got ax = {}",
            acc[1].x,
        );
        // Tangential velocity is in +y; with n̂·v = 0 the v-term vanishes,
        // so the correction is purely radial.
        assert!(
            acc[1].y.abs() < 1e-20,
            "ay should be ~0 at purely-tangential perihelion, got {}",
            acc[1].y,
        );
        assert!(
            acc[1].z.abs() < 1e-20,
            "az should be ~0 in the orbital plane (z = vz = 0), got {}",
            acc[1].z,
        );
    }

    /// The 1PN perturbation must declare that it requires the exact 1/r
    /// Newtonian base (Exactness) and a smooth kernel (Continuity). This
    /// is what the kernel-requirement check inside
    /// `System::add_perturbation` matches against the active kernel's
    /// properties.
    #[test]
    fn kernel_requirements_are_exact_and_smooth() {
        use apsis::physics::gravity::kernel::{Continuity, Exactness, KernelRequirements};
        let req = PostNewtonian1PN::solar_units().kernel_requirements();
        assert_eq!(req, KernelRequirements::exact_and_smooth());
        assert_eq!(req.required_exactness, Some(Exactness::Exact));
        assert_eq!(req.min_continuity, Some(Continuity::Smooth));
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
        let mercury = Body::rocky(M_MERCURY).at(r_peri, 0.0).with_velocity(0.0, v_peri);

        let bodies = vec![sun, mercury];
        let mut acc = vec![Vec3::ZERO; 2];
        PostNewtonian1PN::solar_units().accumulate(&bodies, &mut acc);

        let a_pn = acc[1].length();
        let a_newt = 1.0 / (r_peri * r_peri); // G M / r²
        let ratio = a_pn / a_newt;

        assert!(
            ratio > 1e-8 && ratio < 1e-6,
            "Mercury a_pn/a_newt = {ratio:.3e}, expected band 10⁻⁸–10⁻⁶"
        );
    }
}
