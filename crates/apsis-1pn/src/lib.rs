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
//! sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::solar_units()));
//!
//! sys.integrate_for(200.0 * std::f64::consts::PI);
//! println!("dE/E = {:.3e}", sys.energy_delta());
//! ```

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::{HamiltonianOperator, Operator, ParameterValidationError};
use apsis::units::UnitSystem;

/// Speed of light in m/s — CODATA exact by SI definition.
const C_SI: f64 = 299_792_458.0;

/// Relative tolerance for `from_raw_c_validated`. A bit looser than the
/// f64 round-off floor because the validator derives `c_expected` via
/// floating-point conversion from SI, which itself carries ~3 ULP of
/// noise. 1e-9 catches unit-mismatch errors (typical relative error
/// 1e-2 to 1e+8) without flagging legitimate truncation.
const C_VALIDATION_TOLERANCE: f64 = 1.0e-9;

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
    const AU_SI: f64 = 149_597_870_700.0; // m, IAU 2012
    const YEAR_S: f64 = 365.25 * 86_400.0;
    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
    C_SI * (YEAR_S / TWO_PI) / AU_SI
};

/// First post-Newtonian gravitational correction (Schwarzschild,
/// test-particle form applied pairwise).
///
/// Register via
/// [`System::add_hamiltonian_perturbation`](apsis::core::system::System::add_hamiltonian_perturbation):
///
/// ```ignore
/// sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::solar_units()));
/// ```
///
/// Stateless. Safe to share across threads.
///
/// # Cross-reference to REBOUNDx
///
/// This implementation corresponds to REBOUNDx's `gr` effect (Anderson
/// et al. 1975 test-particle 1PN, velocity-dependent in `(r̂·v)·v` and
/// `v² · r̂` terms). It is **not** the same as REBOUNDx's
/// `gr_potential`, which is a velocity-independent effective potential
/// (Nobili & Roxburgh 1986) that gets pericenter precession right but
/// the mean motion wrong by `O(GM/(a·c²))`. The Nobili–Roxburgh form
/// would be a separate operator with closed-form `potential` and
/// WHFast-symplectic-friendly dispatch; it is not implemented here.
///
/// The full N-body Einstein–Infeld–Hoffmann Hamiltonian — the rigorous
/// form when masses are comparable — is also out of scope. For the
/// Sun–Mercury validation regime (`m_Mercury / m_Sun ≈ 2 × 10⁻⁷`) the
/// test-particle simplification is canonical and gates the 4.4 ppm
/// agreement reported in [`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`].
#[derive(Debug, Clone, Copy)]
pub struct PostNewtonian1PN {
    /// Speed of light in the caller's unit system.
    c: f64,
}

impl PostNewtonian1PN {
    // ── Named-regime constructors (Pattern A) ─────────────────────────────────

    /// Construct for the apsis-1pn historical solar convention
    /// ([`C_SOLAR_UNITS`]): G = 1, length = AU, mass = M☉, time =
    /// year/(2π). This is the unit system the Mercury perihelion gate
    /// and the rest of the 1PN validation portfolio run in.
    ///
    /// **Distinct from [`apsis::units::UnitSystem::solar`]**, which
    /// uses the IAU convention (G ≈ 4π², length = AU, mass = M☉,
    /// time = 1 year). Use [`for_units`](Self::for_units) with
    /// [`UnitSystem::solar`] for the IAU convention; use this
    /// constructor for the G = 1 convention.
    pub const fn solar_units() -> Self {
        Self { c: C_SOLAR_UNITS }
    }

    /// Construct for an arbitrary [`UnitSystem`], deriving `c` from the
    /// SI definition of the speed of light in the chosen L/T scaling.
    /// The recommended path for any non-canonical unit choice — the
    /// user picks the units, `c` is computed exactly so the relativistic
    /// correction stays consistent with the rest of the integration.
    pub fn for_units(units: UnitSystem) -> Self {
        Self { c: C_SI * units.time_scale_si() / units.length_scale_si() }
    }

    // ── Raw escape (with optional validation) ─────────────────────────────────

    /// Construct from a raw `c` value in the caller's unit system. No
    /// validation that the value is consistent with any known unit
    /// system — use [`from_raw_c_validated`](Self::from_raw_c_validated)
    /// for cross-checked construction, or [`for_units`](Self::for_units)
    /// to skip raw input entirely.
    pub const fn from_raw_c(c: f64) -> Self {
        Self { c }
    }

    /// Construct from a raw `c` value, cross-checking against the `c`
    /// derived from the given [`UnitSystem`]. Returns
    /// [`ParameterValidationError`] when the relative error between the
    /// supplied and derived values exceeds `1e-9`. Typical use case:
    /// reading `c` from an external config or another simulator and
    /// confirming it matches before attaching the perturbation.
    pub fn from_raw_c_validated(
        c: f64,
        units: UnitSystem,
    ) -> Result<Self, ParameterValidationError> {
        let expected = C_SI * units.time_scale_si() / units.length_scale_si();
        let rel_err = (c - expected).abs() / expected.abs();
        if rel_err > C_VALIDATION_TOLERANCE {
            return Err(ParameterValidationError {
                operator: "PostNewtonian1PN",
                parameter: "c",
                got: c,
                expected,
                tolerance: C_VALIDATION_TOLERANCE,
                message: format!(
                    "speed of light {c:.6e} disagrees with the value derived from the \
                     supplied UnitSystem ({expected:.6e}); check that the UnitSystem L/T \
                     scaling matches the unit system the `c` value was originally measured in"
                ),
            });
        }
        Ok(Self { c })
    }

    /// Speed of light this instance was configured with.
    pub const fn c(&self) -> f64 {
        self.c
    }
}

impl Operator for PostNewtonian1PN {
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
}

impl HamiltonianOperator for PostNewtonian1PN {
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]) {
        debug_assert_eq!(
            bodies.len(),
            acc.len(),
            "HamiltonianOperator contract: acc must be sized to bodies"
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
            let v_i = Vec3::new(b_i.vel_x, b_i.vel_y, b_i.vel_z);
            let v2_i = v_i.length_squared();

            let mut a = Vec3::ZERO;

            for j in 0..bodies.len() {
                if i == j {
                    continue;
                }
                let b_j = &bodies[j];
                let dx = b_j.pos_x - b_i.pos_x;
                let dy = b_j.pos_y - b_i.pos_y;
                let dz = b_j.pos_z - b_i.pos_z;
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

            acc[i] += a;
        }
    }

    // `potential` inherits the default [`Potential::NotAvailable`].
    //
    // Test-particle pairwise 1PN is derived as a force expansion of the
    // geodesic equation, not from a two-body scalar Hamiltonian. The
    // rigorous form is the full Einstein–Infeld–Hoffmann N-body
    // Hamiltonian — out of scope for this test-particle approximation
    // crate. `System::total_energy` therefore excludes this operator's
    // contribution; `System::conservation_report` surfaces the exclusion
    // by classifying the system as `HamiltonianForceOnly` whenever 1PN
    // is the only registered Hamiltonian operator.
}

/// Federation entry point — a [`HamiltonianOperatorDescriptor`] that
/// consumers register without ever naming [`PostNewtonian1PN`] directly.
///
/// The descriptor delegates `kernel_requirements` to the produced
/// operator; metadata (`name`, `description`) is the single authoritative
/// source for any UI surface that lists this plugin.
///
/// [`HamiltonianOperatorDescriptor`]: apsis::physics::integrator::HamiltonianOperatorDescriptor
pub struct Descriptor;

impl apsis::physics::integrator::HamiltonianOperatorDescriptor for Descriptor {
    fn name(&self) -> &str {
        "General Relativity (1PN)"
    }

    fn description(&self) -> &str {
        "Schwarzschild perihelion advance — Mercury 43 arcsec/century"
    }

    fn kernel_requirements(&self) -> KernelRequirements {
        <PostNewtonian1PN as Operator>::kernel_requirements(&PostNewtonian1PN::solar_units())
    }

    fn build(&self) -> Box<dyn HamiltonianOperator> {
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

    /// In Hénon canonical units (L = 1 m, T = 1 s, mass adjusted for
    /// G = 1) the speed of light numerically equals its SI value. The
    /// derivation chain is `c = c_SI · T_scale / L_scale`; both scales
    /// are 1 here. Cleanest verification of the derivation.
    #[test]
    fn for_units_canonical_gives_si_c() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::canonical());
        assert!(
            (pn.c() - 299_792_458.0).abs() < 1.0,
            "canonical units expect c_SI = 299_792_458, got {}",
            pn.c(),
        );
    }

    /// `for_units(UnitSystem::solar())` returns `c` in the IAU solar
    /// convention (L = 1 AU, T = 1 year, M = 1 M☉ → G ≈ 4π²). This
    /// is **distinct** from `PostNewtonian1PN::solar_units()`, which
    /// pins `c` to [`C_SOLAR_UNITS`] in the apsis-1pn historical
    /// convention (G = 1, L = AU, T = yr/(2π), M = M☉). The two
    /// "solar" names refer to different unit systems by design — IAU
    /// is what [`apsis::units::UnitSystem::solar`] exposes, while
    /// [`C_SOLAR_UNITS`] is the G=1 reference baked into apsis-1pn.
    #[test]
    fn for_units_solar_uses_iau_convention() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar());
        let c_si: f64 = 299_792_458.0;
        let year_s: f64 = 365.25 * 86_400.0;
        let au_m: f64 = 1.495_978_707e11;
        let expected = c_si * year_s / au_m;
        assert!(
            (pn.c() - expected).abs() / expected < 1e-12,
            "for_units(solar) c={} disagrees with derived IAU c={}",
            pn.c(),
            expected,
        );
    }

    /// `solar_units()` returns `C_SOLAR_UNITS` unchanged — locks the
    /// public constant against accidental drift inside the constructor.
    #[test]
    fn solar_units_returns_c_solar_units() {
        assert_eq!(PostNewtonian1PN::solar_units().c(), C_SOLAR_UNITS);
    }

    /// Raw escape constructor accepts any numeric value with no
    /// validation. The whole point of the escape is to make
    /// unvalidated construction syntactically obvious through the
    /// `from_raw_` prefix.
    #[test]
    fn from_raw_c_accepts_arbitrary_value() {
        let pn = PostNewtonian1PN::from_raw_c(1.234e5);
        assert_eq!(pn.c(), 1.234e5);
    }

    /// Validated escape accepts the `c` derived from the same
    /// `UnitSystem` it validates against — self-consistency of the
    /// derivation chain `for_units → from_raw_c_validated`. Mirrors
    /// the case where a user reads `c` from an external source and
    /// wants the simulator to confirm it.
    #[test]
    fn from_raw_c_validated_accepts_matching_c() {
        let units = UnitSystem::canonical();
        let c_derived = PostNewtonian1PN::for_units(units).c();
        let pn = PostNewtonian1PN::from_raw_c_validated(c_derived, units)
            .expect("c derived from for_units must validate against the same UnitSystem");
        assert_eq!(pn.c(), c_derived);
    }

    /// Validated escape rejects a `c` value inconsistent with the
    /// supplied unit system — the protection against unit-mismatch
    /// errors that would silently distort the relativistic correction.
    /// Passing `c = 1.0` (geometric-units value) against
    /// `UnitSystem::canonical()` (which expects `c_SI ≈ 3e8`) must
    /// fail with `parameter == "c"` and `expected ≈ c_SI`.
    #[test]
    fn from_raw_c_validated_rejects_mismatched_c() {
        let err = PostNewtonian1PN::from_raw_c_validated(1.0, UnitSystem::canonical())
            .expect_err("c=1 against canonical units must not validate");
        assert_eq!(err.operator, "PostNewtonian1PN");
        assert_eq!(err.parameter, "c");
        assert_eq!(err.got, 1.0);
        assert!(
            (err.expected - 299_792_458.0).abs() < 1.0,
            "expected ≈ c_SI, got {}",
            err.expected,
        );
    }

    /// Acceleration on a lone body must be zero — no self-interaction.
    #[test]
    fn isolated_body_feels_no_pn_force() {
        let bodies = vec![Body::star(1.0)];
        let mut acc = vec![Vec3::ZERO];
        PostNewtonian1PN::solar_units().accumulate_force(&bodies, &mut acc);
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
        pn.accumulate_force(&bodies, &mut once);

        let mut twice = vec![Vec3::ZERO; 2];
        pn.accumulate_force(&bodies, &mut twice);
        pn.accumulate_force(&bodies, &mut twice);

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
        let pn = PostNewtonian1PN::from_raw_c(1e20);
        let mut acc = vec![Vec3::ZERO; 2];
        pn.accumulate_force(&bodies, &mut acc);
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
        PostNewtonian1PN::solar_units().accumulate_force(&bodies, &mut acc);

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
    /// `System::add_hamiltonian_perturbation` matches against the active
    /// kernel's properties.
    #[test]
    fn kernel_requirements_are_exact_and_smooth() {
        use apsis::physics::gravity::kernel::{Continuity, Exactness, KernelRequirements};
        let req =
            <PostNewtonian1PN as Operator>::kernel_requirements(&PostNewtonian1PN::solar_units());
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
        PostNewtonian1PN::solar_units().accumulate_force(&bodies, &mut acc);

        let a_pn = acc[1].length();
        let a_newt = 1.0 / (r_peri * r_peri); // G M / r²
        let ratio = a_pn / a_newt;

        assert!(
            ratio > 1e-8 && ratio < 1e-6,
            "Mercury a_pn/a_newt = {ratio:.3e}, expected band 10⁻⁸–10⁻⁶"
        );
    }
}
