//! First post-Newtonian gravitational correction (Schwarzschild,
//! test-particle form applied pairwise) per Anderson et al. (1975,
//! *ApJ* 200, 221). Recovers the GR perihelion precession
//! `Δφ = 6π G M / (c² a (1 − e²))` at leading order. The CI gate
//! reproduces Mercury's 43 arcsec/century to within 100 ppm of GR.
//!
//! Force expression on receiver `i` from source `j ≠ i`:
//!
//! ```text
//!   a_1PN(i ← j) = G m_j / (c² r²) · [ (4 G m_j / r − v_i²) · r̂ + 4 (r̂ · v_i) v_i ]
//! ```
//!
//! # Critical precondition
//!
//! Attaching 1PN to a softened kernel (`NewtonKernel::new(ε > 0)`)
//! invalidates the physical model: numerical apsidal precession from a
//! Plummer-style 1/√(r²+ε²) potential is ~5 × 10⁴ larger than the
//! relativistic signal and of opposite sign at Mercury's orbit. A
//! full-potential apsidal-angle quadrature puts it at ϖ̇ ≈ −2.29 × 10⁶
//! arcsec/century (the leading-order −3 n ε² / [2 a² (1 − e²)²] closed
//! form is +2.7 % above that); see paper §3.2. The default
//! `NewtonKernel::exact()` is silent against the kernel-requirement
//! check; opting into ε > 0 emits a structured warning.
//!
//! # Use
//!
//! ```ignore
//! let units   = UnitSystem::solar_canonical();
//! let sun     = Body::star(1.0);
//! let mercury = Body::rocky(1.66e-7).at(0.307, 0.0).with_velocity(0.0, 2.078);
//! let mut sys = System::new(vec![sun, mercury], units)
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-4);
//! sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(units)))?;
//! ```
//!
//! # Reference
//!
//! Anderson, J. D., Esposito, P. B., Martin, W., Thornton, C. L., &
//! Muhleman, D. O. (1975). DOI:
//! [10.1086/153180](https://doi.org/10.1086/153180).

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::regime::{classify_mass_ratio, mass_ratio};
use apsis::physics::integrator::{
    Citation, HamiltonianOperator, Operator, RegimeViolation, Severity,
};
use apsis::units::UnitSystem;

/// Mass ratio at which the test-particle approximation starts losing
/// accuracy. Calibrated so Sun–Jupiter (≈ 9.5 × 10⁻⁴) sits inside.
const PN1_MASS_RATIO_WARN: f64 = 1.0e-2;

/// Mass ratio above which the test-particle approximation no longer
/// applies; the full EIH N-body Hamiltonian is the rigorous form.
const PN1_MASS_RATIO_HARD: f64 = 1.0e-1;

/// Speed of light, m/s.
const C_SI: f64 = 299_792_458.0;

/// Speed of light in solar canonical units (AU per year/(2π), G = 1),
/// derived at compile time from `c_SI · (year_s / 2π) / AU_SI`.
/// Current value ≈ `10_065.130`.
pub const C_SOLAR_UNITS: f64 = {
    const AU_SI: f64 = 149_597_870_700.0; // m, IAU 2012
    const YEAR_S: f64 = 365.25 * 86_400.0;
    const TWO_PI: f64 = 2.0 * std::f64::consts::PI;
    C_SI * (YEAR_S / TWO_PI) / AU_SI
};

/// First post-Newtonian Schwarzschild correction, test-particle form
/// applied pairwise. Carries the [`UnitSystem`] it was constructed for;
/// `System` registration returns `Err(UnitSystemMismatch)` if it
/// disagrees with the system's own unit system.
#[derive(Debug, Clone, Copy)]
pub struct PostNewtonian1PN {
    /// Speed of light in `units`.
    c: f64,
    units: UnitSystem,
}

impl PostNewtonian1PN {
    /// Construct for the supplied [`UnitSystem`], deriving `c` from
    /// `c_SI · T_scale / L_scale`. Recommended constructor — pass the
    /// same units used to build the `System`. For IAU solar (G ≈ 4π²)
    /// use [`UnitSystem::solar`]; for canonical (G = 1) use
    /// [`UnitSystem::solar_canonical`].
    pub fn for_units(units: UnitSystem) -> Self {
        Self { c: C_SI * units.time_scale_si() / units.length_scale_si(), units }
    }

    /// Construct from an explicit `c` value pinned to `units`. Use when
    /// `c` is computed externally; prefer [`for_units`](Self::for_units)
    /// otherwise. `c` is unchecked; the `units` cross-check at
    /// registration still applies.
    pub const fn from_raw_c(c: f64, units: UnitSystem) -> Self {
        Self { c, units }
    }

    pub const fn c(&self) -> f64 {
        self.c
    }

    pub const fn units(&self) -> UnitSystem {
        self.units
    }
}

impl Operator for PostNewtonian1PN {
    fn declared_units(&self) -> Option<UnitSystem> {
        Some(self.units)
    }

    /// 1PN expands around the bit-exact 1/r Hamiltonian and requires a
    /// smooth Hamiltonian flow for symplectic integration. Softening or
    /// force discontinuities break the derivation.
    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::exact_and_smooth()
    }

    /// Test-particle pairwise 1PN assumes `m_secondary ≪ m_primary` for
    /// every secondary; `bodies[0]` is the primary by convention.
    fn check_regime(&self, bodies: &[Body], _t: f64) -> Vec<RegimeViolation> {
        let mut violations = Vec::new();
        if bodies.len() < 2 {
            return violations;
        }
        for i in 1..bodies.len() {
            let Some(ratio) = mass_ratio(bodies, 0, i) else {
                continue;
            };
            let Some((severity, threshold)) =
                classify_mass_ratio(ratio, PN1_MASS_RATIO_WARN, PN1_MASS_RATIO_HARD)
            else {
                continue;
            };
            // One bound key per severity tier: the dedup state in
            // `System` would otherwise suppress an Approaching → Hard
            // escalation if both fire across two cadences.
            let bound = match severity {
                Severity::Approaching => "max_secondary_to_primary_mass_ratio.approaching",
                Severity::Exceeded => "max_secondary_to_primary_mass_ratio.exceeded",
                Severity::Hard => "max_secondary_to_primary_mass_ratio.hard",
                // `Severity` is non_exhaustive; future variants degrade
                // safely into the generic key rather than silently
                // suppressing the warning.
                _ => "max_secondary_to_primary_mass_ratio.unknown_severity",
            };
            violations.push(RegimeViolation {
                operator: self.name(),
                bound,
                value: ratio,
                threshold,
                severity,
                body_index: Some(i),
                message: "test-particle pairwise 1PN derivation assumes \
                          m_secondary / m_primary ≪ 1; the full Einstein–\
                          Infeld–Hoffmann N-body Hamiltonian is the rigorous \
                          form for comparable masses",
            });
        }
        violations
    }

    /// Mass ratio is static; one check per ~Mercury-orbit worth of
    /// steps is sufficient.
    fn regime_check_cadence(&self) -> usize {
        15_000
    }

    /// Anderson et al. (1975) for the Schwarzschild test-particle form;
    /// Will (1993) for the EIH background.
    fn citation(&self) -> Option<Citation> {
        Some(Citation {
            bibtex: PN1_BIBTEX,
            doi: Some("10.1086/153180"),
            crate_name: env!("CARGO_PKG_NAME"),
            crate_version: env!("CARGO_PKG_VERSION"),
            commit_hash: option_env!("APSIS_1PN_GIT_COMMIT").filter(|s| !s.is_empty()),
            description: Some("First-post-Newtonian Schwarzschild correction"),
            url: Some("https://github.com/GabrielEstefanski/apsis"),
            author: Some("Estefanski, G. B."),
        })
    }
}

const PN1_BIBTEX: &str = r#"@article{anderson1975,
  author  = {Anderson, J. D. and Esposito, P. B. and Martin, W. and Thornton, C. L. and Muhleman, D. O.},
  title   = {Experimental test of general relativity using time-delay data from Mariner 6 and Mariner 7},
  journal = {Astrophysical Journal},
  volume  = {200},
  pages   = {221--233},
  year    = {1975},
  doi     = {10.1086/153180}
}
@book{will1993,
  author    = {Will, C. M.},
  title     = {Theory and Experiment in Gravitational Physics},
  publisher = {Cambridge University Press},
  year      = {1993},
  edition   = {2}
}"#;

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

/// [`HamiltonianOperatorDescriptor`] for plugin registries — produces
/// a `PostNewtonian1PN` for the supplied [`UnitSystem`].
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
        // Kernel requirements are unit-system-independent; pick any
        // UnitSystem to satisfy the constructor signature.
        <PostNewtonian1PN as Operator>::kernel_requirements(&PostNewtonian1PN::for_units(
            UnitSystem::solar_canonical(),
        ))
    }

    fn build(&self, units: UnitSystem) -> Box<dyn HamiltonianOperator> {
        Box::new(PostNewtonian1PN::for_units(units))
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
    /// convention (L = 1 AU, T = 1 year, M = 1 M☉ → G ≈ 4π²).
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

    /// `for_units(solar_canonical)` c uses Gaussian time (`sqrt(AU³/GM_sun)`),
    /// the G=1 portfolio baseline; it differs from the IAU `year/(2π)`
    /// literal `C_SOLAR_UNITS` by ~19 ppm (the Gaussian-vs-IAU-year gap).
    #[test]
    fn for_units_solar_canonical_close_to_c_solar_units() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let rel_diff = (pn.c() - C_SOLAR_UNITS).abs() / C_SOLAR_UNITS;
        // ~19 ppm gap; bound guards a GM-primitive regression (the old G·M
        // primitive widened it to ~110 ppm).
        assert!(
            rel_diff < 5e-5,
            "for_units(solar_canonical) c={} vs C_SOLAR_UNITS={} gap {:.3e}, expected ~19 ppm",
            pn.c(),
            C_SOLAR_UNITS,
            rel_diff,
        );
        assert!(
            pn.c() != C_SOLAR_UNITS,
            "for_units should produce Gaussian-time c, not the IAU C_SOLAR_UNITS literal",
        );
    }

    /// `from_raw_c` accepts any value but pins the operator to the
    /// supplied `UnitSystem`. The unit-system binding is what protects
    /// against silent unit-mismatch at registration; the raw `c` value
    /// itself is unchecked at construction.
    #[test]
    fn from_raw_c_accepts_arbitrary_value_with_units() {
        let pn = PostNewtonian1PN::from_raw_c(1.234e5, UnitSystem::canonical());
        assert_eq!(pn.c(), 1.234e5);
        assert_eq!(pn.units(), UnitSystem::canonical());
    }

    /// Test-particle pairwise 1PN regime: Sun + Mercury sits well
    /// inside the regime (m_M / m_S ≈ 1.7e-7); `check_regime` returns
    /// no violations. Locks in that the validation portfolio's
    /// canonical scenario does not spuriously trigger the bound.
    #[test]
    fn check_regime_silent_for_sun_mercury() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let bodies =
            vec![Body::star(1.0), Body::rocky(1.66e-7).at(0.387, 0.0).with_velocity(0.0, 1.61)];
        let violations = pn.check_regime(&bodies, 0.0);
        assert!(
            violations.is_empty(),
            "Sun + Mercury should be inside 1PN regime; got {violations:?}"
        );
    }

    /// Test-particle 1PN regime: equal-mass binary is far outside the
    /// envelope. `check_regime` must return a `Hard` violation
    /// referencing the offending body and the mass-ratio bound.
    #[test]
    fn check_regime_flags_equal_mass_binary_as_hard() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let bodies = vec![Body::star(1.0), Body::star(1.0)];
        let violations = pn.check_regime(&bodies, 0.0);
        assert_eq!(violations.len(), 1, "expected exactly one violation for equal-mass binary");
        let v = &violations[0];
        assert_eq!(v.severity, Severity::Hard);
        assert_eq!(v.body_index, Some(1));
        assert_eq!(v.value, 1.0);
        assert!(v.bound.starts_with("max_secondary_to_primary_mass_ratio"));
    }

    /// Sun + Jupiter (m_J / m_S ≈ 9.5e-4) sits just inside the warn
    /// threshold (1e-2). Documents the calibration: realistic solar-
    /// system 1PN runs do not trigger the bound. If this test fails
    /// because Jupiter's mass-ratio crosses the warn threshold, the
    /// threshold itself is wrong (Jupiter is a known good test case).
    #[test]
    fn check_regime_silent_for_sun_jupiter() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let bodies = vec![Body::star(1.0), Body::gas_giant(9.547e-4)];
        let violations = pn.check_regime(&bodies, 0.0);
        assert!(
            violations.is_empty(),
            "Sun + Jupiter should be inside 1PN regime; got {violations:?}"
        );
    }

    /// `declared_units` returns `Some(units)` matching what the
    /// operator was constructed with. The `System` registration check
    /// reads this to detect unit-system mismatch.
    #[test]
    fn declared_units_returns_constructor_units() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        assert_eq!(pn.declared_units(), Some(UnitSystem::solar_canonical()));

        let pn_iau = PostNewtonian1PN::for_units(UnitSystem::solar());
        assert_eq!(pn_iau.declared_units(), Some(UnitSystem::solar()));
    }

    /// Acceleration on a lone body must be zero — no self-interaction.
    #[test]
    fn isolated_body_feels_no_pn_force() {
        let bodies = vec![Body::star(1.0)];
        let mut acc = vec![Vec3::ZERO];
        PostNewtonian1PN::for_units(UnitSystem::solar_canonical())
            .accumulate_force(&bodies, &mut acc);
        assert_eq!(acc[0], Vec3::ZERO);
    }

    /// 1PN perturbation is additive: register twice and the effect doubles.
    /// Guards the "must add, not overwrite" clause of the trait contract.
    #[test]
    fn accumulate_is_additive() {
        let bodies =
            vec![Body::star(1.0), Body::rocky(1e-6).at(0.5, 0.0).with_velocity(0.0, 1.414)];
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());

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
        let pn = PostNewtonian1PN::from_raw_c(1e20, UnitSystem::solar_canonical());
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
        PostNewtonian1PN::for_units(UnitSystem::solar_canonical())
            .accumulate_force(&bodies, &mut acc);

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
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let req = <PostNewtonian1PN as Operator>::kernel_requirements(&pn);
        assert_eq!(req, KernelRequirements::exact_and_smooth());
        assert_eq!(req.required_exactness, Some(Exactness::Exact));
        assert_eq!(req.min_continuity, Some(Continuity::Smooth));
    }

    /// `citation()` returns Anderson 1975 + Will 1993 with the
    /// implementing crate's name + version pinned at build time. The
    /// commit hash is `Some(...)` when the build came from a git
    /// checkout (the common case in CI and dev) and `None` otherwise.
    /// We only assert the structural shape — name and version are
    /// known at compile time, the commit-hash branch is environment-
    /// dependent and tested by `c_solar_units_matches_si_derivation`-
    /// style sanity rather than direct equality.
    #[test]
    fn citation_pins_anderson_1975_and_will_1993() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let c = pn.citation().expect("PostNewtonian1PN must publish a citation");
        assert_eq!(c.crate_name, "apsis-1pn", "crate_name must match this crate");
        assert_eq!(c.crate_version, env!("CARGO_PKG_VERSION"));
        assert_eq!(c.doi, Some("10.1086/153180"));
        assert!(c.bibtex.contains("anderson1975"), "bibtex missing primary reference");
        assert!(c.bibtex.contains("will1993"), "bibtex missing textbook reference");
        assert!(c.bibtex.contains("Mariner"), "Anderson 1975 abstract phrase missing from bibtex");
        // Commit hash, if present, must look like a SHA — no
        // accidental garbage like "fatal: not a git repository".
        if let Some(h) = c.commit_hash {
            assert!(h.chars().all(|ch| ch.is_ascii_hexdigit()), "commit_hash not a hex SHA: {h}");
            assert!(h.len() >= 7, "commit_hash suspiciously short: {h}");
        }
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
        PostNewtonian1PN::for_units(UnitSystem::solar_canonical())
            .accumulate_force(&bodies, &mut acc);

        let a_pn = acc[1].length();
        let a_newt = 1.0 / (r_peri * r_peri); // G M / r²
        let ratio = a_pn / a_newt;

        assert!(
            ratio > 1e-8 && ratio < 1e-6,
            "Mercury a_pn/a_newt = {ratio:.3e}, expected band 10⁻⁸–10⁻⁶"
        );
    }
}
