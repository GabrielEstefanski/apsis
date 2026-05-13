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
//! use apsis::units::UnitSystem;
//! use apsis_1pn::PostNewtonian1PN;
//!
//! let units   = UnitSystem::solar_canonical();
//! let sun     = Body::star(1.0).unsoftened();
//! let mercury = Body::rocky(1.66e-7).at(0.307, 0.0).with_velocity(0.0, 2.078).unsoftened();
//! let mut sys = System::new(vec![sun, mercury], units)
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-4);
//!
//! // The same `units` flows through System and operator. The
//! // registration check panics on mismatch.
//! sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(units)));
//!
//! sys.integrate_for(200.0 * std::f64::consts::PI);
//! println!("dE/E = {:.3e}", sys.energy_delta());
//! ```

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::regime::{classify_mass_ratio, mass_ratio};
use apsis::physics::integrator::{HamiltonianOperator, Operator, RegimeViolation, Severity};
use apsis::units::UnitSystem;

/// Mass ratio above which the test-particle pairwise 1PN approximation
/// starts losing accuracy noticeably (warn level). Calibrated so that
/// the Sun–Jupiter case (m_J/m_Sun ≈ 9.5e-4) sits inside the regime,
/// while equal-mass binaries are well outside.
const PN1_MASS_RATIO_WARN: f64 = 1.0e-2;

/// Mass ratio above which the test-particle pairwise 1PN approximation
/// is fundamentally inappropriate (hard level). The full Einstein–
/// Infeld–Hoffmann N-body Hamiltonian is the rigorous form for
/// comparable masses.
const PN1_MASS_RATIO_HARD: f64 = 1.0e-1;

/// Speed of light in m/s — CODATA exact by SI definition.
const C_SI: f64 = 299_792_458.0;

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
/// use apsis::units::UnitSystem;
/// let units = UnitSystem::solar_canonical();
/// let mut sys = System::new(bodies, units).with_integrator(IntegratorKind::Ias15);
/// sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(units)));
/// ```
///
/// The operator carries the [`UnitSystem`] it was constructed for; the
/// `System`'s registration check panics on mismatch (see
/// [`Operator::declared_units`](apsis::physics::integrator::Operator::declared_units)).
/// Constructing with `UnitSystem::solar_canonical()` and registering
/// against a `System` built with `UnitSystem::solar()` (IAU convention)
/// will panic — the two are different unit systems despite the shared
/// "solar" label.
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
    /// Speed of light in `units`.
    c: f64,
    /// Unit system this instance was constructed for. Used by
    /// [`Operator::declared_units`](apsis::physics::integrator::Operator::declared_units)
    /// to fail loudly at registration when paired with a `System`
    /// integrating in a different unit system.
    units: UnitSystem,
}

impl PostNewtonian1PN {
    // ── Named-regime constructor (Pattern A) ──────────────────────────────────

    /// Construct for an arbitrary [`UnitSystem`], deriving `c` exactly
    /// from `c_SI · T_scale / L_scale`. The recommended constructor —
    /// the user picks the unit system once (typically the same one
    /// passed to [`System::new`]), `c` is computed so the relativistic
    /// correction stays consistent with the rest of the integration.
    ///
    /// The most common solar-physics choice is
    /// [`UnitSystem::solar_canonical`] (G = 1, AU, year/(2π), M☉),
    /// matching the apsis-1pn validation portfolio; for IAU
    /// compatibility (G ≈ 4π², AU, year, M☉) use [`UnitSystem::solar`].
    ///
    /// [`System::new`]: apsis::core::system::System::new
    pub fn for_units(units: UnitSystem) -> Self {
        Self { c: C_SI * units.time_scale_si() / units.length_scale_si(), units }
    }

    // ── Raw escape ────────────────────────────────────────────────────────────

    /// Construct from a raw `c` value with the operator pinned to the
    /// supplied [`UnitSystem`]. No cross-check between `c` and `units`
    /// — `c` is taken as given. The `System` registration check still
    /// validates that `units` matches the `System`'s own `UnitSystem`,
    /// so the value cannot land silently in the wrong frame.
    ///
    /// Use when `c` is computed by neighbouring code (so cross-checking
    /// against `units` is redundant), or for hypothetical experiments
    /// where `c` is intentionally non-physical (e.g. "what if `c` were
    /// 5 % larger?"). Prefer [`for_units`](Self::for_units) for normal
    /// physics — it derives `c` from the unit system, eliminating the
    /// raw value entirely.
    pub const fn from_raw_c(c: f64, units: UnitSystem) -> Self {
        Self { c, units }
    }

    /// Speed of light this instance was configured with.
    pub const fn c(&self) -> f64 {
        self.c
    }

    /// Unit system this instance was configured with.
    pub const fn units(&self) -> UnitSystem {
        self.units
    }
}

impl Operator for PostNewtonian1PN {
    fn declared_units(&self) -> Option<UnitSystem> {
        Some(self.units)
    }

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

    /// Test-particle pairwise 1PN assumes `m_secondary ≪ m_primary`
    /// for every secondary in the system (treating `bodies[0]` as the
    /// primary). Beyond ~1 % the leading-order error from the
    /// dropped EIH cross-terms becomes comparable to the 1PN signal
    /// itself; beyond ~10 % the operator's output is no longer
    /// physics — the full EIH N-body Hamiltonian is required.
    ///
    /// `bodies[0]` is treated as the primary by convention (consistent
    /// with the rest of the apsis solar-system fixtures: Sun first,
    /// planets following).
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

    /// 1PN's only declared regime bound is the mass ratio, which is
    /// static for any sane simulation (masses do not change). Check
    /// once per Mercury orbit's worth of steps — about 15 000 IAS15
    /// substeps for the 500-orbit gate. Effectively a no-op in the
    /// hot loop.
    fn regime_check_cadence(&self) -> usize {
        15_000
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

    /// `for_units(UnitSystem::solar_canonical())` returns `c` in the
    /// G=1 solar convention (L = AU, T = Gaussian time, M = M☉) — the
    /// apsis-1pn validation portfolio baseline. The Gaussian time
    /// scale (`sqrt(AU³/(G·M))`) is what makes `G_code = 1` exactly;
    /// it differs from the IAU `year/(2π)` by ~0.009 % (the historical
    /// astrodynamics gap between the Gaussian and IAU definitions).
    /// `C_SOLAR_UNITS` is therefore close to but not bit-equal to
    /// `for_units(solar_canonical).c()` — the constant uses the IAU
    /// year for backwards compatibility, the constructor uses the
    /// Gaussian time so the integrator sees `G = 1` exactly.
    #[test]
    fn for_units_solar_canonical_close_to_c_solar_units() {
        let pn = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
        let rel_diff = (pn.c() - C_SOLAR_UNITS).abs() / C_SOLAR_UNITS;
        // ~0.009 % gap between Gaussian and IAU year definitions.
        assert!(
            rel_diff < 1e-3,
            "for_units(solar_canonical) c={} differs from C_SOLAR_UNITS={} by {:.3e}, \
             expected gap < 0.1 %",
            pn.c(),
            C_SOLAR_UNITS,
            rel_diff,
        );
        // But it's NOT bit-equal — the IAU/Gaussian mismatch is the
        // whole point of using Gaussian time. Lock that.
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
        let bodies = vec![
            Body::star(1.0).unsoftened(),
            Body::rocky(1.66e-7).at(0.387, 0.0).with_velocity(0.0, 1.61).unsoftened(),
        ];
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
        let bodies = vec![Body::star(1.0).unsoftened(), Body::star(1.0).unsoftened()];
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
        let bodies = vec![Body::star(1.0).unsoftened(), Body::gas_giant(9.547e-4).unsoftened()];
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
