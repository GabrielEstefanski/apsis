//! Generalized central force `a = A · r^γ` per Tamayo, Rein, Shi &
//! Hernandez (2019, *MNRAS* 491, 2885). First federation-grade
//! exemplar of the **observable-inversion constructor** locked in
//! [`apsis::contract`] and ADR-005 / ADR-006:
//! [`CentralForce::from_apsidal_rate`] takes a measured apsidal
//! precession rate ω̇ and inverts it into the coupling `A` that
//! reproduces it on the supplied orbit.
//!
//! The naming is deliberate: apsis is named after the apsidal axis,
//! and the headline constructor is the one that reproduces an apsidal
//! observable. The library's name and its highest-leverage feature
//! refer to the same object.
//!
//! # Physics
//!
//! For a body of mass m at separation r from a source body, the
//! generalized central force is
//!
//! ```text
//!   a_central(i ← source) = A · r^γ · r̂
//! ```
//!
//! where r̂ points outward from source to receiver. `γ` parameterizes
//! the radial dependence:
//!
//! | γ      | Force law                    | Notable use case                                         |
//! |--------|------------------------------|----------------------------------------------------------|
//! | `-3`   | `A / r³`                     | Effective Schwarzschild precession (Nobili & Roxburgh 1986) |
//! | `-2`   | `A / r²`                     | Degenerate — looks like gravity, no apsidal precession   |
//! | `-1`   | `A / r`                      | Logarithmic potential (galactic halo flat rotation)      |
//! | `+1`   | `A · r`                      | Hooke / harmonic oscillator                              |
//!
//! Conservative — derives from the closed-form potential
//! `V_central = − A · r^(γ+1) / (γ + 1)` for `γ ≠ −1`, and
//! `V_central = − A · ln(r)` for `γ = −1` (the logarithmic singularity).
//! Both forms are published through
//! [`HamiltonianOperator::potential`](apsis::physics::integrator::HamiltonianOperator::potential)
//! so [`System::energy`](apsis::core::system::System::energy) accounts for
//! the radial contribution.
//!
//! # Observable inversion
//!
//! The Tamayo et al. 2019 result is that a near-circular orbit at
//! instantaneous separation `d` and mean motion `n` precesses at
//!
//! ```text
//!   ω̇ = (1 + γ/2) · A · d^(γ + 2) · n / (G · M_source)
//! ```
//!
//! Inverting:
//!
//! ```text
//!   A = G · M_source · ω̇ / [(1 + γ/2) · d^(γ + 2) · n]
//! ```
//!
//! [`CentralForce::from_apsidal_rate`] consumes the inverted form: the
//! caller supplies an observed (or desired) ω̇, the operator computes
//! `A`, the integrator reproduces ω̇ on a circular orbit. The
//! synthetic round-trip gate (`tests/round_trip_gate.rs`) closes the
//! loop: register at ω̇ = X, integrate, measure ω̇, assert X within
//! tolerance.
//!
//! Limitation: the inversion assumes near-circular geometry. For high
//! eccentricity the apsidal rate picks up an `e`-dependent correction
//! that this constructor does not apply.
//!
//! # Use
//!
//! ```ignore
//! use apsis::core::system::System;
//! use apsis::domain::body::Body;
//! use apsis::physics::integrator::IntegratorKind;
//! use apsis::units::UnitSystem;
//! use apsis_central::CentralForce;
//!
//! let units = UnitSystem::solar_canonical();
//! let sun = Body::star(1.0);
//! let mercury =
//!     Body::rocky(1.66e-7).at(0.387, 0.0).with_velocity(0.0, 1.61);
//! let bodies = vec![sun, mercury];
//!
//! // Observable inversion: pick a desired apsidal rate ω̇ and let
//! // the operator compute the coupling that reproduces it.
//! let omega_dot_desired = 5e-9; // rad per Gaussian time unit
//! let force = CentralForce::from_apsidal_rate(
//!     0,                              // source
//!     1,                              // target
//!     omega_dot_desired,
//!     -3.0,                            // γ = −3 (Schwarzschild-effective)
//!     &bodies,
//!     units,
//! )?;
//!
//! let mut sys = System::new(bodies, units)
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-3);
//! sys.add_hamiltonian_perturbation(Box::new(force))?;
//! ```
//!
//! # Reference
//!
//! Tamayo, D., Rein, H., Shi, P., & Hernandez, D. M. (2019). REBOUNDx:
//! a library for adding conservative and dissipative forces to
//! otherwise symplectic N-body integrations. *MNRAS* 491, 2885–2901.
//! DOI: [10.1093/mnras/stz3018](https://doi.org/10.1093/mnras/stz3018).

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use std::fmt;

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::{
    Citation, HamiltonianOperator, Operator, Potential, RegimeViolation, Severity,
};
use apsis::physics::orbital::compute_elements;
use apsis::units::UnitSystem;

/// Generalized central force `a = A · r^γ` from a single source body
/// onto every other body in the system.
///
/// Stateless. Safe to share across threads.
///
/// Conservative — derives from the closed-form potential exposed by
/// [`HamiltonianOperator::potential`]. Composes additively with apsis's
/// Newtonian gravity: the receiver feels base gravity *and* this
/// central force, both applied to the same `(receiver, source)` pair.
#[derive(Debug, Clone, Copy)]
pub struct CentralForce {
    /// Index of the source body.
    source: usize,
    /// Coupling coefficient `A` in `a = A · r^γ`.
    a_central: f64,
    /// Radial power `γ`.
    gamma: f64,
    /// Unit system this operator was built for.
    units: UnitSystem,
}

impl CentralForce {
    // ── Raw escape ────────────────────────────────────────────────────────────

    /// Construct from an explicit `(A, γ)` pair, pinned to the supplied
    /// [`UnitSystem`]. Use when `A` is computed by neighbouring code or
    /// for direct exploration of the parameter space.
    ///
    /// `γ = -2` is allowed by the constructor — the force is then
    /// proportional to `1/r²`, indistinguishable in shape from
    /// gravity, and produces no apsidal precession. Useful as a
    /// counter-test, not as a physics model.
    pub const fn from_raw(source: usize, a_central: f64, gamma: f64, units: UnitSystem) -> Self {
        Self { source, a_central, gamma, units }
    }

    // ── Observable inversion ──────────────────────────────────────────────────

    /// **Observable-inversion exemplar.** Construct from a desired apsidal
    /// precession rate `omega_dot` (radians per time, in `units`) on
    /// the orbit of `target` around `source`. Inverts the Tamayo
    /// et al. 2019 secular result:
    ///
    /// ```text
    ///   A = G · M_source · ω̇ / [(1 + γ/2) · d^(γ+2) · n]
    /// ```
    ///
    /// where `d` is the instantaneous separation and `n` is the mean
    /// motion of the target's current orbit.
    ///
    /// Assumes near-circular geometry — high `e` introduces an
    /// `e`-dependent correction the inversion does not apply.
    ///
    /// # Errors
    ///
    /// - [`ApsidalInversionError::DegenerateGamma`] — `γ ≈ -2`. The
    ///   precession vanishes as `(1 + γ/2) → 0`, so `A` diverges and
    ///   no finite coupling reproduces the requested rate.
    /// - [`ApsidalInversionError::IndexOutOfRange`] — `source` or
    ///   `target` is past the end of `bodies`.
    /// - [`ApsidalInversionError::SourceEqualsTarget`] — the source
    ///   cannot be its own target.
    /// - [`ApsidalInversionError::UnboundOrbit`] — the target is not
    ///   on a bound orbit around the source; mean motion is undefined.
    pub fn from_apsidal_rate(
        source: usize,
        target: usize,
        omega_dot: f64,
        gamma: f64,
        bodies: &[Body],
        units: UnitSystem,
    ) -> Result<Self, ApsidalInversionError> {
        if source >= bodies.len() {
            return Err(ApsidalInversionError::IndexOutOfRange {
                kind: IndexKind::Source,
                idx: source,
                len: bodies.len(),
            });
        }
        if target >= bodies.len() {
            return Err(ApsidalInversionError::IndexOutOfRange {
                kind: IndexKind::Target,
                idx: target,
                len: bodies.len(),
            });
        }
        if source == target {
            return Err(ApsidalInversionError::SourceEqualsTarget { idx: source });
        }
        // γ = -2 is the degenerate case: precession vanishes, so the
        // inversion has no finite solution.
        if (gamma + 2.0).abs() < 1e-12 {
            return Err(ApsidalInversionError::DegenerateGamma { gamma });
        }

        // Need orbital elements of `target` around `source`. apsis's
        // `compute_elements` takes G as g_factor; canonical units have
        // G = 1, IAU solar has G ≈ 4π² — derive from the unit system
        // so the constructor stays unit-agnostic.
        let g_code = units.g();
        let elems = compute_elements(bodies, target, source, g_code)
            .ok_or(ApsidalInversionError::UnboundOrbit { idx: target })?;
        if !elems.is_bound() {
            return Err(ApsidalInversionError::UnboundOrbit { idx: target });
        }
        let n = elems.mean_motion();
        if !n.is_finite() || n <= 0.0 {
            return Err(ApsidalInversionError::UnboundOrbit { idx: target });
        }

        // Instantaneous separation d = |r_target − r_source|.
        let src = &bodies[source];
        let tgt = &bodies[target];
        let dx = tgt.pos_x - src.pos_x;
        let dy = tgt.pos_y - src.pos_y;
        let dz = tgt.pos_z - src.pos_z;
        let d = (dx * dx + dy * dy + dz * dz).sqrt();

        // µ = G · M_source. apsis stores `mass` in code units; G is
        // already 1 in canonical, ≈ 4π² in IAU solar. The combination
        // is what the inversion needs.
        let mu_source = g_code * src.mass;

        // Inversion: A = µ · ω̇ / [(1 + γ/2) · d^(γ+2) · n]
        let denom = (1.0 + gamma / 2.0) * libm::pow(d, gamma + 2.0) * n;
        let a_central = mu_source * omega_dot / denom;

        Ok(Self { source, a_central, gamma, units })
    }

    // ── Accessors ────────────────────────────────────────────────────────────

    pub const fn a_central(&self) -> f64 {
        self.a_central
    }
    pub const fn gamma(&self) -> f64 {
        self.gamma
    }
    pub const fn source(&self) -> usize {
        self.source
    }
    pub const fn units(&self) -> UnitSystem {
        self.units
    }
}

/// Failure modes for [`CentralForce::from_apsidal_rate`]. `Result`
/// rather than panic so a long pipeline (parameter scan, sensitivity
/// analysis) can decide policy per-failure: skip degenerate γ, log
/// out-of-range indices, propagate unbound orbits.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum ApsidalInversionError {
    /// `γ ≈ -2` — precession identically vanishes for `1/r²` forces,
    /// so no finite `A` reproduces a non-zero ω̇.
    DegenerateGamma { gamma: f64 },
    /// One of `source` or `target` is past the end of `bodies`.
    IndexOutOfRange { kind: IndexKind, idx: usize, len: usize },
    /// `source == target`. A body cannot apsidally-precess about itself.
    SourceEqualsTarget { idx: usize },
    /// The target's orbit around the source is hyperbolic, parabolic,
    /// or degenerate. Mean motion is undefined and the inversion has
    /// no solution.
    UnboundOrbit { idx: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Source,
    Target,
}

impl fmt::Display for ApsidalInversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DegenerateGamma { gamma } => write!(
                f,
                "apsidal inversion: γ = {gamma} is degenerate (precession vanishes for 1/r²); \
                 pick γ ≠ -2"
            ),
            Self::IndexOutOfRange { kind, idx, len } => {
                let which = match kind {
                    IndexKind::Source => "source",
                    IndexKind::Target => "target",
                };
                write!(
                    f,
                    "apsidal inversion: {which} index {idx} out of range for body vector of length {len}"
                )
            },
            Self::SourceEqualsTarget { idx } => {
                write!(f, "apsidal inversion: source and target are the same body (index {idx})")
            },
            Self::UnboundOrbit { idx } => write!(
                f,
                "apsidal inversion: body {idx} is not on a bound orbit; mean motion undefined"
            ),
        }
    }
}

impl std::error::Error for ApsidalInversionError {}

// ── Operator impls ──────────────────────────────────────────────────────────

impl Operator for CentralForce {
    fn name(&self) -> &'static str {
        "CentralForce"
    }

    fn declared_units(&self) -> Option<UnitSystem> {
        Some(self.units)
    }

    /// Generalized central force is well-defined under any kernel —
    /// it's a separate force pair, not a modification of the base
    /// gravity. No requirement on Exactness or Continuity.
    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::none()
    }

    /// Surfaces the source-out-of-range case as `Hard`. β-table style
    /// per-body invariants do not apply here (the operator carries a
    /// single `(A, γ)` rather than a per-body table).
    fn check_regime(&self, bodies: &[Body], _t: f64) -> Vec<RegimeViolation> {
        let mut violations = Vec::new();
        if self.source >= bodies.len() {
            violations.push(RegimeViolation {
                operator: self.name(),
                bound: "source_index_in_range",
                value: self.source as f64,
                threshold: bodies.len() as f64,
                severity: Severity::Hard,
                body_index: Some(self.source),
                message: "source index past end of body vector; no force will be applied",
            });
        }
        violations
    }

    fn regime_check_cadence(&self) -> usize {
        15_000
    }

    fn citation(&self) -> Option<Citation> {
        Some(Citation {
            bibtex: TAMAYO_2019_BIBTEX,
            doi: Some("10.1093/mnras/stz3018"),
            crate_name: env!("CARGO_PKG_NAME"),
            crate_version: env!("CARGO_PKG_VERSION"),
            commit_hash: option_env!("APSIS_CENTRAL_GIT_COMMIT").filter(|s| !s.is_empty()),
            description: Some("General central force after Tamayo et al. 2019"),
            url: Some("https://github.com/GabrielEstefanski/apsis"),
        })
    }
}

impl HamiltonianOperator for CentralForce {
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]) {
        debug_assert_eq!(
            bodies.len(),
            acc.len(),
            "HamiltonianOperator contract: acc must be sized to bodies"
        );
        if self.source >= bodies.len() {
            return;
        }
        let src = &bodies[self.source];
        let m_src = src.mass;
        for i in 0..bodies.len() {
            if i == self.source {
                continue;
            }
            let b_i = &bodies[i];
            // r̂ points from source to receiver (force is "outward
            // from central particle" per Tamayo 2019).
            let dx = b_i.pos_x - src.pos_x;
            let dy = b_i.pos_y - src.pos_y;
            let dz = b_i.pos_z - src.pos_z;
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-30 {
                continue;
            }
            // Per-axis: a = A · r^γ · r̂ = A · r^(γ−1) · (Δ).
            // Compute r^(γ−1) via r²^((γ−1)/2) so we avoid one sqrt
            // when γ is integer-friendly.
            let prefac = self.a_central * libm::pow(r2, (self.gamma - 1.0) / 2.0);
            acc[i].x += prefac * dx;
            acc[i].y += prefac * dy;
            acc[i].z += prefac * dz;
            // Newton's third law on the source: equal-and-opposite
            // contribution scaled by m_receiver / m_source. Without
            // this, momentum drifts — the force pair is asymmetric
            // otherwise.
            let recoil = -b_i.mass / m_src;
            acc[self.source].x += recoil * prefac * dx;
            acc[self.source].y += recoil * prefac * dy;
            acc[self.source].z += recoil * prefac * dz;
        }
    }

    /// Closed-form V from integrating `F = −∇V` for the radial
    /// `a = A · r^γ` law:
    ///
    /// - `γ ≠ −1`: `V = − A · r^(γ+1) / (γ + 1)` per receiver
    /// - `γ = −1`: `V = − A · ln(r)` per receiver (logarithmic
    ///   singularity)
    ///
    /// Sign chosen so that `−∂V/∂r = +A · r^γ` (outward), matching
    /// the force expression. Summed over receivers (m_i · V_per_unit).
    fn potential(&self, bodies: &[Body]) -> Potential {
        if self.source >= bodies.len() {
            return Potential::Value(0.0);
        }
        let src = &bodies[self.source];
        let mut v = 0.0_f64;
        let logarithmic = (self.gamma + 1.0).abs() < 1e-12;
        for i in 0..bodies.len() {
            if i == self.source {
                continue;
            }
            let b_i = &bodies[i];
            let dx = b_i.pos_x - src.pos_x;
            let dy = b_i.pos_y - src.pos_y;
            let dz = b_i.pos_z - src.pos_z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-15 {
                continue;
            }
            let v_per_unit = if logarithmic {
                -self.a_central * libm::log(r)
            } else {
                -self.a_central * libm::pow(r, self.gamma + 1.0) / (self.gamma + 1.0)
            };
            v += b_i.mass * v_per_unit;
        }
        Potential::Value(v)
    }
}

const TAMAYO_2019_BIBTEX: &str = r#"@article{tamayo2019,
  author  = {Tamayo, D. and Rein, H. and Shi, P. and Hernandez, D. M.},
  title   = {{REBOUNDx}: a library for adding conservative and dissipative forces to otherwise symplectic {N}-body integrations},
  journal = {Monthly Notices of the Royal Astronomical Society},
  volume  = {491},
  number  = {2},
  pages   = {2885--2901},
  year    = {2019},
  doi     = {10.1093/mnras/stz3018}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    fn solar() -> UnitSystem {
        UnitSystem::solar_canonical()
    }

    #[test]
    fn from_raw_round_trips_inputs() {
        let f = CentralForce::from_raw(0, 1.5e-7, -3.0, solar());
        assert_eq!(f.source(), 0);
        assert_eq!(f.a_central(), 1.5e-7);
        assert_eq!(f.gamma(), -3.0);
        assert_eq!(f.units(), solar());
    }

    /// Source body must not feel its own central force in the
    /// per-receiver loop. The recoil term is the only contribution
    /// applied to the source index.
    #[test]
    fn source_receives_only_recoil_not_self_force() {
        // Single body at origin = source; nothing for it to act on.
        let bodies = vec![Body::star(1.0)];
        let f = CentralForce::from_raw(0, 1.0, -3.0, solar());
        let mut acc = vec![Vec3::ZERO];
        f.accumulate_force(&bodies, &mut acc);
        assert_eq!(acc[0], Vec3::ZERO, "isolated source must feel zero central force");
    }

    /// Force is purely radial. Tangential placement (receiver at
    /// (1, 0, 0) from source) gives a force along +x with no y/z
    /// component.
    #[test]
    fn force_is_purely_radial() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0)];
        let f = CentralForce::from_raw(0, 1.0, -3.0, solar());
        let mut acc = vec![Vec3::ZERO; 2];
        f.accumulate_force(&bodies, &mut acc);
        // Force on receiver = A · r^γ · r̂. r=1, A=1, γ=-3 → magnitude 1, +x dir.
        assert!((acc[1].x - 1.0).abs() < 1e-14);
        assert!(acc[1].y.abs() < 1e-14);
        assert!(acc[1].z.abs() < 1e-14);
    }

    /// Newton's third law on the source: the recoil contribution
    /// equals `-m_recv/m_src` times the receiver's force. Locks the
    /// momentum-conservation invariant across the operator.
    #[test]
    fn third_law_recoil_matches_mass_ratio() {
        let m_src = 1.0;
        let m_recv = 1e-6;
        let bodies = vec![Body::star(m_src), Body::rocky(m_recv).at(1.0, 0.0)];
        let f = CentralForce::from_raw(0, 1.0, -3.0, solar());
        let mut acc = vec![Vec3::ZERO; 2];
        f.accumulate_force(&bodies, &mut acc);
        // Total momentum derivative = m_src·a_src + m_recv·a_recv.
        // Must be zero (the central force is internal).
        let p_dot_x = m_src * acc[0].x + m_recv * acc[1].x;
        assert!(p_dot_x.abs() < 1e-18, "central force violates momentum conservation: {p_dot_x}");
    }

    /// Force is additive: registering twice doubles the contribution.
    #[test]
    fn force_is_additive() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0)];
        let f = CentralForce::from_raw(0, 0.5, -3.0, solar());

        let mut once = vec![Vec3::ZERO; 2];
        f.accumulate_force(&bodies, &mut once);
        let mut twice = vec![Vec3::ZERO; 2];
        f.accumulate_force(&bodies, &mut twice);
        f.accumulate_force(&bodies, &mut twice);

        for i in 0..2 {
            assert!((twice[i].x - 2.0 * once[i].x).abs() < 1e-14);
            assert!((twice[i].y - 2.0 * once[i].y).abs() < 1e-14);
            assert!((twice[i].z - 2.0 * once[i].z).abs() < 1e-14);
        }
    }

    /// Closed-form V matches `−∇V = accumulate_force` for γ ≠ −1
    /// (polynomial branch).
    #[test]
    fn potential_matches_force_polynomial_branch() {
        let make_bodies = |x: f64| vec![Body::star(1.0), Body::rocky(1e-10).at(x, 0.0)];
        let f = CentralForce::from_raw(0, 0.3, -3.0, solar());
        let h = 1e-6;
        let v_plus = match f.potential(&make_bodies(1.0 + h)) {
            Potential::Value(v) => v,
            _ => panic!("must return Value"),
        };
        let v_minus = match f.potential(&make_bodies(1.0 - h)) {
            Potential::Value(v) => v,
            _ => panic!("must return Value"),
        };
        let m_recv = 1e-10;
        let expected_ax = -(v_plus - v_minus) / (2.0 * h * m_recv);

        let mut acc = vec![Vec3::ZERO; 2];
        f.accumulate_force(&make_bodies(1.0), &mut acc);
        let rel = (acc[1].x - expected_ax).abs() / expected_ax.abs().max(1e-30);
        assert!(rel < 1e-6, "force vs −∇V mismatch: ax = {}, expected {}", acc[1].x, expected_ax);
    }

    /// Logarithmic branch (γ = −1) handled separately. Same gradient
    /// check.
    #[test]
    fn potential_matches_force_logarithmic_branch() {
        let make_bodies = |x: f64| vec![Body::star(1.0), Body::rocky(1e-10).at(x, 0.0)];
        let f = CentralForce::from_raw(0, 0.3, -1.0, solar());
        let h = 1e-6;
        let v_plus = match f.potential(&make_bodies(1.0 + h)) {
            Potential::Value(v) => v,
            _ => panic!(),
        };
        let v_minus = match f.potential(&make_bodies(1.0 - h)) {
            Potential::Value(v) => v,
            _ => panic!(),
        };
        let m_recv = 1e-10;
        let expected_ax = -(v_plus - v_minus) / (2.0 * h * m_recv);

        let mut acc = vec![Vec3::ZERO; 2];
        f.accumulate_force(&make_bodies(1.0), &mut acc);
        let rel = (acc[1].x - expected_ax).abs() / expected_ax.abs().max(1e-30);
        assert!(rel < 1e-6, "log branch: ax = {}, expected {}", acc[1].x, expected_ax);
    }

    // ── Observable-inversion error paths ──────────────────────────────────────

    fn circular_pair() -> Vec<Body> {
        vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 1.0)]
    }

    #[test]
    fn from_apsidal_rate_rejects_degenerate_gamma() {
        let bodies = circular_pair();
        let err = CentralForce::from_apsidal_rate(0, 1, 1e-9, -2.0, &bodies, solar()).unwrap_err();
        assert!(matches!(err, ApsidalInversionError::DegenerateGamma { .. }));
    }

    #[test]
    fn from_apsidal_rate_rejects_source_eq_target() {
        let bodies = circular_pair();
        let err = CentralForce::from_apsidal_rate(0, 0, 1e-9, -3.0, &bodies, solar()).unwrap_err();
        assert!(matches!(err, ApsidalInversionError::SourceEqualsTarget { idx: 0 }));
    }

    #[test]
    fn from_apsidal_rate_rejects_index_out_of_range() {
        let bodies = circular_pair();
        let err = CentralForce::from_apsidal_rate(0, 5, 1e-9, -3.0, &bodies, solar()).unwrap_err();
        assert!(matches!(
            err,
            ApsidalInversionError::IndexOutOfRange { kind: IndexKind::Target, idx: 5, len: 2 }
        ));
    }

    /// Hyperbolic flyby: mean motion undefined → UnboundOrbit.
    #[test]
    fn from_apsidal_rate_rejects_unbound_orbit() {
        let bodies = vec![
            Body::star(1.0),
            // Velocity above escape (v_escape at r=1 in canonical = √2);
            // 2.5 > √2 → hyperbolic.
            Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 2.5),
        ];
        let err = CentralForce::from_apsidal_rate(0, 1, 1e-9, -3.0, &bodies, solar()).unwrap_err();
        assert!(matches!(err, ApsidalInversionError::UnboundOrbit { idx: 1 }));
    }

    /// Observable-inversion happy path: round-trip the inversion algebraically.
    /// Compute A from a desired ω̇, then plug A back into the forward
    /// formula and confirm it produces the same ω̇. This locks the
    /// inversion arithmetic without needing a long integration.
    #[test]
    fn from_apsidal_rate_inversion_is_self_consistent() {
        let bodies = circular_pair();
        let omega_dot_in = 5e-9_f64;
        let gamma = -3.0_f64;
        let f = CentralForce::from_apsidal_rate(0, 1, omega_dot_in, gamma, &bodies, solar())
            .expect("circular pair must invert");

        // Forward formula: ω̇ = (1 + γ/2) · A · d^(γ+2) · n / µ
        let units = solar();
        let g_code = units.g();
        let elems = compute_elements(&bodies, 1, 0, g_code).expect("elements");
        let n = elems.mean_motion();
        let dx = bodies[1].pos_x - bodies[0].pos_x;
        let dy = bodies[1].pos_y - bodies[0].pos_y;
        let dz = bodies[1].pos_z - bodies[0].pos_z;
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        let mu_src = g_code * bodies[0].mass;

        let omega_dot_out =
            (1.0 + gamma / 2.0) * f.a_central() * libm::pow(d, gamma + 2.0) * n / mu_src;
        let rel = ((omega_dot_out - omega_dot_in) / omega_dot_in).abs();
        assert!(
            rel < 1e-12,
            "inversion not self-consistent: in {omega_dot_in}, out {omega_dot_out}, rel {rel:.3e}"
        );
    }

    #[test]
    fn citation_pins_tamayo_2019() {
        let f = CentralForce::from_raw(0, 1.0, -3.0, solar());
        let c = f.citation().expect("CentralForce must publish a citation");
        assert_eq!(c.crate_name, "apsis-central");
        assert_eq!(c.crate_version, env!("CARGO_PKG_VERSION"));
        assert_eq!(c.doi, Some("10.1093/mnras/stz3018"));
        assert!(c.bibtex.contains("tamayo2019"));
        if let Some(h) = c.commit_hash {
            assert!(h.chars().all(|ch| ch.is_ascii_hexdigit()), "bad commit_hash: {h}");
            assert!(h.len() >= 7);
        }
    }

    #[test]
    fn declared_units_returns_constructor_units() {
        let f = CentralForce::from_raw(0, 1.0, -3.0, solar());
        assert_eq!(f.declared_units(), Some(solar()));
    }
}
