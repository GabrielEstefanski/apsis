//! Radiation-pressure perturbations per Burns, Lamy & Soter
//! (1979, *Icarus* 40, 1).
//!
//! [`RadiationPressure`] is the radial 1/r² force opposing gravity
//! (Hamiltonian, conservative). [`PoyntingRobertsonDrag`] is the
//! relativistic v/c term causing semi-major-axis decay
//! (non-conservative). Both take per-body
//! `β = F_rad / F_grav` as input.
//!
//! # Conventions
//!
//! - `bodies[source]` is the radiating body. `betas[source]` must be 0.
//! - `r̂` points from receiver to source.
//! - `β > 1` is the unbound "blowout" regime, integrated without warning.
//!
//! # Use
//!
//! ```ignore
//! let units = UnitSystem::solar_canonical();
//! let sun = Body::star(1.0);
//! let dust = Body::rocky(1e-15).at(1.0, 0.0).with_velocity(0.0, 1.0);
//! let mut sys = System::new(vec![sun, dust], units)
//!     .with_integrator(IntegratorKind::Ias15)
//!     .with_dt(1e-3);
//! sys.add_hamiltonian_perturbation(Box::new(
//!     RadiationPressure::from_raw_betas(0, vec![0.0, 0.1], units),
//! ))?;
//! sys.add_non_conservative_perturbation(Box::new(
//!     PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.1], units),
//! ))?;
//! ```
//!
//! # Reference
//!
//! Burns, J. A., Lamy, P. L., & Soter, S. (1979).
//! DOI: [10.1016/0019-1035(79)90050-2](https://doi.org/10.1016/0019-1035(79)90050-2).

#![deny(unsafe_code)]
#![allow(clippy::needless_range_loop)]

use apsis::domain::body::Body;
use apsis::math::Vec3;
use apsis::physics::gravity::kernel::KernelRequirements;
use apsis::physics::integrator::{
    Citation, HamiltonianOperator, NonConservativeOperator, Operator, Potential, RegimeViolation,
    Severity,
};
use apsis::units::UnitSystem;

/// Speed of light in m/s — CODATA exact by SI definition.
const C_SI: f64 = 299_792_458.0;

/// Per-body β bookkeeping shared by [`RadiationPressure`] and
/// [`PoyntingRobertsonDrag`].
#[derive(Debug, Clone)]
struct BetaTable {
    source: usize,
    betas: Vec<f64>,
    /// Speed of light in `units`. Cached at construction.
    c: f64,
    units: UnitSystem,
}

impl BetaTable {
    fn new(source: usize, betas: Vec<f64>, units: UnitSystem) -> Self {
        let c = C_SI * units.time_scale_si() / units.length_scale_si();
        Self { source, betas, c, units }
    }

    /// β for body `i`, or 0 if the table is shorter than `bodies`.
    fn beta_for(&self, i: usize) -> f64 {
        self.betas.get(i).copied().unwrap_or(0.0)
    }
}

// ── RadiationPressure (Hamiltonian) ──────────────────────────────────────────

/// Radial radiation pressure: conservative central force scaled per
/// receiver by `β`. Force expression
/// `F_rad(i ← source) = + β_i · G · M_source · m_i / r² · r̂`;
/// closed-form potential `V_rad = − β · G · M · m / r` published via
/// [`HamiltonianOperator::potential`].
#[derive(Debug, Clone)]
pub struct RadiationPressure {
    table: BetaTable,
}

impl RadiationPressure {
    /// `source` is the index of the radiating body; `betas[i]` is body
    /// `i`'s β. `betas[source]` must be 0 (enforced at registration via
    /// [`Operator::check_regime`]).
    pub fn from_raw_betas(source: usize, betas: Vec<f64>, units: UnitSystem) -> Self {
        Self { table: BetaTable::new(source, betas, units) }
    }

    pub fn beta_for(&self, i: usize) -> f64 {
        self.table.beta_for(i)
    }

    pub fn source(&self) -> usize {
        self.table.source
    }

    pub fn units(&self) -> UnitSystem {
        self.table.units
    }
}

impl Operator for RadiationPressure {
    fn name(&self) -> &'static str {
        "RadiationPressure"
    }

    fn declared_units(&self) -> Option<UnitSystem> {
        Some(self.table.units)
    }

    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::none()
    }

    fn check_regime(&self, bodies: &[Body], _t: f64) -> Vec<RegimeViolation> {
        let mut violations = Vec::new();
        if self.table.betas.len() != bodies.len() {
            violations.push(RegimeViolation {
                operator: self.name(),
                bound: "betas_table_length",
                value: self.table.betas.len() as f64,
                threshold: bodies.len() as f64,
                severity: Severity::Hard,
                body_index: None,
                message: "betas table length differs from body vector length; out-of-range \
                          indices are treated as β = 0, which silently disables radiation \
                          on bodies past the table end",
            });
        }
        if self.table.source < self.table.betas.len() && self.table.betas[self.table.source] != 0.0
        {
            violations.push(RegimeViolation {
                operator: self.name(),
                bound: "self_radiation",
                value: self.table.betas[self.table.source],
                threshold: 0.0,
                severity: Severity::Hard,
                body_index: Some(self.table.source),
                message: "source body must not feel its own radiation; \
                          set betas[source] = 0",
            });
        }
        violations
    }

    /// β table is static — the only thing the dynamic check could catch
    /// is bodies being added/removed mid-run. Cadence matches the 1PN
    /// rate so the per-step cost stays negligible.
    fn regime_check_cadence(&self) -> usize {
        15_000
    }

    fn citation(&self) -> Option<Citation> {
        Some(Citation {
            bibtex: BURNS_LAMY_SOTER_BIBTEX,
            doi: Some("10.1016/0019-1035(79)90050-2"),
            crate_name: env!("CARGO_PKG_NAME"),
            crate_version: env!("CARGO_PKG_VERSION"),
            commit_hash: option_env!("APSIS_RADIATION_GIT_COMMIT").filter(|s| !s.is_empty()),
            description: Some("Radiation pressure and Poynting--Robertson drag after Burns 1979"),
            url: Some("https://github.com/GabrielEstefanski/apsis"),
            author: Some("Estefanski, G. B."),
        })
    }
}

impl HamiltonianOperator for RadiationPressure {
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]) {
        debug_assert_eq!(
            bodies.len(),
            acc.len(),
            "HamiltonianOperator contract: acc must be sized to bodies"
        );
        if self.table.source >= bodies.len() {
            return;
        }
        let src = &bodies[self.table.source];
        for i in 0..bodies.len() {
            if i == self.table.source {
                continue;
            }
            let beta = self.beta_for(i);
            if beta == 0.0 {
                continue;
            }
            let b_i = &bodies[i];
            let dx = src.pos_x - b_i.pos_x;
            let dy = src.pos_y - b_i.pos_y;
            let dz = src.pos_z - b_i.pos_z;
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-30 {
                continue;
            }
            let inv_r = r2.sqrt().recip();
            // r̂ points from receiver i to source. Radiation force is
            // outward from source = away from source = −r̂ direction
            // for the receiver. So a_rad on i is along −r̂.
            //
            //   a_rad = − β · (G · M_src / r²) · r̂      with G = 1
            //
            // The minus sign here is what makes the gravitational
            // acceleration on the receiver appear with effective
            // GM_eff = GM · (1 − β): apsis's Newtonian kernel
            // contributes +GM/r²·r̂ on the receiver (towards source);
            // we add −β·GM/r²·r̂ (away).
            let pref = -beta * src.mass / r2;
            acc[i].x += pref * dx * inv_r;
            acc[i].y += pref * dy * inv_r;
            acc[i].z += pref * dz * inv_r;
        }
    }

    /// Closed-form V from the central radiation potential summed over
    /// receivers: `V_rad = − β_i · G · M_src · m_i / r_i`. Sign chosen
    /// so that ∂V/∂r_i = `−` accumulate_force for body i, matching the
    /// HamiltonianOperator contract (force = −∇V).
    fn potential(&self, bodies: &[Body]) -> Potential {
        if self.table.source >= bodies.len() {
            return Potential::Value(0.0);
        }
        let src = &bodies[self.table.source];
        let mut v = 0.0_f64;
        for i in 0..bodies.len() {
            if i == self.table.source {
                continue;
            }
            let beta = self.beta_for(i);
            if beta == 0.0 {
                continue;
            }
            let b_i = &bodies[i];
            let dx = src.pos_x - b_i.pos_x;
            let dy = src.pos_y - b_i.pos_y;
            let dz = src.pos_z - b_i.pos_z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-15 {
                continue;
            }
            // Gravity contributes V_grav = −G·M·m/r (G=1). Radiation
            // pressure reduces effective gravity by factor (1−β), i.e.
            // it contributes V_rad = +β·M·m/r.
            v += beta * src.mass * b_i.mass / r;
        }
        Potential::Value(v)
    }
}

// ── PoyntingRobertsonDrag (non-conservative) ─────────────────────────────────

/// Relativistic angular-momentum loss from re-emitted radiation.
/// Dissipative — register against IAS15 (energy drift IS the physical
/// signal); symplectic integrators emit a warning at registration.
///
/// Stateless. Safe to share across threads.
///
/// Constructed with the same `(source, betas, units)` triple as
/// [`RadiationPressure`] — the two operators share the per-body β
/// definition (Burns et al. eq. 17) and are typically registered
/// together. β = 0 entries silently skip the receiver.
#[derive(Debug, Clone)]
pub struct PoyntingRobertsonDrag {
    table: BetaTable,
}

impl PoyntingRobertsonDrag {
    /// Construct from a raw per-body β array. Same calling convention
    /// as [`RadiationPressure::from_raw_betas`].
    pub fn from_raw_betas(source: usize, betas: Vec<f64>, units: UnitSystem) -> Self {
        Self { table: BetaTable::new(source, betas, units) }
    }

    /// β value applied to body `i`.
    pub fn beta_for(&self, i: usize) -> f64 {
        self.table.beta_for(i)
    }

    /// Index of the radiating body.
    pub fn source(&self) -> usize {
        self.table.source
    }

    /// Unit system this operator was constructed for.
    pub fn units(&self) -> UnitSystem {
        self.table.units
    }
}

impl Operator for PoyntingRobertsonDrag {
    fn name(&self) -> &'static str {
        "PoyntingRobertsonDrag"
    }

    fn declared_units(&self) -> Option<UnitSystem> {
        Some(self.table.units)
    }

    fn kernel_requirements(&self) -> KernelRequirements {
        KernelRequirements::none()
    }

    /// Same regime checks as [`RadiationPressure`] — the operators
    /// share the β-table invariants. Duplicated rather than factored
    /// out so each operator's regime warnings name its own type.
    fn check_regime(&self, bodies: &[Body], _t: f64) -> Vec<RegimeViolation> {
        let mut violations = Vec::new();
        if self.table.betas.len() != bodies.len() {
            violations.push(RegimeViolation {
                operator: self.name(),
                bound: "betas_table_length",
                value: self.table.betas.len() as f64,
                threshold: bodies.len() as f64,
                severity: Severity::Hard,
                body_index: None,
                message: "betas table length differs from body vector length; out-of-range \
                          indices are treated as β = 0, which silently disables PR drag \
                          on bodies past the table end",
            });
        }
        if self.table.source < self.table.betas.len() && self.table.betas[self.table.source] != 0.0
        {
            violations.push(RegimeViolation {
                operator: self.name(),
                bound: "self_radiation",
                value: self.table.betas[self.table.source],
                threshold: 0.0,
                severity: Severity::Hard,
                body_index: Some(self.table.source),
                message: "source body must not feel its own radiation; \
                          set betas[source] = 0",
            });
        }
        violations
    }

    fn regime_check_cadence(&self) -> usize {
        15_000
    }

    fn citation(&self) -> Option<Citation> {
        // Same paper, different equation — Burns et al. eq. (7) for
        // the PR drag, eq. (4–5) for the radial pressure.
        Some(Citation {
            bibtex: BURNS_LAMY_SOTER_BIBTEX,
            doi: Some("10.1016/0019-1035(79)90050-2"),
            crate_name: env!("CARGO_PKG_NAME"),
            crate_version: env!("CARGO_PKG_VERSION"),
            commit_hash: option_env!("APSIS_RADIATION_GIT_COMMIT").filter(|s| !s.is_empty()),
            description: Some("Radiation pressure and Poynting--Robertson drag after Burns 1979"),
            url: Some("https://github.com/GabrielEstefanski/apsis"),
            author: Some("Estefanski, G. B."),
        })
    }
}

impl NonConservativeOperator for PoyntingRobertsonDrag {
    fn accumulate_force(&self, bodies: &[Body], acc: &mut [Vec3]) {
        debug_assert_eq!(
            bodies.len(),
            acc.len(),
            "NonConservativeOperator contract: acc must be sized to bodies"
        );
        if self.table.source >= bodies.len() {
            return;
        }
        let src = &bodies[self.table.source];
        let inv_c = self.table.c.recip();
        for i in 0..bodies.len() {
            if i == self.table.source {
                continue;
            }
            let beta = self.beta_for(i);
            if beta == 0.0 {
                continue;
            }
            let b_i = &bodies[i];
            // Source-relative geometry. r̂ points receiver → source.
            let dx = src.pos_x - b_i.pos_x;
            let dy = src.pos_y - b_i.pos_y;
            let dz = src.pos_z - b_i.pos_z;
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 < 1e-30 {
                continue;
            }
            let inv_r = r2.sqrt().recip();
            let rhat_x = dx * inv_r;
            let rhat_y = dy * inv_r;
            let rhat_z = dz * inv_r;

            // Receiver velocity in inertial frame. The Burns et al.
            // derivation assumes the source is approximately at rest;
            // for a heliocentric problem with bodies[0] = Sun this is
            // the leading-order approximation. Heavy-source corrections
            // are O(v_src/c) smaller and out of scope for this crate.
            let vx = b_i.vel_x;
            let vy = b_i.vel_y;
            let vz = b_i.vel_z;
            let v_radial = vx * rhat_x + vy * rhat_y + vz * rhat_z;

            let pref = beta * src.mass / r2;
            // Burns et al. eq. (7): F_PR = − (β·GM/r²) · [(2·v_r/c)·r̂ + v/c]
            //
            // The sign is opposite to radiation pressure on the radial
            // component — PR drag opposes motion, radial pressure pushes
            // outward. The (v/c) term decelerates the receiver's
            // tangential motion → angular momentum loss → semi-major
            // axis decay.
            let s_rhat = -2.0 * v_radial * inv_c;
            let s_v = -inv_c;
            acc[i].x += pref * (s_rhat * rhat_x + s_v * vx);
            acc[i].y += pref * (s_rhat * rhat_y + s_v * vy);
            acc[i].z += pref * (s_rhat * rhat_z + s_v * vz);
        }
    }
}

const BURNS_LAMY_SOTER_BIBTEX: &str = r#"@article{burns1979,
  author  = {Burns, J. A. and Lamy, P. L. and Soter, S.},
  title   = {Radiation forces on small particles in the solar system},
  journal = {Icarus},
  volume  = {40},
  number  = {1},
  pages   = {1--48},
  year    = {1979},
  doi     = {10.1016/0019-1035(79)90050-2}
}"#;

#[cfg(test)]
mod tests {
    use super::*;

    fn solar() -> UnitSystem {
        UnitSystem::solar_canonical()
    }

    /// `from_raw_betas` accepts any β table and pins the operator to
    /// the supplied unit system. Reading back via the accessors must
    /// return what was stored.
    #[test]
    fn from_raw_betas_round_trips_inputs() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.5, 0.1], solar());
        assert_eq!(r.source(), 0);
        assert_eq!(r.beta_for(0), 0.0);
        assert_eq!(r.beta_for(1), 0.5);
        assert_eq!(r.beta_for(2), 0.1);
        assert_eq!(r.units(), solar());
    }

    /// Radiation pressure on the source body itself must be zero.
    /// Locks the "skip i == source" branch in accumulate_force.
    #[test]
    fn radiation_pressure_zero_on_source() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 1.0)];
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.5], solar());
        let mut acc = vec![Vec3::ZERO; 2];
        r.accumulate_force(&bodies, &mut acc);
        assert_eq!(acc[0], Vec3::ZERO, "source body must not feel its own radiation");
    }

    /// β = 0 receiver feels nothing. Compact bodies (1PN regime) should
    /// be transparent to the radiation operator.
    #[test]
    fn radiation_pressure_skips_zero_beta_receivers() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-7).at(1.0, 0.0).with_velocity(0.0, 1.0)];
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.0], solar());
        let mut acc = vec![Vec3::ZERO; 2];
        r.accumulate_force(&bodies, &mut acc);
        assert_eq!(acc[1], Vec3::ZERO, "β = 0 receiver must feel no radiation force");
    }

    /// Effective-gravity equivalence: applying RadiationPressure with β
    /// must reduce the total radial acceleration on the receiver by a
    /// factor β relative to pure gravity at the same geometry. This is
    /// the federation pattern's strongest claim — radiation pressure
    /// is gravity with reduced GM.
    #[test]
    fn radiation_pressure_reduces_effective_gravity_by_beta() {
        // Receiver at +1 AU on the +x axis; gravity points in −x.
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0)];
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.3], solar());
        let mut acc = vec![Vec3::ZERO; 2];
        r.accumulate_force(&bodies, &mut acc);
        // Gravity on the receiver is +1 AU away → a_grav points in −x
        // with magnitude G·M/r² = 1 (G=1, M=1, r=1). Radiation
        // pressure adds +β in +x. Net acc on receiver should be +0.3·x̂.
        assert!((acc[1].x - 0.3).abs() < 1e-14, "expected +β x̂ contribution, got {}", acc[1].x);
        assert!(acc[1].y.abs() < 1e-14, "ay must be 0 for purely radial source");
        assert!(acc[1].z.abs() < 1e-14, "az must be 0 in plane");
    }

    /// Closed-form potential matches −∇V = accumulate_force for a
    /// simple two-body geometry. Numerical gradient with central
    /// differences keeps tolerances loose enough to survive ULP noise.
    #[test]
    fn radiation_pressure_potential_matches_force_to_central_diff() {
        let beta = 0.4;
        let make_bodies =
            |x: f64| vec![Body::star(1.0), Body::rocky(1e-10).at(x, 0.0).with_velocity(0.0, 0.0)];
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, beta], solar());

        let h = 1e-6;
        let v_plus = match r.potential(&make_bodies(1.0 + h)) {
            Potential::Value(v) => v,
            _ => panic!("must return Value"),
        };
        let v_minus = match r.potential(&make_bodies(1.0 - h)) {
            Potential::Value(v) => v,
            _ => panic!("must return Value"),
        };
        // Force on body 1 in x-direction: F_x = − ∂V/∂x. Body 1 has
        // mass 1e-10 → divide gradient by mass to get acceleration.
        let m_receiver = 1e-10;
        let expected_ax = -(v_plus - v_minus) / (2.0 * h * m_receiver);

        let mut acc = vec![Vec3::ZERO; 2];
        r.accumulate_force(&make_bodies(1.0), &mut acc);

        let rel = (acc[1].x - expected_ax).abs() / expected_ax.abs().max(1e-30);
        assert!(rel < 1e-6, "force vs −∇V mismatch: ax = {}, expected {}", acc[1].x, expected_ax);
    }

    /// PR drag on the source body itself must be zero.
    #[test]
    fn pr_drag_zero_on_source() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 1.0)];
        let p = PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.5], solar());
        let mut acc = vec![Vec3::ZERO; 2];
        p.accumulate_force(&bodies, &mut acc);
        assert_eq!(acc[0], Vec3::ZERO);
    }

    /// PR drag opposes tangential motion: for a circular orbit
    /// (purely tangential v at radius r), the drag must have a
    /// negative tangential component (same sign as −v).
    #[test]
    fn pr_drag_opposes_tangential_motion() {
        // Receiver at +x axis with +y tangential velocity. Drag should
        // be in −y direction (decelerating).
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 1.0)];
        let p = PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.5], solar());
        let mut acc = vec![Vec3::ZERO; 2];
        p.accumulate_force(&bodies, &mut acc);
        assert!(acc[1].y < 0.0, "PR drag must oppose tangential motion (ay < 0), got {}", acc[1].y);
    }

    /// At c → ∞ the PR drag vanishes. The 1/c prefactor in Burns
    /// eq. (7) governs the scale.
    #[test]
    fn pr_drag_vanishes_at_infinite_c() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0).with_velocity(0.0, 1.0)];
        // Override c by constructing a unit system whose length/time
        // scale push c → very large. Easier: use canonical (c = c_SI
        // ≈ 3e8) and just use a large value via raw_betas with a
        // hand-built unit system.
        let huge = UnitSystem::custom(1.0, 1e20, 1.0).expect("custom units");
        let p = PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.5], huge);
        let mut acc = vec![Vec3::ZERO; 2];
        p.accumulate_force(&bodies, &mut acc);
        assert!(acc[1].length() < 1e-20, "PR drag should vanish at c → ∞, got {}", acc[1].length());
    }

    /// Mis-sized β table is flagged as Hard by check_regime.
    #[test]
    fn check_regime_flags_mismatched_betas_length() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.0], solar()); // 1 entry
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10)]; // 2 bodies
        let v = r.check_regime(&bodies, 0.0);
        assert!(v.iter().any(|x| x.bound == "betas_table_length" && x.severity == Severity::Hard));
    }

    /// Source-self-radiation invariant flagged when betas[source] != 0.
    #[test]
    fn check_regime_flags_self_radiation() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.5, 0.1], solar());
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10)];
        let v = r.check_regime(&bodies, 0.0);
        assert!(v.iter().any(|x| x.bound == "self_radiation" && x.severity == Severity::Hard));
    }

    /// In-regime configuration is silent.
    #[test]
    fn check_regime_silent_for_well_formed_table() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.5], solar());
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10)];
        assert!(r.check_regime(&bodies, 0.0).is_empty());
    }

    /// Both operators publish citations pointing at Burns/Lamy/Soter
    /// 1979 with the apsis-radiation crate's name + version.
    #[test]
    fn both_operators_cite_burns_1979() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.5], solar());
        let p = PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.5], solar());

        for op_name in [
            r.citation().expect("RadiationPressure must publish citation"),
            p.citation().expect("PoyntingRobertsonDrag must publish citation"),
        ] {
            assert_eq!(op_name.crate_name, "apsis-radiation");
            assert_eq!(op_name.crate_version, env!("CARGO_PKG_VERSION"));
            assert_eq!(op_name.doi, Some("10.1016/0019-1035(79)90050-2"));
            assert!(op_name.bibtex.contains("burns1979"), "bibtex missing key reference");
        }
    }

    /// `declared_units` returns `Some(units)` for both operators —
    /// drives the System registration check.
    #[test]
    fn declared_units_returns_constructor_units() {
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.5], solar());
        assert_eq!(r.declared_units(), Some(solar()));
        let p = PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.5], solar());
        assert_eq!(p.declared_units(), Some(solar()));
    }

    /// Force is additive: registering twice doubles the contribution.
    /// Locks the "must add, not overwrite" trait clause.
    #[test]
    fn radiation_pressure_force_is_additive() {
        let bodies = vec![Body::star(1.0), Body::rocky(1e-10).at(1.0, 0.0)];
        let r = RadiationPressure::from_raw_betas(0, vec![0.0, 0.3], solar());

        let mut once = vec![Vec3::ZERO; 2];
        r.accumulate_force(&bodies, &mut once);
        let mut twice = vec![Vec3::ZERO; 2];
        r.accumulate_force(&bodies, &mut twice);
        r.accumulate_force(&bodies, &mut twice);

        for i in 0..2 {
            assert!((twice[i].x - 2.0 * once[i].x).abs() < 1e-14);
            assert!((twice[i].y - 2.0 * once[i].y).abs() < 1e-14);
            assert!((twice[i].z - 2.0 * once[i].z).abs() < 1e-14);
        }
    }
}
