//! WHFast — Wisdom-Holman Fast (Rein & Tamayo 2015, *MNRAS* 452, 376).
//!
//! Symplectic mixed-variable integrator in democratic-heliocentric
//! coordinates. Same KDK structure as [`super::wisdom_holman::WisdomHolman`]
//! with persistent per-body compensators (Neumaier) on the position
//! and velocity accumulators, reducing the round-off envelope on
//! length-N sums from `O(N · ε)` to `O(√N · ε)` (Higham 2002 §4.5).
//!
//! Lab notebook:
//! `docs/experiments/2026-05-13-whfast-integrator.md`.
//!
//! # References
//!
//! - Rein, H. & Tamayo, D. (2015). *MNRAS* 452, 376.
//! - Wisdom, J. & Holman, M. (1991). *AJ* 102, 1528.

use crate::domain::body::Body;
use crate::math::{CompensatedF64, Vec3};
use crate::physics::integrator::dense::{DenseSnapshot, WhDenseData};
use crate::physics::integrator::helpers::{evaluate, scale_acc_and_pe};
use crate::physics::integrator::kepler::kepler_step;
use crate::physics::integrator::operator_dispatch::accumulate_perturbation_forces;
use crate::physics::integrator::traits::{
    HierarchySignal, Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Minimum `M_central / Σ m_i (i > 0)` for the WH-class derivation.
const WHFAST_DOMINANCE_RATIO: f64 = 10.0;

/// WHFast integrator. Carries persistent per-body compensators
/// (`cs_pos`, `cs_vel`) that survive across step calls and accumulate
/// the round-off bits a naive `f64 +=` would drop. Compensators
/// resize lazily on first step; a body-count change zeroes them.
///
/// The slot-0 (central body) compensators are allocated but never
/// written: the Sun's inertial state is reconstructed analytically from
/// barycenter conservation each step, not accumulated. The unused slot
/// keeps `body_idx = i + 1` indexing trivial throughout.
#[derive(Default)]
pub struct WHFast {
    cs_pos: Vec<Vec3>,
    cs_vel: Vec<Vec3>,
}

impl WHFast {
    pub fn new() -> Self {
        Self::default()
    }

    /// `true` when `bodies[0]` dominates the system mass distribution
    /// to the WH-class threshold (`M_central ≥ 10 · Σ m_rest`).
    pub fn is_suitable_for(bodies: &[Body]) -> bool {
        if bodies.len() < 2 {
            return false;
        }
        let m0 = bodies[0].mass;
        let m_rest: f64 = bodies[1..].iter().map(|b| b.mass).sum();
        let max_other = bodies[1..].iter().map(|b| b.mass).fold(0.0_f64, f64::max);
        m0 >= max_other && m0 >= WHFAST_DOMINANCE_RATIO * m_rest
    }

    fn ensure_state_size(&mut self, n: usize) {
        if self.cs_pos.len() != n {
            self.cs_pos.clear();
            self.cs_pos.resize(n, Vec3::ZERO);
            self.cs_vel.clear();
            self.cs_vel.resize(n, Vec3::ZERO);
        }
    }
}

/// Neumaier-compensated `+=` on the three axes of a `(value, comp)`
/// Vec3 pair.
#[inline(always)]
fn compensated_add(
    value_x: &mut f64,
    value_y: &mut f64,
    value_z: &mut f64,
    comp: &mut Vec3,
    delta: Vec3,
) {
    let mut cx = CompensatedF64::new(*value_x, comp.x);
    let mut cy = CompensatedF64::new(*value_y, comp.y);
    let mut cz = CompensatedF64::new(*value_z, comp.z);
    cx += delta.x;
    cy += delta.y;
    cz += delta.z;
    *value_x = cx.value;
    *value_y = cy.value;
    *value_z = cz.value;
    comp.x = cx.comp;
    comp.y = cy.comp;
    comp.z = cz.comp;
}

impl Integrator for WHFast {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        if bodies.len() < 2 {
            return StepResult {
                consumed_dt: dt,
                potential_energy: 0.0,
                used_fallback: false,
                step_snapshot: None,
                degraded: false,
                hierarchy_signal: Some(HierarchySignal::Violated),
            };
        }

        self.ensure_state_size(bodies.len());

        let m0 = bodies[0].mass;
        let m_total: f64 = bodies.iter().map(|b| b.mass).sum();
        let mu = ctx.g_factor * m0;

        let pre_x0_inertial: Vec<Vec3> =
            bodies.iter().map(|b| Vec3::new(b.pos_x, b.pos_y, b.pos_z)).collect();
        let pre_v0_inertial: Vec<Vec3> =
            bodies.iter().map(|b| Vec3::new(b.vel_x, b.vel_y, b.vel_z)).collect();
        let pre_a0_inertial: Vec<Vec3> =
            if acc.len() == bodies.len() { acc.clone() } else { vec![Vec3::ZERO; bodies.len()] };

        // ── Galilean shift to the centre-of-mass rest frame ──────────────
        let p_total = bodies
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z));
        let v_com = p_total / m_total;
        for b in bodies.iter_mut() {
            b.vel_x -= v_com.x;
            b.vel_y -= v_com.y;
            b.vel_z -= v_com.z;
        }

        let r0_in = Vec3::new(bodies[0].pos_x, bodies[0].pos_y, bodies[0].pos_z);
        let m_q_in: Vec3 = bodies[1..]
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * (Vec3::new(b.pos_x, b.pos_y, b.pos_z) - r0_in));

        // ── Translate planets to heliocentric coordinates ────────────────
        for b in bodies[1..].iter_mut() {
            b.pos_x -= r0_in.x;
            b.pos_y -= r0_in.y;
            b.pos_z -= r0_in.z;
        }

        let q0_helio_rest: Vec<Vec3> =
            bodies[1..].iter().map(|b| Vec3::new(b.pos_x, b.pos_y, b.pos_z)).collect();
        let v0_inertial_rest: Vec<Vec3> =
            bodies[1..].iter().map(|b| Vec3::new(b.vel_x, b.vel_y, b.vel_z)).collect();
        let planet_masses: Vec<f64> = bodies[1..].iter().map(|b| b.mass).collect();

        // ── First half-kick (compensated on velocity) ────────────────────
        let pe = whfast_kick(bodies, ctx, 0.5 * dt, acc, mu, &mut self.cs_vel);

        // ── Drift: H_K (Kepler around fixed central, per planet) ──────────
        // Compensated on the position delta. Velocity delta from kepler_step
        // is also compensated (kepler_step changes v by an O(dt) amount per
        // step).
        for (i, b) in bodies[1..].iter_mut().enumerate() {
            let body_idx = i + 1;
            let q = Vec3::new(b.pos_x, b.pos_y, b.pos_z);
            let v = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
            let (q_new, v_new) = kepler_step(q, v, dt, mu);
            let dq = q_new - q;
            let dv = v_new - v;
            compensated_add(
                &mut b.pos_x,
                &mut b.pos_y,
                &mut b.pos_z,
                &mut self.cs_pos[body_idx],
                dq,
            );
            compensated_add(
                &mut b.vel_x,
                &mut b.vel_y,
                &mut b.vel_z,
                &mut self.cs_vel[body_idx],
                dv,
            );
        }

        // ── Drift: H_indirect (uniform shift on heliocentric positions) ──
        let p_planets_post_kepler: Vec3 = bodies[1..]
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z));
        let indirect_shift = (p_planets_post_kepler / m0) * dt;
        for (i, b) in bodies[1..].iter_mut().enumerate() {
            let body_idx = i + 1;
            compensated_add(
                &mut b.pos_x,
                &mut b.pos_y,
                &mut b.pos_z,
                &mut self.cs_pos[body_idx],
                indirect_shift,
            );
        }

        // ── Second half-kick ─────────────────────────────────────────────
        let _ = whfast_kick(bodies, ctx, 0.5 * dt, acc, mu, &mut self.cs_vel);

        // ── Reconstruct central body inertial state via barycenter ───────
        let m_q_out: Vec3 = bodies[1..]
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.pos_x, b.pos_y, b.pos_z));
        let r0_out = r0_in + (m_q_in - m_q_out) / m_total;

        // ── Translate planets back from heliocentric ─────────────────────
        for b in bodies[1..].iter_mut() {
            b.pos_x += r0_out.x;
            b.pos_y += r0_out.y;
            b.pos_z += r0_out.z;
        }

        let p_planets_out: Vec3 = bodies[1..]
            .iter()
            .fold(Vec3::ZERO, |s, b| s + b.mass * Vec3::new(b.vel_x, b.vel_y, b.vel_z));
        let v0_out_rest = -p_planets_out / m0;

        bodies[0].pos_x = r0_out.x;
        bodies[0].pos_y = r0_out.y;
        bodies[0].pos_z = r0_out.z;
        bodies[0].vel_x = v0_out_rest.x;
        bodies[0].vel_y = v0_out_rest.y;
        bodies[0].vel_z = v0_out_rest.z;

        // ── Inverse Galilean shift back to the original frame ────────────
        let dr_com = v_com * dt;
        for b in bodies.iter_mut() {
            b.pos_x += dr_com.x;
            b.pos_y += dr_com.y;
            b.pos_z += dr_com.z;
            b.vel_x += v_com.x;
            b.vel_y += v_com.y;
            b.vel_z += v_com.z;
        }

        // ── Body-aligned acceleration buffer ─────────────────────────────
        let n = bodies.len();
        let r0_inertial = Vec3::new(bodies[0].pos_x, bodies[0].pos_y, bodies[0].pos_z);
        let mut planet_total_acc = Vec::with_capacity(n - 1);
        let mut sun_acc = Vec3::ZERO;
        for (i, b) in bodies[1..].iter().enumerate() {
            let q = Vec3::new(b.pos_x, b.pos_y, b.pos_z) - r0_inertial;
            let r2 = q.length_squared().max(1e-60);
            let inv_r3 = 1.0 / (r2 * r2.sqrt());
            let kepler_pull_on_planet = -mu * q * inv_r3;
            planet_total_acc.push(acc[i] + kepler_pull_on_planet);
            sun_acc += ctx.g_factor * b.mass * q * inv_r3;
        }
        acc.resize(n, Vec3::ZERO);
        acc[0] = sun_acc;
        for (i, &a) in planet_total_acc.iter().enumerate() {
            acc[i + 1] = a;
        }

        let masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
        let signal = HierarchySignal::classify(&masses);

        let wh_data = WhDenseData {
            mu,
            m_sun: m0,
            m_total,
            v_com,
            r0_sun_rest: r0_in,
            m_q_in,
            q0_helio_rest,
            v0_inertial_rest,
            planet_masses,
        };
        let step_snapshot = DenseSnapshot {
            t0: 0.0,
            dt,
            x0: pre_x0_inertial,
            v0: pre_v0_inertial,
            a0: pre_a0_inertial,
            b: Vec::new(),
            kind: IntegratorKind::WHFast,
            wh_data: Some(wh_data),
        };

        StepResult {
            consumed_dt: dt,
            potential_energy: pe,
            used_fallback: false,
            step_snapshot: Some(step_snapshot),
            degraded: false,
            hierarchy_signal: Some(signal),
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::WHFast
    }

    fn resume_state(&self) -> Vec<u8> {
        whfast_resume::encode(&self.cs_pos, &self.cs_vel)
    }

    fn restore_resume_state(
        &mut self,
        bytes: &[u8],
    ) -> Result<(), crate::physics::integrator::traits::ResumeError> {
        let (cs_pos, cs_vel) = whfast_resume::decode(bytes)?;
        self.cs_pos = cs_pos;
        self.cs_vel = cs_vel;
        Ok(())
    }
}

mod whfast_resume {
    use crate::math::Vec3;
    use crate::physics::integrator::traits::ResumeError;

    /// Layout: `magic(b"WHF")` ‖ `version(u8 = 1)` ‖ `n(u32 LE)` ‖
    /// per-body `[cs_pos.{x,y,z}, cs_vel.{x,y,z}]` as f64 LE.
    const MAGIC: &[u8; 3] = b"WHF";
    const VERSION: u8 = 1;

    pub fn encode(cs_pos: &[Vec3], cs_vel: &[Vec3]) -> Vec<u8> {
        debug_assert_eq!(cs_pos.len(), cs_vel.len());
        let n = cs_pos.len();
        let mut out = Vec::with_capacity(3 + 1 + 4 + n * 48);
        out.extend_from_slice(MAGIC);
        out.push(VERSION);
        out.extend_from_slice(&(n as u32).to_le_bytes());
        for (p, v) in cs_pos.iter().zip(cs_vel.iter()) {
            for c in [p.x, p.y, p.z, v.x, v.y, v.z] {
                out.extend_from_slice(&c.to_le_bytes());
            }
        }
        out
    }

    pub fn decode(bytes: &[u8]) -> Result<(Vec<Vec3>, Vec<Vec3>), ResumeError> {
        if bytes.len() < 8 || &bytes[..3] != MAGIC || bytes[3] != VERSION {
            return Err(ResumeError::UnsupportedFormat);
        }
        let n = u32::from_le_bytes(bytes[4..8].try_into().unwrap()) as usize;
        let needed = 8 + n * 48;
        if bytes.len() < needed {
            return Err(ResumeError::Truncated);
        }
        let mut cs_pos = Vec::with_capacity(n);
        let mut cs_vel = Vec::with_capacity(n);
        let mut off = 8;
        for _ in 0..n {
            let read = |off: usize| f64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
            cs_pos.push(Vec3::new(read(off), read(off + 8), read(off + 16)));
            cs_vel.push(Vec3::new(read(off + 24), read(off + 32), read(off + 40)));
            off += 48;
        }
        Ok((cs_pos, cs_vel))
    }
}

/// Half-kick on planet velocities, compensated against `cs_vel`.
fn whfast_kick(
    bodies: &mut [Body],
    ctx: &mut IntegratorContext<'_>,
    dt: f64,
    acc: &mut Vec<Vec3>,
    mu: f64,
    cs_vel: &mut [Vec3],
) -> f64 {
    let m0 = bodies[0].mass;

    let raw_pe = evaluate(&bodies[1..], ctx.force, acc);

    let bary_acc_raw =
        acc.iter().zip(bodies[1..].iter()).fold(Vec3::ZERO, |a, (a_i, b)| a + b.mass * *a_i);
    let indirect_raw = -bary_acc_raw / m0;

    let pe_inter = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
    let pe_central = central_potential(bodies, mu);

    let n = bodies.len();
    let mut pert_acc: Vec<Vec3> = vec![Vec3::ZERO; n];
    accumulate_perturbation_forces(
        bodies,
        &mut pert_acc,
        ctx.hamiltonian_perturbations,
        ctx.non_conservative_perturbations,
    );

    let indirect = indirect_raw * ctx.g_factor;
    for (i, ai) in acc.iter().enumerate() {
        let body_idx = i + 1;
        let kick = (*ai + pert_acc[body_idx] + indirect) * dt;
        let b = &mut bodies[body_idx];
        compensated_add(&mut b.vel_x, &mut b.vel_y, &mut b.vel_z, &mut cs_vel[body_idx], kick);
    }

    pe_inter + pe_central
}

/// Central Keplerian potential `−μ Σ m_i / |q_i|` in the
/// heliocentric frame.
fn central_potential(bodies: &[Body], mu: f64) -> f64 {
    bodies[1..]
        .iter()
        .map(|b| {
            let r = (b.pos_x * b.pos_x + b.pos_y * b.pos_y + b.pos_z * b.pos_z).sqrt().max(1e-30);
            -mu * b.mass / r
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::integrator::WisdomHolman;
    use crate::physics::integrator::force_model::GravityForceModel;

    /// Sun + 2 widely-separated planets on circular Keplerian orbits at
    /// G·M = 1 in canonical units. Mirrors `mercurius::tests`.
    fn quiet_planetary() -> Vec<Body> {
        vec![
            Body::star(1.0),
            Body::rocky(1.0e-6).at(1.0, 0.0).with_velocity(0.0, 1.0),
            Body::rocky(1.0e-6).at(2.0, 0.0).with_velocity(0.0, std::f64::consts::FRAC_1_SQRT_2),
        ]
    }

    /// Sun + Earth-mass test particle on a circular orbit at r = 1.
    /// Period is 2π, so N orbits = N · 2π in time.
    fn two_body_circular() -> Vec<Body> {
        vec![Body::star(1.0), Body::rocky(1.0e-9).at(1.0, 0.0).with_velocity(0.0, 1.0)]
    }

    fn step_via(integrator: &mut dyn Integrator, bodies: &mut [Body], dt: f64, n_steps: usize) {
        let mut force = GravityForceModel::new(0.5, 16);
        let mut acc: Vec<Vec3> = vec![Vec3::ZERO; bodies.len()];
        let hamiltonian: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> = Vec::new();
        let non_conservative: Vec<Box<dyn crate::physics::integrator::NonConservativeOperator>> =
            Vec::new();
        let mut observers: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        for _ in 0..n_steps {
            let mut ctx = IntegratorContext {
                force: &mut force,
                g_factor: 1.0,
                hamiltonian_perturbations: &hamiltonian,
                non_conservative_perturbations: &non_conservative,
                observers: &mut observers,
            };
            integrator.step(bodies, &mut ctx, dt, &mut acc);
        }
    }

    fn total_energy(bs: &[Body]) -> f64 {
        let mut ke = 0.0;
        let mut pe = 0.0;
        for (i, b) in bs.iter().enumerate() {
            ke += 0.5 * b.mass * (b.vel_x.powi(2) + b.vel_y.powi(2) + b.vel_z.powi(2));
            for j in (i + 1)..bs.len() {
                let dx = b.pos_x - bs[j].pos_x;
                let dy = b.pos_y - bs[j].pos_y;
                let dz = b.pos_z - bs[j].pos_z;
                let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1.0e-30);
                pe -= b.mass * bs[j].mass / r;
            }
        }
        ke + pe
    }

    /// `is_suitable_for` mirrors the WH 1991 dominance threshold.
    #[test]
    fn is_suitable_for_matches_wh_dominance() {
        let solar_like = vec![Body::star(1.0), Body::rocky(3e-6), Body::gas_giant(1e-3)];
        assert!(WHFast::is_suitable_for(&solar_like));

        let comparable = vec![Body::star(1.0), Body::star(0.5)];
        assert!(!WHFast::is_suitable_for(&comparable));

        let single = vec![Body::star(1.0)];
        assert!(!WHFast::is_suitable_for(&single));
    }

    /// Compensators resize lazily on first step / N change.
    #[test]
    fn ensure_state_size_resizes_compensators() {
        let mut wh = WHFast::new();
        assert!(wh.cs_pos.is_empty());
        wh.ensure_state_size(3);
        assert_eq!(wh.cs_pos.len(), 3);
        assert_eq!(wh.cs_vel.len(), 3);
        assert!(wh.cs_pos.iter().all(|v| *v == Vec3::ZERO));
        assert!(wh.cs_vel.iter().all(|v| *v == Vec3::ZERO));
        wh.ensure_state_size(5);
        assert_eq!(wh.cs_pos.len(), 5);
    }

    /// Two-body Kepler closure on a circular orbit. After N orbits the
    /// heliocentric position should return to its starting state up to
    /// the integrator's truncation floor. Closure is checked in the
    /// central-body-relative frame because the inertial planet position
    /// carries the COM drift `v_com · t_total ≈ 6e-7` over 100 orbits.
    ///
    /// The DH Kepler kernel propagates each planet around `mu = G·M_central`,
    /// not `G·M_total`; the indirect-drift Hamiltonian compensates over a
    /// full orbit but the residual is `O((dt/T)² · m_p/m_0)` in radial
    /// position and `O((dt/T) · m_p/m_0 · t)` in along-track phase. At
    /// `N = 100`, `steps_per_orbit = 1024`, `m_p/m_0 = 1e-9`, the dominant
    /// drift is the along-track phase: `~ 1 µAU` after 100 orbits.
    /// Eccentricity and semi-major axis stay bound by the symplectic
    /// invariant — only the orbital phase walks.
    ///
    /// Bound `1e-5`: catches a broken Kepler step, missed indirect drift,
    /// or a frame-leak; loose enough to track the published `O(N · dt²)`
    /// phase residual.
    #[test]
    fn two_body_closure_returns_after_n_orbits() {
        let n_orbits = 100;
        let steps_per_orbit = 1024;
        let dt = 2.0 * std::f64::consts::PI / steps_per_orbit as f64;
        let bodies0 = two_body_circular();
        let mut bodies = bodies0.clone();
        let mut wh = WHFast::new();
        step_via(&mut wh, &mut bodies, dt, n_orbits * steps_per_orbit);

        let q0 = Vec3::new(
            bodies0[1].pos_x - bodies0[0].pos_x,
            bodies0[1].pos_y - bodies0[0].pos_y,
            bodies0[1].pos_z - bodies0[0].pos_z,
        );
        let q = Vec3::new(
            bodies[1].pos_x - bodies[0].pos_x,
            bodies[1].pos_y - bodies[0].pos_y,
            bodies[1].pos_z - bodies[0].pos_z,
        );
        let dr = (q - q0).length();
        let dr_radial = (q.length() - q0.length()).abs();
        assert!(
            dr < 1.0e-5,
            "two-body closure after {n_orbits} orbits: |Δr| = {dr:.3e}, expected < 1e-5"
        );
        assert!(
            dr_radial < 1.0e-10,
            "two-body radial closure after {n_orbits} orbits: |Δr_rad| = {dr_radial:.3e}, expected < 1e-10"
        );
    }

    /// Brouwer-style energy conservation on a quiet planetary system.
    /// At `N ≤ 10⁴` both WHFast and WH 1991 sit at the IEEE-754 truncation
    /// floor; the test guards against gross algorithmic regressions
    /// (kick / drift / indirect-shift ordering, COM drift, frame leak)
    /// rather than measuring the compensated-summation advantage. The
    /// `O(√N · ε)` vs `O(N · ε)` separation lives in the cross-implementation
    /// lab notebook at `N ≳ 10⁸`.
    #[test]
    fn quiet_system_energy_conservation() {
        let dt = 1.0e-3;
        let n_steps = 10_000;
        let mut bodies = quiet_planetary();
        let e0 = total_energy(&bodies);

        let mut wh = WHFast::new();
        step_via(&mut wh, &mut bodies, dt, n_steps);

        let e1 = total_energy(&bodies);
        let rel = ((e1 - e0) / e0).abs();
        assert!(
            rel < 1.0e-10,
            "WHFast quiet-system energy drift over {n_steps} steps: |ΔE/E| = {rel:.3e}, expected < 1e-10"
        );
    }

    /// At small step counts WHFast and WH 1991 share the same KDK
    /// structure, drift kernel (`kepler_step`), and DH frame change;
    /// they differ only in the compensated-summation accumulator. At
    /// `N = 200` the per-step compensator updates are below the f64
    /// ULP of the running positions, so trajectories agree to the
    /// truncation floor. Tightening below `1e-10` would catch a
    /// compensator-ordering bug; loosening above `1e-8` would miss a
    /// frame-leak regression.
    #[test]
    fn whfast_matches_wisdom_holman_at_short_horizon() {
        let dt = 1.0e-3;
        let n_steps = 200;
        let mut bodies_whfast = quiet_planetary();
        let mut bodies_wh = quiet_planetary();

        let mut whfast = WHFast::new();
        let mut wh = WisdomHolman::new();
        step_via(&mut whfast, &mut bodies_whfast, dt, n_steps);
        step_via(&mut wh, &mut bodies_wh, dt, n_steps);

        for (bw, bf) in bodies_whfast.iter().zip(bodies_wh.iter()) {
            let dx = bw.pos_x - bf.pos_x;
            let dy = bw.pos_y - bf.pos_y;
            let dz = bw.pos_z - bf.pos_z;
            let r = (bw.pos_x.powi(2) + bw.pos_y.powi(2) + bw.pos_z.powi(2)).sqrt().max(1.0e-30);
            let rel = (dx * dx + dy * dy + dz * dz).sqrt() / r;
            assert!(
                rel < 1.0e-9,
                "WHFast vs WisdomHolman short-horizon parity: |Δr|/r = {rel:.3e}, expected < 1e-9"
            );
        }
    }

    /// Refusing a non-hierarchical scenario is the integrator's contract;
    /// `is_suitable_for` is the gate, but the per-step path also has to
    /// signal violation when called on an inadequate system. Validates
    /// the `HierarchySignal::Violated` arm rather than panic-on-bad-input.
    #[test]
    fn step_signals_hierarchy_on_short_input() {
        let mut bodies = vec![Body::star(1.0)];
        let mut force = GravityForceModel::new(0.5, 16);
        let mut acc: Vec<Vec3> = vec![Vec3::ZERO; bodies.len()];
        let hamiltonian: Vec<Box<dyn crate::physics::integrator::HamiltonianOperator>> = Vec::new();
        let non_conservative: Vec<Box<dyn crate::physics::integrator::NonConservativeOperator>> =
            Vec::new();
        let mut observers: Vec<Box<dyn crate::physics::integrator::Operator>> = Vec::new();
        let mut ctx = IntegratorContext {
            force: &mut force,
            g_factor: 1.0,
            hamiltonian_perturbations: &hamiltonian,
            non_conservative_perturbations: &non_conservative,
            observers: &mut observers,
        };
        let mut wh = WHFast::new();
        let result = wh.step(&mut bodies, &mut ctx, 0.01, &mut acc);
        assert_eq!(result.hierarchy_signal, Some(HierarchySignal::Violated));
    }
}
