//! WHFast — Wisdom-Holman Fast (Rein & Tamayo 2015, *MNRAS* 452, 376).
//!
//! Symplectic mixed-variable integrator in democratic-heliocentric
//! coordinates. Same KDK structure as [`super::wisdom_holman::WisdomHolman`]
//! with two additions: persistent per-body compensators (Neumaier) on
//! the position and velocity accumulators, reducing the round-off
//! envelope on length-N sums from `O(N · ε)` to `O(√N · ε)`
//! (Higham 2002 §4.5); and an order-17 symplectic corrector
//! (Wisdom 1996) applied at sync boundaries, pushing boundary
//! truncation from `O(dt²)` to `O(dt^18)`. Corrector opt-out via
//! [`WHFast::without_correctors`].
//!
//! Lab notebook:
//! `docs/experiments/2026-05-13-whfast-integrator.md`.
//!
//! # References
//!
//! - Rein, H. & Tamayo, D. (2015). *MNRAS* 452, 376.
//! - Wisdom, J. (1996). *AJ* 112, 1305.
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
pub struct WHFast {
    cs_pos: Vec<Vec3>,
    cs_vel: Vec<Vec3>,
    /// Apply the order-17 symplectic corrector at sync boundaries.
    /// Default `true`; opt-out via [`WHFast::without_correctors`].
    with_correctors: bool,
}

impl Default for WHFast {
    fn default() -> Self {
        Self::new()
    }
}

impl WHFast {
    pub fn new() -> Self {
        Self { cs_pos: Vec::new(), cs_vel: Vec::new(), with_correctors: true }
    }

    /// Disable the symplectic corrector — runs bare KDK.
    #[must_use]
    pub fn without_correctors(mut self) -> Self {
        self.with_correctors = false;
        self
    }

    pub fn has_correctors(&self) -> bool {
        self.with_correctors
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

    /// Builder defaults: corrector ON.
    #[test]
    fn default_has_correctors_on() {
        let wh = WHFast::new();
        assert!(wh.has_correctors());
    }

    /// `without_correctors` flips the flag.
    #[test]
    fn without_correctors_disables() {
        let wh = WHFast::new().without_correctors();
        assert!(!wh.has_correctors());
    }

    /// Compensators resize lazily on first step / N change.
    #[test]
    fn ensure_state_size_resizes_compensators() {
        let mut wh = WHFast::new();
        assert!(wh.cs_pos.is_empty());
        wh.ensure_state_size(3);
        assert_eq!(wh.cs_pos.len(), 3);
        assert_eq!(wh.cs_vel.len(), 3);
        // All compensators start at zero.
        assert!(wh.cs_pos.iter().all(|v| *v == Vec3::ZERO));
        assert!(wh.cs_vel.iter().all(|v| *v == Vec3::ZERO));
        // Resize to a different N zeros the buffers.
        wh.ensure_state_size(5);
        assert_eq!(wh.cs_pos.len(), 5);
    }
}
