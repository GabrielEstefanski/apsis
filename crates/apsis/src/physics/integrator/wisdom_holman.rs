//! Wisdom–Holman mixed-variable symplectic integrator (2nd-order).
//!
//! Scheme (heliocentric frame):
//!
//!   kick_pert(½dt)  →  drift_Kepler(dt)  →  kick_pert(½dt)
//!
//! The Hamiltonian is split as H = H_Kepler + H_pert.  H_Kepler is integrated
//! exactly via the analytic universal-variable propagator; H_pert contributes
//! velocity kicks that include the heliocentric indirect term (momentum
//! cross-term) required to preserve symplecticity.
//!
//! # Fallback
//!
//! When `bodies[0]` does not satisfy the dominance criterion (M₀ / Σ mᵢ < 10),
//! the integrator silently falls back to [`Yoshida4`] for that step.  This
//! produces physically correct (if slower) results rather than a silently
//! wrong trajectory.
//!
//! # References
//! - Wisdom & Holman (1991). *Astron. J.* 102, 1528–1538.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::helpers::{
    apply_perturbations_planets, evaluate, scale_acc_and_pe,
};
use crate::physics::integrator::kepler::kepler_step;
use crate::physics::integrator::traits::{
    Integrator, IntegratorContext, IntegratorKind, StepResult,
};
use crate::physics::integrator::yoshida4::Yoshida4;

/// Minimum ratio `M_central / Σ m_i (i > 0)` required for WH validity.
const WH_DOMINANCE_RATIO: f64 = 10.0;

/// Wisdom–Holman mixed-variable symplectic map for Keplerian systems.
pub struct WisdomHolman {
    /// Fallback integrator used when the dominance criterion fails.
    fallback: Yoshida4,
}

impl Default for WisdomHolman {
    fn default() -> Self {
        Self::new()
    }
}

impl WisdomHolman {
    pub fn new() -> Self {
        Self { fallback: Yoshida4 }
    }

    /// Returns `true` if `bodies[0]` dominates the system.
    ///
    /// Public so that `System` can expose this check to the UI without
    /// coupling to the WH internals.
    pub fn is_suitable_for(bodies: &[Body]) -> bool {
        if bodies.len() < 2 {
            return false;
        }
        let m0 = bodies[0].mass;
        let m_rest: f64 = bodies[1..].iter().map(|b| b.mass).sum();
        let max_other = bodies[1..].iter().map(|b| b.mass).fold(0.0_f64, f64::max);
        m0 >= max_other && m0 >= WH_DOMINANCE_RATIO * m_rest
    }

    /// Perturbation kick: evaluate inter-planetary forces, compute the
    /// heliocentric indirect term, and apply a velocity kick of `dt` to
    /// all planets.
    ///
    /// Returns the total potential (inter-planetary + central) scaled by
    /// `g_factor`.
    fn wh_kick(
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
        mu: f64,
    ) -> f64 {
        let total_m0 = bodies[0].mass;

        // Evaluate forces on planets only (bodies[1..]).
        let raw_pe = evaluate(&bodies[1..], ctx.force, acc);

        // Indirect term: a_indirect = −(Σ mⱼ aⱼ) / M₀
        // Computed from raw (pre-scaled) accelerations.
        let (ax_bary_raw, ay_bary_raw, az_bary_raw) =
            acc.iter().zip(bodies[1..].iter()).fold((0.0_f64, 0.0_f64, 0.0_f64), |a, (ai, b)| {
                (a.0 + b.mass * ai.x, a.1 + b.mass * ai.y, a.2 + b.mass * ai.z)
            });

        let indirect_x_raw = -ax_bary_raw / total_m0;
        let indirect_y_raw = -ay_bary_raw / total_m0;
        let indirect_z_raw = -az_bary_raw / total_m0;

        // Scale gravitational accelerations and compute central potential.
        let pe_inter = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
        let pe_central = central_potential(bodies, mu);
        let potential = pe_inter + pe_central;

        // Non-gravitational perturbations on planets.
        apply_perturbations_planets(&bodies[1..], acc, ctx.perturbations);

        // Apply kick with scaled indirect term.
        let indirect_x = indirect_x_raw * ctx.g_factor;
        let indirect_y = indirect_y_raw * ctx.g_factor;
        let indirect_z = indirect_z_raw * ctx.g_factor;
        for (i, ai) in acc.iter().enumerate() {
            bodies[i + 1].vx += (ai.x + indirect_x) * dt;
            bodies[i + 1].vy += (ai.y + indirect_y) * dt;
            bodies[i + 1].vz += (ai.z + indirect_z) * dt;
        }

        potential
    }
}

impl Integrator for WisdomHolman {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        // Fallback to Yoshida4 if the system is not hierarchical.
        if !Self::is_suitable_for(bodies) {
            let mut result = self.fallback.step(bodies, ctx, dt, acc);
            result.used_fallback = true;
            return result;
        }

        let mu = ctx.g_factor * bodies[0].mass;
        let total_m0 = bodies[0].mass;

        // ── To heliocentric frame ────────────────────────────────────────
        let (cx0, cy0, cvx0, cvy0) = (bodies[0].x, bodies[0].y, bodies[0].vx, bodies[0].vy);
        for b in &mut bodies[1..] {
            b.x -= cx0;
            b.y -= cy0;
            b.vx -= cvx0;
            b.vy -= cvy0;
        }

        // ── First half-kick ──────────────────────────────────────────────
        let _ = Self::wh_kick(bodies, ctx, 0.5 * dt, acc, mu);

        // ── Exact Keplerian drift ────────────────────────────────────────
        for b in bodies[1..].iter_mut() {
            let r0 = Vec3::new(b.x, b.y, b.z);
            let v0 = Vec3::new(b.vx, b.vy, b.vz);
            let (r1, v1) = kepler_step(r0, v0, dt, mu);
            b.x = r1.x;
            b.y = r1.y;
            b.z = r1.z;
            b.vx = v1.x;
            b.vy = v1.y;
            b.vz = v1.z;
        }

        // ── Second half-kick ─────────────────────────────────────────────
        let pe = Self::wh_kick(bodies, ctx, 0.5 * dt, acc, mu);

        // ── Back to inertial (barycentric) frame ─────────────────────────
        let (px, py) = bodies[1..]
            .iter()
            .fold((0.0_f64, 0.0_f64), |(px, py), b| (px + b.mass * b.vx, py + b.mass * b.vy));

        bodies[0].vx = -px / total_m0;
        bodies[0].vy = -py / total_m0;
        bodies[0].x += bodies[0].vx * dt;
        bodies[0].y += bodies[0].vy * dt;

        let (cx1, cy1, cvx1, cvy1) = (bodies[0].x, bodies[0].y, bodies[0].vx, bodies[0].vy);
        for b in &mut bodies[1..] {
            b.x += cx1;
            b.y += cy1;
            b.vx += cvx1;
            b.vy += cvy1;
        }

        StepResult {
            consumed_dt: dt,
            potential_energy: pe,
            used_fallback: false,
            step_snapshot: None,
            degraded: false,
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::WisdomHolman
    }
}

/// Central Keplerian potential: −μ Σ mᵢ / rᵢ  (heliocentric).
fn central_potential(bodies: &[Body], mu: f64) -> f64 {
    bodies[1..]
        .iter()
        .map(|b| {
            let r = (b.x * b.x + b.y * b.y).sqrt().max(1e-30);
            -mu * b.mass / r
        })
        .sum()
}
