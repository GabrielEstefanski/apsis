//! Implicit Midpoint — A-stable symplectic 2nd-order integrator.
//!
//! Single-stage Gauss-Legendre. Advances `y_{n+1} = y_n + dt · f((y_n +
//! y_{n+1}) / 2)`, solved iteratively. Symplectic, time-symmetric,
//! A-stable. No central-mass dominance assumption — accepts any topology.
//!
//! [`Solver::Newton`] is reserved in the v1 API; only `Picard` is
//! implemented.
//!
//! A-stable, not L-stable: `R(z) → -1` as `Re(z) → -∞`. Stiff modes
//! oscillate. Dissipation-dominant regimes need Radau IIA / BDF / TR-BDF2.
//!
//! Lab notebook: `docs/experiments/2026-05-14-implicit-midpoint-integrator.md`.
//!
//! # References
//!
//! - Hairer, Lubich, Wanner (2006). *Geometric Numerical Integration*,
//!   2nd ed. Springer. Chapters II.1.4, VI.1.
//! - Hairer, Wanner (1996). *Solving ODE II: Stiff Problems*, 2nd ed.
//!   Springer. §IV.6, §IV.8.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::helpers::{evaluate, scale_acc_and_pe};
use crate::physics::integrator::operator_dispatch::accumulate_perturbation_forces;
use crate::physics::integrator::traits::{
    AdaptiveStats, Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Iteration solver for the implicit midpoint equation. Fixed at
/// construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Solver {
    /// Fixed-point on the midpoint state. Converges for non-stiff and
    /// mildly-stiff conservative dynamics.
    #[default]
    Picard,
    /// Reserved in the v1 API; not yet implemented.
    /// Selecting this variant in v1 panics on first step.
    Newton,
}

/// Outcome of a single iteration loop within an integration step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IterationOutcome {
    Converged { iterations: u32, residual: f64 },
    MaxIterExhausted { residual: f64 },
}

const DEFAULT_MAX_ITERATIONS: u32 = 10;
const DEFAULT_TOLERANCE: f64 = 8.0 * f64::EPSILON;

/// Implicit Midpoint integrator. Carries scratch buffers, the chosen
/// solver / `max_iterations` / `tolerance`, and cumulative diagnostic
/// counters surfaced through [`AdaptiveStats`]. Buffers resize lazily.
pub struct ImplicitMidpoint {
    solver: Solver,
    max_iterations: u32,
    tolerance: f64,
    q0: Vec<Vec3>,
    v0: Vec<Vec3>,
    q_prev: Vec<Vec3>,
    v_prev: Vec<Vec3>,
    v_avg: Vec<Vec3>,
    /// Accepted steps so far.
    cum_steps: u64,
    /// Sum of iteration counts across all step calls. Mean iteration
    /// count = `cum_iterations / cum_steps`.
    cum_iterations: u64,
    /// Number of steps that exhausted `max_iterations` without converging.
    /// Healthy runs see this stay at zero; nonzero is the trigger for
    /// switching to `Solver::Newton` once #117 lands.
    cum_max_iter_hits: u64,
}

impl Default for ImplicitMidpoint {
    fn default() -> Self {
        Self::new()
    }
}

impl ImplicitMidpoint {
    pub fn new() -> Self {
        Self {
            solver: Solver::default(),
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
            q0: Vec::new(),
            v0: Vec::new(),
            q_prev: Vec::new(),
            v_prev: Vec::new(),
            v_avg: Vec::new(),
            cum_steps: 0,
            cum_iterations: 0,
            cum_max_iter_hits: 0,
        }
    }

    pub fn with_solver(mut self, solver: Solver) -> Self {
        self.solver = solver;
        self
    }

    pub fn with_iteration_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance.max(f64::EPSILON);
        self
    }

    pub fn with_max_iterations(mut self, max_iterations: u32) -> Self {
        self.max_iterations = max_iterations.max(1);
        self
    }

    /// IM has no hierarchy assumption — accepts any ≥ 2-body system.
    pub fn is_suitable_for(bodies: &[Body]) -> bool {
        bodies.len() >= 2
    }

    fn ensure_state_size(&mut self, n: usize) {
        if self.q0.len() != n {
            self.q0.resize(n, Vec3::ZERO);
            self.v0.resize(n, Vec3::ZERO);
            self.q_prev.resize(n, Vec3::ZERO);
            self.v_prev.resize(n, Vec3::ZERO);
            self.v_avg.resize(n, Vec3::ZERO);
        }
    }
}

impl Integrator for ImplicitMidpoint {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        let n = bodies.len();
        if n < 2 {
            return StepResult {
                consumed_dt: dt,
                potential_energy: 0.0,
                used_fallback: false,
                step_snapshot: None,
                degraded: false,
                hierarchy_signal: None,
            };
        }

        self.ensure_state_size(n);

        for (i, b) in bodies.iter().enumerate() {
            self.q0[i] = Vec3::new(b.pos_x, b.pos_y, b.pos_z);
            self.v0[i] = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
        }

        let outcome = iterate_to_convergence(
            self.solver,
            self.max_iterations,
            self.tolerance,
            bodies,
            ctx,
            dt,
            acc,
            &self.q0,
            &self.v0,
            &mut self.q_prev,
            &mut self.v_prev,
            &mut self.v_avg,
        );

        let (iterations, degraded) = match outcome {
            IterationOutcome::Converged { iterations, .. } => (iterations as u64, false),
            IterationOutcome::MaxIterExhausted { residual } => {
                crate::warn_diag!(
                    crate::core::log::Source::Integrator,
                    "ImplicitMidpoint iteration did not converge within max_iterations",
                    max_iterations = self.max_iterations,
                    final_residual = residual,
                    tolerance = self.tolerance,
                    bodies = n,
                    dt = dt,
                );
                (self.max_iterations as u64, true)
            },
        };
        self.cum_steps = self.cum_steps.saturating_add(1);
        self.cum_iterations = self.cum_iterations.saturating_add(iterations);
        if degraded {
            self.cum_max_iter_hits = self.cum_max_iter_hits.saturating_add(1);
        }

        // Mid-iteration force evaluation lands at the midpoint state, not
        // the end state. Diagnostics (`System::step` → `Diagnostics::compute`)
        // and the next outer step both expect end-state acc. One refresh.
        let raw_pe = evaluate(bodies, ctx.force, acc);
        let pe = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
        accumulate_perturbation_forces(
            bodies,
            acc,
            ctx.hamiltonian_perturbations,
            ctx.non_conservative_perturbations,
        );

        StepResult {
            consumed_dt: dt,
            potential_energy: pe,
            used_fallback: false,
            step_snapshot: None,
            degraded,
            hierarchy_signal: None,
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::ImplicitMidpoint
    }

    /// Surface the iteration counters via the existing IAS15-flavoured
    /// fields: `substeps` = step count, `picard_iters` = total Picard
    /// iterations across all steps, `degraded` = max-iter exhaustions.
    /// IAS15-only fields stay at their `Default` zeros.
    fn adaptive_stats(&self) -> Option<AdaptiveStats> {
        Some(AdaptiveStats {
            substeps: self.cum_steps,
            picard_iters: self.cum_iterations,
            degraded: self.cum_max_iter_hits,
            ..AdaptiveStats::default()
        })
    }
}

/// Single iteration loop shared between Picard and (future #117) Newton.
/// Per-iteration update rule is the only branch.
#[allow(clippy::too_many_arguments)]
fn iterate_to_convergence(
    solver: Solver,
    max_iterations: u32,
    tolerance: f64,
    bodies: &mut [Body],
    ctx: &mut IntegratorContext<'_>,
    dt: f64,
    acc: &mut Vec<Vec3>,
    q0: &[Vec3],
    v0: &[Vec3],
    q_prev: &mut [Vec3],
    v_prev: &mut [Vec3],
    v_avg: &mut [Vec3],
) -> IterationOutcome {
    let mut residual = f64::INFINITY;
    for k in 0..max_iterations {
        match solver {
            Solver::Picard => picard_step(bodies, ctx, dt, acc, q0, v0, q_prev, v_prev, v_avg),
            Solver::Newton => unimplemented!(
                "ImplicitMidpoint Solver::Newton is reserved in the v1 API \
                 and not yet implemented"
            ),
        }
        residual = relative_state_delta(bodies, q_prev, v_prev);
        if residual < tolerance {
            return IterationOutcome::Converged { iterations: k + 1, residual };
        }
    }
    IterationOutcome::MaxIterExhausted { residual }
}

/// One Picard iteration: snapshot iterate, write midpoint to bodies,
/// evaluate, write next iterate.
#[allow(clippy::too_many_arguments)]
fn picard_step(
    bodies: &mut [Body],
    ctx: &mut IntegratorContext<'_>,
    dt: f64,
    acc: &mut Vec<Vec3>,
    q0: &[Vec3],
    v0: &[Vec3],
    q_prev: &mut [Vec3],
    v_prev: &mut [Vec3],
    v_avg: &mut [Vec3],
) {
    for (i, b) in bodies.iter().enumerate() {
        q_prev[i] = Vec3::new(b.pos_x, b.pos_y, b.pos_z);
        v_prev[i] = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
    }

    for (i, b) in bodies.iter_mut().enumerate() {
        b.pos_x = 0.5 * (q0[i].x + q_prev[i].x);
        b.pos_y = 0.5 * (q0[i].y + q_prev[i].y);
        b.pos_z = 0.5 * (q0[i].z + q_prev[i].z);
        b.vel_x = 0.5 * (v0[i].x + v_prev[i].x);
        b.vel_y = 0.5 * (v0[i].y + v_prev[i].y);
        b.vel_z = 0.5 * (v0[i].z + v_prev[i].z);
        v_avg[i] = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
    }

    let raw_pe = evaluate(bodies, ctx.force, acc);
    let _ = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
    accumulate_perturbation_forces(
        bodies,
        acc,
        ctx.hamiltonian_perturbations,
        ctx.non_conservative_perturbations,
    );

    for (i, b) in bodies.iter_mut().enumerate() {
        b.pos_x = q0[i].x + dt * v_avg[i].x;
        b.pos_y = q0[i].y + dt * v_avg[i].y;
        b.pos_z = q0[i].z + dt * v_avg[i].z;
        b.vel_x = v0[i].x + dt * acc[i].x;
        b.vel_y = v0[i].y + dt * acc[i].y;
        b.vel_z = v0[i].z + dt * acc[i].z;
    }
}

/// `‖y_k − y_{k-1}‖ / ‖y_k‖`, two-norm over positions + velocities.
fn relative_state_delta(bodies: &[Body], q_prev: &[Vec3], v_prev: &[Vec3]) -> f64 {
    let mut delta_sq = 0.0_f64;
    let mut state_sq = 0.0_f64;
    for (i, b) in bodies.iter().enumerate() {
        let q = Vec3::new(b.pos_x, b.pos_y, b.pos_z);
        let v = Vec3::new(b.vel_x, b.vel_y, b.vel_z);
        let dq = q - q_prev[i];
        let dv = v - v_prev[i];
        delta_sq += dq.length_squared() + dv.length_squared();
        state_sq += q.length_squared() + v.length_squared();
    }
    (delta_sq / state_sq.max(f64::EPSILON)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::integrator::WisdomHolman;
    use crate::physics::integrator::force_model::GravityForceModel;

    /// Sun + 2 widely-separated planets on circular Keplerian orbits at
    /// G·M = 1. Mirrors the `mercurius` / `whfast` test fixtures.
    fn quiet_planetary() -> Vec<Body> {
        vec![
            Body::star(1.0),
            Body::rocky(1.0e-6).at(1.0, 0.0).with_velocity(0.0, 1.0),
            Body::rocky(1.0e-6).at(2.0, 0.0).with_velocity(0.0, std::f64::consts::FRAC_1_SQRT_2),
        ]
    }

    /// Sun + Earth-mass test particle on a circular orbit at r = 1.
    fn two_body_circular() -> Vec<Body> {
        vec![Body::star(1.0), Body::rocky(1.0e-9).at(1.0, 0.0).with_velocity(0.0, 1.0)]
    }

    /// Equal-mass binary at separation `2`. Differentiator scenario for
    /// IM: WH/WHFast/Mercurius refuse this on hierarchy grounds; IM
    /// integrates it without ceremony.
    fn equal_mass_binary() -> Vec<Body> {
        vec![
            Body::rocky(1.0).at(-1.0, 0.0).with_velocity(0.0, -0.5),
            Body::rocky(1.0).at(1.0, 0.0).with_velocity(0.0, 0.5),
        ]
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

    /// IM accepts any ≥ 2-body topology, including configurations
    /// WH-class integrators reject.
    #[test]
    fn is_suitable_for_accepts_any_topology() {
        assert!(ImplicitMidpoint::is_suitable_for(&quiet_planetary()));
        assert!(ImplicitMidpoint::is_suitable_for(&equal_mass_binary()));
        assert!(!ImplicitMidpoint::is_suitable_for(&[Body::star(1.0)]));
    }

    #[test]
    fn builders_set_solver_tolerance_and_max_iter() {
        let im = ImplicitMidpoint::new()
            .with_solver(Solver::Picard)
            .with_iteration_tolerance(1.0e-12)
            .with_max_iterations(20);
        assert_eq!(im.solver, Solver::Picard);
        assert!((im.tolerance - 1.0e-12).abs() < 1.0e-30);
        assert_eq!(im.max_iterations, 20);
    }

    /// Sub-`f64::EPSILON` tolerance is unreachable; clamp it.
    #[test]
    fn tolerance_clamps_to_epsilon() {
        let im = ImplicitMidpoint::new().with_iteration_tolerance(0.0);
        assert!(im.tolerance >= f64::EPSILON);
    }

    /// Closure on a circular orbit. Along-track drifts at the IM2
    /// truncation floor `O(N · dt² · ω) ≈ 1.5e-2` for `(N, dt) =
    /// (100·1024, 2π/1024)`; radial stays tight because symplecticity
    /// preserves `(a, e)`.
    #[test]
    fn two_body_closure_returns_after_n_orbits() {
        let n_orbits = 100;
        let steps_per_orbit = 1024;
        let dt = 2.0 * std::f64::consts::PI / steps_per_orbit as f64;
        let bodies0 = two_body_circular();
        let mut bodies = bodies0.clone();
        let mut im = ImplicitMidpoint::new();
        step_via(&mut im, &mut bodies, dt, n_orbits * steps_per_orbit);

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
            dr < 5.0e-2,
            "two-body closure after {n_orbits} orbits: |Δr| = {dr:.3e}, expected < 5e-2"
        );
        assert!(
            dr_radial < 1.0e-8,
            "two-body radial closure after {n_orbits} orbits: |Δr_rad| = {dr_radial:.3e}, expected < 1e-8"
        );
    }

    /// Brouwer energy conservation on a quiet hierarchical system at
    /// `N = 10⁴`. Inline bound `1e-9` is looser than the lab-notebook
    /// Tier 2 bound (`1e-10` over `10⁶` steps) because at `N = 10⁴` the
    /// floor is iteration-tolerance-accumulated noise (`~ √N · ε_iter
    /// · |E|` ≈ `1e-13`) plus per-step truncation, not the secular
    /// scaling that dominates at `N = 10⁶`.
    #[test]
    fn quiet_system_energy_conservation() {
        let dt = 1.0e-3;
        let n_steps = 10_000;
        let mut bodies = quiet_planetary();
        let e0 = total_energy(&bodies);

        let mut im = ImplicitMidpoint::new();
        step_via(&mut im, &mut bodies, dt, n_steps);

        let e1 = total_energy(&bodies);
        let rel = ((e1 - e0) / e0).abs();
        assert!(
            rel < 1.0e-9,
            "ImplicitMidpoint quiet-system energy drift over {n_steps} steps: |ΔE/E| = {rel:.3e}, expected < 1e-9"
        );
    }

    /// Forward N + backward N closes to iteration tolerance. Inline
    /// bound `1e-10` is looser than the Tier 3 protocol bound (`1e-12`)
    /// because the round-trip residual scales `~ √(2N) · ε_iter · |q|`
    /// — at `N = 10³` this gives `~ 7e-15 · |q|`, well inside `1e-10`.
    /// Tighter bounds are reachable but require longer round-trips
    /// (Tier 3 uses `N = 10³` with cleaner controlled state).
    #[test]
    fn time_symmetry_round_trip_returns_to_initial() {
        let dt = 1.0e-3;
        let n_steps = 1_000;
        let bodies0 = quiet_planetary();
        let mut bodies = bodies0.clone();

        let mut im = ImplicitMidpoint::new();
        step_via(&mut im, &mut bodies, dt, n_steps);
        step_via(&mut im, &mut bodies, -dt, n_steps);

        for (i, (b, b0)) in bodies.iter().zip(bodies0.iter()).enumerate() {
            let dq = ((b.pos_x - b0.pos_x).powi(2)
                + (b.pos_y - b0.pos_y).powi(2)
                + (b.pos_z - b0.pos_z).powi(2))
            .sqrt();
            let dv = ((b.vel_x - b0.vel_x).powi(2)
                + (b.vel_y - b0.vel_y).powi(2)
                + (b.vel_z - b0.vel_z).powi(2))
            .sqrt();
            let q0_norm = (b0.pos_x.powi(2) + b0.pos_y.powi(2) + b0.pos_z.powi(2)).sqrt().max(1.0);
            let v0_norm = (b0.vel_x.powi(2) + b0.vel_y.powi(2) + b0.vel_z.powi(2)).sqrt().max(1.0);
            assert!(
                dq / q0_norm < 1.0e-10,
                "body {i}: round-trip position drift {:.3e} exceeds 1e-10",
                dq / q0_norm,
            );
            assert!(
                dv / v0_norm < 1.0e-10,
                "body {i}: round-trip velocity drift {:.3e} exceeds 1e-10",
                dv / v0_norm,
            );
        }
    }

    /// Equal-mass binary — refused by WH-class integrators, integrated
    /// by IM. Energy stays bound over 100 orbits (period `4π`).
    #[test]
    fn equal_mass_binary_stays_bound() {
        let period = 4.0 * std::f64::consts::PI;
        let dt = period / 1024.0;
        let n_steps = 1024 * 100;
        let mut bodies = equal_mass_binary();
        let e0 = total_energy(&bodies);

        let mut im = ImplicitMidpoint::new();
        step_via(&mut im, &mut bodies, dt, n_steps);

        let e1 = total_energy(&bodies);
        let rel = ((e1 - e0) / e0).abs();
        assert!(
            rel < 1.0e-8,
            "equal-mass binary energy drift over 100 orbits: |ΔE/E| = {rel:.3e}, expected < 1e-8"
        );
        // Liveness check — energy bound above is the actual gate; this
        // catches NaN/blow-up that the energy diff would mask.
        for (i, b) in bodies.iter().enumerate() {
            let r = (b.pos_x.powi(2) + b.pos_y.powi(2) + b.pos_z.powi(2)).sqrt();
            assert!(
                r.is_finite() && r < 10.0,
                "body {i}: |r| = {r} after binary integration — unbounded or non-finite"
            );
        }
    }

    /// Figure-eight choreography (Chenciner & Montgomery 2000) — three
    /// equal masses on a planar figure-8, no close encounters, no
    /// dominant primary. WH-class integrators refuse this via the
    /// hierarchy gate; IM integrates it. Period `T ≈ 6.3259`; bound
    /// `1e-6` over 10 periods catches algorithmic regressions in the
    /// non-hierarchy code path while staying inside IM2's truncation
    /// floor for the smooth choreography.
    #[test]
    fn figure_eight_three_body_stays_bound() {
        let r1 = (0.970_004_36, -0.243_087_53);
        let v3 = (-0.932_407_37, -0.864_731_46);
        let v_outer = (-v3.0 * 0.5, -v3.1 * 0.5);
        let bodies0 = vec![
            Body::rocky(1.0).at(r1.0, r1.1).with_velocity(v_outer.0, v_outer.1),
            Body::rocky(1.0).at(-r1.0, -r1.1).with_velocity(v_outer.0, v_outer.1),
            Body::rocky(1.0).at(0.0, 0.0).with_velocity(v3.0, v3.1),
        ];
        let mut bodies = bodies0.clone();
        let e0 = total_energy(&bodies);

        let period = 6.325_9;
        let n_periods = 10;
        let steps_per_period = 1024;
        let dt = period / steps_per_period as f64;
        let n_steps = n_periods * steps_per_period;
        let mut im = ImplicitMidpoint::new();
        step_via(&mut im, &mut bodies, dt, n_steps);

        let e1 = total_energy(&bodies);
        let rel = ((e1 - e0) / e0).abs();
        assert!(
            rel < 1.0e-6,
            "figure-8 energy drift over {n_periods} periods: |ΔE/E| = {rel:.3e}, expected < 1e-6"
        );
        for (i, b) in bodies.iter().enumerate() {
            assert!(
                b.pos_x.is_finite() && b.pos_y.is_finite() && b.pos_z.is_finite(),
                "body {i}: non-finite kinematics after figure-8 integration"
            );
        }
    }

    /// Cross-integrator drift IM-vs-WH: WH uses analytical Kepler, IM
    /// iterates numerically — they differ at `O(N · dt² · ω)`. Bound
    /// `1e-5` is a regression guard; observed `~7.5e-7`.
    #[test]
    fn matches_wisdom_holman_on_quiet_hierarchical_system() {
        let dt = 1.0e-3;
        let n_steps = 200;
        let mut bodies_im = quiet_planetary();
        let mut bodies_wh = quiet_planetary();

        let mut im = ImplicitMidpoint::new();
        let mut wh = WisdomHolman::new();
        step_via(&mut im, &mut bodies_im, dt, n_steps);
        step_via(&mut wh, &mut bodies_wh, dt, n_steps);

        for (bi, bw) in bodies_im.iter().zip(bodies_wh.iter()) {
            let dx = bi.pos_x - bw.pos_x;
            let dy = bi.pos_y - bw.pos_y;
            let dz = bi.pos_z - bw.pos_z;
            let r = (bi.pos_x.powi(2) + bi.pos_y.powi(2) + bi.pos_z.powi(2)).sqrt().max(1.0e-30);
            let rel = (dx * dx + dy * dy + dz * dz).sqrt() / r;
            assert!(
                rel < 1.0e-5,
                "IM vs WisdomHolman short-horizon parity: |Δr|/r = {rel:.3e}, expected < 1e-5"
            );
        }
    }

    /// Selecting `Solver::Newton` in v1 panics.
    #[test]
    #[should_panic(expected = "reserved in the v1 API")]
    fn newton_solver_is_unimplemented_in_v1() {
        let mut bodies = two_body_circular();
        let mut im = ImplicitMidpoint::new().with_solver(Solver::Newton);
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
        im.step(&mut bodies, &mut ctx, 1e-3, &mut acc);
    }

    /// Single body — no force pairs, no-op step.
    #[test]
    fn single_body_input_is_no_op() {
        let mut bodies = vec![Body::star(1.0)];
        let mut im = ImplicitMidpoint::new();
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
        let result = im.step(&mut bodies, &mut ctx, 1e-3, &mut acc);
        assert_eq!(result.consumed_dt, 1e-3);
        assert!(!result.degraded);
        assert!(result.hierarchy_signal.is_none());
    }
}
