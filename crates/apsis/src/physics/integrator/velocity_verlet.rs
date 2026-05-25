//! Velocity Verlet (leapfrog KDK) — 2nd-order symplectic integrator.
//!
//! Scheme: F(t) → kick(½dt) → drift(dt) → F(t+dt) → kick(½dt)
//!
//! The two half-kicks bracketing the drift are equivalent to a single
//! full kick at the midpoint, giving 2nd-order accuracy with one force
//! evaluation per amortised step.
//!
//! # References
//! - Verlet (1967). *Phys. Rev.* 159, 98.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::helpers::{evaluate, scale_acc_and_pe};
use crate::physics::integrator::operator_dispatch::accumulate_perturbation_forces;
use crate::physics::integrator::primitives::{drift, kick};
use crate::physics::integrator::traits::{
    Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Velocity Verlet (leapfrog KDK) — 2nd-order symplectic integrator.
pub struct VelocityVerlet;

impl Integrator for VelocityVerlet {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        let raw_pe = evaluate(bodies, ctx.force, acc);
        scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
        accumulate_perturbation_forces(
            bodies,
            acc,
            ctx.hamiltonian_perturbations,
            ctx.non_conservative_perturbations,
        );

        kick(bodies, acc, 0.5 * dt);
        drift(bodies, dt);

        let raw_pe = evaluate(bodies, ctx.force, acc);
        let pe = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
        accumulate_perturbation_forces(
            bodies,
            acc,
            ctx.hamiltonian_perturbations,
            ctx.non_conservative_perturbations,
        );

        kick(bodies, acc, 0.5 * dt);

        StepResult {
            consumed_dt: dt,
            potential_energy: pe,
            used_fallback: false,
            step_snapshot: None,
            degraded: false,
            hierarchy_signal: None,
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::VelocityVerlet
    }
}
