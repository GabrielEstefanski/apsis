//! Yoshida 4th-order (Forest–Ruth DKD composition) — 4th-order symplectic.
//!
//! Scheme: drift(c₀) → F → kick(d₀) → drift(c₁) → F → kick(d₁)
//!       → drift(c₂) → F → kick(d₂) → drift(c₃) → F (consistency eval)
//!
//! The middle kick coefficient d₁ = w₀ ≈ −1.70 is negative, meaning the
//! second sub-step is a backward kick in time.  This is not a bug — it is
//! the mechanism by which leading error terms cancel to achieve 4th-order.
//!
//! # References
//! - Forest & Ruth (1990). *Nucl. Instrum. Methods Phys. Res.* A 290, 395–400.
//! - Yoshida (1990). *Phys. Lett. A* 150, 262–268.

use crate::domain::body::Body;
use crate::physics::integrator::coefficients::{Y4_C, Y4_D};
use crate::physics::integrator::helpers::{apply_perturbations, evaluate, scale_acc_and_pe};
use crate::physics::integrator::primitives::{drift, kick};
use crate::physics::integrator::traits::{
    Integrator, IntegratorContext, IntegratorKind, StepResult,
};

/// Yoshida / Forest–Ruth 4th-order symplectic composition.
pub struct Yoshida4;

impl Integrator for Yoshida4 {
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<(f64, f64)>,
    ) -> StepResult {
        // Three DKD sub-steps with Yoshida coefficients.
        for i in 0..3 {
            drift(bodies, Y4_C[i] * dt);

            let raw_pe = evaluate(bodies, ctx.force, acc);
            scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
            apply_perturbations(bodies, acc, ctx.perturbations);

            kick(bodies, acc, Y4_D[i] * dt);
        }

        // Final drift to complete the step.
        drift(bodies, Y4_C[3] * dt);

        // ── Consistent energy snapshot ────────────────────────────────────
        // After the final drift the phase-space point is (q(t+dt), v(t+dt)).
        // `last_potential` still holds PE(q‴) — the potential at the positions
        // BEFORE the drift — which is inconsistent with the current body state.
        // Without this correction, energy diagnostics report O(dt) error
        // instead of O(dt⁴), making Y4 appear worse than VV in the UI.
        //
        // Re-evaluating the potential at q(t+dt) costs one additional force
        // call per step (3 → 4 total).  The accelerations are also updated
        // so that `acc` is consistent with the final positions.
        let raw_pe = evaluate(bodies, ctx.force, acc);
        let pe = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);

        StepResult { potential_energy: pe, used_fallback: false }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Yoshida4
    }
}
