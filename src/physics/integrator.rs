//! Symplectic N-body integrators.
//!
//! # Available integrators
//!
//! | Variant | Order | Force evals/step | Notes |
//! |---------|-------|-----------------|-------|
//! | [`Integrator::VelocityVerlet`] | 2nd | 2 (amortised 1) | Standard leapfrog KDK |
//! | [`Integrator::Yoshida4`]       | 4th | 3               | Forest–Ruth / Yoshida (1990) |
//!
//! # Choosing an integrator
//!
//! **Velocity Verlet** — phase error ∝ dt². Good for exploratory runs and
//! real-time visualisation. The symplectic property prevents secular energy
//! drift but phase accuracy degrades for long integrations (> ~100 orbits).
//!
//! **Yoshida 4th-order** — phase error ∝ dt⁴. Three times more force
//! evaluations per step, but you can use a 5–10× larger dt for the same energy
//! conservation. Recommended for publication-quality runs and long-term
//! stability studies.
//!
//! # Yoshida-4 coefficients (Forest–Ruth)
//!
//! ```text
//! cbrt2  = 2^(1/3)
//! w₁     = 1 / (2 − cbrt2)  ≈  1.3512071919596578
//! w₀     = 1 − 2 w₁          ≈ −1.7024143839193156
//!
//! Drift coefficients:  c = [w₁/2,  (w₁+w₀)/2,  (w₀+w₁)/2,  w₁/2]
//! Kick  coefficients:  d = [w₁,     w₀,          w₁        ]
//! ```
//!
//! The middle kick coefficient w₀ is negative — the second sub-step runs
//! "backwards" in time. This is correct and essential to the 4th-order
//! cancellation.
//!
//! # References
//! - Verlet (1967). *Phys. Rev.* 159, 98.
//! - Forest & Ruth (1990). *Nucl. Instrum. Methods Phys. Res.* A 290, 395–400.
//! - Yoshida (1990). *Phys. Lett. A* 150, 262–268.

use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;

// ── Integrator enum ───────────────────────────────────────────────────────────

/// Selects the symplectic integration algorithm used by [`System::step`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Integrator {
    /// Velocity Verlet (leapfrog KDK) — 2nd-order symplectic.
    VelocityVerlet,

    /// Yoshida / Forest–Ruth composition — 4th-order symplectic.
    Yoshida4,
}

impl Integrator {
    /// Short human-readable label used in the UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "Velocity Verlet (2nd)",
            Self::Yoshida4 => "Yoshida 4th-order",
        }
    }

    /// Formal convergence order in the time step.
    pub fn order(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
        }
    }

    /// Number of force evaluations consumed per full time step.
    ///
    /// VV's two evaluations can be amortised to one by sharing the endpoint
    /// across consecutive steps, but this implementation keeps them explicit
    /// for simplicity and restart-safety.
    pub fn force_evals_per_step(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 3,
        }
    }

    /// One-line description shown in the UI tooltip.
    pub fn description(self) -> &'static str {
        match self {
            Self::VelocityVerlet =>
                "2nd-order symplectic leapfrog. Fast; energy oscillates around \
                 the initial value. Phase error ∝ dt². Good for real-time \
                 visualisation and short integrations.",
            Self::Yoshida4 =>
                "4th-order symplectic composition (Forest–Ruth). 3 force evals \
                 per step but phase error ∝ dt⁴ — allows 5–10× larger dt for \
                 the same energy conservation. Required for publication-quality \
                 long-term runs.",
        }
    }
}

// ── Yoshida-4 (Forest–Ruth) coefficients ─────────────────────────────────────

/// 2^(1/3)
const CBRT2: f64 = 1.2599210498948732_f64;

/// w₁ = 1 / (2 − 2^(1/3))
pub const Y4_W1: f64 = 1.0_f64 / (2.0_f64 - CBRT2);

/// w₀ = 1 − 2 w₁  (negative — middle sub-step goes backwards)
pub const Y4_W0: f64 = 1.0_f64 - 2.0_f64 * Y4_W1;

/// Drift (position) coefficients: c[i] applied before force eval i, plus final drift.
pub const Y4_C: [f64; 4] = [
    Y4_W1 * 0.5,
    (Y4_W1 + Y4_W0) * 0.5,
    (Y4_W0 + Y4_W1) * 0.5, // == Y4_C[1] by symmetry
    Y4_W1 * 0.5,            // == Y4_C[0] by symmetry
];

/// Kick (velocity) coefficients: d[i] applied after force eval i.
pub const Y4_D: [f64; 3] = [Y4_W1, Y4_W0, Y4_W1];

// ── Primitive kernels ─────────────────────────────────────────────────────────

/// Rebuild the gravity structure and fill `scratch_acc` with accelerations.
///
/// Returns the raw (unscaled) gravitational potential energy.
pub fn evaluate_accelerations(
    bodies: &[Body],
    theta: f64,
    engine: &mut BarnesHutEngine,
    scratch_acc: &mut Vec<(f64, f64)>,
) -> f64 {
    if scratch_acc.len() != bodies.len() {
        scratch_acc.resize(bodies.len(), (0.0, 0.0));
    }
    engine.build(bodies);
    engine.evaluate(bodies, theta, scratch_acc)
}

/// Apply a velocity kick: `v += a · dt`.
///
/// Pass `0.5 * dt` for a half-kick (VV), or any scaled `w · dt` for
/// Yoshida sub-steps (including negative w for the middle sub-step).
pub fn kick(bodies: &mut [Body], acc: &[(f64, f64)], dt: f64) {
    for (body, &(ax, ay)) in bodies.iter_mut().zip(acc.iter()) {
        body.vx += ax * dt;
        body.vy += ay * dt;
    }
}

/// Thin alias kept for backward compatibility.
#[inline(always)]
pub fn half_kick(bodies: &mut [Body], acc: &[(f64, f64)], dt: f64) {
    kick(bodies, acc, dt);
}

/// Advance all positions using the current velocities: `x += v · dt`.
pub fn drift(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}
