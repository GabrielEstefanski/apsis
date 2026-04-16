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

use crate::core::body::Body;
use crate::physics::gravity::BarnesHutEngine;

// ── Integrator enum ───────────────────────────────────────────────────────────

/// Selects the symplectic integration algorithm used by [`System::step`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Integrator {
    /// Velocity Verlet (leapfrog KDK) — 2nd-order symplectic.
    VelocityVerlet,

    /// Yoshida / Forest–Ruth composition — 4th-order symplectic.
    Yoshida4,

    /// Wisdom & Holman (1991) symplectic map for Keplerian systems with small perturbations.
    WisdomHolman,
}

impl Integrator {
    /// Short human-readable label used in the UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "Velocity Verlet (2nd)",
            Self::Yoshida4 => "Yoshida 4th-order",
            Self::WisdomHolman => "Wisdom–Holman (2nd, Keplerian)",
        }
    }

    /// Formal convergence order in the time step.
    pub fn order(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 2,
        }
    }

    /// Number of force evaluations consumed per full time step.
    ///
    /// VV's two evaluations can be amortised to one by sharing the endpoint
    /// across consecutive steps, but this implementation keeps them explicit
    /// for simplicity and restart-safety.
    ///
    /// Yoshida4 performs 4 evaluations: 3 for the Forest–Ruth composition plus
    /// 1 final evaluation at q(t+dt) after the closing drift.  Without the 4th
    /// evaluation `last_potential` is at q‴ (pre-drift), making the energy
    /// diagnostic O(dt) instead of O(dt⁴).  The extra call is cheap relative
    /// to the 5–10× larger dt that Yoshida4 tolerates.
    pub fn force_evals_per_step(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 1,
        }
    }

    /// One-line description shown in the UI tooltip.
    pub fn description(self) -> &'static str {
        match self {
            Self::VelocityVerlet => {
                "2nd-order symplectic leapfrog. Fast; energy oscillates around \
                 the initial value. Phase error ∝ dt². Good for real-time \
                 visualisation and short integrations."
            },
            Self::Yoshida4 => {
                "4th-order symplectic composition (Forest–Ruth). 4 force evals \
                 per step but phase error ∝ dt⁴ — allows 5–10× larger dt for \
                 the same energy conservation. Required for publication-quality \
                 long-term runs."
            },
            Self::WisdomHolman => {
                "Symplectic map for Keplerian systems with small perturbations. \
                 Not implemented yet."
            },
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
    Y4_W1 * 0.5,           // == Y4_C[0] by symmetry
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

/// Advance all positions using the current velocities: `x += v · dt`.
pub fn drift(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}

pub trait PerturbationForce: Send + Sync {
    /// Accumulates non-gravitational accelerations into `scratch_acc`.
    ///
    /// `scratch_acc[i]` corresponds to `bodies[i]`. Implementations must
    /// **add** to existing values, not overwrite, so multiple perturbations
    /// compose correctly.
    fn accumulate(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)]);

    /// Accumulates accelerations for a sub-slice of bodies starting at
    /// global index `offset` within `System::bodies`.
    ///
    /// Used by [`System::apply_perturbations_planets`] during the
    /// Wisdom–Holman sub-step, where `scratch_acc` covers only `bodies[1..]`
    /// and the global index of each entry is `local_index + offset`.
    ///
    /// The default implementation ignores `offset` and delegates to
    /// [`accumulate`] — correct for perturbations that derive params
    /// dynamically from the body slice rather than from a pre-indexed vec.
    fn accumulate_offset(&self, bodies: &[Body], scratch_acc: &mut [(f64, f64)], offset: usize) {
        let _ = offset;
        self.accumulate(bodies, scratch_acc);
    }
}
