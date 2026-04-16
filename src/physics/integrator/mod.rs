//! Symplectic N-body integrators.
//!
//! # Available integrators
//!
//! | Variant | Order | Force evals/step | Notes |
//! |---------|-------|-----------------|-------|
//! | [`Integrator::VelocityVerlet`] | 2nd | 2 (amortised 1) | Standard leapfrog KDK |
//! | [`Integrator::Yoshida4`]       | 4th | 3               | Forest–Ruth / Yoshida (1990) |
//! | [`Integrator::WisdomHolman`]   | 2nd | 1               | Keplerian + perturbation split |
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
//! **Wisdom–Holman** — mixed-variable symplectic map. Keplerian two-body
//! propagation is analytic (see [`kepler`]); only perturbations are stepped
//! numerically. Intended for hierarchical planetary systems.
//!
//! # Module layout
//!
//! - [`coefficients`] — Yoshida-4 composition constants.
//! - [`primitives`]   — `evaluate_accelerations`, `kick`, `drift` kernels.
//! - [`perturbation`] — public [`PerturbationForce`] extension trait.
//! - [`kepler`]       — universal-variable two-body propagator (WH core).
//!
//! # References
//! - Verlet (1967). *Phys. Rev.* 159, 98.
//! - Forest & Ruth (1990). *Nucl. Instrum. Methods Phys. Res.* A 290, 395–400.
//! - Yoshida (1990). *Phys. Lett. A* 150, 262–268.
//! - Wisdom & Holman (1991). *Astron. J.* 102, 1528–1538.

pub mod coefficients;
pub mod kepler;
pub mod perturbation;
pub mod primitives;

pub use coefficients::{Y4_C, Y4_D, Y4_W0, Y4_W1};
pub use kepler::kepler_step;
pub use perturbation::PerturbationForce;
pub use primitives::{drift, evaluate_accelerations, kick};

// ── Integrator enum ───────────────────────────────────────────────────────────

/// Selects the symplectic integration algorithm used by [`crate::core::system::System::step`].
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
                "Mixed-variable symplectic map. Keplerian two-body motion is \
                 solved analytically; perturbations are stepped numerically. \
                 Designed for hierarchical planetary systems."
            },
        }
    }
}
