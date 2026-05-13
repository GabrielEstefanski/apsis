//! Symplectic N-body integrators.
//!
//! # Available integrators
//!
//! | Struct | Order | Force evals/step | Notes |
//! |--------|-------|-----------------|-------|
//! | [`VelocityVerlet`]  | 2nd  | 2 (amortised 1) | Standard leapfrog KDK |
//! | [`Yoshida4`]        | 4th  | 4               | Forest–Ruth / Yoshida (1990) |
//! | [`WisdomHolman`]    | 2nd  | 1               | Keplerian + perturbation split |
//! | [`Ias15`]           | 15th | ~14 (adaptive)  | Gauss-Radau, Rein & Spiegel (2015) — **default** |
//!
//! # Architecture
//!
//! Integrators implement the [`Integrator`] trait and receive an
//! [`IntegratorContext`] carrying the force model and physical parameters.
//! This decouples the integration algorithm from both the force engine
//! (Barnes-Hut, direct O(N²), GPU, …) and the simulation orchestrator.
//!
//! [`IntegratorKind`] is a plain enum used for UI display, snapshot
//! serialisation, and `Metrics` — it carries no stepping logic.
//!
//! # Module layout
//!
//! - [`coefficients`]       — Yoshida-4 composition constants.
//! - [`primitives`]         — `kick`, `drift` kernels.
//! - [`operator`]           — public [`Operator`] / [`HamiltonianOperator`] / [`NonConservativeOperator`] extension traits.
//! - [`operator_dispatch`]  — helpers integrators call to apply operators at canonical positions.
//! - [`kepler`]             — universal-variable two-body propagator (WH core).
//! - [`force_model`]        — [`ForceModel`] trait + [`GravityForceModel`] wrapper.
//! - [`helpers`]            — shared `evaluate`, `scale_acc_and_pe` helpers.
//! - [`traits`]             — [`Integrator`] trait, [`IntegratorContext`], [`StepResult`], [`IntegratorKind`].
//! - [`velocity_verlet`], [`yoshida4`], [`wisdom_holman`], [`ias15`] — integrator implementations.
//!
//! # References
//! - Verlet (1967). *Phys. Rev.* 159, 98.
//! - Forest & Ruth (1990). *Nucl. Instrum. Methods Phys. Res.* A 290, 395–400.
//! - Yoshida (1990). *Phys. Lett. A* 150, 262–268.
//! - Wisdom & Holman (1991). *Astron. J.* 102, 1528–1538.
//! - Everhart (1985). *Dyn. Comets* 115, 185–202 (original RADAU15).
//! - Rein & Spiegel (2015). *MNRAS* 446, 1424–1437 (IAS15).

pub mod citation;
pub mod coefficients;
pub mod conservation;
pub mod dense;
pub mod force_model;
pub mod helpers;
pub mod ias15;
pub mod kepler;
pub mod mercurius;
pub mod operator;
pub mod operator_dispatch;
pub mod primitives;
pub mod regime;
pub mod traits;
pub mod velocity_verlet;
pub mod wisdom_holman;
pub mod yoshida4;

// ── Re-exports ────────────────────────────────────────────────────────────────

pub use citation::{Citation, render_provenance};
pub use coefficients::{Y4_C, Y4_D, Y4_W0, Y4_W1};
pub use conservation::{
    ConservationClass, ConservationContribution, ConservationReport, EnergyImpact, OperatorRole,
    PotentialStatus,
};
pub use dense::{DenseCoeffs, DenseSnapshot};
pub use force_model::{ForceModel, GravityForceModel};
pub use helpers::{evaluate, scale_acc_and_pe};
pub use ias15::Ias15;
pub use kepler::kepler_step;
pub use mercurius::Mercurius;
pub use operator::{
    HamiltonianOperator, HamiltonianOperatorDescriptor, NonConservativeOperator, Operator,
    Potential, UnitSystemMismatch,
};
pub use primitives::{drift, kick};
pub use regime::{RegimeViolation, Severity};
pub use traits::{AdaptiveStats, Integrator, IntegratorContext, IntegratorKind, StepResult};
pub use velocity_verlet::VelocityVerlet;
pub use wisdom_holman::WisdomHolman;
pub use yoshida4::Yoshida4;

// ── Factory ───────────────────────────────────────────────────────────────────

/// Create a boxed integrator from a [`IntegratorKind`] discriminant.
pub fn make_integrator(kind: IntegratorKind) -> Box<dyn Integrator> {
    match kind {
        IntegratorKind::VelocityVerlet => Box::new(VelocityVerlet),
        IntegratorKind::Yoshida4 => Box::new(Yoshida4),
        IntegratorKind::WisdomHolman => Box::new(WisdomHolman::new()),
        IntegratorKind::Ias15 => Box::new(Ias15::new()),
        IntegratorKind::Mercurius => Box::new(Mercurius::new()),
    }
}
