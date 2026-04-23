//! Core integrator abstraction: trait, context, result, and kind enum.
//!
//! # Architecture
//!
//! ```text
//! ┌────────────┐     ┌──────────────────┐     ┌───────────┐
//! │ Integrator │────▶│ IntegratorContext │────▶│ ForceModel│
//! │  (trait)   │     │  (force+params)  │     │  (trait)  │
//! └────────────┘     └──────────────────┘     └───────────┘
//!       │
//!       ▼
//! ┌────────────┐
//! │ StepResult │  ← returned after each integration step
//! └────────────┘
//! ```
//!
//! The [`Integrator`] trait replaces the old `Integrator` enum, enabling
//! new integration schemes to be added without touching the core.
//!
//! [`IntegratorKind`] is a plain enum retained for UI display, snapshot
//! serialisation, and `Metrics`.  It is **not** used for dispatch.

use crate::domain::body::Body;
use crate::physics::integrator::force_model::ForceModel;
use crate::physics::integrator::perturbation::PerturbationForce;

// ── IntegratorKind (serialisable enum) ────────────────────────────────────────

/// Identifies an integration algorithm without carrying behaviour.
///
/// Used for snapshot serialisation, UI combo-boxes, and `Metrics`.
/// The actual stepping logic lives in structs that implement [`Integrator`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegratorKind {
    VelocityVerlet,
    Yoshida4,
    WisdomHolman,
    Ias15,
}

impl IntegratorKind {
    /// Short human-readable label used in the UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "Velocity Verlet (2nd)",
            Self::Yoshida4 => "Yoshida 4th-order",
            Self::WisdomHolman => "Wisdom–Holman (2nd, Keplerian)",
            Self::Ias15 => "IAS15 (15th, adaptive)",
        }
    }

    /// Formal convergence order in the time step.
    pub fn order(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 2,
            Self::Ias15 => 15,
        }
    }

    /// Nominal number of force evaluations consumed per full time step.
    ///
    /// IAS15 is adaptive: each accepted substep uses ~7 evals with ~2 Picard
    /// iterations; the quoted number is an amortised average — the true
    /// count varies per step and is recorded in `Metrics`.
    pub fn force_evals_per_step(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 1,
            Self::Ias15 => 14,
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
            Self::Ias15 => {
                "15th-order adaptive Gauss-Radau integrator (Rein & Spiegel \
                 2015). Non-symplectic but conserves energy to machine \
                 precision via step-size control; handles close encounters \
                 and high eccentricities without artefacts. Default choice \
                 for long-term, publication-quality integration."
            },
        }
    }

    /// All known variants, in the order shown in the UI combo-box.
    pub const ALL: [IntegratorKind; 4] = [
        IntegratorKind::Ias15,
        IntegratorKind::Yoshida4,
        IntegratorKind::VelocityVerlet,
        IntegratorKind::WisdomHolman,
    ];

    /// Canonical string slug used in `run.toml` and CLI arguments.
    pub fn slug(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "velocity_verlet",
            Self::Yoshida4 => "yoshida4",
            Self::WisdomHolman => "wisdom_holman",
            Self::Ias15 => "ias15",
        }
    }
}

impl std::str::FromStr for IntegratorKind {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "velocity_verlet" => Ok(Self::VelocityVerlet),
            "yoshida4" => Ok(Self::Yoshida4),
            "wisdom_holman" => Ok(Self::WisdomHolman),
            "ias15" => Ok(Self::Ias15),
            _ => Err(format!(
                "unknown integrator {:?}; valid slugs: {}",
                s,
                Self::ALL.iter().map(|k| k.slug()).collect::<Vec<_>>().join(", "),
            )),
        }
    }
}

// ── IntegratorContext ─────────────────────────────────────────────────────────

/// Everything an integrator needs from the simulation besides bodies and dt.
///
/// Passed as `&mut` so the integrator can call `force.compute()` (which
/// requires `&mut self` for tree rebuilds, etc.).
///
/// # Design
///
/// This struct exists to **avoid coupling integrators to `System`**.
/// An integrator never sees the full `System`; it only sees this narrow
/// interface of force model + physical parameters.
pub struct IntegratorContext<'a> {
    /// The force model (e.g. Barnes-Hut gravity).
    pub force: &'a mut dyn ForceModel,

    /// Gravitational scaling factor: `G_eff = G₀ · g_factor`.
    pub g_factor: f64,

    /// Non-gravitational perturbation forces (radiation, drag, …).
    pub perturbations: &'a [Box<dyn PerturbationForce>],

    /// Optional cooperative wall-clock deadline. Adaptive integrators
    /// (IAS15) check this after each rejection in the retry loop; when
    /// the deadline is passed they accept the current attempt rather
    /// than spending more wall time shrinking `dt`, and mark the step
    /// as [`StepResult::degraded`].
    ///
    /// This is a courtesy knob so the physics-thread batch loop stays
    /// responsive to the UI even in a pathological scene; it does not
    /// strictly bound the step because the current attempt is still
    /// allowed to run to completion. Fixed-step integrators ignore it.
    pub deadline: Option<std::time::Instant>,
}

// ── StepResult ────────────────────────────────────────────────────────────────

/// Output produced by a single integration step.
///
/// Centralises the physical diagnostics that `System` needs after each step,
/// so no integrator-specific logic leaks into the orchestrator.
pub struct StepResult {
    /// Simulated time actually advanced by this call, in the same units as
    /// the `dt` input. For fixed-step integrators (VV, Y4, WH) this is
    /// always equal to the requested `dt`. For adaptive integrators (IAS15)
    /// the step is a single sub-step whose size is chosen by the error
    /// controller; `consumed_dt` reports the accepted size — the caller is
    /// responsible for re-invoking `step` until the desired budget is met.
    ///
    /// `System::step` advances `System::t` by **this value**, never by the
    /// requested `dt`. This keeps `System::t` consistent with the physical
    /// state of the bodies when an adaptive integrator accepts a sub-step
    /// smaller than the caller's budget (Rein & Spiegel 2015, §2.3).
    pub consumed_dt: f64,

    /// Gravitational potential energy at the end-of-step positions,
    /// **already scaled** by `g_factor`.
    pub potential_energy: f64,

    /// `true` if the integrator fell back to a different algorithm this step
    /// (e.g. Wisdom–Holman → Yoshida4 when the dominance criterion fails).
    pub used_fallback: bool,

    /// Dense-output snapshot for sub-step interpolation.
    ///
    /// IAS15 fills this with the accepted sub-step's b-coefficients, valid
    /// over `[t − consumed_dt, t]`. Other integrators leave it `None`;
    /// [`System::step`] supplies an Order-2 fallback using the pre-step
    /// kinematics.
    pub step_snapshot: Option<super::dense::DenseSnapshot>,

    /// `true` when the integrator accepted a sub-step under duress rather
    /// than on merit. For IAS15 this means the error controller wanted to
    /// shrink `dt` further but hit the `DT_MIN` floor; the step was taken
    /// anyway to avoid stalling the simulation, but the local truncation
    /// bound `ε` was not actually satisfied. Fixed-step integrators (VV,
    /// Y4, WH) always report `false`. Callers that care about energy-budget
    /// quality can use this to surface a warning or log a degraded-step
    /// event (REBOUND logs an equivalent `ias15.min_dt` warning).
    pub degraded: bool,
}

// ── ExecutionProfile ──────────────────────────────────────────────────────────

/// How an integrator expects to be driven.
///
/// This is part of the *contract* an integrator advertises — not a hint.
/// Consumers (physics thread, UI, benchmark) read it to choose their
/// execution discipline and to adapt their feedback to the user.
///
/// Two profiles are recognised today:
///
/// * [`Realtime`](ExecutionProfile::Realtime) — per-step wall time is
///   bounded by the force evaluation cost (O(N) or O(N log N) per step)
///   and the user-facing `dt`. Safe to drive from a render loop at
///   60 Hz. Fixed-step integrators (Velocity Verlet, Yoshida 4,
///   Wisdom–Holman) fall here.
///
/// * [`Precision`](ExecutionProfile::Precision) — per-step wall time is
///   unbounded in the adversarial case. IAS15's adaptive controller can
///   shrink `dt` arbitrarily in a stiff regime, spending seconds to
///   minutes on a single visible frame. These integrators must be run
///   off the render thread to completion ("precision run" UI), with a
///   progress indicator rather than a frame budget.
///
/// New integrators should default to `Realtime` unless their
/// algorithmic structure makes per-step wall time unbounded.
///
/// # Future evolution
///
/// Two values is enough today. If a future integrator sits between
/// these (e.g. FMM with amortised per-step cost that spikes on
/// reorganisation), extend this enum rather than adding parallel flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionProfile {
    /// Bounded per-step wall time; drive from a real-time loop.
    Realtime,
    /// Unbounded per-step wall time; run to completion off-thread.
    Precision,
}

// ── Integrator trait ──────────────────────────────────────────────────────────

/// A symplectic (or general) N-body integrator.
///
/// # Contract
///
/// - `step()` advances `bodies` by **at most** the requested `dt`.
///   - Fixed-step integrators (VV, Y4, WH) always consume exactly `dt`.
///   - Adaptive integrators (IAS15) take **one** internal sub-step whose
///     size is chosen by the error controller; the accepted size is
///     reported via [`StepResult::consumed_dt`]. The caller loops until
///     the desired simulation time is reached, mirroring the REBOUND API
///     (`reb_integrator_ias15_part1/2` — Rein & Spiegel 2015).
/// - `step()` may call `ctx.force.compute()` one or more times.
/// - `step()` must leave `acc` consistent with the final body positions
///   (so that diagnostics can read it).
/// - `step()` must apply `ctx.g_factor` scaling and `ctx.perturbations`
///   at the appropriate points in the scheme.
///
/// # Mutability
///
/// `&mut self` is required because some integrators carry internal state
/// across steps (e.g. Wisdom–Holman's fallback integrator, IAS15's
/// predictor–corrector history).
pub trait Integrator: Send {
    /// Advance the system by one time step.
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<(f64, f64)>,
    ) -> StepResult;

    /// Which algorithm this integrator represents.
    fn kind(&self) -> IntegratorKind;

    /// Set the error tolerance for adaptive integrators (IAS15).
    /// No-op for fixed-step integrators.
    fn set_epsilon(&mut self, _eps: f64) {}

    /// Return the current error tolerance, if applicable.
    fn epsilon(&self) -> Option<f64> {
        None
    }

    /// Cumulative adaptive-integrator counters. `None` for fixed-step
    /// integrators (they have no sub-step / rejection / Picard notion).
    /// See [`AdaptiveStats`] for field semantics.
    fn adaptive_stats(&self) -> Option<AdaptiveStats> {
        None
    }

    /// Execution discipline this integrator requires of its caller.
    /// Default is [`Realtime`](ExecutionProfile::Realtime); adaptive /
    /// implicit schemes whose per-step wall time is unbounded should
    /// return [`Precision`](ExecutionProfile::Precision).
    fn execution_profile(&self) -> ExecutionProfile {
        ExecutionProfile::Realtime
    }

    /// Whether this integrator requires the force model to be a
    /// deterministic function of state — i.e. `f(x, v, t)` bit-reproducible
    /// across calls with the same `(x, v, t)` to within f64 ULP.
    ///
    /// # Why this matters — Picard continuity
    ///
    /// High-order implicit methods (IAS15 in particular) solve an
    /// implicit system by **Picard predictor–corrector iteration**
    /// within each adaptive sub-step. The iteration converges iff its
    /// operator is a contraction, which requires the force function
    /// to be *continuous and deterministic in state across
    /// iterations*:
    ///
    /// 1. Between iterations, body positions drift by a small amount —
    ///    that is the point of the corrector.
    /// 2. The corrector consults `ForceModel::compute` on the
    ///    perturbed state.
    /// 3. If the force function has position-dependent *topological*
    ///    discontinuities (e.g. Barnes-Hut: a body near a cell
    ///    boundary crosses leaves in response to sub-ULP drift → the
    ///    multipole approximation for that body's far-field changes
    ///    discretely), the Picard operator is not a contraction. The
    ///    iteration oscillates at the discontinuity scale; the outer
    ///    controller reads the oscillation as truncation error,
    ///    rejects the step, shrinks `dt`, and cascades toward
    ///    `DT_MIN` — regardless of how physically benign the scenario
    ///    actually is.
    ///
    /// This is why REBOUND pairs IAS15 exclusively with direct O(N²)
    /// summation. The pairing is a mathematical prerequisite of the
    /// method, not an implementation shortcut.
    ///
    /// # Why other integrators do not need this
    ///
    /// Low-order explicit / symplectic schemes (Verlet, Yoshida,
    /// WHFast) do not solve an implicit system. Their per-step error
    /// bound absorbs force-evaluation noise at O(dt²) or better, so
    /// Barnes-Hut's tree-rebuild variation is invisible at the
    /// trajectory level. They return `false` and may be paired with
    /// any force model.
    ///
    /// # Enforcement
    ///
    /// `System::set_integrator` reads this together with the force
    /// model's [`is_deterministic`](crate::physics::integrator::force_model::ForceModel::is_deterministic)
    /// and auto-reconfigures the force model (raises the exact
    /// threshold so BH is bypassed) when they conflict. The
    /// auto-correction emits a structured
    /// [`warn_diag!`](crate::warn_diag) event so the user sees
    /// exactly what changed and why.
    ///
    /// # Future evolution
    ///
    /// Returning a boolean here will eventually be upgraded to a
    /// required `DeterminismLevel` (`Strict` / `Approximate { bound }`
    /// / `Nondeterministic`), once a second non-trivial force model
    /// (FMM with bounded multipole error, GPU kernel with reduction
    /// noise) makes the distinction load-bearing. Until then the
    /// boolean is sufficient and does not encode spurious precision.
    fn requires_deterministic_force(&self) -> bool {
        false
    }
}

/// Per-integrator lifetime counters exposed by [`Integrator::adaptive_stats`].
///
/// Counts are **cumulative** from integrator construction. Compute rates in
/// the caller (e.g. `rejections / substeps` as an acceptance-efficiency
/// indicator, or `picard_iters / attempts` for mean inner-loop work).
#[derive(Debug, Clone, Copy, Default)]
pub struct AdaptiveStats {
    /// Accepted sub-steps.
    pub substeps: u64,
    /// Total rejected attempts (controller shrank `dt` and retried).
    /// Sum of `rejections_picard` + `rejections_truncation`.
    pub rejections: u64,
    /// Rejections caused by Picard predictor–corrector non-convergence.
    /// A high ratio vs. `rejections_truncation` means the step size
    /// routinely exceeds the local Lipschitz bound (stiff / high-e regime).
    pub rejections_picard: u64,
    /// Rejections where Picard converged but truncation error exceeded
    /// `ε`. This is the "well-behaved" rejection class handled by the
    /// standard `(ε/err)^(1/7)` controller.
    pub rejections_truncation: u64,
    /// Total Picard iterations across all attempts (accepted + rejected).
    pub picard_iters: u64,
    /// Accepted sub-steps that hit the `DT_MIN` escape or deadline
    /// without meeting tolerance. Should be zero in healthy scenes.
    pub degraded: u64,
}
