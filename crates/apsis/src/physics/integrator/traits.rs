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
use crate::math::Vec3;
use crate::physics::integrator::force_model::ForceModel;
use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Operator,
};

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
    WHFast,
    Ias15,
    Mercurius,
    ImplicitMidpoint,
}

impl IntegratorKind {
    /// Short human-readable label used in the UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::VelocityVerlet => "Velocity Verlet (2nd)",
            Self::Yoshida4 => "Yoshida 4th-order",
            Self::WisdomHolman => "Wisdom–Holman (2nd, Keplerian)",
            Self::WHFast => "WHFast (2nd, Keplerian, compensated)",
            Self::Ias15 => "IAS15 (15th, adaptive)",
            Self::Mercurius => "Mercurius (hybrid, WH + IAS15)",
            Self::ImplicitMidpoint => "Implicit Midpoint (2nd, A-stable)",
        }
    }

    /// Execution profile without constructing an integrator instance.
    ///
    /// Provides the UI with a cheap way to ask "is this a precision
    /// integrator?" before deciding whether to show the Precision
    /// Run panel, the setup modal, or the confirmation dialog.
    /// Mirrors the value that [`Integrator::execution_profile`]
    /// returns on the constructed instance.
    pub fn execution_profile(self) -> ExecutionProfile {
        match self {
            Self::Ias15 => ExecutionProfile::Precision,
            Self::VelocityVerlet
            | Self::Yoshida4
            | Self::WisdomHolman
            | Self::WHFast
            | Self::Mercurius
            | Self::ImplicitMidpoint => ExecutionProfile::Realtime,
        }
    }

    /// Formal convergence order in the time step.
    pub fn order(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 2,
            Self::WHFast => 2,
            Self::Ias15 => 15,
            Self::Mercurius => 2,
            Self::ImplicitMidpoint => 2,
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
            Self::WHFast => 1,
            Self::Ias15 => 14,
            // Mercurius: 2 K-weighted half-kicks + analytical Kepler drift +
            // an IAS15 sub-integration whose cost is data-dependent. The
            // quoted number is an amortised lower bound assuming no pair
            // is in close encounter; engaging encounters add IAS15 substeps.
            Self::Mercurius => 2,
            // ImplicitMidpoint: one force eval per fixed-point iteration
            // plus one final eval at the end-state. Mean iteration count
            // is 3-6 for non-stiff conservative gravity; quoted figure is
            // an upper bound assuming `max_iterations = 10`.
            Self::ImplicitMidpoint => 11,
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
            Self::WHFast => {
                "Wisdom-Holman split with compensated summation on per-step \
                 position and velocity accumulators (Rein & Tamayo 2015). \
                 Round-off envelope reduced from O(N · ε) to O(√N · ε), \
                 unlocking long-horizon planetary integration. Same \
                 hierarchical-mass requirement as WH."
            },
            Self::Ias15 => {
                "15th-order adaptive Gauss-Radau integrator (Rein & Spiegel \
                 2015). Non-symplectic but conserves energy to machine \
                 precision via step-size control; handles close encounters \
                 and high eccentricities without artefacts. Default choice \
                 for long-term, publication-quality integration."
            },
            Self::Mercurius => {
                "Hybrid symplectic integrator (Rein et al. 2019). Wisdom-Holman \
                 outer step with K-weighted planet-planet kicks; IAS15 \
                 sub-integrates the (1-K)-weighted close-encounter residual \
                 over the same outer interval. Localises encounter cost to \
                 the encountering pair while preserving secular stability \
                 elsewhere. Requires a hierarchical mass distribution."
            },
            Self::ImplicitMidpoint => {
                "Single-stage Gauss-Legendre symplectic method (Hairer-Lubich-Wanner \
                 2006, Chapter II.1.4). A-stable on the entire left half-plane and \
                 time-symmetric. Iterates per step on the implicit midpoint state \
                 (Picard, max 10 iterations by default; Newton-Krylov reserved in \
                 the API). Makes no central-mass dominance assumption — accepts \
                 BH binaries, equal-mass triples, particle clouds. Not L-stable; \
                 dissipation-dominant extreme regimes need Radau IIA or BDF."
            },
        }
    }

    /// All known variants, in the order shown in the UI combo-box.
    pub const ALL: [IntegratorKind; 7] = [
        IntegratorKind::Ias15,
        IntegratorKind::Mercurius,
        IntegratorKind::WHFast,
        IntegratorKind::ImplicitMidpoint,
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
            Self::WHFast => "whfast",
            Self::Ias15 => "ias15",
            Self::Mercurius => "mercurius",
            Self::ImplicitMidpoint => "implicit_midpoint",
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
            "whfast" => Ok(Self::WHFast),
            "ias15" => Ok(Self::Ias15),
            "mercurius" => Ok(Self::Mercurius),
            "implicit_midpoint" => Ok(Self::ImplicitMidpoint),
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

    /// Hamiltonian operators registered on the system. Each contributes
    /// a force (via `accumulate_force`) and an energy term (via
    /// `energy_contribution`). 1PN GR correction is the canonical
    /// example. Symplectic integrators preserve their conservation
    /// invariants when *only* operators of this class are registered.
    pub hamiltonian_perturbations: &'a [Box<dyn HamiltonianOperator>],

    /// Non-conservative operators registered on the system. Each
    /// contributes a force but no Hamiltonian (drag, radiation
    /// reaction, dissipative coupling). Symplectic integrators degrade
    /// silently with these registered; the system emits a `warn_diag`
    /// at registration time so the broken invariant is documented.
    pub non_conservative_perturbations: &'a [Box<dyn NonConservativeOperator>],

    /// Pure observers registered on the system. No force, no energy,
    /// just step-boundary `observe` calls (Shadow Hamiltonian tracker,
    /// audit trail emitters, etc.). Dispatched at synchronized state
    /// after the integrator has fully resolved the outer step.
    pub observers: &'a mut [Box<dyn Operator>],

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
    /// event.
    pub degraded: bool,

    /// Mass-distribution-based hierarchy regime classification, populated
    /// only by integrators whose derivation depends on a hierarchical mass
    /// distribution. `None` when the concept does not apply (VV, Y4, IAS15
    /// have no hierarchy assumption); `Some(_)` for Wisdom-Holman with
    /// values graded according to the WH derivation's small-parameter
    /// regime. Observability only — no integrator branches on the value.
    pub hierarchy_signal: Option<HierarchySignal>,
}

// ── HierarchySignal ───────────────────────────────────────────────────────────

/// Classification of a system's mass distribution against the dominance
/// criterion underlying the Wisdom-Holman perturbation expansion.
///
/// The variants grade where the system sits relative to the WH derivation's
/// small-parameter regime. Two ratios drive the classification:
///
/// 1. The central body must be the most massive (no other body more massive
///    than `bodies[0]`).
/// 2. The central-to-rest mass ratio `m_0 / Σ_{i≥1} m_i` must clear a
///    dominance threshold for the perturbation series to converge in the
///    asymptotic sense WH 1991 §III assumes.
///
/// The signal is observability-only. The integrator does not change
/// behaviour based on it; the value is surfaced through [`StepResult`] so
/// callers (UI, logging, validation) can detect when the system has left
/// the validated regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HierarchySignal {
    /// Central body is dominant; the WH derivation operates inside its
    /// validated regime. `m_0 ≥ max(m_i)` and `m_0 / Σ m_i ≥ 10`.
    Hierarchical,
    /// Central body is the most massive but its dominance over the rest of
    /// the system is approaching the threshold. The perturbation expansion
    /// is at the edge of its small-parameter regime; observed energy drift
    /// may exceed the WH 1991 published floor without indicating a defect.
    /// `m_0 ≥ max(m_i)` and `5 ≤ m_0 / Σ m_i < 10`.
    Borderline,
    /// Central body fails the dominance criterion. The WH derivation does
    /// not apply; observed conservation should not be expected to match the
    /// validated regime. `m_0 < max(m_i)` or `m_0 / Σ m_i < 5`.
    Violated,
}

impl HierarchySignal {
    /// Classify a mass distribution by the Wisdom-Holman dominance criterion.
    ///
    /// Returns [`HierarchySignal::Violated`] for trivially degenerate inputs
    /// (fewer than two bodies, all zero mass).
    pub fn classify(masses: &[f64]) -> Self {
        if masses.len() < 2 {
            return Self::Violated;
        }
        let m0 = masses[0];
        let m_rest: f64 = masses[1..].iter().sum();
        let max_other = masses[1..].iter().copied().fold(0.0_f64, f64::max);

        if m0 < max_other || m_rest <= 0.0 {
            return Self::Violated;
        }

        let ratio = m0 / m_rest;
        if ratio >= 10.0 {
            Self::Hierarchical
        } else if ratio >= 5.0 {
            Self::Borderline
        } else {
            Self::Violated
        }
    }

    /// Short human-readable label, suitable for diagnostic output and UI.
    pub fn label(self) -> &'static str {
        match self {
            Self::Hierarchical => "hierarchical",
            Self::Borderline => "borderline",
            Self::Violated => "violated",
        }
    }
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
/// - `step()` advances `bodies` by **at most** the requested `dt_hint`.
///   - Fixed-step integrators (VV, Y4, WH) always consume exactly
///     `dt_hint`. They report
///     [`controls_own_step_size`](Self::controls_own_step_size) as
///     `false` and leave [`proposed_next_dt`](Self::proposed_next_dt)
///     at its `None` default.
///   - Self-adaptive integrators (IAS15) treat `dt_hint` as a
///     *first-call seed* for their internal step controller; on
///     subsequent calls the controller-chosen step is used instead.
///     The accepted size of the single sub-step is reported via
///     [`StepResult::consumed_dt`]. The caller loops until the desired
///     simulation time is reached — the substep-granularity contract
///     specified in Rein & Spiegel (2015) §2.3. These integrators
///     report [`controls_own_step_size`](Self::controls_own_step_size)
///     as `true` so the orchestrator can read their controller's
///     next-step proposal via
///     [`proposed_next_dt`](Self::proposed_next_dt) and surface it to
///     the user / UI as the simulation's "current dt".
/// - `step()` may call `ctx.force.compute()` one or more times.
/// - `step()` must leave `acc` consistent with the final body positions
///   (so that diagnostics can read it).
/// - `step()` must apply `ctx.g_factor` scaling and `ctx.perturbations`
///   at the appropriate points in the scheme.
///
/// # Why the hint-vs-cap distinction matters
///
/// Treating the orchestrator's `dt_hint` as a hard per-call cap silently
/// pins a self-adaptive integrator to whatever step the user supplied as
/// an initial guess: the controller can shrink below the cap on stiff
/// regions but can never grow above it, even when the local truncation
/// error would permit it. This was the IAS15 floor-cascade bug
/// investigated in
/// `docs/experiments/2026-04-26-ias15-warmstart-bug.md` — the integrator
/// behaved like a fixed-step scheme with internal sub-stepping below the
/// cap, producing correct trajectories at orders-of-magnitude unnecessary
/// substep counts.
///
/// The two-method contract here makes the distinction explicit at the
/// trait boundary: a future SABA, Hermite, or MERCURIUS implementation
/// declares its discipline up front, and the orchestrator routes the
/// hint vs cap accordingly.
///
/// # Mutability
///
/// `&mut self` is required because some integrators carry internal state
/// across steps (e.g. Wisdom–Holman's fallback integrator, IAS15's
/// predictor–corrector history).
pub trait Integrator: Send {
    /// Advance the system by one time step.
    ///
    /// `dt_hint` is interpreted per the contract advertised by
    /// [`controls_own_step_size`](Self::controls_own_step_size):
    /// fixed-step schemes consume exactly `dt_hint`; self-adaptive
    /// schemes use it as a first-call seed and otherwise let their
    /// controller pick.
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt_hint: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult;

    /// Which algorithm this integrator represents.
    fn kind(&self) -> IntegratorKind;

    /// Whether the integrator's internal controller decides the actual
    /// step size, treating the orchestrator's `dt_hint` as a *hint*
    /// rather than a hard cap. `false` means the integrator consumes
    /// exactly `dt_hint`; `true` means it picks its own step (possibly
    /// after seeding from `dt_hint` on the first call).
    ///
    /// The orchestrator uses this to decide whether to reset its own
    /// `current_dt` field after each step (fixed-step → reset to
    /// `user_dt`; self-adaptive → adopt
    /// [`proposed_next_dt`](Self::proposed_next_dt)) and whether to
    /// surface a controller-driven step size in UI panels.
    fn controls_own_step_size(&self) -> bool {
        false
    }

    /// The step size the integrator's controller proposes for the next
    /// call, if applicable. Returns `None` for fixed-step schemes; for
    /// self-adaptive schemes it tracks the controller's most-recent
    /// `dt_next` recommendation.
    ///
    /// The orchestrator may use this to update its `current_dt` field
    /// so external observers (UI dt readout, headless CSV columns)
    /// reflect the controller's actual cadence rather than the user's
    /// initial guess.
    fn proposed_next_dt(&self) -> Option<f64> {
        None
    }

    /// Apply a uniform translation `(-dx, -dy)` to every body, routing the
    /// shift through whatever compensated-summation accumulators the
    /// integrator maintains for body position.
    ///
    /// # Why this matters — bit-reproducibility
    ///
    /// IAS15 carries a Neumaier-style compensation buffer (`csx`) that
    /// pairs with each body's stored position. The pair `(x, csx)`
    /// represents an extended-precision running sum: every accepted
    /// substep updates `x` via `add_cs(p = body.pos_x, csp = csx, inp =
    /// position increment)`, which preserves low-order bits across the
    /// long integration horizons IAS15 advertises (~10⁹ orbits at f64
    /// machine precision, per Rein & Spiegel 2015 §3).
    ///
    /// External translations of body position — most commonly the
    /// periodic COM-recentering applied by `System::step` — disrupt this
    /// invariant when written as a bare `body.pos_x -= dx`. The compensation
    /// `csx` then references the rounding history of the *pre-shift*
    /// running sum, but `body.pos_x` has been arbitrarily perturbed; the next
    /// `add_cs` call wipes the prior compensation rather than continuing
    /// to track it. For a single sub-ULP shift this loss is negligible,
    /// but it accumulates into a bit-reproducibility gap on long runs and
    /// undermines snapshot-replay determinism — exactly the property the
    /// reproducibility ADR commits to.
    ///
    /// Routing the shift through `add_cs` (or its arithmetic equivalent)
    /// preserves the invariant: the new `(x, csx)` pair is the
    /// compensated representation of `(x_old, csx_old) - (dx, dy)`.
    /// Default impl performs an uncompensated subtraction (correct for
    /// integrators with no per-body compensation buffer); IAS15
    /// overrides to use its own buffers.
    fn recenter_bodies(&mut self, bodies: &mut [Body], dx: f64, dy: f64) {
        for b in bodies.iter_mut() {
            b.pos_x -= dx;
            b.pos_y -= dy;
        }
    }

    /// Set the error tolerance for adaptive integrators (IAS15).
    /// No-op for fixed-step integrators.
    fn set_epsilon(&mut self, _eps: f64) {}

    /// Return the current error tolerance, if applicable.
    fn epsilon(&self) -> Option<f64> {
        None
    }

    /// Set the Hill-radius multiplier for hybrid close-encounter
    /// integrators (Mercurius). No-op for integrators without a
    /// changeover-based encounter trigger.
    fn set_hill_factor(&mut self, _alpha: f64) {}

    /// Return the active Hill-radius multiplier, if applicable.
    fn hill_factor(&self) -> Option<f64> {
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
    /// The pairing of IAS15 with direct O(N²) summation is therefore
    /// a mathematical prerequisite of the method (Rein & Spiegel 2015
    /// §2.1 explicitly assumes a deterministic force law for the
    /// predictor–corrector convergence proof), not an implementation
    /// shortcut.
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
///
/// # Instrumentation policy
///
/// Every counter here is incremented unconditionally at runtime — the
/// cost is one `saturating_add` per accept/reject path, well below the
/// Picard arithmetic that surrounds it. Enabling them by default makes
/// the counters available to *any* caller (benchmarks, headless CSV,
/// the lab notebooks under `validation/`) without recompilation, which
/// is the only way to catch the slow-onset cumulative failures that
/// motivated the figure-8 cascade investigation
/// (`docs/experiments/2026-04-26-ias15-warmstart-bug.md`): a regression
/// in those counters can be diagnosed from a single CSV column rather
/// than from a recompile-and-rerun loop.
///
/// More expensive per-step trace data — b-coefficient norms after
/// warmstart, cross-term-vs-diagonal contribution ratios, the
/// step-size ratio history — sit behind the optional `ias15-diag`
/// Cargo feature, since their per-step overhead is non-trivial and
/// they are useful only when actively investigating a controller
/// regression.
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
    /// Picard predictor–corrector early-exits via the stagnation guard
    /// (residual stopped decreasing for two consecutive iterations,
    /// counted as success-by-saturation). Healthy scenes should see
    /// `picard_stagnations / substeps ≪ 1`; a sustained high ratio
    /// indicates the warmstart is biasing `b` outside Picard's basin
    /// of attraction — which on the figure-8 cascade was the
    /// *symptom* of the missing-cross-terms warmstart bug.
    pub picard_stagnations: u64,
    /// Number of "shrink → grow" reversals in the controller's `dt_next`.
    /// A reversal is counted whenever `dt_next` increases after an
    /// accept whose previous accept's `dt_next` had decreased. The
    /// frequency reveals controller chatter: a healthy run on a smooth
    /// regime sees `shrink_grow_cycles / substeps ≈ 0`; a run plagued
    /// by wrong-warmstart bias oscillates as the controller alternately
    /// over- and under-shrinks. Together with `picard_stagnations` this
    /// is the cheapest "is something off?" signal the orchestrator can
    /// surface without enabling the `ias15-diag` feature.
    pub shrink_grow_cycles: u64,
}

#[cfg(test)]
mod hierarchy_signal_tests {
    use super::HierarchySignal;

    #[test]
    fn solar_mercury_is_hierarchical() {
        let masses = [1.0, 1.66e-7];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Hierarchical);
    }

    #[test]
    fn dominance_ratio_at_threshold_is_hierarchical() {
        // m_0 / Σ m_i = 1 / 0.1 = 10 — exactly at the threshold, accepted.
        let masses = [1.0, 0.1];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Hierarchical);
    }

    #[test]
    fn dominance_ratio_in_borderline_band_is_borderline() {
        // m_0 / Σ m_i = 1 / 0.15 ≈ 6.67 — in [5, 10) band, borderline.
        let masses = [1.0, 0.15];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Borderline);
    }

    #[test]
    fn dominance_ratio_below_5x_is_violated() {
        // m_0 / Σ m_i = 1 / 0.25 = 4 — below 5, perturbation expansion fails.
        let masses = [1.0, 0.25];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }

    #[test]
    fn equal_mass_binary_is_violated() {
        let masses = [1.0, 1.0];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }

    #[test]
    fn another_body_more_massive_than_central_is_violated() {
        let masses = [1.0, 10.0, 1.0e-6];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }

    #[test]
    fn single_body_is_violated() {
        let masses = [1.0];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }

    #[test]
    fn empty_mass_distribution_is_violated() {
        let masses: [f64; 0] = [];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }

    #[test]
    fn all_planet_masses_zero_is_violated() {
        // Edge case: avoid division by zero or false-positive Hierarchical.
        let masses = [1.0, 0.0, 0.0];
        assert_eq!(HierarchySignal::classify(&masses), HierarchySignal::Violated);
    }
}
