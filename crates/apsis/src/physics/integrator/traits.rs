//! [`Integrator`] trait, per-step [`IntegratorContext`], [`StepResult`],
//! and the [`IntegratorKind`] enum used for serialisation and metrics.

use crate::domain::body::Body;
use crate::math::Vec3;
use crate::physics::integrator::force_model::ForceModel;
use crate::physics::integrator::operator::{
    HamiltonianOperator, NonConservativeOperator, Operator,
};

// ── IntegratorKind (serialisable enum) ────────────────────────────────────────

/// Identifies an integration algorithm without carrying behaviour.
///
/// Used for snapshot serialisation and `Metrics`. The actual stepping
/// logic lives in structs that implement [`Integrator`].
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
    /// Short human-readable label, suitable for diagnostic output.
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

    /// Nominal force evaluations per full step. Amortised averages for
    /// adaptive schemes (IAS15 ~7 evals × ~2 Picard iters per substep)
    /// and amortised lower bounds for hybrid schemes (Mercurius, no
    /// active encounter; ImplicitMidpoint, `max_iterations = 10`). The
    /// true per-step count is recorded in `Metrics`.
    pub fn force_evals_per_step(self) -> u32 {
        match self {
            Self::VelocityVerlet => 2,
            Self::Yoshida4 => 4,
            Self::WisdomHolman => 1,
            Self::WHFast => 1,
            Self::Ias15 => 14,
            Self::Mercurius => 2,
            Self::ImplicitMidpoint => 11,
        }
    }

    /// All known variants, ordered roughly by typical use frequency.
    pub const ALL: [IntegratorKind; 7] = [
        IntegratorKind::Ias15,
        IntegratorKind::Mercurius,
        IntegratorKind::WHFast,
        IntegratorKind::ImplicitMidpoint,
        IntegratorKind::Yoshida4,
        IntegratorKind::VelocityVerlet,
        IntegratorKind::WisdomHolman,
    ];

    /// Canonical string slug used in serialisation and CLI arguments.
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
pub struct IntegratorContext<'a> {
    /// The force model (e.g. Barnes-Hut gravity).
    pub force: &'a mut dyn ForceModel,

    /// Gravitational scaling factor: `G_eff = G₀ · g_factor`.
    pub g_factor: f64,

    /// Hamiltonian operators: contribute a force and an energy term.
    /// Symplectic integrators preserve their invariants under these.
    pub hamiltonian_perturbations: &'a [Box<dyn HamiltonianOperator>],

    /// Non-conservative operators: force only, no Hamiltonian (drag,
    /// radiation reaction). Break symplecticity by construction.
    pub non_conservative_perturbations: &'a [Box<dyn NonConservativeOperator>],

    /// Step-boundary observers, dispatched at synchronised state after
    /// the integrator has resolved the outer step.
    pub observers: &'a mut [Box<dyn Operator>],
}

// ── StepResult ────────────────────────────────────────────────────────────────

/// Output produced by a single integration step.
pub struct StepResult {
    /// Simulated time actually advanced by this call. Fixed-step
    /// integrators always return the requested `dt`; adaptive
    /// integrators return the accepted sub-step size. `System::step`
    /// advances `System::t` by this value (Rein & Spiegel 2015, §2.3).
    pub consumed_dt: f64,

    /// Gravitational potential energy at the end-of-step positions,
    /// **already scaled** by `g_factor`.
    pub potential_energy: f64,

    /// `true` if the integrator fell back to a different algorithm this step
    /// (e.g. Wisdom–Holman → Yoshida4 when the dominance criterion fails).
    pub used_fallback: bool,

    /// Dense-output snapshot for sub-step interpolation, valid over
    /// `[t − consumed_dt, t]`. Integrators with dense output fill this;
    /// others leave it `None` and [`System::step`] supplies an Order-2
    /// fallback using the pre-step kinematics.
    pub step_snapshot: Option<super::dense::DenseSnapshot>,

    /// `true` when an adaptive integrator hit the `DT_MIN` floor and
    /// accepted a sub-step without actually meeting its local error
    /// bound. Fixed-step integrators always report `false`.
    pub degraded: bool,

    /// Hierarchy regime classification — `Some(_)` only for integrators
    /// whose derivation depends on a hierarchical mass distribution
    /// (Wisdom-Holman family); `None` otherwise.
    pub hierarchy_signal: Option<HierarchySignal>,
}

// ── HierarchySignal ───────────────────────────────────────────────────────────

/// Classification of a system's mass distribution against the dominance
/// criterion underlying the Wisdom-Holman perturbation expansion (WH
/// 1991 §III). Two predicates: the central body must be the most
/// massive, and `m_0 / Σ_{i≥1} m_i` must clear a dominance threshold.
///
/// Observability only — surfaced through [`StepResult`] so callers can
/// detect when the system has left the validated regime; no integrator
/// branches on the value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HierarchySignal {
    /// `m_0 ≥ max(m_i)` and `m_0 / Σ m_i ≥ 10`. Inside the validated
    /// WH regime.
    Hierarchical,
    /// `m_0 ≥ max(m_i)` and `5 ≤ m_0 / Σ m_i < 10`. At the edge of the
    /// small-parameter regime; expect WH energy drift above the 1991
    /// published floor.
    Borderline,
    /// `m_0 < max(m_i)` or `m_0 / Σ m_i < 5`. WH derivation does not
    /// apply.
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

    /// Short human-readable label, suitable for diagnostic output.
    pub fn label(self) -> &'static str {
        match self {
            Self::Hierarchical => "hierarchical",
            Self::Borderline => "borderline",
            Self::Violated => "violated",
        }
    }
}

// ── Integrator trait ──────────────────────────────────────────────────────────

/// A symplectic (or general) N-body integrator.
///
/// # Contract
///
/// - `step()` advances `bodies` by **at most** the requested `dt_hint`.
///   Fixed-step schemes consume exactly `dt_hint`; self-adaptive
///   schemes treat it as a first-call seed and report the accepted
///   sub-step via [`StepResult::consumed_dt`] (Rein & Spiegel 2015 §2.3).
/// - `step()` may call `ctx.force.compute()` one or more times.
/// - `step()` leaves `acc` consistent with the final body positions.
/// - `step()` applies `ctx.g_factor` and the operators in
///   `ctx.{hamiltonian,non_conservative}_perturbations` at the
///   appropriate points in the scheme.
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

    /// `true` if per-step wall time is unbounded in the adversarial
    /// case (controller may shrink `dt` arbitrarily). Default `false`;
    /// override on schemes whose per-step cost is not bounded by the
    /// force evaluation alone.
    fn is_adaptive(&self) -> bool {
        false
    }

    /// Whether the integrator's internal controller decides the actual
    /// step size, treating `dt_hint` as a *hint* rather than a hard cap.
    /// The orchestrator reads this to decide whether to reset
    /// `current_dt` to `user_dt` after each step (fixed-step) or to
    /// adopt [`proposed_next_dt`](Self::proposed_next_dt) (self-adaptive).
    fn controls_own_step_size(&self) -> bool {
        false
    }

    /// The step size the integrator's controller proposes for the next
    /// call. `None` for fixed-step schemes; self-adaptive schemes
    /// return their most-recent `dt_next` recommendation.
    fn proposed_next_dt(&self) -> Option<f64> {
        None
    }

    /// Apply a uniform translation `(-dx, -dy)` to every body, routing
    /// the shift through whatever compensated-summation accumulators
    /// the integrator maintains. Default impl performs an uncompensated
    /// subtraction; integrators carrying per-body compensation buffers
    /// override to preserve them (necessary for bit-reproducibility
    /// across COM-recentering shifts).
    fn recenter_bodies(&mut self, bodies: &mut [Body], dx: f64, dy: f64) {
        for b in bodies.iter_mut() {
            b.pos_x -= dx;
            b.pos_y -= dy;
        }
    }

    /// One-shot cap on the next `step()` call's time advance. Set by
    /// `System::integrate_until` to clip the final step onto `t_end`
    /// (exact finish time). Self-adaptive integrators must honour the
    /// cap for exactly one call and leave their controller rhythm
    /// unchanged by the capped step. No-op for fixed-step schemes —
    /// the orchestrator caps their `dt_hint` directly.
    fn cap_next_step(&mut self, _max_dt: f64) {}

    /// Set the error tolerance for adaptive integrators. No-op
    /// otherwise.
    fn set_epsilon(&mut self, _eps: f64) {}

    /// Return the current error tolerance, if applicable.
    fn epsilon(&self) -> Option<f64> {
        None
    }

    /// Set the Hill-radius multiplier for hybrid close-encounter
    /// integrators. No-op for integrators without a changeover-based
    /// encounter trigger.
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

    /// `true` if the integrator's Picard fixed-point iteration requires
    /// `f(x, v, t)` to be bit-reproducible across calls. Enforced at
    /// `System::set_integrator` against
    /// [`ForceModel::is_deterministic`](crate::physics::integrator::force_model::ForceModel::is_deterministic);
    /// see [ADR-003](../../docs/adr/003-integrator-execution-profile.md).
    fn requires_deterministic_force(&self) -> bool {
        false
    }

    /// Serialise the integrator's per-step scratch (predictor history,
    /// Kepler accumulators, hybrid-mode flags, etc.) for inclusion in a
    /// mid-run snapshot. Empty by default for stateless integrators.
    fn resume_state(&self) -> Vec<u8> {
        Vec::new()
    }

    /// Restore previously serialised scratch state. Default no-op accepts
    /// only an empty payload — stateless integrators have nothing to
    /// restore. Override to validate and deserialise.
    fn restore_resume_state(&mut self, bytes: &[u8]) -> Result<(), ResumeError> {
        if bytes.is_empty() { Ok(()) } else { Err(ResumeError::UnsupportedFormat) }
    }
}

/// Errors returned by [`Integrator::restore_resume_state`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResumeError {
    /// Payload version, magic, or schema is unrecognised.
    UnsupportedFormat,
    /// Payload ended before all expected fields had been read.
    Truncated,
    /// Body count encoded in the payload disagrees with the System's
    /// current `bodies.len()`.
    BodyCountMismatch { expected: usize, found: usize },
}

impl std::fmt::Display for ResumeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedFormat => write!(f, "unsupported resume-state format"),
            Self::Truncated => write!(f, "resume-state payload truncated"),
            Self::BodyCountMismatch { expected, found } => {
                write!(f, "resume-state body count mismatch: expected {expected}, found {found}")
            },
        }
    }
}

impl std::error::Error for ResumeError {}

/// Per-integrator lifetime counters exposed by [`Integrator::adaptive_stats`].
///
/// Counts are **cumulative** from integrator construction. Always-on:
/// one `saturating_add` per accept/reject path. See
/// [`docs/integrator.md`](../../docs/integrator.md) §"Diagnostic
/// counters in AdaptiveStats" for the healthy regime each counter
/// implies and how to interpret elevated values.
#[derive(Debug, Clone, Copy, Default)]
pub struct AdaptiveStats {
    /// Accepted sub-steps.
    pub substeps: u64,
    /// Total rejected attempts (controller shrank `dt` and retried).
    /// Sum of `rejections_picard` + `rejections_truncation`.
    pub rejections: u64,
    /// Rejections caused by Picard predictor–corrector non-convergence.
    pub rejections_picard: u64,
    /// Rejections where Picard converged but truncation error exceeded `ε`.
    pub rejections_truncation: u64,
    /// Total Picard iterations across all attempts (accepted + rejected).
    pub picard_iters: u64,
    /// Accepted sub-steps that hit the `DT_MIN` floor without meeting
    /// tolerance.
    pub degraded: u64,
    /// Picard early-exits via the stagnation guard (residual stopped
    /// decreasing for two consecutive iterations, accepted as
    /// success-by-saturation).
    pub picard_stagnations: u64,
    /// "Shrink → grow" reversals in the controller's `dt_next`: a
    /// reversal is counted whenever `dt_next` increases after an
    /// accept whose previous accept's `dt_next` had decreased.
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
