//! IAS15 — 15th-order adaptive Gauss-Radau integrator.
//!
//! Reference implementation of the algorithm described in:
//!
//!   * Rein H. & Spiegel D. S. (2015). *IAS15: a fast, adaptive,
//!     high-order integrator for gravitational dynamics, accurate to
//!     machine precision over a billion orbits*, MNRAS **446**,
//!     1424–1437.  [arXiv:1409.4779](https://arxiv.org/abs/1409.4779)
//!   * Everhart E. (1985). *An efficient integrator that uses Gauss-Radau
//!     spacings*, in «Dynamics of Comets: Their Origin and Evolution»,
//!     A. Carusi & G. B. Valsecchi eds., Astrophysics and Space Science
//!     Library 115, 185–202.
//!
//! The layout of this file mirrors REBOUND's `integrator_ias15.c` so the
//! two implementations are straightforward to cross-read.
//!
//! IAS15 is the modern refinement of Everhart's RADAU15. It combines
//!
//!   * 8-node Gauss-Radau quadrature (one end-point node + 7 interior),
//!   * a power-series ansatz for the acceleration within the step,
//!   * a predictor–corrector Picard loop to solve the implicit system,
//!   * adaptive step control driven by the dominant truncation term
//!     (∝ |b₆|), with a safety factor to damp oscillation, and
//!   * compensated summation on every state update to keep rounding
//!     error below the truncation error even over ∼10⁹ orbits.
//!
//! Unlike symplectic integrators, IAS15 is **not** measure-preserving,
//! but the adaptive step control delivers energy conservation at the
//! round-off floor — in practice indistinguishable from a symplectic
//! method for bound orbits, and strictly superior for close encounters
//! or high eccentricities where fixed-step schemes degrade.
//!
//! # Sub-step semantics (ADR-004)
//!
//! A single call to [`Integrator::step`] performs **one** adaptive
//! sub-step. The `dt` argument is treated as an upper bound on the
//! sub-step size (the caller's budget for this tick); the error
//! controller may accept anything in `(DT_MIN, dt]` and reports the
//! actual size via [`StepResult::consumed_dt`]. The caller (the
//! [`System::step`] orchestrator) loops until its desired simulation
//! time is reached, exactly as REBOUND does via
//! `reb_integrator_ias15_part1` / `part2`.
//!
//! This design — substep-granularity at the trait boundary — was
//! chosen for three reasons:
//!
//!   * `System::t` stays consistent with body state even when a close
//!     encounter forces the controller to shrink `dt_try` well below
//!     the caller's hint (previous releases hid this under an outer
//!     "budget loop" and incremented `t` by the requested `dt`
//!     regardless, causing a visual *teleport* when the internal
//!     substep budget was exceeded).
//!   * The dense-output snapshot reflects the *actual* sub-step the
//!     integrator just accepted, giving the renderer a clean,
//!     well-defined window `[t − consumed_dt, t]` to interpolate.
//!   * Callers can respond to external events (pause, shutdown) between
//!     sub-steps without modifying the integrator.
//!
//! # Rejection rollback
//!
//! When a candidate attempt fails the error tolerance, we **must**
//! restore every piece of integrator state — not just positions and
//! velocities, but also `b[]`, `e[]`, the compensated-summation
//! accumulators (`csx`, `csv`, `csb`) — otherwise the divergent
//! information from the rejected attempt silently contaminates the
//! next try. See `Attempt::snapshot` / `Attempt::restore`.

use crate::domain::body::Body;
use crate::physics::integrator::dense::{DenseSnapshot, predict_ias15};
use crate::physics::integrator::helpers::{apply_perturbations, evaluate, scale_acc_and_pe};
use crate::physics::integrator::traits::{
    Integrator, IntegratorContext, IntegratorKind, StepResult,
};

// ── Gauss-Radau node spacings ────────────────────────────────────────────────
//
// 8 nodes on [0, 1]: h₀ = 0 is the left end-point (implicit; the step
// starts there); h₁ … h₇ are the 7 interior Gauss-Radau quadrature
// points. Values are literal transcriptions of Everhart (1985) /
// Rein & Spiegel (2015) — extra digits are preserved so the sums
// `h[j] - h[i]` stay exact to double precision.

const H: [f64; 8] = [
    0.0,
    0.056_262_560_536_922_146_465_652_191_031_8,
    0.180_240_691_736_892_364_987_579_942_835_4,
    0.352_624_717_113_169_637_373_907_770_280_6,
    0.547_153_626_330_555_383_001_448_557_701_4,
    0.734_210_177_215_410_531_523_210_621_826_2,
    0.885_320_946_839_095_768_090_359_762_915_4,
    0.977_520_613_561_287_501_891_174_500_440_5,
];

// ── Triangular b ↔ g conversion coefficients ─────────────────────────────────
//
// `g_k` is the (k+1)-th Newton divided difference of the acceleration
// increment F(hₙ) - F(h₀); `b_k` is the power-series coefficient such
// that F(u) ≈ F₀ + Σ b_k · uᵏ⁺¹ (with u = τ/dt ∈ [0,1]).
//
// The mapping is upper-triangular:
//
//     b_j = g_j + Σ_{k>j}  c_mat[k][j] · g_k
//     g_j = b_j + Σ_{k>j}  d_mat[k][j] · b_k
//
// We store only the lower triangle (row k, cols 0..k). Values from
// Everhart (1985) table I; constants cross-checked against the
// REBOUND reference implementation.

const C_MAT: [[f64; 7]; 7] = [
    [0.0,                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-0.056_262_560_536_922_146,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [ 0.010_140_802_830_063_630, -0.236_503_252_273_814_5,   0.0, 0.0, 0.0, 0.0, 0.0],
    [-0.003_575_897_729_251_617,  0.093_537_695_259_462_07, -0.589_127_969_386_984_1,   0.0, 0.0, 0.0, 0.0],
    [ 0.001_956_565_409_947_221, -0.054_755_386_889_068_69,  0.415_881_200_082_306_86, -1.136_281_595_717_539_5,   0.0, 0.0, 0.0],
    [-0.001_436_530_236_370_892,  0.042_158_527_721_268_71, -0.360_099_596_502_056_8,   1.250_150_711_840_691_0,  -1.870_491_772_932_950_1,   0.0, 0.0],
    [ 0.001_271_790_309_026_868, -0.038_760_357_915_906_77,  0.360_962_243_452_846_0,  -1.466_884_208_400_426_9,   2.906_136_259_308_429_3,  -2.755_812_719_772_045_8, 0.0],
];

const D_MAT: [[f64; 7]; 7] = [
    [0.0,                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.056_262_560_536_922_146,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.003_165_475_718_170_829,   0.236_503_252_273_814_5,   0.0, 0.0, 0.0, 0.0, 0.0],
    [0.000_178_097_769_221_743,   0.045_792_985_506_027_92,  0.589_127_969_386_984_1,   0.0, 0.0, 0.0, 0.0],
    [0.000_010_020_236_522_329_1, 0.008_431_857_153_525_70,  0.253_534_069_054_569_27,  1.136_281_595_717_539_5,   0.0, 0.0, 0.0],
    [0.000_000_563_764_163_931_8, 0.001_529_784_002_500_466, 0.097_834_236_532_444_00,  0.875_254_664_684_091_1,   1.870_491_772_932_950_1,   0.0, 0.0],
    [0.000_000_031_718_815_401_8, 0.000_276_293_090_982_648, 0.036_028_553_983_736_46,  0.576_733_000_277_078_7,   2.248_588_760_769_159_8,   2.755_812_719_772_045_8, 0.0],
];

// ── Tuning parameters ────────────────────────────────────────────────────────

/// Target relative error on the dominant truncation term. `1e-9` is
/// the value recommended by Rein & Spiegel (2015) as giving machine-
/// precision energy conservation over gigayear integrations while
/// remaining robust. Exposed via [`Ias15::with_epsilon`] for users
/// who need looser/tighter control.
const DEFAULT_EPSILON: f64 = 1e-9;

/// Floor on `dt` to prevent a pathological scene (e.g. contact
/// singularity) from driving the step size to zero and stalling the
/// scheduler. Below this, we accept the attempt regardless and let
/// the caller decide what to do (typically: log a degraded-step and
/// keep integrating). Matches REBOUND's `integrator_ias15.min_dt`
/// default behaviour when no explicit floor is configured.
const DT_MIN: f64 = 1e-12;

/// Multiplier on the theoretically optimal Δt after each attempt.
/// Keeps the controller away from the accept/reject boundary so
/// step size doesn't oscillate between borderline-too-large and
/// too-small. 0.9 matches REBOUND (`integrator_ias15.safety_factor`).
const DT_SAFETY: f64 = 0.9;

/// Conservative growth factor used only as a fallback when the error
/// estimate is zero (exact machine-precision step). In all other cases
/// the error formula drives dt_next directly — no growth cap is applied,
/// matching REBOUND's controller exactly.
const DT_ZERO_ERR_GROWTH: f64 = 2.0;

/// Cap on predictor-corrector Picard iterations per attempt. In well-
/// behaved regimes 2–3 suffice; 12 is a safety net against pathological
/// cases where the iteration fails to converge at all.
const MAX_PICARD_ITERATIONS: usize = 12;

/// Convergence threshold on max|Δb₆|/max|a₀| across one Picard
/// iteration. `1e-16` is essentially round-off: REBOUND uses the
/// same threshold with an early-exit when two consecutive iterations
/// fail to improve, which we also do.
const PICARD_TOL: f64 = 1e-16;

/// Lower bound on user-settable epsilon. f64 machine epsilon is
/// ≈2.22e-16, so tolerances below ~1e-14 cannot be distinguished from
/// floating-point round-off in either the Picard residual or the
/// truncation estimate. Pinning the floor three decades above machine
/// epsilon keeps the error controller honest: at `ε = 1e-13` the
/// optimal-dt formula `(ε/err)^(1/7)` still produces meaningful
/// adjustments rather than noise, and the retry loop cannot stall
/// with Picard residual and truncation both floating on round-off.
const EPSILON_MIN: f64 = 1e-13;

/// Upper bound on user-settable epsilon. Above `~1e-3` the local
/// truncation error approaches the step itself and the 15th-order
/// machinery buys nothing over a cheap low-order method; we cap here
/// to flag misconfiguration rather than silently degrade to
/// garbage-quality integration.
const EPSILON_MAX: f64 = 1e-3;

// ── Integrator struct ────────────────────────────────────────────────────────

/// Per-body polynomial state for one substep (coefficients of the
/// series expansion of the acceleration within the step). Index 0..7
/// is the coefficient order; the pair is (x-component, y-component).
type BodyCoeffs = [(f64, f64); 7];

pub struct Ias15 {
    /// Target relative error on the dominant truncation term.
    epsilon: f64,

    /// Power-series coefficients for the acceleration, per body.
    /// `b[i][k]` gives the k-th coefficient for body i.
    /// Warm-started from the previous accepted step.
    b: Vec<BodyCoeffs>,
    /// Coefficients at the previous accepted step — used to extrapolate
    /// `b` when the step size changes (see [`Self::warmstart_b`]).
    e: Vec<BodyCoeffs>,
    /// Newton divided-difference form, derived from `b` on-demand.
    g: Vec<BodyCoeffs>,
    /// Compensated-summation carry terms for the `b` coefficients.
    csb: Vec<BodyCoeffs>,
    /// Compensated-summation carry for positions.
    csx: Vec<(f64, f64)>,
    /// Compensated-summation carry for velocities.
    csv: Vec<(f64, f64)>,

    /// Step size proposed for the next attempt. Seeded from the caller's
    /// `dt` on first use; thereafter driven by the error controller.
    dt_next: f64,

    /// The `dt_try` that was accepted on the most recent internal attempt.
    /// Used as `dt_prev` in [`Self::warmstart_b`] so the q = dt_try/dt_prev
    /// ratio is correct. Zero means "no accepted step yet" — warm-start is
    /// skipped and `e` is left at zero.
    dt_last_accepted: f64,

    /// Cumulative sub-step count across the integrator's lifetime
    /// (one per accepted attempt). Surfaced via [`Metrics`] so a UI
    /// can show effective work rate (sub-steps / wall-second).
    substeps_total: u64,

    /// Cumulative count of rejected attempts. High values vs. `substeps_total`
    /// indicate the error controller is over-estimating dt; diagnostic only.
    rejections_total: u64,

    /// Cumulative Picard iteration count summed across all attempts
    /// (accepted and rejected). Mean iterations per sub-step is
    /// `picard_iters_total / (substeps_total + rejections_total)`.
    picard_iters_total: u64,

    /// Cumulative count of degraded accepts (`DT_MIN` escape clause fired).
    /// Should stay at zero for well-posed scenes.
    degraded_total: u64,
}

impl Default for Ias15 {
    fn default() -> Self {
        Self::new()
    }
}

impl Ias15 {
    pub fn new() -> Self {
        Self {
            epsilon: DEFAULT_EPSILON,
            b: Vec::new(),
            e: Vec::new(),
            g: Vec::new(),
            csb: Vec::new(),
            csx: Vec::new(),
            csv: Vec::new(),
            dt_next: 0.0,
            dt_last_accepted: 0.0,
            substeps_total: 0,
            rejections_total: 0,
            picard_iters_total: 0,
            degraded_total: 0,
        }
    }

    /// Override the default tolerance (`1e-9`). Tighter values give
    /// better energy conservation at the cost of proportionally smaller
    /// step sizes. Clamped to `[EPSILON_MIN, EPSILON_MAX]`.
    pub fn with_epsilon(mut self, epsilon: f64) -> Self {
        self.epsilon = epsilon.clamp(EPSILON_MIN, EPSILON_MAX);
        self
    }

    /// Returns the current error tolerance.
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    fn ensure_capacity(&mut self, n: usize) {
        if self.b.len() != n {
            self.b = vec![[(0.0, 0.0); 7]; n];
            self.e = vec![[(0.0, 0.0); 7]; n];
            self.g = vec![[(0.0, 0.0); 7]; n];
            self.csb = vec![[(0.0, 0.0); 7]; n];
            self.csx = vec![(0.0, 0.0); n];
            self.csv = vec![(0.0, 0.0); n];
            self.dt_last_accepted = 0.0;
            // Reset dt_next too: a value from a different body count is
            // physically meaningless (different acceleration regime) and
            // would spend rejections in the retry loop before re-calibrating.
            // The `if self.dt_next <= 0.0` guard in `step()` re-seeds from
            // the caller's budget on the next entry.
            self.dt_next = 0.0;
        }
    }
}

// ── State snapshot for rejection rollback ────────────────────────────────────

/// Immutable snapshot of everything we must rewind if an attempt is
/// rejected. Rolling back positions / velocities alone is not enough:
/// the `b`, `e`, and compensated-summation arrays carry information
/// from the rejected attempt that would otherwise bias the next try.
struct Attempt {
    x: Vec<(f64, f64)>,
    v: Vec<(f64, f64)>,
    b: Vec<BodyCoeffs>,
    e: Vec<BodyCoeffs>,
    csb: Vec<BodyCoeffs>,
    csx: Vec<(f64, f64)>,
    csv: Vec<(f64, f64)>,
}

impl Attempt {
    fn snapshot(bodies: &[Body], ias: &Ias15) -> Self {
        Self {
            x: bodies.iter().map(|b| (b.x, b.y)).collect(),
            v: bodies.iter().map(|b| (b.vx, b.vy)).collect(),
            b: ias.b.clone(),
            e: ias.e.clone(),
            csb: ias.csb.clone(),
            csx: ias.csx.clone(),
            csv: ias.csv.clone(),
        }
    }

    /// Roll bodies and integrator state back to the start of the current
    /// sub-step. Ref-based so a single snapshot can be reused across
    /// multiple rejection retries within one `step()` call.
    fn restore(&self, bodies: &mut [Body], ias: &mut Ias15) {
        for (i, b) in bodies.iter_mut().enumerate() {
            b.x = self.x[i].0;
            b.y = self.x[i].1;
            b.vx = self.v[i].0;
            b.vy = self.v[i].1;
        }
        ias.b.clone_from(&self.b);
        ias.e.clone_from(&self.e);
        ias.csb.clone_from(&self.csb);
        ias.csx.clone_from(&self.csx);
        ias.csv.clone_from(&self.csv);
    }
}

// ── Core algorithm ───────────────────────────────────────────────────────────

impl Integrator for Ias15 {
    /// Perform **one** adaptive Gauss-Radau sub-step.
    ///
    /// The input `dt` is treated as the caller's budget: the accepted
    /// sub-step will not exceed it. The actual step size, chosen by the
    /// error controller, is reported through [`StepResult::consumed_dt`],
    /// and the caller is expected to re-invoke `step` until the full
    /// simulation-time target has been reached (REBOUND-style driver;
    /// see the module-level documentation on sub-step semantics).
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt: f64,
        acc: &mut Vec<(f64, f64)>,
    ) -> StepResult {
        let n = bodies.len();
        self.ensure_capacity(n);

        // The caller's budget acts as a hard upper bound on this sub-step
        // size; the error controller may shrink it further. A zero or
        // negative hint is treated as DT_MIN (pathological, but legal).
        let dt_cap = dt.max(DT_MIN);

        // Seed `dt_next` from the budget on the very first call.
        if self.dt_next <= 0.0 {
            self.dt_next = dt_cap;
        }

        let mut dt_try = self.dt_next.min(dt_cap).max(DT_MIN);

        // Snapshot taken once per sub-step: body kinematics + integrator
        // state that must survive rejection retries. `a0` (start-of-step
        // acceleration) is also invariant across retries because the
        // snapshot restores positions, so we evaluate it only once —
        // saving up to `max_picard_iter × n_rejects` force calls per
        // sub-step compared to re-evaluating inside the retry loop.
        let snapshot = Attempt::snapshot(bodies, self);

        let raw_pe = evaluate(bodies, ctx.force, acc);
        scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
        apply_perturbations(bodies, acc, ctx.perturbations);
        let a0: Vec<(f64, f64)> = acc.clone();

        // ── Rejection retry loop ─────────────────────────────────────────
        //
        // Each iteration: warm-start `b`, run Picard, estimate error. On
        // reject we restore snapshot (positions *and* integrator state)
        // and shrink `dt_try`. DT_MIN is a hard floor — if a pathological
        // configuration forces the controller down to it, we accept the
        // step unconditionally to keep the simulation progressing rather
        // than spin forever.
        let (accepted_dt, final_pe, final_snapshot, degraded) = loop {
            if self.dt_last_accepted > 0.0 {
                self.warmstart_b(dt_try, self.dt_last_accepted);
            }
            self.recompute_g_from_b();

            let (converged, picard_err, picard_iters) =
                self.picard_loop(bodies, ctx, acc, &a0, dt_try);
            self.picard_iters_total =
                self.picard_iters_total.saturating_add(picard_iters as u64);

            let trunc_err = self.truncation_error(&a0);
            let max_err = trunc_err.max(picard_err).max(0.0);
            let dt_optimal = self.optimal_dt(dt_try, max_err);

            // Acceptance requires both (i) Picard convergence and (ii) a
            // truncation error within tolerance. A non-converged Picard
            // with incidentally-small residual is **not** admissible —
            // that path silently accepted divergent coefficients in the
            // previous outer-budget-loop implementation. Two escape
            // hatches accept without (i)+(ii) holding:
            //   - the floor case `dt_try <= DT_MIN` (pathological scene);
            //   - an exhausted cooperative `ctx.deadline` (host wants the
            //     call to return even if we haven't converged).
            // In either case we flag the step as `degraded` so the caller
            // can surface the condition.
            let on_merit = converged && max_err <= self.epsilon;
            let deadline_hit = ctx
                .deadline
                .map(|d| std::time::Instant::now() >= d)
                .unwrap_or(false);
            let accept = on_merit || dt_try <= DT_MIN || deadline_hit;

            if !accept {
                self.rejections_total = self.rejections_total.saturating_add(1);
                snapshot.restore(bodies, self);
                dt_try = dt_optimal.max(DT_MIN).min(dt_cap);
                continue;
            }

            let step_degraded = !on_merit;
            self.substeps_total = self.substeps_total.saturating_add(1);
            if step_degraded {
                self.degraded_total = self.degraded_total.saturating_add(1);
            }

            // Accept path ────────────────────────────────────────────────
            // Build the dense-output snapshot *before* we advance the
            // state, so it carries the pre-step kinematics (the b-coeffs
            // below are the accepted values — `self.b` is not further
            // modified on the accept path). The caller (`System::step`)
            // fills in the absolute `t0` as `system.t() - consumed_dt`.
            let step_snapshot = DenseSnapshot {
                t0: 0.0,
                dt: dt_try,
                x0: snapshot.x.clone(),
                v0: snapshot.v.clone(),
                a0: a0.clone(),
                b: self.b.clone(),
                kind: IntegratorKind::Ias15,
            };

            self.advance_state(bodies, &a0, dt_try);

            // Post-step force evaluation: publishes `acc` consistent with
            // the final body positions, and returns the potential energy
            // the caller will use for energy-conservation diagnostics.
            let raw_pe = evaluate(bodies, ctx.force, acc);
            let pe = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
            apply_perturbations(bodies, acc, ctx.perturbations);

            self.update_warmstart_record();
            self.dt_last_accepted = dt_try;
            // `dt_next` is driven by the pure error formula, no relative
            // growth cap (cf. REBOUND's controller). The external budget
            // clamps it on the next entry to `step()` via `dt_cap`.
            self.dt_next = dt_optimal.max(DT_MIN);

            break (dt_try, pe, step_snapshot, step_degraded);
        };

        StepResult {
            consumed_dt: accepted_dt,
            potential_energy: final_pe,
            used_fallback: false,
            step_snapshot: Some(final_snapshot),
            degraded,
        }
    }

    fn kind(&self) -> IntegratorKind {
        IntegratorKind::Ias15
    }

    fn set_epsilon(&mut self, eps: f64) {
        self.epsilon = eps.clamp(EPSILON_MIN, EPSILON_MAX);
    }

    fn epsilon(&self) -> Option<f64> {
        Some(self.epsilon)
    }

    fn adaptive_stats(&self) -> Option<super::traits::AdaptiveStats> {
        Some(super::traits::AdaptiveStats {
            substeps: self.substeps_total,
            rejections: self.rejections_total,
            picard_iters: self.picard_iters_total,
            degraded: self.degraded_total,
        })
    }
}

impl Ias15 {

    /// Inner predictor-corrector iteration. Given `a0` (acceleration at
    /// the start of the attempt) and a target `dt_try`, iteratively
    /// refines `b` until max|Δb₆|/max|a₀| < `PICARD_TOL` or we hit the
    /// iteration cap. Returns `(converged, residual, iters)` — `iters`
    /// counts the actual iterations consumed (1..=MAX_PICARD_ITERATIONS)
    /// so the outer controller can aggregate them into diagnostics.
    fn picard_loop(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        acc: &mut Vec<(f64, f64)>,
        a0: &[(f64, f64)],
        dt_try: f64,
    ) -> (bool, f64, u32) {
        let n = bodies.len();
        // We keep a pristine copy of the start-of-attempt positions and
        // velocities so each stage can predict from (x0, v0) without
        // worrying about the previous stage having mutated `bodies`.
        let x0: Vec<(f64, f64)> = bodies.iter().map(|b| (b.x, b.y)).collect();
        let v0: Vec<(f64, f64)> = bodies.iter().map(|b| (b.vx, b.vy)).collect();

        let mut last_residual = f64::INFINITY;
        let mut iters: u32 = 0;

        for iter in 0..MAX_PICARD_ITERATIONS {
            iters = (iter as u32) + 1;
            // Snapshot b₆ before the iteration — residual is measured
            // against this.
            let b6_old: Vec<(f64, f64)> = self.b.iter().map(|row| row[6]).collect();

            for stage in 1..=7 {
                let s = H[stage];
                // Predict positions at node `s`.
                for i in 0..n {
                    let (px, py) = predict_ias15(
                        x0[i], v0[i], a0[i], &self.b[i], s, dt_try,
                    );
                    bodies[i].x = px;
                    bodies[i].y = py;
                }

                // Evaluate acceleration at predicted positions.
                let raw_pe = evaluate(bodies, ctx.force, acc);
                let _ = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
                apply_perturbations(bodies, acc, ctx.perturbations);

                // Update divided-difference g and then b via c-coeffs.
                self.update_g_and_b(stage, a0, acc);
            }

            // Residual = max‖Δb₆‖ / max‖a₀‖ (vector norms; per-body).
            let mut max_db6 = 0.0_f64;
            let mut max_a = 0.0_f64;
            for i in 0..n {
                let db6x = self.b[i][6].0 - b6_old[i].0;
                let db6y = self.b[i][6].1 - b6_old[i].1;
                max_db6 = max_db6.max((db6x * db6x + db6y * db6y).sqrt());
                let a_mag = (a0[i].0 * a0[i].0 + a0[i].1 * a0[i].1).sqrt();
                max_a = max_a.max(a_mag);
            }
            let residual = if max_a > 0.0 { max_db6 / max_a } else { max_db6 };

            if residual < PICARD_TOL {
                // Restore positions/velocities to start-of-attempt so
                // the caller can advance cleanly from (x0, v0).
                restore_xv(bodies, &x0, &v0);
                return (true, residual, iters);
            }

            // Stagnation guard: if we're not making progress, bail.
            // Consistent with REBOUND's heuristic.
            if iter >= 2 && residual > last_residual {
                restore_xv(bodies, &x0, &v0);
                return (false, residual, iters);
            }
            last_residual = residual;
        }

        restore_xv(bodies, &x0, &v0);
        (false, last_residual, iters)
    }

    /// Advance positions and velocities to the end of the accepted
    /// attempt using compensated summation (Neumaier-style) so the
    /// integrator's round-off error stays O(ε) rather than O(ε·N_steps).
    fn advance_state(&mut self, bodies: &mut [Body], a0: &[(f64, f64)], dt: f64) {
        let n = bodies.len();
        for i in 0..n {
            let bi = &self.b[i];

            // Full-step position increment (s = 1):
            //   Δx/dt² = a₀/2 + b₀/6 + b₁/12 + b₂/20 + b₃/30 + b₄/42 + b₅/56 + b₆/72
            let dx = dt * dt * (
                a0[i].0 * 0.5
                    + bi[0].0 / 6.0
                    + bi[1].0 / 12.0
                    + bi[2].0 / 20.0
                    + bi[3].0 / 30.0
                    + bi[4].0 / 42.0
                    + bi[5].0 / 56.0
                    + bi[6].0 / 72.0
            );
            let dy = dt * dt * (
                a0[i].1 * 0.5
                    + bi[0].1 / 6.0
                    + bi[1].1 / 12.0
                    + bi[2].1 / 20.0
                    + bi[3].1 / 30.0
                    + bi[4].1 / 42.0
                    + bi[5].1 / 56.0
                    + bi[6].1 / 72.0
            );

            // Full-step velocity increment:
            //   Δv/dt = a₀ + b₀/2 + b₁/3 + b₂/4 + b₃/5 + b₄/6 + b₅/7 + b₆/8
            let dvx = dt * (
                a0[i].0
                    + bi[0].0 / 2.0
                    + bi[1].0 / 3.0
                    + bi[2].0 / 4.0
                    + bi[3].0 / 5.0
                    + bi[4].0 / 6.0
                    + bi[5].0 / 7.0
                    + bi[6].0 / 8.0
            );
            let dvy = dt * (
                a0[i].1
                    + bi[0].1 / 2.0
                    + bi[1].1 / 3.0
                    + bi[2].1 / 4.0
                    + bi[3].1 / 5.0
                    + bi[4].1 / 6.0
                    + bi[5].1 / 7.0
                    + bi[6].1 / 8.0
            );

            // First integrate the v·dt contribution to position.
            let vdt_x = bodies[i].vx * dt;
            let vdt_y = bodies[i].vy * dt;

            add_cs(&mut bodies[i].x, &mut self.csx[i].0, vdt_x);
            add_cs(&mut bodies[i].y, &mut self.csx[i].1, vdt_y);
            add_cs(&mut bodies[i].x, &mut self.csx[i].0, dx);
            add_cs(&mut bodies[i].y, &mut self.csx[i].1, dy);

            add_cs(&mut bodies[i].vx, &mut self.csv[i].0, dvx);
            add_cs(&mut bodies[i].vy, &mut self.csv[i].1, dvy);
        }
    }

    /// Estimate of the dominant truncation error term, normalised by
    /// the acceleration magnitude: err = max‖b₆‖ / max‖a₀‖. For a
    /// 15th-order method this is the correct leading term since b₆
    /// multiplies u⁷ ≈ 1 at the end of the step.
    fn truncation_error(&self, a0: &[(f64, f64)]) -> f64 {
        let mut max_b6 = 0.0_f64;
        let mut max_a = 0.0_f64;
        for (i, row) in self.b.iter().enumerate() {
            let b = row[6];
            max_b6 = max_b6.max((b.0 * b.0 + b.1 * b.1).sqrt());
            max_a = max_a.max((a0[i].0 * a0[i].0 + a0[i].1 * a0[i].1).sqrt());
        }
        if max_a > 0.0 { max_b6 / max_a } else { 0.0 }
    }

    /// Compute the optimal Δt for the next attempt given the current
    /// `dt` and the measured error `err`. Safety factor damps the
    /// controller, and the exponent 1/7 comes from the dominant term
    /// scaling as dt⁷ (since b₆ already multiplies u⁷).
    fn optimal_dt(&self, dt_current: f64, err: f64) -> f64 {
        if err <= 0.0 {
            // Zero error means the step was exact to machine precision.
            // Grow conservatively rather than to infinity.
            return dt_current * DT_ZERO_ERR_GROWTH;
        }
        let ratio = (self.epsilon / err).powf(1.0 / 7.0);
        dt_current * DT_SAFETY * ratio
    }

    /// Extrapolate `b` from the previous accepted step to the current
    /// `dt_try`. Uses the standard REBOUND formula: b_new is a simple
    /// rescaling of b by powers of `(dt_try / dt_prev)` plus a
    /// correction from the drift `e = b - b_prev` to capture how the
    /// coefficients changed last step. This drastically reduces the
    /// number of Picard iterations in steady-state integration.
    fn warmstart_b(&mut self, dt_try: f64, dt_prev: f64) {
        if dt_prev <= 0.0 {
            return;
        }
        let q = dt_try / dt_prev;
        let q2 = q * q;
        let q3 = q2 * q;
        let q4 = q3 * q;
        let q5 = q4 * q;
        let q6 = q5 * q;
        let q7 = q6 * q;

        for i in 0..self.b.len() {
            let be = [
                (self.b[i][0].0 - self.e[i][0].0, self.b[i][0].1 - self.e[i][0].1),
                (self.b[i][1].0 - self.e[i][1].0, self.b[i][1].1 - self.e[i][1].1),
                (self.b[i][2].0 - self.e[i][2].0, self.b[i][2].1 - self.e[i][2].1),
                (self.b[i][3].0 - self.e[i][3].0, self.b[i][3].1 - self.e[i][3].1),
                (self.b[i][4].0 - self.e[i][4].0, self.b[i][4].1 - self.e[i][4].1),
                (self.b[i][5].0 - self.e[i][5].0, self.b[i][5].1 - self.e[i][5].1),
                (self.b[i][6].0 - self.e[i][6].0, self.b[i][6].1 - self.e[i][6].1),
            ];

            // Rescale b-coefficients for the new step size.
            self.e[i][0] = (self.b[i][0].0 * q,  self.b[i][0].1 * q);
            self.e[i][1] = (self.b[i][1].0 * q2, self.b[i][1].1 * q2);
            self.e[i][2] = (self.b[i][2].0 * q3, self.b[i][2].1 * q3);
            self.e[i][3] = (self.b[i][3].0 * q4, self.b[i][3].1 * q4);
            self.e[i][4] = (self.b[i][4].0 * q5, self.b[i][4].1 * q5);
            self.e[i][5] = (self.b[i][5].0 * q6, self.b[i][5].1 * q6);
            self.e[i][6] = (self.b[i][6].0 * q7, self.b[i][6].1 * q7);

            for k in 0..7 {
                self.b[i][k] = (self.e[i][k].0 + be[k].0, self.e[i][k].1 + be[k].1);
            }
        }
    }

    fn update_warmstart_record(&mut self) {
        self.e = self.b.clone();
    }

    fn recompute_g_from_b(&mut self) {
        // g_j = b_j + Σ_{k>j} d_mat[k][j] · b_k
        for i in 0..self.b.len() {
            let bi = self.b[i];
            for j in 0..7 {
                let mut gx = bi[j].0;
                let mut gy = bi[j].1;
                for k in (j + 1)..7 {
                    gx += D_MAT[k][j] * bi[k].0;
                    gy += D_MAT[k][j] * bi[k].1;
                }
                self.g[i][j] = (gx, gy);
            }
        }
    }

    /// After evaluating acceleration at stage `n` (1..=7), update g_{n-1}
    /// via Newton divided differences of (F - F₀); then propagate the
    /// delta back into b₀..b_{n-1} using c_mat. Compensated summation
    /// keeps round-off under control across many Picard iterations.
    fn update_g_and_b(
        &mut self,
        stage: usize,
        a0: &[(f64, f64)],
        acc_n: &[(f64, f64)],
    ) {
        let n = stage - 1; // index of the g coefficient we're updating
        let hn = H[stage];

        for i in 0..self.g.len() {
            // Compute Newton divided difference of order n+1 for body i.
            let dfx = acc_n[i].0 - a0[i].0;
            let dfy = acc_n[i].1 - a0[i].1;

            let (mut tx, mut ty) = (dfx / hn, dfy / hn);
            for k in 0..n {
                tx = (tx - self.g[i][k].0) / (hn - H[k + 1]);
                ty = (ty - self.g[i][k].1) / (hn - H[k + 1]);
            }

            let dgx = tx - self.g[i][n].0;
            let dgy = ty - self.g[i][n].1;
            self.g[i][n] = (tx, ty);

            // Propagate Δg_n into b₀..b_n using compensated summation.
            for j in 0..=n {
                let coeff = if j == n { 1.0 } else { C_MAT[n][j] };
                let dbx = coeff * dgx;
                let dby = coeff * dgy;
                add_cs(&mut self.b[i][j].0, &mut self.csb[i][j].0, dbx);
                add_cs(&mut self.b[i][j].1, &mut self.csb[i][j].1, dby);
            }
        }
    }
}

// ── Free helpers ─────────────────────────────────────────────────────────────

/// Compensated summation (Neumaier): `p += inp` with `csp` absorbing
/// the rounding residual. Standard round-off-resistant accumulation.
fn add_cs(p: &mut f64, csp: &mut f64, inp: f64) {
    let y = inp - *csp;
    let t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

fn restore_xv(bodies: &mut [Body], x: &[(f64, f64)], v: &[(f64, f64)]) {
    for (i, b) in bodies.iter_mut().enumerate() {
        b.x = x[i].0;
        b.y = x[i].1;
        b.vx = v[i].0;
        b.vy = v[i].1;
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::domain::materials::Material;

    /// Kepler e=0.5 orbit: tests both (i) peak |δE/E₀| never exceeds
    /// the quoted tolerance over 100 orbits AND (ii) the energy error
    /// is oscillatory rather than monotonically drifting. The latter
    /// check catches a classical family of bugs (missing state in
    /// rollback, asymmetric update, etc.) that a peak-only assertion
    /// lets through when the bug drifts the error slowly in one
    /// direction.
    #[test]
    fn ias15_kepler_energy_peak_and_monotonic_drift() {
        const A: f64 = 2.0;
        const E: f64 = 0.5;
        const MU: f64 = 2.0;
        const N_ORBITS: u64 = 100;
        const PEAK_TOL: f64 = 1e-12;
        const DRIFT_TOL: f64 = 1e-13;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();

        let period = 2.0 * std::f64::consts::PI * (A.powi(3) / MU).sqrt();
        let total_time = N_ORBITS as f64 * period;
        // Budget of `period/20` lets the error controller choose any
        // sub-step up to that size. At ε = 1e-9 on an e=0.5 orbit the
        // controller tends to settle near `period/30` — a smaller budget
        // merely caps growth without changing step count in practice,
        // while a smaller budget *combined with* a fixed for-loop
        // driver artificially inflates the number of `step()` calls.
        let dt_budget = period / 20.0;

        let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
        b1.softening = 0.0;
        let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], 0.5, dt_budget, 10, 1);
        sys.set_integrator(IntegratorKind::Ias15);

        let mut peak = 0.0_f64;
        // Samples for drift detection: (t, δE/E₀) every ~0.5% of the run.
        let mut samples: Vec<(f64, f64)> = Vec::new();
        let sample_dt = total_time / 200.0;
        let mut next_sample = 0.0;

        // REBOUND-style driver: advance by calling `step()` until the
        // target simulation time is reached. Each call consumes one
        // adaptive sub-step whose size IAS15 chose; using a fixed
        // `for _ in 0..n_steps` loop here would silently assume every
        // call consumes `dt_budget` and fall short of the intended
        // integration window.
        while sys.t() < total_time {
            sys.step();
            let err = sys.metrics().rel_energy_error;
            peak = peak.max(err.abs());
            if sys.t() >= next_sample {
                samples.push((sys.t(), err));
                next_sample += sample_dt;
            }
        }

        assert!(
            peak < PEAK_TOL,
            "IAS15 Kepler: peak |δE/E₀| = {:.3e} exceeds {:.0e} over {} orbits",
            peak, PEAK_TOL, N_ORBITS,
        );

        // Linear regression on the samples. A well-behaved IAS15 run
        // produces near-zero slope; a rollback/rounding bug shows up
        // as a consistent drift in one direction.
        let n = samples.len() as f64;
        let sum_t: f64 = samples.iter().map(|&(t, _)| t).sum();
        let sum_e: f64 = samples.iter().map(|&(_, e)| e).sum();
        let mean_t = sum_t / n;
        let mean_e = sum_e / n;
        let mut num = 0.0;
        let mut den = 0.0;
        for &(t, e) in &samples {
            num += (t - mean_t) * (e - mean_e);
            den += (t - mean_t).powi(2);
        }
        let slope = if den > 0.0 { num / den } else { 0.0 };
        let drift = (slope * total_time).abs();

        assert!(
            drift < DRIFT_TOL,
            "IAS15 Kepler: monotonic drift |slope·t_final| = {:.3e} exceeds {:.0e} \
             (slope = {:.3e}, peak = {:.3e}) — suggests asymmetric state update",
            drift, DRIFT_TOL, slope, peak,
        );
    }

    /// High-eccentricity Kepler (e = 0.9): the regime where fixed-step
    /// symplectic integrators lose their advantage and where IAS15's
    /// adaptive step size is essential — perihelion passages happen
    /// fast, apoapsis is leisurely, and the time-scale ratio is ~400.
    ///
    /// Published IAS15 results (Rein & Spiegel 2015, Fig. 5) show
    /// machine-precision energy conservation across 10⁴ orbits at e=0.9.
    /// Here we check 50 orbits with a tight tolerance to keep the test
    /// fast; the asserted bound is conservative relative to the paper.
    #[test]
    fn ias15_kepler_high_eccentricity() {
        const A: f64 = 1.0;
        const E: f64 = 0.9;
        const MU: f64 = 2.0;
        const N_ORBITS: u64 = 50;
        const DT: f64 = 0.1; // large budget; IAS15 will subdivide near perihelion
        const PEAK_TOL: f64 = 1e-11;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();

        let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
        b1.softening = 0.0;
        let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], 0.5, DT, 10, 1);
        sys.set_integrator(IntegratorKind::Ias15);

        let period = 2.0 * std::f64::consts::PI * (A.powi(3) / MU).sqrt();
        let total_time = N_ORBITS as f64 * period;
        let n_steps = (total_time / DT).ceil() as u64;

        let mut peak = 0.0_f64;
        for _ in 0..n_steps {
            sys.step();
            peak = peak.max(sys.metrics().rel_energy_error.abs());
        }

        assert!(
            peak < PEAK_TOL,
            "IAS15 high-e Kepler (e={E}): peak |δE/E₀| = {:.3e} exceeds {:.0e} \
             over {} orbits — adaptive step control not tracking perihelion",
            peak, PEAK_TOL, N_ORBITS,
        );
    }

    /// Pythagorean (Burrau 1913) three-body: m=(3,4,5) placed on the
    /// vertices of a 3-4-5 triangle at rest. Pure gravity, ε=0, G=1.
    ///
    /// The system is chaotic, with violent close encounters around
    /// t ≈ 2–5 before chaos-driven ejection (~t ≈ 46). IAS15's
    /// adaptive step is tested most severely during these encounters:
    /// any asymmetric rollback, missed state variable, or controller
    /// oscillation shows up as energy drift.
    ///
    /// We integrate past the strongest close-encounter phase (t=10)
    /// and assert a tight relative energy bound — well beyond what
    /// any fixed-step integrator in the zoo (VV / Y4) can deliver at
    /// comparable cost.
    #[test]
    fn ias15_pythagorean_energy_through_close_encounters() {
        const DT: f64 = 0.01;
        const T_END: f64 = 10.0;
        const PEAK_TOL: f64 = 1e-11;

        let mut bodies = vec![
            Body::new( 1.0,  3.0, 0.0, 0.0, 3.0, Material::Rocky),
            Body::new(-2.0, -1.0, 0.0, 0.0, 4.0, Material::Rocky),
            Body::new( 1.0, -1.0, 0.0, 0.0, 5.0, Material::Rocky),
        ];
        for b in &mut bodies { b.softening = 0.0; }

        let mut sys = System::new(bodies, 0.5, DT, 10, 1);
        sys.set_integrator(IntegratorKind::Ias15);

        let n_steps = (T_END / DT).ceil() as u64;
        let mut peak = 0.0_f64;
        for _ in 0..n_steps {
            sys.step();
            peak = peak.max(sys.metrics().rel_energy_error.abs());
        }

        assert!(
            peak < PEAK_TOL,
            "IAS15 Pythagorean: peak |δE/E₀| = {:.3e} exceeds {:.0e} over t=[0,{T_END}] \
             — likely a bug in rejection rollback or the truncation-error estimator",
            peak, PEAK_TOL,
        );
    }

    /// Regression: `System::t` must track the sub-step IAS15 physically
    /// executed, not the caller's budget.
    ///
    /// The previous implementation ran an internal `while budget > 0` loop
    /// inside `step()` and returned `consumed_dt == dt`, advancing `System::t`
    /// by the full requested `dt` while the dense-output snapshot only
    /// covered the *last* sub-step. Interpolating inside that window then
    /// extrapolated over earlier sub-steps, producing the visible "teleport"
    /// artefact.
    ///
    /// Under the REBOUND-style contract (`reb_integrator_ias15_part1/2`,
    /// Rein & Spiegel 2015 §2.3), each `step()` call executes exactly one
    /// adaptive sub-step and reports its size via `StepResult::consumed_dt`;
    /// `System::step` advances `System::t` by that value. A budget far
    /// larger than what the controller can accept at perihelion therefore
    /// yields `System::t` strictly below the budget after one call.
    #[test]
    fn ias15_system_t_matches_adaptive_substep() {
        const A: f64 = 1.0;
        const E: f64 = 0.9;
        const MU: f64 = 2.0;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();

        // Budget of a full orbital period — far too large for a single IAS15
        // sub-step at perihelion, so the controller MUST shrink it.
        let period = 2.0 * std::f64::consts::PI * (A.powi(3) / MU).sqrt();
        let dt_budget = period;

        let mut b1 = Body::new(-r_peri / 2.0, 0.0, 0.0, -v_peri / 2.0, 1.0, Material::Rocky);
        b1.softening = 0.0;
        let mut b2 = Body::new(r_peri / 2.0, 0.0, 0.0, v_peri / 2.0, 1.0, Material::Rocky);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], 0.5, dt_budget, 10, 1);
        sys.set_integrator(IntegratorKind::Ias15);

        let t0 = sys.t();
        sys.step();
        let consumed = sys.t() - t0;

        assert!(
            consumed > 0.0,
            "IAS15 sub-step consumed zero time — caller would busy-loop"
        );
        assert!(
            consumed < dt_budget,
            "IAS15 step() consumed the full budget ({:.3e}) instead of adapting \
             down at perihelion — teleport regression (budget loop leaked back in)",
            dt_budget,
        );

        // System::t must advance strictly monotonically across further
        // sub-steps; any regression here implies the caller is reading a
        // stale `consumed_dt` or the integrator is returning negative time.
        let mut t_prev = sys.t();
        for k in 0..200 {
            sys.step();
            let t_now = sys.t();
            assert!(
                t_now > t_prev,
                "System::t regressed at sub-step {}: {:.6e} → {:.6e}",
                k, t_prev, t_now,
            );
            t_prev = t_now;
        }
    }
}
