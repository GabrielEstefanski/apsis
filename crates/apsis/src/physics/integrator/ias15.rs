//! IAS15 — 15th-order adaptive Gauss-Radau integrator.
//!
//! Implementation of the algorithm specified in:
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
//! Section references throughout the module text point at the Rein &
//! Spiegel paper (the algorithmic specification); divergences from the
//! specification — and the empirical analyses motivating each choice —
//! are noted at the relevant call site, with cross-references to the
//! `docs/experiments/` lab notebooks. The independent implementation in
//! REBOUND's `reb_integrator_ias15` provides the cross-implementation
//! parity reference used by the validation suite under
//! [`validation/rebound-parity/`](../../../../validation/rebound-parity/);
//! it is **not** the source of any code in this file.
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
//! # Sub-step semantics (ADR-004, revised 2026-04-26)
//!
//! A single call to [`Integrator::step`] performs **one** adaptive
//! sub-step. The `dt` argument is the *first-call seed* for the
//! controller's `dt_next`, not a per-call hard cap. After the first
//! call the controller is fully responsible for choosing the next
//! step size — the canonical IAS15 semantics from Rein & Spiegel 2015
//! §2.3, where the user's `dt` is treated as an initial guess that
//! the controller then mutates and reuses on subsequent calls.
//!
//! Treating `dt` as a per-call hard cap was the previous behaviour
//! and was a defect: when the orchestrator passed the user's initial
//! guess `user_dt` on every call (the natural pattern for fixed-step
//! integrators in [`crate::core::system::System::step`]), IAS15 was
//! pinned to that guess and unable to grow to the natural step size
//! that local truncation error allowed. On the figure-8 choreography
//! this manifested as a substep cascade where the controller
//! repeatedly shrunk to `DT_MIN` after close approaches and recovered
//! only ~7× per accepted sub-step — never reaching the natural step
//! because every recovery hit the artificial `user_dt` ceiling.
//!
//! Substep-granularity at the trait boundary is preserved for the
//! same reasons that motivated the original design:
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
use crate::math::Vec3;
use crate::physics::integrator::dense::{DenseSnapshot, predict_ias15, predict_v_ias15};
use crate::physics::integrator::helpers::{evaluate, scale_acc_and_pe};
use crate::physics::integrator::operator_dispatch::accumulate_perturbation_forces;
use crate::physics::integrator::traits::{
    Integrator, IntegratorContext, IntegratorKind, StepResult,
};

// ── Phase-timing instrumentation (feature-gated) ─────────────────────────────
//
// When compiled with `--features ias15-profile`, every wrapped phase accumulates
// wall time and call count into a thread-local [`profile::PhaseTimings`]. When
// the feature is off, [`time_phase!`] expands to its block expression unchanged
// — zero call overhead, zero codegen footprint.
//
// The instrumentation is exposed only through [`profile::snapshot`] and
// [`profile::reset`] (feature-gated). The benchmark harness is the only
// consumer; production builds are unaffected.

#[cfg(feature = "ias15-profile")]
pub mod profile {
    use std::cell::RefCell;
    use std::time::Duration;

    #[derive(Default, Debug, Clone, Copy)]
    pub struct PhaseEntry {
        pub total: Duration,
        pub count: u64,
    }

    #[derive(Default, Debug, Clone)]
    pub struct PhaseTimings {
        pub evaluate: PhaseEntry,
        pub update_g_and_b: PhaseEntry,
        pub warmstart_b: PhaseEntry,
        pub recompute_g_from_b: PhaseEntry,
        pub advance_state: PhaseEntry,
        pub residual_compute: PhaseEntry,
        pub snapshot_capture: PhaseEntry,
        pub snapshot_restore: PhaseEntry,
        /// Cost of `let a0 = acc.clone();` at the top of `step()`.
        /// Allocates a fresh `Vec<(f64, f64)>` of length N per sub-step
        /// — an alloc path entirely independent of the persistent
        /// rollback snapshot and therefore not caught by
        /// `snapshot_capture`. Called out because at large N it
        /// becomes a non-trivial fraction of per-sub-step work, and
        /// is an obvious candidate for persistent-buffer reuse.
        pub a0_clone: PhaseEntry,
        /// Cost of constructing the `DenseSnapshot` on the accept
        /// path — the 4 `Vec::clone()` calls that copy x/v/a0/b into
        /// a fresh owned snapshot for downstream consumers (renderer
        /// interpolation). At N=641 this is ~100 KB of alloc+memcpy
        /// per accepted sub-step; at the accept rate the IAS15
        /// controller typically runs, this is the dominant source
        /// of allocator pressure visible as render-thread stutter.
        pub dense_snapshot_build: PhaseEntry,

        /// Wall time of Barnes-Hut **tree construction** specifically.
        /// Set from inside `GravityForceModel::compute` when the
        /// `ias15-profile` feature is compiled into the force model.
        /// Paired with `tree_traverse` it decomposes the `evaluate`
        /// phase into its two structural halves — answering whether
        /// the per-call cost at large N is dominated by rebuilding
        /// the tree (fixable by caching across Picard iterations)
        /// or by traversing it (requires traversal-level SIMD).
        pub tree_build: PhaseEntry,

        /// Wall time of Barnes-Hut **tree traversal** (the force
        /// accumulation itself). Complements `tree_build`; the two
        /// together approximately reconstruct `evaluate` — the gap
        /// between `evaluate - (build + traverse)` is the per-call
        /// function-dispatch / bookkeeping overhead and should be
        /// small.
        pub tree_traverse: PhaseEntry,
    }

    thread_local! {
        pub(super) static TIMINGS: RefCell<PhaseTimings> =
            RefCell::new(PhaseTimings::default());
    }

    /// Snapshot of the current accumulated timings. Cheap clone — the
    /// struct is a few Durations and counters.
    pub fn snapshot() -> PhaseTimings {
        TIMINGS.with(|t| t.borrow().clone())
    }

    /// Zero the accumulator. Called by the bench harness between scenarios
    /// so each scenario's breakdown is reported independently.
    pub fn reset() {
        TIMINGS.with(|t| *t.borrow_mut() = PhaseTimings::default());
    }

    /// Record a Barnes-Hut tree-build sample. Invoked from inside
    /// `GravityForceModel::compute` (i.e. crossing out of
    /// `ias15.rs`), so it needs a free-function entry point rather
    /// than the `time_phase!` macro which is scoped to this file.
    pub fn record_tree_build(elapsed: std::time::Duration) {
        TIMINGS.with(|t| {
            let mut s = t.borrow_mut();
            s.tree_build.total += elapsed;
            s.tree_build.count += 1;
        });
    }

    /// Record a Barnes-Hut tree-traversal sample. Paired with
    /// [`record_tree_build`] across `GravityForceModel::compute`
    /// to split the `evaluate` phase.
    pub fn record_tree_traverse(elapsed: std::time::Duration) {
        TIMINGS.with(|t| {
            let mut s = t.borrow_mut();
            s.tree_traverse.total += elapsed;
            s.tree_traverse.count += 1;
        });
    }
}

#[cfg(feature = "ias15-profile")]
macro_rules! time_phase {
    ($field:ident, $block:block) => {{
        let __profile_start = std::time::Instant::now();
        let __profile_result = $block;
        let __profile_elapsed = __profile_start.elapsed();
        $crate::physics::integrator::ias15::profile::TIMINGS.with(|t| {
            let mut s = t.borrow_mut();
            s.$field.total += __profile_elapsed;
            s.$field.count += 1;
        });
        __profile_result
    }};
}

#[cfg(not(feature = "ias15-profile"))]
macro_rules! time_phase {
    ($field:ident, $block:block) => {{ $block }};
}

// ── Per-step diagnostic trace (feature-gated) ────────────────────────────────
//
// `ias15-diag` compiles in detailed per-step trace emission useful for
// investigating slow-onset cumulative failures of the IAS15 controller —
// the figure-8 cascade documented in
// `docs/experiments/2026-04-26-ias15-warmstart-bug.md` is the
// motivating example. The feature is off by default (zero codegen
// footprint); when on, emission is *additionally* gated on the env var
// `APSIS_IAS15_TRACE=1` so a feature-enabled binary stays silent until
// the operator explicitly opts in.
//
// Each trace line is tab-separated with a stable column schema so
// downstream analysis (Python/pandas, awk pipelines) can parse without
// regex. The columns are documented at the emission site.

#[cfg(feature = "ias15-diag")]
pub mod diag {
    use std::cell::Cell;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    static TRACE_ENABLED: AtomicBool = AtomicBool::new(false);
    static TRACE_INITIALISED: AtomicBool = AtomicBool::new(false);
    static EVENT_CAP: AtomicUsize = AtomicUsize::new(2000);

    thread_local! {
        /// Monotonic event counter so the trace is grep-friendly even
        /// when several runs interleave their output.
        pub(super) static EVENT_COUNTER: Cell<u64> = const { Cell::new(0) };
    }

    /// Read the env var lazily on first call, then memoise. The env var
    /// is the *runtime* gate; the Cargo feature is the *compile-time*
    /// gate. Both must be set for any output.
    ///
    /// Also reads `APSIS_IAS15_TRACE_CAP` if present and overrides the
    /// default per-thread emission cap of 2 000 events. Keeping a cap
    /// matters because cascade scenarios emit ≥10⁸ attempts; without
    /// throttling the trace itself becomes the cost.
    pub(super) fn trace_enabled() -> bool {
        if !TRACE_INITIALISED.load(Ordering::Relaxed) {
            let on = std::env::var("APSIS_IAS15_TRACE").map(|v| v == "1").unwrap_or(false);
            if let Ok(cap) = std::env::var("APSIS_IAS15_TRACE_CAP") {
                if let Ok(n) = cap.parse::<usize>() {
                    EVENT_CAP.store(n, Ordering::Relaxed);
                }
            }
            TRACE_ENABLED.store(on, Ordering::Relaxed);
            TRACE_INITIALISED.store(true, Ordering::Relaxed);
        }
        TRACE_ENABLED.load(Ordering::Relaxed)
    }

    /// Bump and return the next event id, or `None` if the cap is reached.
    /// Callers should skip emission when `None` is returned.
    pub(super) fn next_event() -> Option<u64> {
        EVENT_COUNTER.with(|c| {
            let n = c.get();
            if (n as usize) >= EVENT_CAP.load(Ordering::Relaxed) {
                return None;
            }
            c.set(n + 1);
            Some(n)
        })
    }

    /// Reset the event counter. Useful between scenarios in a benchmark
    /// harness so each scenario's trace has independent ids.
    pub fn reset_events() {
        EVENT_COUNTER.with(|c| c.set(0));
    }
}

/// Emit a warmstart diagnostic line. Invoked *after* `warmstart_b` has
/// run on the IAS15 instance. Reports two L₂ norms over all bodies and
/// components:
///
/// * `b_norm` — the warmstart's full prediction (Pascal cross-terms
///   plus the previous-step Picard residual `be`); this is the b that
///   Picard will refine on the upcoming attempt.
/// * `e_norm` — the pure Pascal-extrapolated piece, before the `be`
///   correction; lets the trace separate "what the polynomial-basis
///   transform alone would predict" from "what we actually feed Picard".
///
/// On smooth motion at `q ≈ 1` both norms track each other and the
/// ratio is near 1; under aggressive `dt` changes (close encounters)
/// the cross-terms shift the prediction substantially and the two
/// diverge. A persistent divergence after Picard refinement (i.e. b
/// norms still inflated post-converge) is the post-bug-fix signature
/// we track to make sure no equivalent regression slips back in.
#[cfg(feature = "ias15-diag")]
fn diag_emit_warmstart(ias: &Ias15, q: f64, dt_try: f64) {
    if !diag::trace_enabled() {
        return;
    }
    let id = match diag::next_event() {
        Some(n) => n,
        None => return,
    };
    let mut b_norm_sq = 0.0_f64;
    for row in &ias.b {
        for c in row {
            b_norm_sq += c.0 * c.0 + c.1 * c.1;
        }
    }
    let mut e_norm_sq = 0.0_f64;
    for row in &ias.e {
        for c in row {
            e_norm_sq += c.0 * c.0 + c.1 * c.1;
        }
    }
    eprintln!(
        "[ias15-diag]\twarmstart\tev={}\tdt_try={:.6e}\tq={:.6e}\tb_norm={:.6e}\te_norm={:.6e}",
        id,
        dt_try,
        q,
        b_norm_sq.sqrt(),
        e_norm_sq.sqrt(),
    );
}

/// Emit a per-attempt diagnostic line *after* `decide_dt`. Captures the
/// actual error signal the controller saw and the action it chose, so a
/// post-mortem trace can correlate stagnation events with truncation
/// rejections and `dt_next` proposals.
#[cfg(feature = "ias15-diag")]
fn diag_emit_attempt(
    ias: &Ias15,
    dt_try: f64,
    dt_next_after: f64,
    trunc_err: f64,
    picard_converged: bool,
    picard_iters: u32,
    decision_label: &'static str,
) {
    if !diag::trace_enabled() {
        return;
    }
    let id = match diag::next_event() {
        Some(n) => n,
        None => return,
    };
    eprintln!(
        "[ias15-diag]\tattempt\tev={}\tsubstep={}\tdt_try={:.6e}\tdt_next={:.6e}\ttrunc_err={:.6e}\tpicard_conv={}\tpicard_iters={}\tstagnations={}\tcycles={}\tdecision={}",
        id,
        ias.substeps_total,
        dt_try,
        dt_next_after,
        trunc_err,
        picard_converged as u8,
        picard_iters,
        ias.picard_stagnations_total,
        ias.shrink_grow_cycles_total,
        decision_label,
    );
}

// No-op shims when the feature is off — the calls compile out entirely.
#[cfg(not(feature = "ias15-diag"))]
fn diag_emit_warmstart(_ias: &Ias15, _q: f64, _dt_try: f64) {}

#[cfg(not(feature = "ias15-diag"))]
fn diag_emit_attempt(
    _ias: &Ias15,
    _dt_try: f64,
    _dt_next_after: f64,
    _trunc_err: f64,
    _picard_converged: bool,
    _picard_iters: u32,
    _decision_label: &'static str,
) {
}

// ── Gauss-Radau node spacings ────────────────────────────────────────────────
//
// 8 nodes on [0, 1]: h₀ = 0 is the left end-point (implicit; the step
// starts there); h₁ … h₇ are the 7 interior Gauss-Radau quadrature
// points. Values are literal transcriptions of Everhart (1985) /
// Rein & Spiegel (2015) — extra digits are preserved so the sums
// `h[j] - h[i]` stay exact to double precision.

const H: [f64; 8] = [
    0.0,
    0.056_262_560_536_922_15,
    0.180_240_691_736_892_36,
    0.352_624_717_113_169_6,
    0.547_153_626_330_555_4,
    0.734_210_177_215_410_5,
    0.885_320_946_839_095_8,
    0.977_520_613_561_287_5,
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
// Everhart (1985) table I; the same constants appear in Rein &
// Spiegel (2015) §2 and in any specification-correct IAS15
// implementation, including the independent C implementation in
// REBOUND used as the parity reference.

const C_MAT: [[f64; 7]; 7] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [-0.056_262_560_536_922_15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.010_140_802_830_063_63, -0.236_503_252_273_814_5, 0.0, 0.0, 0.0, 0.0, 0.0],
    [
        -0.003_575_897_729_251_617,
        0.093_537_695_259_462_07,
        -0.589_127_969_386_984_1,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.001_956_565_409_947_221,
        -0.054_755_386_889_068_69,
        0.415_881_200_082_306_86,
        -1.136_281_595_717_539_5,
        0.0,
        0.0,
        0.0,
    ],
    [
        -0.001_436_530_236_370_892,
        0.042_158_527_721_268_71,
        -0.360_099_596_502_056_8,
        1.250_150_711_840_691,
        -1.870_491_772_932_95,
        0.0,
        0.0,
    ],
    [
        0.001_271_790_309_026_868,
        -0.038_760_357_915_906_77,
        0.360_962_243_452_846,
        -1.466_884_208_400_427,
        2.906_136_259_308_429_3,
        -2.755_812_719_772_045_8,
        0.0,
    ],
];

const D_MAT: [[f64; 7]; 7] = [
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.056_262_560_536_922_15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.003_165_475_718_170_829, 0.236_503_252_273_814_5, 0.0, 0.0, 0.0, 0.0, 0.0],
    [
        0.000_178_097_769_221_743,
        0.045_792_985_506_027_92,
        0.589_127_969_386_984_1,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.000_010_020_236_522_329_1,
        0.008_431_857_153_525_70,
        0.253_534_069_054_569_27,
        1.136_281_595_717_539_5,
        0.0,
        0.0,
        0.0,
    ],
    [
        0.000_000_563_764_163_931_8,
        0.001_529_784_002_500_466,
        0.097_834_236_532_444,
        0.875_254_664_684_091_1,
        1.870_491_772_932_95,
        0.0,
        0.0,
    ],
    [
        0.000_000_031_718_815_401_8,
        0.000_276_293_090_982_648,
        0.036_028_553_983_736_46,
        0.576_733_000_277_078_7,
        2.248_588_760_769_16,
        2.755_812_719_772_045_8,
        0.0,
    ],
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
/// keep integrating). The IAS15 specification (Rein & Spiegel 2015
/// §3.4) leaves this floor as an implementation choice; `1e-12` is
/// three decades above f64 machine epsilon and matches the value
/// used by the reference C implementation when no explicit floor is
/// configured by the caller.
const DT_MIN: f64 = 1e-12;

/// Multiplier on the theoretically optimal Δt after each attempt.
/// Keeps the controller away from the accept/reject boundary so
/// step size doesn't oscillate between borderline-too-large and
/// too-small. `0.9` is the value specified in Rein & Spiegel (2015)
/// §3.4.
const DT_SAFETY: f64 = 0.9;

/// Conservative growth factor used only as a fallback when the error
/// estimate is zero (exact machine-precision step). In all other cases
/// the error formula drives `dt_next` directly, capped above by
/// [`DT_GROWTH_LIMIT`].
const DT_ZERO_ERR_GROWTH: f64 = 2.0;

/// Maximum step-size growth ratio per accepted sub-step.
///
/// The IAS15 specification (Rein & Spiegel 2015 §3.4) caps `dt_new`
/// at `7 × dt_current` to prevent the controller from overshooting
/// the locally-stable step when the previous attempt's truncation
/// error sits well below `ε`. Without this cap the formula
/// `dt_next = dt · 0.9 · (ε/err)^{1/7}` proposes unbounded growth in
/// smooth regions of the trajectory; the next sub-step then attempts
/// a step many times larger than physically warranted, the
/// truncation gate rejects it, and the controller spends a full
/// shrink cascade recovering — a substep pattern that compounds
/// across close encounters in scenarios like the figure-8
/// choreography (see
/// `docs/experiments/2026-04-26-ias15-warmstart-bug.md`,
/// §"controller substep blowup", for the empirical analysis on
/// `apsis`).
///
/// `7.0` is the specification value; tightening it would slow
/// recovery from over-shrinks (post-close-encounter), while
/// loosening it would re-introduce the overshoot pathology.
const DT_GROWTH_LIMIT: f64 = 7.0;

/// Cap on predictor-corrector Picard iterations per attempt. In well-
/// behaved regimes 2–3 suffice; 12 is a safety net against pathological
/// cases where the iteration fails to converge at all.
const MAX_PICARD_ITERATIONS: usize = 12;

/// Convergence threshold on the Picard residual across one
/// predictor–corrector iteration. `1e-16` is essentially f64 round-off
/// — Rein & Spiegel (2015) §3.3 specify exactly this floor and pair
/// it with an early-exit on two consecutive non-improving iterations,
/// which we also do (see [`Ias15::picard_loop_inner`]).
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

/// Shrink factor applied when the Picard predictor–corrector fails to
/// converge. Divergence is a Lipschitz-regime problem (the step is
/// simply too large for the local dynamics); a fixed halving is the
/// canonical IAS15 response (Rein & Spiegel 2015 §3.4) and converges
/// faster than the truncation formula `(ε/err)^{1/7}` when the error
/// comes from non-convergence rather than from dt⁷-scaled truncation.
const PICARD_SHRINK: f64 = 0.5;

// ── Controller decision type ─────────────────────────────────────────────────

/// Outcome of one attempt at an IAS15 sub-step, as decided by
/// [`decide_dt`] from the two independent error signals.
///
/// Keeping this as a small, `PartialEq`-able enum makes the branching
/// trivial to unit-test in isolation, and was the whole point of
/// factoring it out of the retry loop in [`Integrator::step`] — prior
/// revisions collapsed the two signals through `max(picard, trunc)`
/// and so lost the information of *which* signal was failing, which
/// drives the right shrink strategy (TD-004).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DtDecision {
    /// Accept this attempt. `degraded = true` when the acceptance was
    /// forced by an escape hatch (`DT_MIN` floor or cooperative deadline)
    /// rather than convergence + tolerance actually being met; the caller
    /// should then surface [`StepResult::degraded`].
    Accept { degraded: bool },
    /// Picard predictor–corrector did not converge. Apply [`PICARD_SHRINK`]
    /// (a fixed halving); the truncation formula would under-estimate the
    /// needed shrink in this regime because it assumes the dt⁷ scaling
    /// that only holds when Picard has actually produced valid `b`
    /// coefficients.
    RejectPicard,
    /// Picard converged but truncation bound `|b₆|/|a₀|` is above `ε`.
    /// Shrink using the standard controller formula
    /// `dt · safety · (ε / trunc_err)^(1/7)` — this is the well-posed
    /// signal the formula was derived for.
    RejectTruncation,
}

/// Pure decision function for the IAS15 adaptive controller.
///
/// # Why pure
///
/// The retry loop around this function holds `&mut self` on a bunch of
/// coefficient buffers and runs a non-trivial sequence of force
/// evaluations. By isolating the *decision* (which is first-order
/// logic on two floats + two bools) from the *mutation* (which is
/// second-order through the force model), we get cheap unit tests on
/// the control behaviour without standing up a full `System`.
///
/// # Arguments
///
/// - `picard_converged` — whether the predictor–corrector satisfied
///   [`PICARD_TOL`] before hitting [`MAX_PICARD_ITERATIONS`]. Note this
///   is **independent** of `trunc_err`: Picard can diverge while
///   `trunc_err` is incidentally small, and that case must still be
///   rejected rather than silently accepted.
/// - `trunc_err` — the truncation-error estimate `max|b₆|/max|a₀|`,
///   scaling as `dt⁷`.
/// - `dt_try` — step size of the attempt being judged. Used only for
///   the `DT_MIN` escape check.
/// - `eps` — user's target tolerance; clamped to `[EPSILON_MIN, EPSILON_MAX]`
///   on the way in by [`Ias15::set_epsilon`].
/// - `deadline_hit` — cooperative wall-clock budget has been exceeded;
///   used to short-circuit retry spins in pathological scenes.
///
/// # Decision table
///
/// | `converged` | `trunc ≤ ε` | `dt ≤ DT_MIN` | `deadline` | → |
/// |---|---|---|---|---|
/// | T | T | — | — | `Accept { degraded: false }` |
/// | F | — | T | — | `Accept { degraded: true }` |
/// | T | F | T | — | `Accept { degraded: true }` |
/// | — | — | — | T | `Accept { degraded: true }` |
/// | F | — | F | F | `RejectPicard` |
/// | T | F | F | F | `RejectTruncation` |
fn decide_dt(
    picard_converged: bool,
    trunc_err: f64,
    dt_try: f64,
    eps: f64,
    deadline_hit: bool,
) -> DtDecision {
    let on_merit = picard_converged && trunc_err <= eps;
    if on_merit {
        return DtDecision::Accept { degraded: false };
    }
    if dt_try <= DT_MIN || deadline_hit {
        return DtDecision::Accept { degraded: true };
    }
    if !picard_converged { DtDecision::RejectPicard } else { DtDecision::RejectTruncation }
}

// ── Integrator struct ────────────────────────────────────────────────────────

/// Per-body polynomial state for one substep (coefficients of the
/// series expansion of the acceleration within the step). Index 0..7
/// is the coefficient order; the value is the 3D acceleration
/// coefficient.
type BodyCoeffs = [Vec3; 7];

// Layout guard: the snapshot path uses `copy_from_slice` for tight
// memcpy semantics and would silently copy padding if a future
// refactor broke the tight packing. `Vec3` is `#[repr(C)]` with three
// `f64` fields and no padding, so `BodyCoeffs` packs as 7 × 24 = 168
// bytes. Caught at compile time — zero runtime cost.
const _: () = {
    assert!(
        std::mem::size_of::<BodyCoeffs>() == 7 * 24,
        "BodyCoeffs layout changed — verify snapshot copy_from_slice still hits \
         the intended byte range"
    );
    assert!(std::mem::align_of::<BodyCoeffs>() == 8);
};

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
    csx: Vec<Vec3>,
    /// Compensated-summation carry for velocities.
    csv: Vec<Vec3>,

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

    /// Cumulative rejections caused by Picard predictor–corrector not
    /// converging within [`MAX_PICARD_ITERATIONS`]. High values here
    /// (relative to `rejections_truncation_total`) indicate a stiff
    /// or high-eccentricity regime where the step size exceeds the
    /// Lipschitz bound — the error controller cannot see this signal
    /// through the truncation estimate alone (cf. TD-004).
    rejections_picard_total: u64,

    /// Cumulative rejections where Picard converged but the truncation
    /// bound `|b₆|/|a₀|` exceeded `ε`. This is the "well-posed"
    /// rejection class — the controller formula applies and the next
    /// attempt uses `(ε/err)^(1/7)` scaling.
    rejections_truncation_total: u64,

    /// Cumulative Picard iteration count summed across all attempts
    /// (accepted and rejected). Mean iterations per attempt is
    /// `picard_iters_total / (substeps_total + rejections_picard_total
    /// + rejections_truncation_total)`.
    picard_iters_total: u64,

    /// Cumulative count of degraded accepts (`DT_MIN` escape clause
    /// or deadline fired). Should stay at zero for well-posed scenes.
    degraded_total: u64,

    /// Cumulative count of Picard early-exits via the stagnation guard
    /// (residual stopped decreasing for two consecutive iterations). On
    /// smooth motion this counter stays well below `substeps_total`; a
    /// sustained high ratio is the cheapest signal that the warmstart
    /// is biasing `b` outside Picard's basin of attraction. Surfaced
    /// through [`AdaptiveStats::picard_stagnations`].
    picard_stagnations_total: u64,

    /// Cumulative count of "shrink → grow" reversals in the controller.
    /// Detected when `dt_next` increases relative to `dt_last_accepted`
    /// after the previous accept's `dt_next` had decreased relative to
    /// *its* `dt_last_accepted`. Healthy adaptive runs in a smooth
    /// regime register `shrink_grow_cycles_total / substeps_total ≈ 0`;
    /// chatter indicates controller-warmstart oscillation and is the
    /// fingerprint the figure-8 cascade left behind. Surfaced through
    /// [`AdaptiveStats::shrink_grow_cycles`].
    shrink_grow_cycles_total: u64,

    /// Sign of the last `(dt_next - dt_last_accepted)` direction:
    /// `-1` means the controller just shrunk, `+1` grew, `0` means no
    /// previous accept. Updated on every accept in lockstep with
    /// `dt_next`. Used solely to detect reversals — the direction
    /// itself is not exposed.
    dt_dir_prev: i8,

    /// First-Same-As-Last (FSAL) cache flag. `true` iff the caller's
    /// `acc` buffer holds the gravitational acceleration evaluated at
    /// the *current* body positions, i.e. iff the previous accepted
    /// sub-step ended with a force evaluation that produced exactly the
    /// `a₀` the next sub-step needs.
    ///
    /// The flag is set after every accept-path force evaluation and
    /// cleared whenever the integrator's state or the body positions
    /// change in a way that would invalidate the cached `acc`:
    /// `ensure_capacity` resize (body count changed), `recenter_bodies`
    /// (uniform translation applied externally), and the very first
    /// call (no prior accept has occurred).
    ///
    /// External mutation of body positions or of `acc` between
    /// `step()` calls is not detectable from inside the integrator;
    /// callers that perform such mutation are expected to invalidate
    /// the cache explicitly. In `apsis` the only path that mutates
    /// body positions outside the integrator is the periodic COM
    /// recentering, which is now routed through
    /// [`Integrator::recenter_bodies`] precisely so the integrator
    /// can invalidate this flag in the same call.
    has_valid_post_acc: bool,

    /// `ctx.g_factor` from the most recent accept-path force
    /// evaluation. Compared against the incoming `ctx.g_factor` on the
    /// next call: a mismatch invalidates [`Self::has_valid_post_acc`]
    /// because the cached `acc` has been scaled by the old value.
    cached_g_factor: f64,

    /// Length of `ctx.perturbations` from the most recent accept-path
    /// force evaluation. Compared against the incoming length on the
    /// next call: a mismatch invalidates the FSAL cache because
    /// perturbations have been added or cleared since the cache was
    /// populated. Replacement of a perturbation in-place at unchanged
    /// length is *not* detected here — callers that swap perturbations
    /// without changing the length must invalidate the integrator
    /// state explicitly (typically by recreating the integrator).
    cached_perturbation_count: usize,

    // ── Picard scratch buffers ───────────────────────────────────────
    //
    // Start-of-attempt positions and velocities, and the previous
    // iteration's `b₆` snapshot. These are logically local to
    // [`Self::picard_loop`] — moving them into the struct swaps a
    // per-retry `Vec` allocation for a `clear() + extend()` reuse of
    // the existing heap buffer, which steady-state is zero-alloc.
    //
    // Left in a possibly-stale state between calls: the Picard
    // implementation always re-fills them via `clear() + extend`
    // before reading, so the previous run's contents cannot leak.
    pic_x0: Vec<Vec3>,
    pic_v0: Vec<Vec3>,
    pic_b6_old: Vec<Vec3>,

    // ── Rejection-rollback snapshot buffers ──────────────────────────
    //
    // Before the rollback refactor this was a stack-allocated `Attempt`
    // struct that cloned seven `Vec`s per sub-step. Phase profiling
    // (see docs/experiments/2026-04-22-ias15-phase-profile.md) showed
    // 4–8% of total wall time spent in those clones while >99% of
    // sub-steps never reject, making the allocation/memcpy/drop cycle
    // pure waste.
    //
    // Each buffer here is sized to `n` exactly once per body-count
    // change in [`Self::ensure_capacity`]; `capture_snapshot` and
    // `restore_snapshot` then move data via `copy_from_slice` /
    // explicit fills with no further allocation. The
    // [`Self::snapshot_valid`] flag guards the lifecycle: a
    // `restore_snapshot` call without a prior `capture_snapshot`
    // (e.g. after `ensure_capacity` resets buffers) is a programmer
    // error and panics in debug builds rather than silently restoring
    // stale data.
    //
    // `snap_csx` / `snap_csv` are retained despite the invariant that
    // rejections never touch the live `csx`/`csv` (only `advance_state`
    // on the accept path writes them). Keeping them in the snapshot
    // means future code that violates that invariant — e.g. a hook
    // that runs during a rejected attempt — still gets correct
    // rollback semantics without silently corrupting compensated-
    // summation carries. Cost is 32 bytes of memcpy per sub-step at
    // N=2, well inside the wash of the bigger `b`/`e`/`csb` copies.
    snap_x: Vec<Vec3>,
    snap_v: Vec<Vec3>,
    snap_b: Vec<BodyCoeffs>,
    snap_e: Vec<BodyCoeffs>,
    snap_csb: Vec<BodyCoeffs>,
    snap_csx: Vec<Vec3>,
    snap_csv: Vec<Vec3>,
    snapshot_valid: bool,
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
            rejections_picard_total: 0,
            rejections_truncation_total: 0,
            picard_iters_total: 0,
            degraded_total: 0,
            picard_stagnations_total: 0,
            shrink_grow_cycles_total: 0,
            dt_dir_prev: 0,
            has_valid_post_acc: false,
            cached_g_factor: 1.0,
            cached_perturbation_count: 0,
            pic_x0: Vec::new(),
            pic_v0: Vec::new(),
            pic_b6_old: Vec::new(),
            snap_x: Vec::new(),
            snap_v: Vec::new(),
            snap_b: Vec::new(),
            snap_e: Vec::new(),
            snap_csb: Vec::new(),
            snap_csx: Vec::new(),
            snap_csv: Vec::new(),
            snapshot_valid: false,
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
            self.b = vec![[Vec3::ZERO; 7]; n];
            self.e = vec![[Vec3::ZERO; 7]; n];
            self.g = vec![[Vec3::ZERO; 7]; n];
            self.csb = vec![[Vec3::ZERO; 7]; n];
            self.csx = vec![Vec3::ZERO; n];
            self.csv = vec![Vec3::ZERO; n];
            self.dt_last_accepted = 0.0;
            // Reset dt_next too: a value from a different body count is
            // physically meaningless (different acceleration regime) and
            // would spend rejections in the retry loop before re-calibrating.
            // The `if self.dt_next <= 0.0` guard in `step()` re-seeds from
            // the caller's budget on the next entry.
            self.dt_next = 0.0;

            // Snapshot buffers: size to exactly `n`. From this point on,
            // `capture_snapshot` and `restore_snapshot` use length-stable
            // `copy_from_slice` and element-wise fills without reallocating.
            // Invalidate until the first capture so a premature
            // `restore_snapshot` is caught by the debug assertion.
            self.snap_x = vec![Vec3::ZERO; n];
            self.snap_v = vec![Vec3::ZERO; n];
            self.snap_b = vec![[Vec3::ZERO; 7]; n];
            self.snap_e = vec![[Vec3::ZERO; 7]; n];
            self.snap_csb = vec![[Vec3::ZERO; 7]; n];
            self.snap_csx = vec![Vec3::ZERO; n];
            self.snap_csv = vec![Vec3::ZERO; n];
            self.snapshot_valid = false;

            // FSAL cache is tied to the body count: a resize means the
            // caller's `acc` is either still sized for the old N (will
            // be re-evaluated and resized inside `evaluate`) or a fresh
            // buffer with no relation to the current state. Either way
            // we cannot reuse it.
            self.has_valid_post_acc = false;
        }
    }

    /// Populate the rollback snapshot buffers with the current state.
    /// Called once per sub-step, before any `b`/`e`/`csb` modification.
    /// After this returns, [`Self::restore_snapshot`] is allowed.
    ///
    /// Uses per-slot assignment for body-derived fields (which need a
    /// tuple transformation from the `Body` struct) and `copy_from_slice`
    /// for Ias15-internal fields (which are already `Vec<BodyCoeffs>` or
    /// `Vec<(f64, f64)>`). Both paths are pure memcpy on the happy path
    /// where the destination is already sized to `n` by
    /// [`Self::ensure_capacity`].
    fn capture_snapshot(&mut self, bodies: &[Body]) {
        debug_assert_eq!(self.snap_x.len(), bodies.len(), "snapshot buffer size mismatch");

        for (dst, src) in self.snap_x.iter_mut().zip(bodies.iter()) {
            *dst = Vec3::new(src.pos_x, src.pos_y, src.pos_z);
        }
        for (dst, src) in self.snap_v.iter_mut().zip(bodies.iter()) {
            *dst = Vec3::new(src.vel_x, src.vel_y, src.vel_z);
        }
        self.snap_b.copy_from_slice(&self.b);
        self.snap_e.copy_from_slice(&self.e);
        self.snap_csb.copy_from_slice(&self.csb);
        self.snap_csx.copy_from_slice(&self.csx);
        self.snap_csv.copy_from_slice(&self.csv);

        self.snapshot_valid = true;
    }

    /// Roll `bodies` and the integrator state back to the last
    /// [`Self::capture_snapshot`]. Called on rejection branches.
    ///
    /// Debug-asserts the snapshot has actually been captured —
    /// calling `restore` without a prior `capture` (e.g. immediately
    /// after `ensure_capacity`'s buffer reset) is a programmer error
    /// that previously would have silently restored all-zero state;
    /// this fails loud in debug builds and is a single branch-predicted
    /// check in release.
    fn restore_snapshot(&mut self, bodies: &mut [Body]) {
        debug_assert!(
            self.snapshot_valid,
            "restore_snapshot called without a prior capture_snapshot"
        );
        debug_assert_eq!(self.snap_x.len(), bodies.len(), "snapshot buffer size mismatch");

        for (b, src) in bodies.iter_mut().zip(self.snap_x.iter()) {
            b.pos_x = src.x;
            b.pos_y = src.y;
            b.pos_z = src.z;
        }
        for (b, src) in bodies.iter_mut().zip(self.snap_v.iter()) {
            b.vel_x = src.x;
            b.vel_y = src.y;
            b.vel_z = src.z;
        }
        self.b.copy_from_slice(&self.snap_b);
        self.e.copy_from_slice(&self.snap_e);
        self.csb.copy_from_slice(&self.snap_csb);
        self.csx.copy_from_slice(&self.snap_csx);
        self.csv.copy_from_slice(&self.snap_csv);
    }
}

// ── Core algorithm ───────────────────────────────────────────────────────────

impl Integrator for Ias15 {
    /// Perform **one** adaptive Gauss-Radau sub-step.
    ///
    /// The input `dt_hint` is the controller's first-call seed (see the
    /// module-level documentation on sub-step semantics); the actual
    /// step size, chosen by the error controller per the IAS15
    /// specification (Rein & Spiegel 2015 §2.3), is reported through
    /// [`StepResult::consumed_dt`]. The caller is expected to re-invoke
    /// `step` until the full simulation-time target has been reached.
    fn step(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        dt_hint: f64,
        acc: &mut Vec<Vec3>,
    ) -> StepResult {
        let n = bodies.len();
        self.ensure_capacity(n);

        // Defensive contract: every input kinematic component must be
        // finite. NaN/inf in `body.{x, y, z, vx, vy, vz}` would propagate
        // through the Picard substep predictor and the Gauss–Radau b/e/g
        // buffers without producing a usable signal at any controller
        // level. Cheap to assert in debug, free in release.
        debug_assert!(
            bodies.iter().all(|b| b.pos_x.is_finite()
                && b.pos_y.is_finite()
                && b.pos_z.is_finite()
                && b.vel_x.is_finite()
                && b.vel_y.is_finite()
                && b.vel_z.is_finite()),
            "IAS15: non-finite input state — NaN/inf in body kinematics"
        );

        // Caller's `dt_hint` is the *first-call seed* for the controller,
        // not a per-call hard cap. Capping `dt_try` at `dt_hint` on every
        // call would pin the integrator to the user's initial guess and
        // prevent it from ever growing to the natural step that the local
        // truncation error permits — the IAS15 controller would behave
        // like a fixed-step integrator with adaptive sub-stepping below
        // the cap. The canonical IAS15 contract (Rein & Spiegel 2015
        // §2.3) treats the user's `dt` as a hint that the controller
        // mutates and reuses on subsequent calls; we follow the same
        // contract by storing the first-call seed in `dt_next` and
        // letting the controller drive `dt_try` from there onward.
        //
        // A zero or negative hint is treated as DT_MIN (pathological, but
        // legal — e.g. when the caller has no opinion on the step size).
        if self.dt_next <= 0.0 {
            self.dt_next = dt_hint.max(DT_MIN);
        }

        let mut dt_try = self.dt_next.max(DT_MIN);

        // Snapshot taken once per sub-step: body kinematics + integrator
        // state that must survive rejection retries. `a0` (start-of-step
        // acceleration) is also invariant across retries because the
        // snapshot restores positions, so we evaluate it only once —
        // saving up to `max_picard_iter × n_rejects` force calls per
        // sub-step compared to re-evaluating inside the retry loop.
        //
        // The snapshot lives in persistent `snap_*` fields (sized in
        // `ensure_capacity`), so this call path is pure memcpy — no
        // per-sub-step allocation. See the field comments on `snap_*`
        // for the lifecycle contract that `snapshot_valid` enforces.
        time_phase!(snapshot_capture, {
            self.capture_snapshot(bodies);
        });

        // First-Same-As-Last: the previous accept's end-of-step force
        // evaluation produced `acc` at exactly the body positions this
        // sub-step starts from. Reuse it as `a₀` instead of paying a
        // second force evaluation at the same point.
        //
        // Validity of the cache requires four invariants to hold since
        // the post-accept evaluation that populated it:
        //
        //   1. body positions unchanged outside the integrator
        //      (invalidated by `ensure_capacity` resize and by
        //      `recenter_bodies` translation),
        //   2. `acc.len() == n` (the caller has not resized the buffer),
        //   3. `ctx.g_factor` matches the cached value (the cached `acc`
        //      has been scaled by `cached_g_factor`; a new value would
        //      need a re-scale that is cheaper to do via re-evaluation
        //      than to track),
        //   4. `ctx.perturbations.len()` matches the cached count (added
        //      or cleared perturbations would change the cached
        //      contributions). In-place perturbation replacement at
        //      unchanged length is documented as caller responsibility.
        //
        // When all four hold, the fast path simply clones the existing
        // `acc` into `a₀`, saving one force evaluation per accepted
        // sub-step — the canonical FSAL property of any explicit /
        // implicit method whose stage-0 evaluation coincides with the
        // previous step's stage-end.
        let pert_count =
            ctx.hamiltonian_perturbations.len() + ctx.non_conservative_perturbations.len();
        let fsal_valid = self.has_valid_post_acc
            && acc.len() == n
            && self.cached_g_factor == ctx.g_factor
            && self.cached_perturbation_count == pert_count;
        let a0: Vec<Vec3> = if fsal_valid {
            time_phase!(a0_clone, { acc.clone() })
        } else {
            let raw_pe = time_phase!(evaluate, { evaluate(bodies, ctx.force, acc) });
            scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
            accumulate_perturbation_forces(
                bodies,
                acc,
                ctx.hamiltonian_perturbations,
                ctx.non_conservative_perturbations,
            );
            time_phase!(a0_clone, { acc.clone() })
        };

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
                time_phase!(warmstart_b, {
                    self.warmstart_b(dt_try, self.dt_last_accepted);
                });
                diag_emit_warmstart(self, dt_try / self.dt_last_accepted, dt_try);
            }
            time_phase!(recompute_g_from_b, {
                self.recompute_g_from_b();
            });

            let (converged, _picard_err, picard_iters) =
                self.picard_loop(bodies, ctx, acc, &a0, dt_try);
            self.picard_iters_total = self.picard_iters_total.saturating_add(picard_iters as u64);

            let trunc_err = self.truncation_error(&a0);

            // Delegate acceptance/rejection to the pure [`decide_dt`]
            // function. It treats Picard convergence and truncation error
            // as independent signals (TD-004) so the shrink strategy can
            // be picked correctly for each failure class; prior revisions
            // collapsed them with `max(…)` and underfed the controller.
            let deadline_hit =
                ctx.deadline.map(|d| std::time::Instant::now() >= d).unwrap_or(false);

            match decide_dt(converged, trunc_err, dt_try, self.epsilon, deadline_hit) {
                DtDecision::RejectPicard => {
                    self.rejections_picard_total = self.rejections_picard_total.saturating_add(1);
                    time_phase!(snapshot_restore, {
                        self.restore_snapshot(bodies);
                    });
                    // Fixed halving per Rein & Spiegel (2015) §3.4:
                    // Picard divergence means the step exceeds the
                    // local Lipschitz bound, and the (ε/err)^{1/7}
                    // formula — which assumes the dt⁷-scaled
                    // truncation regime — would under-shrink here.
                    let dt_next_attempt = (dt_try * PICARD_SHRINK).max(DT_MIN);
                    diag_emit_attempt(
                        self,
                        dt_try,
                        dt_next_attempt,
                        trunc_err,
                        converged,
                        picard_iters,
                        "reject_picard",
                    );
                    dt_try = dt_next_attempt;
                    continue;
                },
                DtDecision::RejectTruncation => {
                    self.rejections_truncation_total =
                        self.rejections_truncation_total.saturating_add(1);
                    time_phase!(snapshot_restore, {
                        self.restore_snapshot(bodies);
                    });
                    // Halving on truncation rejection per the IAS15
                    // specification (Rein & Spiegel 2015 §3.4).
                    //
                    // Earlier revisions of this file used the
                    // `optimal_dt` formula (`dt · 0.9 · (ε/err)^{1/7}`)
                    // here, mirroring the controller's accept-path
                    // proposal. The trade-off showed up cleanly on the
                    // figure-8 cascade trace: when `err` sits just
                    // above `ε` (the noise-floor regime that emerges
                    // past close-encounter onset), the formula shrinks
                    // by only 10–20 % per attempt; 25–30 retries are
                    // needed to drop `dt` by the factor-of-1000+ that
                    // the local geometry actually demands, and each
                    // retry pays a full force-eval sweep. The
                    // specification's halving converges to an
                    // acceptable `dt` in ~10 retries; we follow it.
                    let dt_next_attempt = (dt_try * 0.5).max(DT_MIN);
                    diag_emit_attempt(
                        self,
                        dt_try,
                        dt_next_attempt,
                        trunc_err,
                        converged,
                        picard_iters,
                        "reject_trunc",
                    );
                    dt_try = dt_next_attempt;
                    continue;
                },
                DtDecision::Accept { degraded: step_degraded } => {
                    self.substeps_total = self.substeps_total.saturating_add(1);
                    if step_degraded {
                        self.degraded_total = self.degraded_total.saturating_add(1);
                        // Distinguish the two causes reported by `decide_dt`:
                        // `dt_try <= DT_MIN` means the adaptive controller
                        // wanted to shrink further but saturated the floor,
                        // which is a **scenario stiffness signal** — the
                        // close-encounter geometry is beyond what IAS15 can
                        // resolve at f64 precision. The deadline branch
                        // (cooperative budget exhausted) is expected in
                        // interactive precision runs and is not a scenario
                        // indictment — silenced here; the cumulative counter
                        // in `AdaptiveStats` still tracks it.
                        //
                        // Log rate: first three occurrences verbatim, then
                        // every power of two (4, 8, 16, 32, ...). Exponentially
                        // thins the emission rate while keeping a running
                        // `floor_hit_count` in every event. Avoids drowning
                        // stderr when a pathological scene hits the floor
                        // thousands of times, without losing the initial
                        // signal.
                        if dt_try <= DT_MIN {
                            let c = self.degraded_total;
                            if c <= 3 || c.is_power_of_two() {
                                crate::warn_diag!(
                                    crate::core::log::Source::Integrator,
                                    "IAS15 dt floor reached; controller accepted degraded step",
                                    dt = dt_try,
                                    floor = DT_MIN,
                                    floor_hit_count = c,
                                    substep = self.substeps_total,
                                    hint = "scenario may be stiff — consider increasing softening, reducing N, or relaxing epsilon",
                                );
                            }
                        }
                    }

                    // Accept path ──────────────────────────────────────
                    // Build the dense-output snapshot *before* we advance
                    // the state, so it carries the pre-step kinematics
                    // (the b-coeffs below are the accepted values —
                    // `self.b` is not further modified on the accept
                    // path). The caller (`System::step`) fills in the
                    // absolute `t0` as `system.t() - consumed_dt`.
                    let step_snapshot = time_phase!(dense_snapshot_build, {
                        DenseSnapshot {
                            t0: 0.0,
                            dt: dt_try,
                            x0: self.snap_x.clone(),
                            v0: self.snap_v.clone(),
                            a0: a0.clone(),
                            b: self.b.clone(),
                            kind: IntegratorKind::Ias15,
                            wh_data: None,
                        }
                    });

                    time_phase!(advance_state, {
                        self.advance_state(bodies, &a0, dt_try);
                    });

                    // Post-step force evaluation: publishes `acc`
                    // consistent with the final body positions, and
                    // returns the potential energy the caller will use
                    // for energy-conservation diagnostics. The same
                    // `acc` is also the next sub-step's start-of-step
                    // `a₀`; we record the parameters used here so the
                    // FSAL fast path on the next call can verify the
                    // cache is still valid (see the FSAL guard at the
                    // top of `step()`).
                    let raw_pe = time_phase!(evaluate, { evaluate(bodies, ctx.force, acc) });
                    let pe = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
                    accumulate_perturbation_forces(
                        bodies,
                        acc,
                        ctx.hamiltonian_perturbations,
                        ctx.non_conservative_perturbations,
                    );

                    self.has_valid_post_acc = true;
                    self.cached_g_factor = ctx.g_factor;
                    self.cached_perturbation_count = ctx.hamiltonian_perturbations.len()
                        + ctx.non_conservative_perturbations.len();

                    self.update_warmstart_record();
                    self.dt_last_accepted = dt_try;
                    // Propose next dt from the clean truncation signal
                    // only. Prior revisions used `max(trunc, picard)`
                    // here, but Picard convergence is binary and
                    // `picard_err` is a residual from a root-finding
                    // iteration — not a truncation-error surrogate.
                    // Feeding it into `optimal_dt` distorted the
                    // controller (cf. TD-004).
                    //
                    // The proposal is then capped above at
                    // `dt_try · DT_GROWTH_LIMIT` per the IAS15
                    // specification's 7×-per-step growth bound (Rein
                    // & Spiegel 2015 §3.4); without the cap,
                    // smooth-region steps with `err ≪ ε` propose
                    // wild overshoots that the next attempt's
                    // truncation gate immediately rejects, costing a
                    // full shrink cascade per close encounter.
                    let raw = self.optimal_dt(dt_try, trunc_err);
                    let new_dt_next = raw.min(dt_try * DT_GROWTH_LIMIT).max(DT_MIN);

                    // Shrink-grow chatter detection: a *reversal* fires
                    // when the current step proposed growth (`dt_next >
                    // dt_try`) immediately after a step that proposed a
                    // shrink (`dt_dir_prev == -1`). On smooth motion the
                    // controller settles on a near-constant `dt_next`
                    // and reversals are rare; persistent chatter
                    // signals warmstart-controller oscillation, which
                    // is the cumulative-failure fingerprint surfaced
                    // through `AdaptiveStats::shrink_grow_cycles`.
                    let dt_dir_now: i8 = if new_dt_next > dt_try {
                        1
                    } else if new_dt_next < dt_try {
                        -1
                    } else {
                        0
                    };
                    if self.dt_dir_prev == -1 && dt_dir_now == 1 {
                        self.shrink_grow_cycles_total =
                            self.shrink_grow_cycles_total.saturating_add(1);
                    }
                    self.dt_dir_prev = dt_dir_now;
                    self.dt_next = new_dt_next;

                    let label = if step_degraded {
                        if deadline_hit { "accept_deadline" } else { "accept_floor" }
                    } else {
                        "accept"
                    };
                    diag_emit_attempt(
                        self,
                        dt_try,
                        new_dt_next,
                        trunc_err,
                        converged,
                        picard_iters,
                        label,
                    );

                    break (dt_try, pe, step_snapshot, step_degraded);
                },
            }
        };

        StepResult {
            consumed_dt: accepted_dt,
            potential_energy: final_pe,
            used_fallback: false,
            step_snapshot: Some(final_snapshot),
            degraded,
            hierarchy_signal: None,
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

    /// IAS15 owns its step size: the caller's `dt_hint` is a first-call
    /// seed, then the Gauss-Radau error controller drives `dt_next`.
    fn controls_own_step_size(&self) -> bool {
        true
    }

    /// The controller's recommended size for the next sub-step. Equals
    /// the caller's hint before the very first call (when `dt_next` is
    /// still 0); otherwise the value computed by `optimal_dt` after the
    /// most-recent accept.
    fn proposed_next_dt(&self) -> Option<f64> {
        if self.dt_next > 0.0 { Some(self.dt_next) } else { None }
    }

    /// Apply a uniform translation `(-dx, -dy)` while keeping IAS15's
    /// per-body Neumaier compensation buffers (`csx`) consistent with
    /// the post-shift body positions. We route the translation through
    /// the same `add_cs` primitive that `advance_state` uses for the
    /// per-step position update, so the compensated representation
    /// `(body.x_new, csx_new)` is the extended-precision continuation
    /// of `(body.x_old, csx_old) − (dx, dy)`.
    ///
    /// A bare subtraction (`body.pos_x -= dx`) without touching `csx`
    /// would silently drop the compensation history accumulated up to
    /// the recentering call: the next `add_cs` invocation
    /// (`y = inp − csx_old`) consumes the stale `csx` entry and the
    /// rounding errors it tracked are gone. For one sub-ULP shift
    /// this loss is below the trajectory-level signal, but on long
    /// integrations under periodic COM recentering it adds up and
    /// breaks the bit-reproducibility property the IAS15 module
    /// header advertises (round-off `O(ε)` rather than `O(ε · N_steps)`).
    fn recenter_bodies(&mut self, bodies: &mut [Body], dx: f64, dy: f64) {
        // The compensation buffers `csx` are sized lazily by
        // `ensure_capacity` on the first `step()` call. If the caller
        // recenters before any step has run (e.g. an explicit
        // `System::recenter_com` invocation right after construction),
        // there is nothing to preserve — fall back to the bare
        // subtraction that the trait default would have done.
        if self.csx.len() != bodies.len() {
            for b in bodies.iter_mut() {
                b.pos_x -= dx;
                b.pos_y -= dy;
            }
            // FSAL cache invalidation: `acc` was at the pre-shift
            // body positions, which are no longer where `bodies` are.
            self.has_valid_post_acc = false;
            return;
        }
        for (i, b) in bodies.iter_mut().enumerate() {
            add_cs(&mut b.pos_x, &mut self.csx[i].x, -dx);
            add_cs(&mut b.pos_y, &mut self.csx[i].y, -dy);
        }
        // FSAL cache invalidation (same reason as the early-return
        // branch above): even with the compensation-aware shift, body
        // positions have moved by `(-dx, -dy)` from where the cached
        // `acc` was evaluated.
        self.has_valid_post_acc = false;
    }

    fn adaptive_stats(&self) -> Option<super::traits::AdaptiveStats> {
        let rejections_picard = self.rejections_picard_total;
        let rejections_truncation = self.rejections_truncation_total;
        Some(super::traits::AdaptiveStats {
            substeps: self.substeps_total,
            rejections: rejections_picard.saturating_add(rejections_truncation),
            rejections_picard,
            rejections_truncation,
            picard_iters: self.picard_iters_total,
            degraded: self.degraded_total,
            picard_stagnations: self.picard_stagnations_total,
            shrink_grow_cycles: self.shrink_grow_cycles_total,
        })
    }

    fn execution_profile(&self) -> super::traits::ExecutionProfile {
        // Adaptive Gauss-Radau with unbounded shrinking toward DT_MIN in
        // stiff regimes; per-step wall time is not bounded by N alone.
        super::traits::ExecutionProfile::Precision
    }

    fn requires_deterministic_force(&self) -> bool {
        // Picard predictor-corrector reaches its fixed point only when
        // f(x, v, t) is bit-reproducible across iterations. BH tree
        // rebuilds are position-dependent and break this; direct O(N²)
        // satisfies it.
        true
    }
}

impl Ias15 {
    /// Drop the FSAL cache. The next [`Self::step`] call will re-evaluate
    /// the start-of-sub-step acceleration from scratch instead of cloning
    /// the previous accept's end-of-step `acc`.
    ///
    /// Required when an external operator (e.g. a hybrid integrator
    /// driving IAS15 as a sub-integrator after applying its own drifts /
    /// kicks) has mutated body positions or velocities since the last
    /// IAS15 call. Without this, the cached `acc` references the wrong
    /// configuration and contaminates the next sub-step's `a₀`.
    pub fn invalidate_force_cache(&mut self) {
        self.has_valid_post_acc = false;
    }

    /// Cap the controller's proposed next sub-step at `cap`. No-op when
    /// the current proposal is already at or below `cap`; otherwise
    /// clips `dt_next` so the next [`Self::step`] call cannot consume
    /// more than `cap` time units.
    ///
    /// Used by hybrid integrators (Mercurius) that drive IAS15 over a
    /// fixed outer window: the controller's natural growth between
    /// sub-steps would otherwise overshoot the window boundary on the
    /// last sub-step and break the outer Hamiltonian split.
    pub fn cap_proposed_dt(&mut self, cap: f64) {
        if cap > 0.0 && self.dt_next > cap {
            self.dt_next = cap;
        }
    }

    /// Inner predictor-corrector iteration. Given `a0` (acceleration at
    /// the start of the attempt) and a target `dt_try`, iteratively
    /// refines `b` until max|Δb₆|/max|a₀| < `PICARD_TOL` or we hit the
    /// iteration cap. Returns `(converged, residual, iters)` — `iters`
    /// counts the actual iterations consumed (1..=MAX_PICARD_ITERATIONS)
    /// so the outer controller can aggregate them into diagnostics.
    ///
    /// Thin wrapper over [`Self::picard_loop_inner`]: moves the
    /// persistent scratch buffers out of `self` for the duration of
    /// the call so the inner function can hold `&mut` on them
    /// simultaneously with `&mut self` (needed to call
    /// [`Self::update_g_and_b`]). The Vec instances are returned to
    /// their fields on every exit path; `mem::take` leaves the fields
    /// as empty `Vec`s during the call, so a panic mid-iteration
    /// would leave the integrator state internally consistent
    /// (capacity lost, but length correct) — not that IAS15 is
    /// expected to survive panics in any meaningful way.
    fn picard_loop(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        acc: &mut Vec<Vec3>,
        a0: &[Vec3],
        dt_try: f64,
    ) -> (bool, f64, u32) {
        let mut x0 = std::mem::take(&mut self.pic_x0);
        let mut v0 = std::mem::take(&mut self.pic_v0);
        let mut b6_old = std::mem::take(&mut self.pic_b6_old);

        let result =
            self.picard_loop_inner(bodies, ctx, acc, a0, dt_try, &mut x0, &mut v0, &mut b6_old);

        self.pic_x0 = x0;
        self.pic_v0 = v0;
        self.pic_b6_old = b6_old;
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn picard_loop_inner(
        &mut self,
        bodies: &mut [Body],
        ctx: &mut IntegratorContext<'_>,
        acc: &mut Vec<Vec3>,
        a0: &[Vec3],
        dt_try: f64,
        x0: &mut Vec<Vec3>,
        v0: &mut Vec<Vec3>,
        b6_old: &mut Vec<Vec3>,
    ) -> (bool, f64, u32) {
        let n = bodies.len();

        // Populate the scratch buffers for this call. `clear() + extend`
        // reuses the existing allocation when `len` fits into `capacity`;
        // when `n` changed since the last call the first push realloc's
        // once and the rest reuse. In steady-state (constant body count)
        // these loops are zero-alloc.
        x0.clear();
        x0.extend(bodies.iter().map(|b| Vec3::new(b.pos_x, b.pos_y, b.pos_z)));
        v0.clear();
        v0.extend(bodies.iter().map(|b| Vec3::new(b.vel_x, b.vel_y, b.vel_z)));

        let mut last_residual = f64::INFINITY;
        let mut no_improve: u32 = 0;
        let mut iters: u32 = 0;

        for iter in 0..MAX_PICARD_ITERATIONS {
            iters = (iter as u32) + 1;
            // Snapshot b₆ before the iteration — residual is measured
            // against this.
            b6_old.clear();
            b6_old.extend(self.b.iter().map(|row| row[6]));

            for stage in 1..=7 {
                let s = H[stage];
                // Predict position AND velocity at node `s` per R&S 2015
                // §2 / eq. 11. Velocity prediction is load-bearing for
                // velocity-dependent perturbations (1PN, drag, radiation,
                // Poynting–Robertson): omitting it leaves `body.(vx, vy)`
                // stale across all seven Gauss–Radau nodes and biases
                // each evaluation by `O(a · dt)`, accumulating to ~10⁻³
                // relative precession error on Mercury 1PN at 500 orbits.
                for i in 0..n {
                    let p = predict_ias15(x0[i], v0[i], a0[i], &self.b[i], s, dt_try);
                    let v = predict_v_ias15(v0[i], a0[i], &self.b[i], s, dt_try);
                    bodies[i].pos_x = p.x;
                    bodies[i].pos_y = p.y;
                    bodies[i].pos_z = p.z;
                    bodies[i].vel_x = v.x;
                    bodies[i].vel_y = v.y;
                    bodies[i].vel_z = v.z;
                }

                // Evaluate acceleration at predicted (x, v).
                let raw_pe = time_phase!(evaluate, { evaluate(bodies, ctx.force, acc) });
                let _ = scale_acc_and_pe(acc, ctx.g_factor, raw_pe);
                accumulate_perturbation_forces(
                    bodies,
                    acc,
                    ctx.hamiltonian_perturbations,
                    ctx.non_conservative_perturbations,
                );

                // Update divided-difference g and then b via c-coeffs.
                time_phase!(update_g_and_b, {
                    self.update_g_and_b(stage, a0, acc);
                });
            }

            // Residual: RMS over per-body relative convergence ratios
            // `||Δb₆[i]|| / ||a₀[i]||`, rather than the ratio of the
            // two maxes across bodies.
            //
            // ## Why per-body, why RMS
            //
            // The previous formulation `max||Δb₆|| / max||a₀||` mixed
            // two maxes over potentially different bodies. At N≈2
            // this is harmless (a single body dominates both). At
            // N≈641 the numerator picks up one body's convergence
            // noise outlier while the denominator picks up a
            // completely different body's acceleration magnitude —
            // the ratio is *not* a convergence criterion any more,
            // it is a noise-to-signal floor that grows with body
            // count. Solar_system-class scenarios then cascade into
            // truncation rejections (measured: 3878 rejections over
            // 2001 accepted sub-steps at N=641 before this change)
            // because the reported residual is artificially high.
            //
            // Per-body ratio keeps each body's convergence
            // self-referential (numerator and denominator are the
            // same body). RMS aggregation is O(N) like max but
            // scales gracefully: a single noisy body contributes
            // 1/√N to the total, so the floor shrinks as N grows
            // rather than staying pinned to the worst outlier.
            // The same per-degree-of-freedom convergence criterion is
            // used in spirit by Rein & Spiegel (2015) §3.3 — they
            // formulate it as a max-max ratio rather than RMS, and
            // the two metrics are numerically close at small N. We
            // pick RMS to stay well-behaved at large N (see the lab
            // notebook on solar-system rejection rates under max-max).
            //
            // Bodies with `||a₀[i]|| == 0` are degenerate and
            // skipped (they do not constrain Picard convergence).
            // If every body has zero acceleration the system is
            // gravity-free and any `b` satisfies the ansatz — we
            // return zero residual rather than an undefined
            // quantity.
            let residual = time_phase!(residual_compute, {
                let mut sum_sq = 0.0_f64;
                let mut count: usize = 0;
                for i in 0..n {
                    // `(x² + y²) + z²` left to right. Re-associating
                    // shifts ULPs and is observable on the Picard
                    // residual; keep the explicit form. See
                    // docs/experiments/2026-04-29-3d-port-baseline.md.
                    let db6x = self.b[i][6].x - b6_old[i].x;
                    let db6y = self.b[i][6].y - b6_old[i].y;
                    let db6z = self.b[i][6].z - b6_old[i].z;
                    let db6 = (db6x * db6x + db6y * db6y + db6z * db6z).sqrt();
                    let a_mag = (a0[i].x * a0[i].x + a0[i].y * a0[i].y + a0[i].z * a0[i].z).sqrt();
                    if a_mag > 0.0 {
                        let rel = db6 / a_mag;
                        sum_sq += rel * rel;
                        count += 1;
                    }
                }
                if count > 0 { (sum_sq / count as f64).sqrt() } else { 0.0 }
            });

            if residual < PICARD_TOL {
                // Restore positions/velocities to start-of-attempt so
                // the caller can advance cleanly from (x0, v0).
                restore_xv(bodies, x0, v0);
                return (true, residual, iters);
            }

            // Stagnation = convergence-by-saturation. Near the f64
            // round-off floor the Picard residual stops decreasing and
            // begins oscillating on ULP noise; in that regime the b
            // coefficients are as accurate as f64 allows, and the
            // remaining accept/reject decision belongs to the
            // truncation-error gate, not to a spurious "Picard didn't
            // converge to PICARD_TOL" reject + dt-halving cascade.
            //
            // The previous formulation returned `(false, …)` on two
            // consecutive non-improving iterations, which forced
            // RejectPicard → dt /= 2. On scenarios with sharp local
            // acceleration (figure-8 central crossings, hard close
            // approaches) the residual plateaus above PICARD_TOL but
            // well below any physically meaningful threshold; halving
            // dt cascades all the way to DT_MIN and the controller
            // gets stuck (see `validation/rebound-parity/figure8/`).
            //
            // The IAS15 specification (Rein & Spiegel 2015 §3.3)
            // treats stagnation as `break` from the predictor–
            // corrector loop: the current b is accepted as best
            // effort and the truncation-error gate decides whether
            // dt was OK. We follow the specification here — return
            // `(true, …)` so `decide_dt` routes through the
            // truncation branch, shrinking dt only when the
            // *physical* error estimate demands it.
            if iter >= 2 && residual > last_residual {
                no_improve += 1;
                if no_improve >= 2 {
                    restore_xv(bodies, x0, v0);
                    self.picard_stagnations_total = self.picard_stagnations_total.saturating_add(1);
                    return (true, residual, iters);
                }
            } else {
                no_improve = 0;
            }
            last_residual = residual;
        }

        // Hit MAX_PICARD_ITERATIONS without convergence and without
        // stagnation-triggered early exit. This is genuine
        // non-convergence (residual still strictly decreasing but
        // hasn't crossed PICARD_TOL within the iteration budget) and
        // is the only path that should drive a Picard reject.
        restore_xv(bodies, x0, v0);
        (false, last_residual, iters)
    }

    /// Advance positions and velocities to the end of the accepted
    /// attempt using compensated summation (Neumaier-style) so the
    /// integrator's round-off error stays O(ε) rather than O(ε·N_steps).
    fn advance_state(&mut self, bodies: &mut [Body], a0: &[Vec3], dt: f64) {
        let n = bodies.len();
        for i in 0..n {
            let bi = &self.b[i];

            // Full-step position increment (s = 1):
            //   Δx/dt² = a₀/2 + b₀/6 + b₁/12 + b₂/20 + b₃/30 + b₄/42 + b₅/56 + b₆/72
            //
            // Summation order: smallest-magnitude term first (b₆/72) up
            // to largest (a₀/2). Natural/left-to-right order accumulates
            // into a growing partial sum, which loses low bits of each
            // subsequent smaller term. Reverse order preserves 1–2 extra
            // bits of precision per step — free, and material over the
            // 10⁹-orbit round-off budget the module advertises.
            let dx = dt
                * dt
                * (bi[6].x / 72.0
                    + bi[5].x / 56.0
                    + bi[4].x / 42.0
                    + bi[3].x / 30.0
                    + bi[2].x / 20.0
                    + bi[1].x / 12.0
                    + bi[0].x / 6.0
                    + a0[i].x * 0.5);
            let dy = dt
                * dt
                * (bi[6].y / 72.0
                    + bi[5].y / 56.0
                    + bi[4].y / 42.0
                    + bi[3].y / 30.0
                    + bi[2].y / 20.0
                    + bi[1].y / 12.0
                    + bi[0].y / 6.0
                    + a0[i].y * 0.5);
            let dz = dt
                * dt
                * (bi[6].z / 72.0
                    + bi[5].z / 56.0
                    + bi[4].z / 42.0
                    + bi[3].z / 30.0
                    + bi[2].z / 20.0
                    + bi[1].z / 12.0
                    + bi[0].z / 6.0
                    + a0[i].z * 0.5);

            // Full-step velocity increment (same ascending-magnitude rule):
            //   Δv/dt = a₀ + b₀/2 + b₁/3 + b₂/4 + b₃/5 + b₄/6 + b₅/7 + b₆/8
            let dvx = dt
                * (bi[6].x / 8.0
                    + bi[5].x / 7.0
                    + bi[4].x / 6.0
                    + bi[3].x / 5.0
                    + bi[2].x / 4.0
                    + bi[1].x / 3.0
                    + bi[0].x / 2.0
                    + a0[i].x);
            let dvy = dt
                * (bi[6].y / 8.0
                    + bi[5].y / 7.0
                    + bi[4].y / 6.0
                    + bi[3].y / 5.0
                    + bi[2].y / 4.0
                    + bi[1].y / 3.0
                    + bi[0].y / 2.0
                    + a0[i].y);
            let dvz = dt
                * (bi[6].z / 8.0
                    + bi[5].z / 7.0
                    + bi[4].z / 6.0
                    + bi[3].z / 5.0
                    + bi[2].z / 4.0
                    + bi[1].z / 3.0
                    + bi[0].z / 2.0
                    + a0[i].z);

            // First integrate the v·dt contribution to position.
            let vdt_x = bodies[i].vel_x * dt;
            let vdt_y = bodies[i].vel_y * dt;
            let vdt_z = bodies[i].vel_z * dt;

            add_cs(&mut bodies[i].pos_x, &mut self.csx[i].x, vdt_x);
            add_cs(&mut bodies[i].pos_y, &mut self.csx[i].y, vdt_y);
            add_cs(&mut bodies[i].pos_z, &mut self.csx[i].z, vdt_z);
            add_cs(&mut bodies[i].pos_x, &mut self.csx[i].x, dx);
            add_cs(&mut bodies[i].pos_y, &mut self.csx[i].y, dy);
            add_cs(&mut bodies[i].pos_z, &mut self.csx[i].z, dz);

            add_cs(&mut bodies[i].vel_x, &mut self.csv[i].x, dvx);
            add_cs(&mut bodies[i].vel_y, &mut self.csv[i].y, dvy);
            add_cs(&mut bodies[i].vel_z, &mut self.csv[i].z, dvz);
        }
    }

    /// Estimate of the dominant truncation error term, normalised by
    /// the acceleration magnitude: per-body ‖b₆[i]‖ / ‖a₀[i]‖
    /// aggregated as an RMS across bodies. For a 15th-order method
    /// this is the correct leading term since b₆ multiplies u⁷ ≈ 1
    /// at the end of the step.
    ///
    /// ## Why per-body RMS and not max-max
    ///
    /// This must use the *same* norm as the Picard convergence check
    /// in `picard_loop_inner`. The two measurements interact through
    /// the outer controller: convergence decides whether b was
    /// computed correctly, truncation decides whether b's magnitude
    /// is acceptable physics. If they use different norms, a step
    /// can "converge" under one definition while "failing truncation"
    /// under the other — producing a cascade of rejections with no
    /// physical cause.
    ///
    /// Empirically this showed up at solar_system-class N (≈641
    /// bodies) as a 194% rejection rate (3878 rejections over 2001
    /// accepted sub-steps), nearly all via `RejectTruncation`. The
    /// original max-max formula treated outliers as if they were
    /// representative of the whole system: `max||b₆||` picked up the
    /// one body whose `b₆` had the largest round-off noise,
    /// `max||a₀||` picked up the Sun's dominant acceleration,
    /// producing a ratio that was a noise-to-signal measurement
    /// rather than a truncation estimate. See the diagnostic write-
    /// up referenced in `picard_loop_inner`.
    fn truncation_error(&self, a0: &[Vec3]) -> f64 {
        let mut sum_sq = 0.0_f64;
        let mut count: usize = 0;
        for (i, row) in self.b.iter().enumerate() {
            let b = row[6];
            // Same `(x² + y²) + z²` reduction as the Picard residual.
            let b6 = (b.x * b.x + b.y * b.y + b.z * b.z).sqrt();
            let a_mag = (a0[i].x * a0[i].x + a0[i].y * a0[i].y + a0[i].z * a0[i].z).sqrt();
            if a_mag > 0.0 {
                let rel = b6 / a_mag;
                sum_sq += rel * rel;
                count += 1;
            }
        }
        if count > 0 { (sum_sq / count as f64).sqrt() } else { 0.0 }
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
    /// `dt_try`. Implements the polynomial-basis transformation
    /// derived in Everhart (1985, eq. III.12) and used in IAS15 (Rein
    /// & Spiegel 2015 §3.2): `b_new` is the Pascal-triangle (binomial)
    /// rescaling of `b` by powers of `(dt_try / dt_prev)`, plus a
    /// correction from the drift `b - e` that carries forward last
    /// step's predictor–corrector residual. This drastically reduces
    /// the number of Picard iterations in steady-state integration.
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

        // Pascal-triangle (binomial) coefficients for the polynomial-basis
        // transformation that maps the previous step's `b` to the next
        // step's `e` after the time variable rescales by `q = dt_try / dt_prev`.
        //
        // The acceleration ansatz is `a(u) = a_0 + b_0 * u + b_1 * u^2 + ...
        // + b_6 * u^7` with `u ∈ [0, 1]` the dimensionless intra-step time.
        // When `dt` changes, the *substantive* polynomial in physical time
        // is unchanged; what changes is the variable substitution
        // `u_new = (dt_prev / dt_try) * u_old = u_old / q`, which when
        // re-expanded in the new `u_new` produces a triangular system:
        //
        //     b_new[k] = q^{k+1} * Σ_{j ≥ k} C(j+1, k+1) * b_old[j]
        //
        // where `C(n, m)` is the binomial coefficient. The diagonal term
        // `j = k` gives the simple `q^{k+1} * b[k]` rescaling; the
        // off-diagonal `j > k` cross-terms are what couple higher-order
        // coefficients into the lower-order prediction. The previous
        // implementation kept only the diagonal term, which is correct
        // for the lone `b_6 → b_6` mapping (single column) but loses
        // every cross-term for `k < 6`. The visible symptom: on stiff
        // scenarios with frequent close-encounter dt changes (figure-8,
        // Mercury 1PN), the warm-started `b` accumulates a per-step
        // bias against the true coefficients, which Picard partially
        // corrects but `truncation_error` reports as elevated; the
        // controller responds by collapsing `dt` to the floor and
        // recovery becomes asymptotically slow as the bias compounds
        // at each rescale. The full Pascal expansion is the unique
        // basis transformation that exactly preserves the polynomial
        // under variable-rescaling, derived in Everhart (1985, §III).
        //
        //          b[0]  b[1]  b[2]  b[3]  b[4]  b[5]  b[6]
        //   e[0]:  1     2     3     4     5     6     7      × q
        //   e[1]:  -     1     3     6     10    15    21     × q^2
        //   e[2]:  -     -     1     4     10    20    35     × q^3
        //   e[3]:  -     -     -     1     5     15    35     × q^4
        //   e[4]:  -     -     -     -     1     6     21     × q^5
        //   e[5]:  -     -     -     -     -     1     7      × q^6
        //   e[6]:  -     -     -     -     -     -     1      × q^7

        for i in 0..self.b.len() {
            let be = [
                Vec3::new(
                    self.b[i][0].x - self.e[i][0].x,
                    self.b[i][0].y - self.e[i][0].y,
                    self.b[i][0].z - self.e[i][0].z,
                ),
                Vec3::new(
                    self.b[i][1].x - self.e[i][1].x,
                    self.b[i][1].y - self.e[i][1].y,
                    self.b[i][1].z - self.e[i][1].z,
                ),
                Vec3::new(
                    self.b[i][2].x - self.e[i][2].x,
                    self.b[i][2].y - self.e[i][2].y,
                    self.b[i][2].z - self.e[i][2].z,
                ),
                Vec3::new(
                    self.b[i][3].x - self.e[i][3].x,
                    self.b[i][3].y - self.e[i][3].y,
                    self.b[i][3].z - self.e[i][3].z,
                ),
                Vec3::new(
                    self.b[i][4].x - self.e[i][4].x,
                    self.b[i][4].y - self.e[i][4].y,
                    self.b[i][4].z - self.e[i][4].z,
                ),
                Vec3::new(
                    self.b[i][5].x - self.e[i][5].x,
                    self.b[i][5].y - self.e[i][5].y,
                    self.b[i][5].z - self.e[i][5].z,
                ),
                Vec3::new(
                    self.b[i][6].x - self.e[i][6].x,
                    self.b[i][6].y - self.e[i][6].y,
                    self.b[i][6].z - self.e[i][6].z,
                ),
            ];

            let b = self.b[i];

            // e[0] = q · (b0 + 2 b1 + 3 b2 + 4 b3 + 5 b4 + 6 b5 + 7 b6)
            let e0_x = q
                * (b[0].x
                    + 2.0 * b[1].x
                    + 3.0 * b[2].x
                    + 4.0 * b[3].x
                    + 5.0 * b[4].x
                    + 6.0 * b[5].x
                    + 7.0 * b[6].x);
            let e0_y = q
                * (b[0].y
                    + 2.0 * b[1].y
                    + 3.0 * b[2].y
                    + 4.0 * b[3].y
                    + 5.0 * b[4].y
                    + 6.0 * b[5].y
                    + 7.0 * b[6].y);
            let e0_z = q
                * (b[0].z
                    + 2.0 * b[1].z
                    + 3.0 * b[2].z
                    + 4.0 * b[3].z
                    + 5.0 * b[4].z
                    + 6.0 * b[5].z
                    + 7.0 * b[6].z);

            // e[1] = q² · (b1 + 3 b2 + 6 b3 + 10 b4 + 15 b5 + 21 b6)
            let e1_x = q2
                * (b[1].x
                    + 3.0 * b[2].x
                    + 6.0 * b[3].x
                    + 10.0 * b[4].x
                    + 15.0 * b[5].x
                    + 21.0 * b[6].x);
            let e1_y = q2
                * (b[1].y
                    + 3.0 * b[2].y
                    + 6.0 * b[3].y
                    + 10.0 * b[4].y
                    + 15.0 * b[5].y
                    + 21.0 * b[6].y);
            let e1_z = q2
                * (b[1].z
                    + 3.0 * b[2].z
                    + 6.0 * b[3].z
                    + 10.0 * b[4].z
                    + 15.0 * b[5].z
                    + 21.0 * b[6].z);

            // e[2] = q³ · (b2 + 4 b3 + 10 b4 + 20 b5 + 35 b6)
            let e2_x = q3 * (b[2].x + 4.0 * b[3].x + 10.0 * b[4].x + 20.0 * b[5].x + 35.0 * b[6].x);
            let e2_y = q3 * (b[2].y + 4.0 * b[3].y + 10.0 * b[4].y + 20.0 * b[5].y + 35.0 * b[6].y);
            let e2_z = q3 * (b[2].z + 4.0 * b[3].z + 10.0 * b[4].z + 20.0 * b[5].z + 35.0 * b[6].z);

            // e[3] = q⁴ · (b3 + 5 b4 + 15 b5 + 35 b6)
            let e3_x = q4 * (b[3].x + 5.0 * b[4].x + 15.0 * b[5].x + 35.0 * b[6].x);
            let e3_y = q4 * (b[3].y + 5.0 * b[4].y + 15.0 * b[5].y + 35.0 * b[6].y);
            let e3_z = q4 * (b[3].z + 5.0 * b[4].z + 15.0 * b[5].z + 35.0 * b[6].z);

            // e[4] = q⁵ · (b4 + 6 b5 + 21 b6)
            let e4_x = q5 * (b[4].x + 6.0 * b[5].x + 21.0 * b[6].x);
            let e4_y = q5 * (b[4].y + 6.0 * b[5].y + 21.0 * b[6].y);
            let e4_z = q5 * (b[4].z + 6.0 * b[5].z + 21.0 * b[6].z);

            // e[5] = q⁶ · (b5 + 7 b6)
            let e5_x = q6 * (b[5].x + 7.0 * b[6].x);
            let e5_y = q6 * (b[5].y + 7.0 * b[6].y);
            let e5_z = q6 * (b[5].z + 7.0 * b[6].z);

            // e[6] = q⁷ · b6   (only column where the diagonal is the full transform)
            let e6_x = q7 * b[6].x;
            let e6_y = q7 * b[6].y;
            let e6_z = q7 * b[6].z;

            self.e[i][0] = Vec3::new(e0_x, e0_y, e0_z);
            self.e[i][1] = Vec3::new(e1_x, e1_y, e1_z);
            self.e[i][2] = Vec3::new(e2_x, e2_y, e2_z);
            self.e[i][3] = Vec3::new(e3_x, e3_y, e3_z);
            self.e[i][4] = Vec3::new(e4_x, e4_y, e4_z);
            self.e[i][5] = Vec3::new(e5_x, e5_y, e5_z);
            self.e[i][6] = Vec3::new(e6_x, e6_y, e6_z);

            for k in 0..7 {
                self.b[i][k] = Vec3::new(
                    self.e[i][k].x + be[k].x,
                    self.e[i][k].y + be[k].y,
                    self.e[i][k].z + be[k].z,
                );
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
                let mut gx = bi[j].x;
                let mut gy = bi[j].y;
                let mut gz = bi[j].z;
                for k in (j + 1)..7 {
                    gx += D_MAT[k][j] * bi[k].x;
                    gy += D_MAT[k][j] * bi[k].y;
                    gz += D_MAT[k][j] * bi[k].z;
                }
                self.g[i][j] = Vec3::new(gx, gy, gz);
            }
        }
    }

    /// After evaluating acceleration at stage `n` (1..=7), update g_{n-1}
    /// via Newton divided differences of (F - F₀); then propagate the
    /// delta back into b₀..b_{n-1} using c_mat. Compensated summation
    /// keeps round-off under control across many Picard iterations.
    fn update_g_and_b(&mut self, stage: usize, a0: &[Vec3], acc_n: &[Vec3]) {
        let n = stage - 1; // index of the g coefficient we're updating
        let hn = H[stage];

        for i in 0..self.g.len() {
            // Compute Newton divided difference of order n+1 for body i.
            let dfx = acc_n[i].x - a0[i].x;
            let dfy = acc_n[i].y - a0[i].y;
            let dfz = acc_n[i].z - a0[i].z;

            let (mut tx, mut ty, mut tz) = (dfx / hn, dfy / hn, dfz / hn);
            for k in 0..n {
                tx = (tx - self.g[i][k].x) / (hn - H[k + 1]);
                ty = (ty - self.g[i][k].y) / (hn - H[k + 1]);
                tz = (tz - self.g[i][k].z) / (hn - H[k + 1]);
            }

            let dgx = tx - self.g[i][n].x;
            let dgy = ty - self.g[i][n].y;
            let dgz = tz - self.g[i][n].z;
            self.g[i][n] = Vec3::new(tx, ty, tz);

            // Propagate Δg_n into b₀..b_n using compensated summation.
            for j in 0..=n {
                let coeff = if j == n { 1.0 } else { C_MAT[n][j] };
                let dbx = coeff * dgx;
                let dby = coeff * dgy;
                let dbz = coeff * dgz;
                add_cs(&mut self.b[i][j].x, &mut self.csb[i][j].x, dbx);
                add_cs(&mut self.b[i][j].y, &mut self.csb[i][j].y, dby);
                add_cs(&mut self.b[i][j].z, &mut self.csb[i][j].z, dbz);
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

fn restore_xv(bodies: &mut [Body], x: &[Vec3], v: &[Vec3]) {
    for (i, b) in bodies.iter_mut().enumerate() {
        b.pos_x = x[i].x;
        b.pos_y = x[i].y;
        b.pos_z = x[i].z;
        b.vel_x = v[i].x;
        b.vel_y = v[i].y;
        b.vel_z = v[i].z;
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::system::System;
    use crate::domain::body::Body;
    use crate::units::UnitSystem;

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

        let mut b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
        b1.softening = 0.0;
        let mut b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(dt_budget)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::Ias15);

        let mut peak = 0.0_f64;
        // Samples for drift detection: (t, δE/E₀) every ~0.5% of the run.
        let mut samples: Vec<(f64, f64)> = Vec::new();
        let sample_dt = total_time / 200.0;
        let mut next_sample = 0.0;

        // Substep-granularity driver per the IAS15 sub-step contract
        // (Rein & Spiegel 2015 §2.3): advance by calling `step()`
        // until the target simulation time is reached. Each call
        // consumes one adaptive sub-step whose size the controller
        // chose; using a fixed `for _ in 0..n_steps` loop here would
        // silently assume every call consumes `dt_budget` and fall
        // short of the intended integration window.
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
            peak,
            PEAK_TOL,
            N_ORBITS,
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
            drift,
            DRIFT_TOL,
            slope,
            peak,
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

        let mut b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
        b1.softening = 0.0;
        let mut b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(DT)
            .with_max_depth(10);
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
            peak,
            PEAK_TOL,
            N_ORBITS,
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
            Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0),
            Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0),
            Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0),
        ];
        for b in &mut bodies {
            b.softening = 0.0;
        }

        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(DT)
            .with_max_depth(10);
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
            peak,
            PEAK_TOL,
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
    /// Under the IAS15 substep-granularity contract (Rein & Spiegel
    /// 2015 §2.3), each `step()` call executes exactly one adaptive
    /// sub-step and reports its size via `StepResult::consumed_dt`;
    /// `System::step` advances `System::t` by that value. A budget far
    /// larger than what the controller can accept at perihelion
    /// therefore yields `System::t` strictly below the budget after
    /// one call.
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

        let mut b1 = Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0);
        b1.softening = 0.0;
        let mut b2 = Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0);
        b2.softening = 0.0;

        let mut sys = System::new(vec![b1, b2], UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(dt_budget)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::Ias15);

        let t0 = sys.t();
        sys.step();
        let consumed = sys.t() - t0;

        assert!(consumed > 0.0, "IAS15 sub-step consumed zero time — caller would busy-loop");
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
                k,
                t_prev,
                t_now,
            );
            t_prev = t_now;
        }
    }

    // ── decide_dt pure-function tests (TD-004) ────────────────────────────
    //
    // The controller's decision logic is factored out as a pure function
    // on two floats + two bools so it can be exhaustively tested without
    // standing up a `System`. Each case covers one row of the decision
    // table documented on [`decide_dt`]; flipping an input and checking
    // the output changes is how we'll catch regressions in future
    // tuning work.

    #[test]
    fn decide_dt_accepts_on_merit() {
        // Picard converged AND truncation within tolerance → clean accept.
        let d = decide_dt(true, 5e-10, 1e-3, 1e-9, false);
        assert_eq!(d, DtDecision::Accept { degraded: false });
    }

    #[test]
    fn decide_dt_rejects_picard_when_not_converged() {
        // Non-convergence dominates: even an incidentally-small trunc
        // must not let us accept divergent `b` coefficients.
        let d = decide_dt(false, 1e-12, 1e-3, 1e-9, false);
        assert_eq!(d, DtDecision::RejectPicard);
    }

    #[test]
    fn decide_dt_rejects_truncation_when_converged_but_over_tol() {
        // Picard fine, but trunc_err above ε → standard controller path.
        let d = decide_dt(true, 1e-6, 1e-3, 1e-9, false);
        assert_eq!(d, DtDecision::RejectTruncation);
    }

    #[test]
    fn decide_dt_dt_min_escape_degrades() {
        // At the floor, we accept regardless of error state so the
        // simulation progresses — but flagged degraded for the caller.
        let d = decide_dt(false, 1.0, DT_MIN, 1e-9, false);
        assert_eq!(d, DtDecision::Accept { degraded: true });

        // Same floor, different failure class (trunc) — still degraded.
        let d = decide_dt(true, 1.0, DT_MIN, 1e-9, false);
        assert_eq!(d, DtDecision::Accept { degraded: true });
    }

    #[test]
    fn decide_dt_deadline_forces_degraded_accept() {
        // Cooperative deadline passed: accept current attempt rather
        // than spend more wall time shrinking.
        let d = decide_dt(false, 1.0, 1e-3, 1e-9, true);
        assert_eq!(d, DtDecision::Accept { degraded: true });
    }

    #[test]
    fn decide_dt_deadline_does_not_demote_clean_accept() {
        // On-merit result takes precedence over deadline — deadline is
        // an escape hatch for *stuck* attempts, not a degrade-poisoner
        // for attempts that converged within tolerance.
        let d = decide_dt(true, 5e-10, 1e-3, 1e-9, true);
        assert_eq!(d, DtDecision::Accept { degraded: false });
    }

    #[test]
    fn decide_dt_trunc_exactly_at_epsilon_is_merit() {
        // Boundary: trunc_err == ε should be accepted (≤, not <),
        // following the threshold convention specified in Rein &
        // Spiegel (2015) §3.4. Flipping this would silently change
        // step-size distributions in benchmarks.
        let d = decide_dt(true, 1e-9, 1e-3, 1e-9, false);
        assert_eq!(d, DtDecision::Accept { degraded: false });
    }

    // ── warmstart_b — direct polynomial-transformation tests ────────────────
    //
    // The acceleration ansatz over one IAS15 step is
    //
    //     a(u) = a₀ + b[0]·u + b[1]·u² + b[2]·u³ + … + b[6]·u⁷
    //
    // with `u ∈ [0, 1]` the dimensionless intra-step time. When the
    // controller picks a different `dt_new = q · dt_prev` for the next
    // step, the *physical* polynomial in absolute time is unchanged; only
    // the parametrisation of `u` changes. Substituting `u_new = (u_old −
    // 1) / q` (i.e. continuing the polynomial past the previous step's
    // end into the next step's start) and re-collecting powers yields
    // the polynomial-basis transformation
    //
    //     b_new[m] = q^{m+1} · Σ_{k ≥ m} C(k+1, m+1) · b_old[k]
    //
    // where `C(n, k)` is the binomial coefficient. The previous
    // implementation kept only the diagonal term (`k = m`) and silently
    // dropped every cross-term — see the design comment on
    // `warmstart_b` for the cascade-failure mode that introduced.
    //
    // These tests pin the transformation in three ways:
    //   1. By formula match against the explicit Pascal-coefficient
    //      table — catches arithmetic typos.
    //   2. By identity at q = 1 — no rescaling should mean no change.
    //   3. By polynomial-equivalence at multiple sample points — the
    //      strongest of the three: even if both the implementation and
    //      the formula table were wrong in the same way, the
    //      polynomial sampled in the new coordinates must agree with
    //      the same polynomial sampled in the old coordinates at the
    //      corresponding `u_old = 1 + q · u_new`.
    //
    // Together they make it nearly impossible to reintroduce the
    // missing-cross-terms regression without one of the three failing.

    /// Helper: construct an `Ias15` with `n` body slots and pre-populate
    /// `b[0]` and `e[0]` to user-supplied coefficients. Uses
    /// `ensure_capacity` (private) to allocate the buffers; sets `e = b`
    /// so `be = 0` in `warmstart_b` and the warmstart contribution is
    /// isolated from the predictor-corrector residual term.
    fn ias15_with_b(b_x: [f64; 7], b_y: [f64; 7]) -> Ias15 {
        let mut ias = Ias15::new();
        ias.ensure_capacity(1);
        for k in 0..7 {
            ias.b[0][k] = Vec3::new(b_x[k], b_y[k], 0.0);
            ias.e[0][k] = Vec3::new(b_x[k], b_y[k], 0.0);
        }
        ias
    }

    /// Reference Pascal-coefficient transformation, written out by hand
    /// from the binomial-expansion derivation. Kept independent of the
    /// implementation under test — a refactor of `warmstart_b` cannot
    /// also "refactor" this reference without the test still catching
    /// the divergence.
    fn pascal_warmstart_reference(b: [f64; 7], q: f64) -> [f64; 7] {
        let q2 = q * q;
        let q3 = q2 * q;
        let q4 = q3 * q;
        let q5 = q4 * q;
        let q6 = q5 * q;
        let q7 = q6 * q;
        [
            q * (b[0]
                + 2.0 * b[1]
                + 3.0 * b[2]
                + 4.0 * b[3]
                + 5.0 * b[4]
                + 6.0 * b[5]
                + 7.0 * b[6]),
            q2 * (b[1] + 3.0 * b[2] + 6.0 * b[3] + 10.0 * b[4] + 15.0 * b[5] + 21.0 * b[6]),
            q3 * (b[2] + 4.0 * b[3] + 10.0 * b[4] + 20.0 * b[5] + 35.0 * b[6]),
            q4 * (b[3] + 5.0 * b[4] + 15.0 * b[5] + 35.0 * b[6]),
            q5 * (b[4] + 6.0 * b[5] + 21.0 * b[6]),
            q6 * (b[5] + 7.0 * b[6]),
            q7 * b[6],
        ]
    }

    /// Evaluate the polynomial `a(u) = Σ_k b[k] · u^{k+1}` at a given
    /// `u`. The constant (`a_0`) term is excluded so the test can isolate
    /// the b-driven part — `warmstart_b` does not touch `a_0` (it lives
    /// in the body state, not in the integrator's coefficient buffer).
    fn poly_b_eval(b: [f64; 7], u: f64) -> f64 {
        let u2 = u * u;
        let u3 = u2 * u;
        let u4 = u3 * u;
        let u5 = u4 * u;
        let u6 = u5 * u;
        let u7 = u6 * u;
        b[0] * u + b[1] * u2 + b[2] * u3 + b[3] * u4 + b[4] * u5 + b[5] * u6 + b[6] * u7
    }

    #[test]
    fn warmstart_b_zero_b_in_zero_b_out() {
        // Cold-start sanity: when both `b` and `e` are zero (as on the
        // very first call, or after `ensure_capacity` zero-allocates
        // the buffers), warmstart must produce zero — there is no
        // signal to extrapolate from. A regression that adds a
        // round-off noise floor here would silently inject a non-zero
        // initial guess into Picard on cold starts and inflate the
        // first few `truncation_error` measurements; the controller
        // would react by shrinking `dt` for no physical reason.
        let mut ias = ias15_with_b([0.0; 7], [0.0; 7]);
        for &q in &[0.1_f64, 0.5, 1.0, 2.0, 10.0] {
            // Reset state for each q in the loop (keep be = 0, b = 0).
            for k in 0..7 {
                ias.b[0][k] = Vec3::ZERO;
                ias.e[0][k] = Vec3::ZERO;
            }
            ias.warmstart_b(q, 1.0);
            for k in 0..7 {
                assert_eq!(
                    ias.b[0][k].x, 0.0,
                    "b[{}].x non-zero at q={} from zero input — round-off injection regression",
                    k, q,
                );
                assert_eq!(
                    ias.b[0][k].y, 0.0,
                    "b[{}].y non-zero at q={} from zero input — round-off injection regression",
                    k, q,
                );
                assert_eq!(
                    ias.b[0][k].z, 0.0,
                    "b[{}].z non-zero at q={} from zero input — round-off injection regression",
                    k, q,
                );
            }
        }
    }

    #[test]
    fn warmstart_b_q_eq_one_reproduces_polynomial_continuation() {
        // At q = 1 (constant dt), the Pascal transformation re-expresses
        // the *same* polynomial in coordinates centred at the next
        // step's start — which corresponds to `u_old = 1`. The new
        // `b` therefore encodes the polynomial's derivatives at
        // `u_old = 1` (in the appropriate scaling), NOT the identity
        // mapping `b_new = b_old`. Naïvely "preserving b at q=1"
        // (the previous diagonal-only formula's behaviour) is the
        // root of the cascade investigated in
        // `docs/experiments/2026-04-26-ias15-warmstart-bug.md`.
        //
        // This test verifies the polynomial-continuation property by
        // checking that for every sampled `u_new ∈ [0, 1]` the new
        // polynomial's b-driven part equals the old polynomial
        // sampled at `u_old = 1 + u_new` minus the old polynomial at
        // `u_old = 1` (the latter goes into `a_0`, not `b`).
        let b_x = [0.7, -0.4, 0.2, -0.15, 0.1, -0.05, 0.03];
        let q = 1.0;

        let mut ias = ias15_with_b(b_x, [0.0; 7]);
        ias.warmstart_b(q, 1.0);

        let new_x: [f64; 7] = std::array::from_fn(|k| ias.b[0][k].x);

        for &u_new in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let u_old = 1.0 + q * u_new;
            let expected = poly_b_eval(b_x, u_old) - poly_b_eval(b_x, 1.0);
            let got = poly_b_eval(new_x, u_new);
            let diff = (got - expected).abs();
            let scale = expected.abs().max(1.0);
            assert!(
                diff <= 1e-12 * scale,
                "q=1 polynomial continuation failed at u_new={}: got {:.18e}, expected {:.18e}",
                u_new,
                got,
                expected,
            );
        }
    }

    #[test]
    fn warmstart_b_matches_pascal_reference_at_q_eq_two() {
        // q = 2 is the typical upper end of step-growth ratio
        // permitted by [`DT_GROWTH_LIMIT`] under the IAS15 controller.
        // With the ICs below every cross-term contribution is O(1),
        // so any missing column of the Pascal table fails the
        // assertion by orders of magnitude rather than by ULPs.
        let b_x = [0.5, -1.0, 0.7, -0.4, 0.3, -0.2, 0.1];
        let b_y = [-0.1, 0.3, -0.5, 0.2, 0.6, -0.4, 0.8];
        let q = 2.0;

        let mut ias = ias15_with_b(b_x, b_y);
        ias.warmstart_b(q, 1.0);

        let ref_x = pascal_warmstart_reference(b_x, q);
        let ref_y = pascal_warmstart_reference(b_y, q);

        for k in 0..7 {
            let diff_x = (ias.b[0][k].x - ref_x[k]).abs();
            let diff_y = (ias.b[0][k].y - ref_y[k]).abs();
            // Tolerance: ~50× f64 ULP scaled by the largest summand
            // for that coefficient (so the bound is meaningful for
            // both the q⁷·b₆ column with magnitude ~1 and the q¹·…
            // column whose accumulator can reach ~50 by the Pascal
            // arithmetic). We do not relax this beyond round-off.
            assert!(
                diff_x <= 1e-13_f64 * ref_x[k].abs().max(1.0),
                "b[{}].x at q={}: got {:.18e}, expected {:.18e}, diff {:.3e}",
                k,
                q,
                ias.b[0][k].x,
                ref_x[k],
                diff_x,
            );
            assert!(
                diff_y <= 1e-13_f64 * ref_y[k].abs().max(1.0),
                "b[{}].y at q={}: got {:.18e}, expected {:.18e}, diff {:.3e}",
                k,
                q,
                ias.b[0][k].y,
                ref_y[k],
                diff_y,
            );
        }
    }

    #[test]
    fn warmstart_b_matches_pascal_reference_at_q_lt_one() {
        // q < 1 is the close-encounter case: dt shrinking after a
        // truncation rejection, or recovery from a degraded floor
        // accept. The diagonal-only formula loses the most accuracy
        // here (small q⁷ exposes the missing cross-terms most
        // visibly), so this is the regime the cascade actually lived
        // in.
        let b_x = [1.2, -0.8, 0.5, -0.3, 0.2, -0.1, 0.05];
        let b_y = [-0.7, 0.4, -0.2, 0.1, -0.05, 0.025, -0.0125];
        let q = 0.3;

        let mut ias = ias15_with_b(b_x, b_y);
        ias.warmstart_b(q, 1.0);

        let ref_x = pascal_warmstart_reference(b_x, q);
        let ref_y = pascal_warmstart_reference(b_y, q);

        for k in 0..7 {
            let diff_x = (ias.b[0][k].x - ref_x[k]).abs();
            let diff_y = (ias.b[0][k].y - ref_y[k]).abs();
            assert!(
                diff_x <= 1e-13_f64 * ref_x[k].abs().max(1.0),
                "b[{}].x at q={}: got {:.18e}, expected {:.18e}, diff {:.3e}",
                k,
                q,
                ias.b[0][k].x,
                ref_x[k],
                diff_x,
            );
            assert!(
                diff_y <= 1e-13_f64 * ref_y[k].abs().max(1.0),
                "b[{}].y at q={}: got {:.18e}, expected {:.18e}, diff {:.3e}",
                k,
                q,
                ias.b[0][k].y,
                ref_y[k],
                diff_y,
            );
        }
    }

    #[test]
    fn warmstart_b_preserves_polynomial_at_sample_points() {
        // The strongest test: independent of the Pascal-coefficient
        // table, the *polynomial* the new `b` represents (in the new
        // step's `u_new ∈ [0, 1]` parametrisation) must agree with the
        // *same* polynomial that the old `b` represents (in the old
        // step's `u_old = 1 + q · u_new` parametrisation) at every
        // sampled `u_new`. If any cross-term is missed by the
        // implementation, this assertion fails by an arbitrarily large
        // margin — far more than ULP — because the missing
        // contributions are O(b_k) for k > m, not O(ULP · b).
        let b_x = [0.7, -0.4, 0.2, -0.15, 0.1, -0.05, 0.03];
        let b_y = [-0.6, 0.5, -0.4, 0.3, -0.2, 0.1, -0.05];

        for &q in &[0.1_f64, 0.5, 0.9, 1.0, 1.5, 2.0, 5.0] {
            let mut ias = ias15_with_b(b_x, b_y);
            ias.warmstart_b(q, 1.0);

            let new_x: [f64; 7] = std::array::from_fn(|k| ias.b[0][k].x);
            let new_y: [f64; 7] = std::array::from_fn(|k| ias.b[0][k].y);

            for &u_new in &[0.0_f64, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
                // Old coordinate corresponding to this u_new.
                let u_old = 1.0 + q * u_new;
                let expected_x = poly_b_eval(b_x, u_old) - poly_b_eval(b_x, 1.0);
                let expected_y = poly_b_eval(b_y, u_old) - poly_b_eval(b_y, 1.0);
                // The new polynomial in u_new is shifted so its
                // constant term equals the old polynomial at u_old=1
                // (which lives in `a_0`, not in `b`). So we compare
                // the b-driven part — the *increment* of the
                // acceleration past the new step's start.
                let got_x = poly_b_eval(new_x, u_new);
                let got_y = poly_b_eval(new_y, u_new);

                let scale_x = expected_x.abs().max(1.0);
                let scale_y = expected_y.abs().max(1.0);
                let diff_x = (got_x - expected_x).abs();
                let diff_y = (got_y - expected_y).abs();
                assert!(
                    diff_x <= 1e-12 * scale_x,
                    "polynomial mismatch at q={}, u_new={}: x got {:.18e}, expected {:.18e}, diff {:.3e}",
                    q,
                    u_new,
                    got_x,
                    expected_x,
                    diff_x,
                );
                assert!(
                    diff_y <= 1e-12 * scale_y,
                    "polynomial mismatch at q={}, u_new={}: y got {:.18e}, expected {:.18e}, diff {:.3e}",
                    q,
                    u_new,
                    got_y,
                    expected_y,
                    diff_y,
                );
            }
        }
    }

    #[test]
    fn warmstart_b_preserves_picard_residual_when_be_nonzero() {
        // `warmstart_b` carries an additive correction `be = b - e`
        // that represents the Picard residual from the previous step
        // (i.e. how much `b` moved beyond the previous warmstart's
        // prediction). Under the rescaling, that residual must be
        // preserved verbatim — it is the part of the polynomial the
        // *previous* step's controller already accounted for, and
        // dropping it would re-introduce the prediction error on
        // every retry.
        //
        // Construction: pick `b` and `e` independently, so `be ≠ 0`.
        // The expected output is `pascal(b) + (b - e)`.
        let b_x = [0.3, -0.2, 0.15, -0.1, 0.07, -0.05, 0.03];
        let e_x = [0.25, -0.15, 0.10, -0.05, 0.03, -0.02, 0.01];
        let q = 0.7;

        let mut ias = Ias15::new();
        ias.ensure_capacity(1);
        for k in 0..7 {
            ias.b[0][k] = Vec3::new(b_x[k], 0.0, 0.0);
            ias.e[0][k] = Vec3::new(e_x[k], 0.0, 0.0);
        }
        ias.warmstart_b(q, 1.0);

        let pascal = pascal_warmstart_reference(b_x, q);
        for k in 0..7 {
            let expected = pascal[k] + (b_x[k] - e_x[k]);
            let diff = (ias.b[0][k].x - expected).abs();
            assert!(
                diff <= 1e-13_f64 * expected.abs().max(1.0),
                "be-correction lost at b[{}]: got {:.18e}, expected {:.18e}, diff {:.3e}",
                k,
                ias.b[0][k].x,
                expected,
                diff,
            );
        }
    }

    // ── 3D validation portfolio (dynamic) ─────────────────────────────────────
    //
    // Where the planar IAS15 unit tests above confirm that `z = vz = 0`
    // input still produces the same energy / drift / peak numbers it did
    // before the 3D port, the tests below confirm that `z != 0` /
    // `vz != 0` motion is integrated *correctly*: the controller responds
    // to close encounters out of the orbital plane the same way it does
    // in the plane, and angular momentum vector conservation holds across
    // a long horizon.
    //
    // Each of the three tests below runs the planar configuration AND the
    // inclined configuration in the same #[test] body, then asserts a
    // **relative** bound: `metric_inclined <= metric_planar · (1 + δ)`.
    // This is stronger than a hardcoded threshold — the planar measurement
    // floats with platform / compiler / cpu-feature variance, but the
    // *ratio* is structurally invariant under rotation. A regression that
    // worsens both planar and inclined silently within their absolute
    // thresholds is caught by the ratio.

    /// Rotate a `(y, z)` plane vector by `angle` around the `x̂` axis.
    /// Inline trigonometry — no dependency on any other rotation helper
    /// in the codebase, so the inclined configuration is constructed
    /// independently of `OrbitalElements::sample_orbit` or any other
    /// path the test exercises.
    fn rotate_around_x(v: crate::math::Vec3, angle: f64) -> crate::math::Vec3 {
        let (s, c) = angle.sin_cos();
        crate::math::Vec3::new(v.x, v.y * c - v.z * s, v.y * s + v.z * c)
    }

    /// Pythagorean three-body, planar vs inclined: peak energy error of
    /// the inclined run must not exceed the planar run by more than 50 %.
    ///
    /// The Pythagorean (Burrau 1913) configuration is the canonical
    /// stress test for the IAS15 controller: violent close encounters
    /// at t ≈ 2–5 force the adaptive step to shrink by orders of
    /// magnitude before recovering. Rotating the entire scenario 30°
    /// out of the plane changes nothing physically — the same close
    /// encounters occur, the same energy must be conserved — so the
    /// ratio of peak errors is a clean signal: any value > 1 + δ
    /// indicates the controller responds differently to in-plane vs
    /// out-of-plane geometry, which is a 3D-specific bug the planar
    /// gate cannot see.
    #[test]
    fn ias15_pythagorean_inclined_matches_planar_within_relative_bound() {
        use crate::math::Vec3;

        const DT: f64 = 0.01;
        const T_END: f64 = 10.0;
        // Allowance for the additional per-axis arithmetic in 3D: ~1 ULP
        // accumulating over the full integration. δ = 0.5 (50 %) is
        // generous against the f64 noise floor while still bounding any
        // systematic regression to under one order of magnitude.
        const RELATIVE_SLACK: f64 = 0.5;
        const INCLINATION: f64 = std::f64::consts::PI / 6.0; // 30°

        // Planar Pythagorean (canonical Burrau initial conditions).
        let planar_bodies = vec![
            Body::rocky(3.0).at(1.0, 3.0).with_velocity(0.0, 0.0).unsoftened(),
            Body::rocky(4.0).at(-2.0, -1.0).with_velocity(0.0, 0.0).unsoftened(),
            Body::rocky(5.0).at(1.0, -1.0).with_velocity(0.0, 0.0).unsoftened(),
        ];

        // Inclined Pythagorean: each body's position rotated 30° around `x̂`.
        // Velocities are zero in the original setup so they need no rotation,
        // but go through the helper anyway so a non-zero IC variant trivially
        // generalises later.
        let inclined_bodies: Vec<Body> = planar_bodies
            .iter()
            .map(|b| {
                let pos = rotate_around_x(Vec3::new(b.pos_x, b.pos_y, b.pos_z), INCLINATION);
                let vel = rotate_around_x(Vec3::new(b.vel_x, b.vel_y, b.vel_z), INCLINATION);
                Body::rocky(b.mass)
                    .at_3d(pos.x, pos.y, pos.z)
                    .with_velocity_3d(vel.x, vel.y, vel.z)
                    .unsoftened()
            })
            .collect();

        let peak = |bodies: Vec<Body>| -> f64 {
            let mut sys = System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(DT)
                .with_max_depth(10);
            sys.set_integrator(IntegratorKind::Ias15);
            let n_steps = (T_END / DT).ceil() as u64;
            let mut p = 0.0_f64;
            for _ in 0..n_steps {
                sys.step();
                p = p.max(sys.metrics().rel_energy_error.abs());
            }
            p
        };

        let peak_planar = peak(planar_bodies);
        let peak_inclined = peak(inclined_bodies);

        let bound = peak_planar * (1.0 + RELATIVE_SLACK);
        assert!(
            peak_inclined <= bound,
            "Pythagorean 3D regression: inclined peak |δE/E₀| = {:.3e} \
             exceeds planar peak {:.3e} × (1 + {RELATIVE_SLACK}) = {:.3e}",
            peak_inclined,
            peak_planar,
            bound,
        );
    }

    /// Kepler `e = 0.9` high-eccentricity, planar vs inclined: same
    /// relative bound as the Pythagorean test above.
    ///
    /// The high-e regime is what stresses the IAS15 sub-step velocity
    /// prediction (`predict_v_ias15`) and the Picard convergence under
    /// rapid acceleration variation near perihelion. A 3D-specific bug
    /// in the velocity-dependent path of any perturbation would surface
    /// here even though no perturbation is registered: the integrator's
    /// own inner loop reads body velocity at substep nodes, and a 3D
    /// regression in the substep predictor would manifest as elevated
    /// inclined peak compared to the planar gauge.
    #[test]
    fn ias15_kepler_high_e_inclined_matches_planar_within_relative_bound() {
        use crate::math::Vec3;

        const A: f64 = 1.0;
        const E: f64 = 0.9;
        const MU: f64 = 2.0;
        const N_ORBITS: u64 = 50;
        const DT: f64 = 0.1;
        const RELATIVE_SLACK: f64 = 0.5;
        const INCLINATION: f64 = std::f64::consts::PI / 6.0;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();
        let period = 2.0 * std::f64::consts::PI * (A.powi(3) / MU).sqrt();
        let total_time = N_ORBITS as f64 * period;

        // Planar two-body high-e Kepler.
        let planar_bodies = vec![
            Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0).unsoftened(),
            Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0).unsoftened(),
        ];

        let inclined_bodies: Vec<Body> = planar_bodies
            .iter()
            .map(|b| {
                let pos = rotate_around_x(Vec3::new(b.pos_x, b.pos_y, b.pos_z), INCLINATION);
                let vel = rotate_around_x(Vec3::new(b.vel_x, b.vel_y, b.vel_z), INCLINATION);
                Body::rocky(b.mass)
                    .at_3d(pos.x, pos.y, pos.z)
                    .with_velocity_3d(vel.x, vel.y, vel.z)
                    .unsoftened()
            })
            .collect();

        let peak = |bodies: Vec<Body>| -> f64 {
            let mut sys = System::new(bodies, UnitSystem::canonical())
                .with_theta(0.5)
                .with_dt(DT)
                .with_max_depth(10);
            sys.set_integrator(IntegratorKind::Ias15);
            let n_steps = (total_time / DT).ceil() as u64;
            let mut p = 0.0_f64;
            for _ in 0..n_steps {
                sys.step();
                p = p.max(sys.metrics().rel_energy_error.abs());
            }
            p
        };

        let peak_planar = peak(planar_bodies);
        let peak_inclined = peak(inclined_bodies);

        let bound = peak_planar * (1.0 + RELATIVE_SLACK);
        assert!(
            peak_inclined <= bound,
            "high-e Kepler 3D regression: inclined peak |δE/E₀| = {:.3e} \
             exceeds planar peak {:.3e} × (1 + {RELATIVE_SLACK}) = {:.3e}",
            peak_inclined,
            peak_planar,
            bound,
        );
    }

    /// Long-run conservation of the angular-momentum vector `h_vec` for
    /// an inclined Kepler orbit, with magnitude and direction asserted
    /// independently.
    ///
    /// `h_vec = r × v` is a vector conserved by two-body dynamics. The
    /// test integrates an inclined `e = 0.5` Kepler orbit over 200
    /// orbits and asserts:
    ///
    ///   - **magnitude drift** `||h(t)| − |h(0)|| / |h(0)|` stays below
    ///     `1e-12` — energetic conservation;
    ///   - **direction drift** `1 − cos(angle(h(t), h(0)))` stays below
    ///     `1e-12` — geometric conservation (orbital plane is fixed).
    ///
    /// Splitting the assertion is what makes a failure diagnosable: a
    /// `mag_drift` failure with `dir_drift = 0` indicates an
    /// energy-leaking integrator step; `dir_drift` failure with
    /// `mag_drift = 0` indicates spurious cross-axis torque (a kernel
    /// asymmetry that the planar tests cannot see by construction). A
    /// joint test like `(h_end - h_start).length()` blends both signals
    /// and tells you only that something drifted.
    #[test]
    fn ias15_inclined_kepler_conserves_h_vec_magnitude_and_direction() {
        use crate::math::Vec3;
        use crate::physics::orbital::compute_invariants;

        const A: f64 = 2.0;
        const E: f64 = 0.5;
        const MU: f64 = 2.0;
        const N_ORBITS: u64 = 200;
        const INCLINATION: f64 = std::f64::consts::PI / 6.0;

        const MAG_TOL: f64 = 1e-12;
        const DIR_TOL: f64 = 1e-12;

        let r_peri = A * (1.0 - E);
        let v_peri = (MU * (1.0 + E) / (A * (1.0 - E))).sqrt();
        let period = 2.0 * std::f64::consts::PI * (A.powi(3) / MU).sqrt();
        let total_time = N_ORBITS as f64 * period;
        let dt_budget = period / 20.0;

        let planar_bodies = [
            Body::rocky(1.0).at(-r_peri / 2.0, 0.0).with_velocity(0.0, -v_peri / 2.0).unsoftened(),
            Body::rocky(1.0).at(r_peri / 2.0, 0.0).with_velocity(0.0, v_peri / 2.0).unsoftened(),
        ];

        let bodies: Vec<Body> = planar_bodies
            .iter()
            .map(|b| {
                let pos = rotate_around_x(Vec3::new(b.pos_x, b.pos_y, b.pos_z), INCLINATION);
                let vel = rotate_around_x(Vec3::new(b.vel_x, b.vel_y, b.vel_z), INCLINATION);
                Body::rocky(b.mass)
                    .at_3d(pos.x, pos.y, pos.z)
                    .with_velocity_3d(vel.x, vel.y, vel.z)
                    .unsoftened()
            })
            .collect();

        // Capture initial h_vec via the orbital-element pipeline rather
        // than reading body fields directly — the test is about the
        // full pipeline holding up over a long horizon.
        let inv0 = compute_invariants(&bodies, 1, 0, 1.0).unwrap();
        let h0 = inv0.h_vec;
        let h0_mag = h0.length();
        assert!(h0_mag > 0.0, "test setup: initial |h_vec| must be non-zero");

        let mut sys = System::new(bodies, UnitSystem::canonical())
            .with_theta(0.5)
            .with_dt(dt_budget)
            .with_max_depth(10);
        sys.set_integrator(IntegratorKind::Ias15);

        while sys.t() < total_time {
            sys.step();
        }

        let inv_end = compute_invariants(sys.bodies(), 1, 0, 1.0).unwrap();
        let h_end = inv_end.h_vec;
        let h_end_mag = h_end.length();

        // Magnitude drift: relative change in |h|.
        let mag_drift = (h_end_mag - h0_mag).abs() / h0_mag;
        assert!(
            mag_drift < MAG_TOL,
            "h_vec magnitude drift {:.3e} exceeds {MAG_TOL:.0e} over {N_ORBITS} orbits — \
             energetic conservation regression",
            mag_drift,
        );

        // Direction drift: 1 − cos(angle between h_end and h0).
        // Stable near 0, well-behaved for arbitrarily small angles.
        let cos_angle = h_end.dot(h0) / (h_end_mag * h0_mag);
        let dir_drift = 1.0 - cos_angle;
        assert!(
            dir_drift < DIR_TOL,
            "h_vec direction drift {:.3e} (1 − cos angle) exceeds {DIR_TOL:.0e} over \
             {N_ORBITS} orbits — orbital plane is precessing spuriously (cross-axis \
             torque leak)",
            dir_drift,
        );
    }
}
