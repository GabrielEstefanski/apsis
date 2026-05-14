//! Neumaier-compensated `f64` accumulator for long-horizon scalar sums.
//!
//! `CompensatedF64` is a drop-in scalar replacement for `f64` whose
//! `+=` and `-=` operators apply Neumaier's variant of Kahan
//! summation, recovering the bits a naive `f64 += dt` loop loses to
//! magnitude-mismatch cancellation. Used wherever a scalar accumulator
//! must stay accurate across many small additions: simulation time
//! `t += dt` over 10⁶ steps, energy baselines, angular-momentum drift
//! tracking, any long-horizon state-of-charge accumulator.
//!
//! # Why
//!
//! Naive `f64` summation `s = s + x` loses the low bits of `x` to
//! `s`'s exponent every time `|s| ≫ |x|`. After `N` steps the
//! accumulated error grows as `~N · ε_machine · max|s|` (worst case)
//! or `~√N · ε_machine · max|s|` (random-sign case). For
//! `dt = 1e-3` over 10⁶ steps, naive `f64` drifts the simulation
//! clock by `~10⁻¹⁰` of the elapsed time — small in absolute terms,
//! large compared to the integrator's own truncation floor (IAS15
//! at `ε = 1e-9` per substep).
//!
//! Neumaier's algorithm carries a running compensator `c` that
//! captures the bits the naive add discards, and folds them back
//! into the next addition. The accumulated error of the compensated
//! sum is `O(ε_machine · max|s|)` independently of `N`. The drift
//! gate (`tests::million_step_drift_beats_naive`) demonstrates the
//! N-independence empirically.
//!
//! # Why not Kahan
//!
//! Plain Kahan compensates only when the new term is smaller than
//! the running sum (`|x| < |s|`). When the order reverses (a large
//! term arrives after many small ones), Kahan loses bits exactly
//! once. Neumaier (1974) handles both directions with one extra
//! `abs` comparison per add and is the variant REBOUND uses in
//! `integrator_whfast.c`. The cost is identical in practice because
//! the comparison branch predicts trivially in either direction.
//!
//! # Reference
//!
//! Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur
//! Summation endlicher Summen. *ZAMM* 54, 39–51.
//!
//! Higham, N. J. (2002). *Accuracy and Stability of Numerical
//! Algorithms*, 2nd ed., §4.3 ("Compensated summation"). SIAM.

use std::ops::{Add, AddAssign, Sub, SubAssign};

/// 16-byte scalar that accumulates additions and subtractions with
/// Neumaier's compensated summation.
///
/// Layout: `value` carries the bits a naive `f64` would also carry;
/// `comp` carries the compensator, the bits that would otherwise be
/// lost to magnitude mismatch. `total()` returns `value + comp` —
/// the settled view, used when reading the accumulator into a plain
/// `f64`.
///
/// `#[repr(C)]` pins the layout so a `&[CompensatedF64]` is bitwise
/// compatible with a `&[[f64; 2]]` for FFI / serialization use cases.
/// All operators are `#[inline(always)]` so LLVM elides the type
/// boundary entirely; the generated assembly for
/// `let mut s: CompensatedF64 = 0.0.into(); s += dt;` is the same
/// instruction sequence REBOUND emits in C.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[repr(C)]
pub struct CompensatedF64 {
    /// Running sum at full precision lost — the bits a naive `f64`
    /// accumulator would carry.
    pub value: f64,

    /// Compensator: the bits the naive add discarded last step,
    /// folded back in on the next add. Zero on a freshly-default
    /// value.
    pub comp: f64,
}

impl CompensatedF64 {
    /// Zero accumulator. Equivalent to `CompensatedF64::default()`.
    pub const ZERO: Self = Self { value: 0.0, comp: 0.0 };

    /// Construct directly from `(value, compensator)` pair. Use only
    /// when materializing from serialized state; for fresh
    /// accumulators prefer [`from_value`](Self::from_value) or
    /// `0.0.into()`.
    #[inline(always)]
    pub const fn new(value: f64, comp: f64) -> Self {
        Self { value, comp }
    }

    /// Construct from a single `f64` with zero compensator. The
    /// settled value is exactly `v`.
    #[inline(always)]
    pub const fn from_value(v: f64) -> Self {
        Self { value: v, comp: 0.0 }
    }

    /// Settled total `value + comp`. Use when reading the accumulator
    /// into a plain `f64` (display, serialization, comparison against
    /// an analytic target).
    ///
    /// Note: `total()` is NOT cumulative — it does not consume the
    /// compensator. Calling it any number of times leaves the
    /// accumulator unchanged.
    #[inline(always)]
    pub fn total(self) -> f64 {
        self.value + self.comp
    }

    /// Apply Neumaier's compensated `+=`. Internal helper used by the
    /// `AddAssign` / `SubAssign` impls; exposed for sites that want
    /// the operation without the syntactic sugar.
    ///
    /// Algorithm (Neumaier 1974):
    /// 1. `t = value + x`
    /// 2. If `|value| ≥ |x|`, the low bits of `x` were dropped during
    ///    the add: recover them via `(value − t) + x`.
    /// 3. Otherwise the low bits of `value` were dropped: recover
    ///    them via `(x − t) + value`.
    /// 4. Either way, accumulate the recovered bits into `comp`.
    #[inline(always)]
    pub fn neumaier_add(&mut self, x: f64) {
        let t = self.value + x;
        let recovered =
            if self.value.abs() >= x.abs() { (self.value - t) + x } else { (x - t) + self.value };
        self.comp += recovered;
        self.value = t;
    }
}

// ── Conversions ─────────────────────────────────────────────────────────────

impl From<f64> for CompensatedF64 {
    #[inline(always)]
    fn from(v: f64) -> Self {
        Self::from_value(v)
    }
}

impl From<CompensatedF64> for f64 {
    /// Settle into a plain `f64`. Equivalent to [`total`](CompensatedF64::total).
    #[inline(always)]
    fn from(c: CompensatedF64) -> Self {
        c.total()
    }
}

// ── Arithmetic with f64 ─────────────────────────────────────────────────────

impl AddAssign<f64> for CompensatedF64 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: f64) {
        self.neumaier_add(rhs);
    }
}

impl SubAssign<f64> for CompensatedF64 {
    #[inline(always)]
    fn sub_assign(&mut self, rhs: f64) {
        self.neumaier_add(-rhs);
    }
}

impl Add<f64> for CompensatedF64 {
    type Output = Self;
    #[inline(always)]
    fn add(mut self, rhs: f64) -> Self {
        self.neumaier_add(rhs);
        self
    }
}

impl Sub<f64> for CompensatedF64 {
    type Output = Self;
    #[inline(always)]
    fn sub(mut self, rhs: f64) -> Self {
        self.neumaier_add(-rhs);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Constructing from a single value and reading back returns the
    /// exact bit pattern. The compensator stays zero.
    #[test]
    fn from_value_and_total_round_trip_exactly() {
        for &v in &[0.0_f64, 1.0, -1.0, 4.567_890_123_456_789, 1e-100, 1e100] {
            let c = CompensatedF64::from_value(v);
            assert_eq!(c.total(), v, "round-trip failed for {v}");
            assert_eq!(c.comp, 0.0, "fresh accumulator must have zero compensator");
        }
    }

    /// Adding zero is exactly idempotent — neither `value` nor `comp`
    /// changes. Failing this means the algorithm is leaking round-off
    /// even on null operations.
    #[test]
    fn adding_zero_is_pointwise_identity() {
        let mut c = CompensatedF64::from_value(1.234_567_e8);
        let before = c;
        c += 0.0;
        assert_eq!(c, before, "+=0 must be exact identity");
    }

    /// Single-add agreement with naive `f64` for arguments where no
    /// cancellation can happen. Locks the algorithm against silently
    /// distorting well-conditioned single sums.
    #[test]
    fn single_add_well_conditioned_matches_naive() {
        for &(a, b) in &[(1.0, 2.0), (0.5, 0.25), (1e6, 3.7), (-2.5, 4.0), (1e-3, 1e-3)] {
            let mut c: CompensatedF64 = a.into();
            c += b;
            assert_eq!(c.total(), a + b, "single add diverged for ({a}, {b})");
        }
    }

    /// **Headline gate.** Naive `f64` accumulator drifts when summing
    /// many small `dt` into a growing total because the low bits of
    /// `dt` fall off the exponent of `t` every step. Neumaier-
    /// compensated keeps every recovered low bit and stays at the
    /// machine epsilon floor independently of the step count.
    ///
    /// At `dt = 0.1`, sum after `10⁶` steps must equal `1e5` exactly
    /// in arithmetic. Naive `f64` accumulates `~N · ε_machine · |t|`
    /// of drift; compensated stays at `O(ε_machine)`.
    ///
    /// The gate asserts the difference is at least 100×, demonstrating
    /// that compensated *is* materially more accurate — not just "no
    /// worse than".
    #[test]
    fn million_step_drift_beats_naive() {
        let dt = 0.1_f64;
        let n_steps = 1_000_000_usize;
        let expected = (n_steps as f64) * dt;

        let mut t_naive = 0.0_f64;
        let mut t_comp = CompensatedF64::ZERO;
        for _ in 0..n_steps {
            t_naive += dt;
            t_comp += dt;
        }

        let err_naive = (t_naive - expected).abs();
        let err_comp = (t_comp.total() - expected).abs();

        // Naive drift is O(N · ε · |t|) ~ 10⁶ · 2.2e-16 · 1e5 ~ 2e-5
        // in the worst case; in practice ~1.4e-5 on this scenario.
        // Lower bound 1e-7 is many orders below the typical value but
        // many orders above the compensated floor — the gate fires
        // only on a regression that brings naive *down to* the
        // compensated level (i.e. a bug that made compensation
        // unnecessary, which is itself a contradiction).
        assert!(
            err_naive > 1e-7,
            "naive f64 should drift after 10⁶ adds; got {err_naive:.3e} (gate: > 1e-7)"
        );
        // Compensated must stay at machine precision ~ 2.2e-16 ·
        // |total| ~ 2.2e-11 worst case. Bound 1e-9 covers two orders
        // of slack for platform variance.
        assert!(
            err_comp < 1e-9,
            "compensated must stay at ε-machine; got {err_comp:.3e} (gate: < 1e-9)"
        );
        assert!(
            err_naive / err_comp.max(f64::MIN_POSITIVE) > 100.0,
            "compensated must beat naive by at least 100×; ratio = {:.1}",
            err_naive / err_comp.max(f64::MIN_POSITIVE)
        );
    }

    /// Compensated summation should be **independent of addition
    /// order** for sums with mixed magnitudes (within ε-machine).
    /// Naive `f64` is order-dependent: summing many small terms
    /// before a big one preserves their contribution; summing the
    /// big one first crushes them under the exponent. This test
    /// exercises both orders and asserts compensated stays aligned.
    #[test]
    fn order_invariance_under_mixed_magnitudes() {
        let big = 1e10_f64;
        let small = 1e-6_f64;
        let n_small = 1_000_usize;

        // Order A: many small first, then big
        let mut a = CompensatedF64::ZERO;
        for _ in 0..n_small {
            a += small;
        }
        a += big;

        // Order B: big first, then many small
        let mut b = CompensatedF64::ZERO;
        b += big;
        for _ in 0..n_small {
            b += small;
        }

        let diff = (a.total() - b.total()).abs();
        assert!(diff < 1e-9, "compensated order-dependence above ε floor: |A − B| = {diff:.3e}");

        // Sanity: both totals are close to the analytic sum.
        let expected = big + (n_small as f64) * small;
        let err_a = (a.total() - expected).abs();
        let err_b = (b.total() - expected).abs();
        assert!(err_a < 1e-9 && err_b < 1e-9, "compensated drifted from analytic");
    }

    /// `From<CompensatedF64> for f64` settles via `total()`. The
    /// type's primary use case is "drop-in scalar that you can
    /// `.into::<f64>()` at the end" — locks that path.
    #[test]
    fn into_f64_settles_via_total() {
        let mut c = CompensatedF64::ZERO;
        for _ in 0..1000 {
            c += 0.001;
        }
        let f: f64 = c.into();
        assert!((f - 1.0).abs() < 1e-12, "into f64 settled to {f}");
    }

    /// Subtraction is an alias for `+= -rhs`. Symmetric round-trip:
    /// `(s += x; s -= x)` must return to the starting value within
    /// ε-machine.
    #[test]
    fn add_then_subtract_round_trips() {
        let mut c = CompensatedF64::from_value(42.0);
        let before = c.total();
        c += 1e-7;
        c -= 1e-7;
        let after = c.total();
        assert!((after - before).abs() < 1e-15, "round-trip drift: {} → {after}", before);
    }
}
