//! Neumaier-compensated `f64` accumulator.
//!
//! [`CompensatedF64`] is an `f64` whose `+=` and `-=` operators apply
//! the Neumaier (1974) variant of compensated summation. The
//! accumulated round-off error is `O(ε_machine · max|s|)` independent
//! of the number of additions, where naive `f64` accumulates
//! `O(N · ε_machine · max|s|)`.
//!
//! Use anywhere a scalar must stay accurate across many small
//! additions (simulation time, energy baselines, drift counters).
//!
//! # References
//!
//! Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur
//! Summation endlicher Summen. *ZAMM* 54, 39–51.
//!
//! Higham, N. J. (2002). *Accuracy and Stability of Numerical
//! Algorithms*, 2nd ed., §4.3. SIAM.

use std::ops::{Add, AddAssign, Sub, SubAssign};

/// Neumaier-compensated scalar.
///
/// `value + comp` is the settled total; [`total`](Self::total) returns
/// it. The compensator captures bits a naive `f64 += x` would lose
/// when `|value| ≫ |x|` (or when the magnitudes invert), and folds
/// them into the next addition.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[repr(C)]
pub struct CompensatedF64 {
    pub value: f64,
    pub comp: f64,
}

impl CompensatedF64 {
    pub const ZERO: Self = Self { value: 0.0, comp: 0.0 };

    #[inline(always)]
    pub const fn new(value: f64, comp: f64) -> Self {
        Self { value, comp }
    }

    #[inline(always)]
    pub const fn from_value(v: f64) -> Self {
        Self { value: v, comp: 0.0 }
    }

    /// Settled value `value + comp`. Non-mutating.
    #[inline(always)]
    pub fn total(self) -> f64 {
        self.value + self.comp
    }

    /// Add `x` and accumulate the lost low bits into `comp`. Used by
    /// the [`AddAssign<f64>`] / [`SubAssign<f64>`] impls; exposed for
    /// call sites that prefer the explicit form.
    #[inline(always)]
    pub fn neumaier_add(&mut self, x: f64) {
        let t = self.value + x;
        let recovered =
            if self.value.abs() >= x.abs() { (self.value - t) + x } else { (x - t) + self.value };
        self.comp += recovered;
        self.value = t;
    }
}

impl From<f64> for CompensatedF64 {
    #[inline(always)]
    fn from(v: f64) -> Self {
        Self::from_value(v)
    }
}

impl From<CompensatedF64> for f64 {
    #[inline(always)]
    fn from(c: CompensatedF64) -> Self {
        c.total()
    }
}

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

    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::{ToPrimitive, Zero};

    /// Settled value round-trips a single `f64` exactly.
    #[test]
    fn from_value_total_round_trip() {
        for &v in &[0.0_f64, 1.0, -1.0, 4.567_890_123_456_789, 1e-100, 1e100] {
            let c = CompensatedF64::from_value(v);
            assert_eq!(c.total(), v);
            assert_eq!(c.comp, 0.0);
        }
    }

    /// Adding zero leaves both `value` and `comp` unchanged.
    #[test]
    fn add_zero_is_identity() {
        let mut c = CompensatedF64::from_value(1.234_567_e8);
        let before = c;
        c += 0.0;
        assert_eq!(c, before);
    }

    /// Single add of well-conditioned arguments matches naive `f64`.
    #[test]
    fn single_add_well_conditioned_matches_naive() {
        for &(a, b) in &[(1.0, 2.0), (0.5, 0.25), (1e6, 3.7), (-2.5, 4.0), (1e-3, 1e-3)] {
            let mut c: CompensatedF64 = a.into();
            c += b;
            assert_eq!(c.total(), a + b);
        }
    }

    /// Add then subtract returns to the starting value within ε.
    #[test]
    fn add_then_subtract_round_trips() {
        let mut c = CompensatedF64::from_value(42.0);
        let before = c.total();
        c += 1e-7;
        c -= 1e-7;
        assert!((c.total() - before).abs() < 1e-15);
    }

    /// Settling via `Into<f64>` matches [`total`].
    #[test]
    fn into_f64_settles_via_total() {
        let mut c = CompensatedF64::ZERO;
        for _ in 0..1000 {
            c += 0.001;
        }
        let f: f64 = c.into();
        assert!((f - 1.0).abs() < 1e-12);
    }

    // ── Drift gates ─────────────────────────────────────────────────────────

    /// Convert an `f64` to its exact `BigRational` representation.
    /// `f64` is `m · 2^e` for integers `m`, `e`; `BigRational` carries
    /// both directions losslessly. Used as the oracle for drift gates
    /// below.
    fn f64_to_rational(x: f64) -> BigRational {
        if x == 0.0 {
            return BigRational::zero();
        }
        let bits = x.to_bits();
        let sign = if bits >> 63 == 1 { -1 } else { 1 };
        let exponent = ((bits >> 52) & 0x7ff) as i64;
        let mantissa = bits & 0x000f_ffff_ffff_ffff;
        let (m, e) = if exponent == 0 {
            (BigInt::from(mantissa as i64 * sign), -1074_i64)
        } else {
            let m = ((1_u64 << 52) | mantissa) as i64;
            (BigInt::from(m * sign), exponent - 1023 - 52)
        };
        if e >= 0 {
            BigRational::from(m << (e as usize))
        } else {
            BigRational::new(m, BigInt::from(1) << ((-e) as usize))
        }
    }

    /// Sum a slice of `f64` exactly via `BigRational`, return the
    /// `f64` closest to the analytic sum. This is the oracle: every
    /// term is converted losslessly, the sum is exact, only the final
    /// projection back to `f64` rounds.
    fn exact_sum_as_f64(xs: &[f64]) -> f64 {
        let mut acc = BigRational::zero();
        for &x in xs {
            acc += f64_to_rational(x);
        }
        acc.to_f64().expect("sum representable as f64")
    }

    /// Constant-`dt` drift over 10⁶ steps. Naive `f64` drifts by
    /// `O(N · ε · |total|)`; compensated stays at the IEEE-754 floor
    /// for the analytic sum.
    ///
    /// Oracle: arbitrary-precision rational sum, projected back to
    /// the nearest `f64`. Independent of the algorithm under test.
    #[test]
    fn constant_dt_million_step_drift_vs_oracle() {
        let dt = 0.1_f64;
        let n_steps = 1_000_000_usize;
        let dts = vec![dt; n_steps];
        let oracle = exact_sum_as_f64(&dts);

        let mut t_naive = 0.0_f64;
        let mut t_comp = CompensatedF64::ZERO;
        for &d in &dts {
            t_naive += d;
            t_comp += d;
        }
        let drift_naive = (t_naive - oracle).abs();
        let drift_comp = (t_comp.total() - oracle).abs();
        assert!(drift_naive > 1e-7, "naive drift below regression bound: {drift_naive:.3e}");
        assert!(drift_comp < 1e-9, "compensated drift above ε floor: {drift_comp:.3e}");
        assert!(
            drift_naive / drift_comp.max(f64::MIN_POSITIVE) > 100.0,
            "compensated must beat naive by ≥ 100×",
        );
    }

    /// Adaptive-cadence drift: `dt` varies pseudo-randomly in
    /// `[0.05, 0.15]` over 10⁶ steps. Same oracle and bounds as the
    /// constant-`dt` gate.
    #[test]
    fn adaptive_dt_million_step_drift_vs_oracle() {
        let n_steps = 1_000_000_usize;
        let mut state: u64 = 0xDEAD_BEEF_CAFE_BABE;
        let mut dts = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let r = (state as f64) / (u64::MAX as f64);
            dts.push(0.05 + 0.10 * r);
        }

        let oracle = exact_sum_as_f64(&dts);

        let mut t_naive = 0.0_f64;
        let mut t_comp = CompensatedF64::ZERO;
        for &d in &dts {
            t_naive += d;
            t_comp += d;
        }
        let drift_naive = (t_naive - oracle).abs();
        let drift_comp = (t_comp.total() - oracle).abs();
        assert!(drift_naive > 1e-9, "naive drift below regression bound: {drift_naive:.3e}");
        assert!(drift_comp < 1e-9, "compensated drift above ε floor: {drift_comp:.3e}");
        assert!(
            drift_naive / drift_comp.max(f64::MIN_POSITIVE) > 100.0,
            "compensated must beat naive by ≥ 100×",
        );
    }

    /// Order independence under mixed magnitudes: summing 1000 ε-scale
    /// terms before a `1e10` term and the reverse order produce the
    /// same compensated total within ε.
    #[test]
    fn order_invariance_under_mixed_magnitudes() {
        let big = 1e10_f64;
        let small = 1e-6_f64;
        let n_small = 1_000_usize;

        let mut a = CompensatedF64::ZERO;
        for _ in 0..n_small {
            a += small;
        }
        a += big;

        let mut b = CompensatedF64::ZERO;
        b += big;
        for _ in 0..n_small {
            b += small;
        }

        let mut terms = vec![small; n_small];
        terms.push(big);
        let oracle = exact_sum_as_f64(&terms);

        assert!((a.total() - b.total()).abs() < 1e-9);
        assert!((a.total() - oracle).abs() < 1e-9);
        assert!((b.total() - oracle).abs() < 1e-9);
    }

    /// Sanity for the rational oracle: an exact integer sum survives
    /// the `f64 → BigRational → f64` round trip with zero drift.
    #[test]
    fn rational_oracle_is_exact_for_representable_sums() {
        let xs = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(exact_sum_as_f64(&xs), 15.0);

        let xs2 = vec![0.5_f64; 100];
        assert_eq!(exact_sum_as_f64(&xs2), 50.0);
    }
}
