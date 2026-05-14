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

    // ── Oracle helpers ──────────────────────────────────────────────────────

    /// Convert an `f64` to its exact `BigRational` (`m · 2^e`).
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

    fn rational_to_f64(r: &BigRational) -> f64 {
        r.to_f64().expect("rational projects to representable f64")
    }

    // ── Drift gates ─────────────────────────────────────────────────────────
    //
    // The oracle is the exact rational sum of the *intent values* the dt
    // sequence is meant to represent (e.g. 1/10, not the f64 closest to
    // 1/10). Comparing summer outputs against the intent measures both
    // representation error (irreducible in f64) and accumulation error
    // (what compensated removes).

    /// Constant `dt = 1/10` over 10⁶ steps. Intent oracle:
    /// `BigRational::new(N, 10)`. Compensated drift bound is one ULP
    /// of the total; naive drift grows worst-case as `O(N · ε · |total|)`.
    ///
    /// Note on the compensated bound: `f64` guarantees `|f64(x) − x| <
    /// 0.5 ULP(x)`, so the per-term representation error is bounded
    /// below the per-term ULP. Compensated reaches the nearest `f64`
    /// to the analytic sum; when that distance is sub-ULP it rounds
    /// to the same `f64` as the projected oracle and the measured
    /// drift is `0.0`. This is the floor in `f64`, not a missing test.
    #[test]
    fn constant_dt_drift_vs_intent_oracle() {
        let n_steps = 1_000_000_usize;
        let dt_f64 = 0.1_f64;
        let dt_intent = BigRational::new(BigInt::from(1), BigInt::from(10));
        let oracle_rat: BigRational = dt_intent * BigRational::from_integer(BigInt::from(n_steps));
        let oracle = rational_to_f64(&oracle_rat);

        let mut t_naive = 0.0_f64;
        let mut t_comp = CompensatedF64::ZERO;
        for _ in 0..n_steps {
            t_naive += dt_f64;
            t_comp += dt_f64;
        }
        let drift_naive = (t_naive - oracle).abs();
        let drift_comp = (t_comp.total() - oracle).abs();
        let ulp_oracle = oracle.abs() * f64::EPSILON;

        assert!(
            drift_naive > 1e-7,
            "naive should drift well above representation floor; got {drift_naive:.3e}"
        );
        assert!(
            drift_comp <= ulp_oracle,
            "compensated drift exceeds 1 ULP of oracle; got {drift_comp:.3e} > {ulp_oracle:.3e}"
        );
        assert!(
            drift_naive > 100.0 * ulp_oracle,
            "naive must drift by ≥ 100 ULPs; got {drift_naive:.3e} vs ULP {ulp_oracle:.3e}"
        );
    }

    /// Adaptive cadence: `dt = 1/20 + (state / u64::MAX) · 1/10`,
    /// `state` from a 64-bit xorshift, over 10⁶ steps. Snapshots every
    /// `N / 10` steps capture the drift trajectory; both the final
    /// drift and the maximum drift across snapshots must clear the
    /// gates. This rules out a clean cancellation at the endpoint
    /// hiding a large excursion mid-run.
    #[test]
    fn adaptive_dt_drift_vs_intent_oracle_with_snapshots() {
        const N_STEPS: usize = 1_000_000;
        const N_SNAPSHOTS: usize = 10;

        // Generate dts in two parallel sequences:
        //   dts_f64    — what the summers consume
        //   dts_intent — what the oracle sums exactly
        let mut state: u64 = 0xDEAD_BEEF_CAFE_BABE;
        let mut dts_f64 = Vec::with_capacity(N_STEPS);
        let mut dts_intent: Vec<BigRational> = Vec::with_capacity(N_STEPS);
        let one_twentieth = BigRational::new(BigInt::from(1), BigInt::from(20));
        let one_tenth = BigRational::new(BigInt::from(1), BigInt::from(10));
        let u64_max_rat = BigRational::from_integer(BigInt::from(u64::MAX));
        for _ in 0..N_STEPS {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let r_f64 = (state as f64) / (u64::MAX as f64);
            dts_f64.push(0.05 + 0.10 * r_f64);
            let r_rat = BigRational::from_integer(BigInt::from(state)) / u64_max_rat.clone();
            dts_intent.push(one_twentieth.clone() + one_tenth.clone() * r_rat);
        }

        let snapshot_stride = N_STEPS / N_SNAPSHOTS;
        let mut oracle_running = BigRational::zero();
        let mut t_naive = 0.0_f64;
        let mut t_comp = CompensatedF64::ZERO;
        let mut max_drift_naive = 0.0_f64;
        let mut max_drift_comp = 0.0_f64;

        for i in 0..N_STEPS {
            t_naive += dts_f64[i];
            t_comp += dts_f64[i];
            oracle_running += &dts_intent[i];

            // Snapshot at every stride boundary (and the final step).
            if (i + 1) % snapshot_stride == 0 || i == N_STEPS - 1 {
                let oracle_now = rational_to_f64(&oracle_running);
                let d_n = (t_naive - oracle_now).abs();
                let d_c = (t_comp.total() - oracle_now).abs();
                if d_n > max_drift_naive {
                    max_drift_naive = d_n;
                }
                if d_c > max_drift_comp {
                    max_drift_comp = d_c;
                }
            }
        }

        let oracle_final = rational_to_f64(&oracle_running);
        let drift_final_naive = (t_naive - oracle_final).abs();
        let drift_final_comp = (t_comp.total() - oracle_final).abs();
        let ulp_final = oracle_final.abs() * f64::EPSILON;

        // Random-sign cancellation in adaptive cadence reduces naive
        // drift to ~√N · ε · |total| ≈ 30 ULPs at N=10⁶. Bound 10
        // ULPs is the floor that proves divergence above pure round-off
        // without overclaiming the random-walk envelope.
        assert!(
            drift_final_naive > 10.0 * ulp_final,
            "naive should drift past 10 ULPs at endpoint; got {drift_final_naive:.3e} vs ULP {ulp_final:.3e}"
        );
        assert!(
            drift_final_comp <= ulp_final,
            "compensated final drift exceeds 1 ULP of oracle; got {drift_final_comp:.3e} > {ulp_final:.3e}"
        );

        assert!(
            max_drift_comp <= 2.0 * ulp_final,
            "compensated max-snapshot drift exceeds 2 ULPs; got {max_drift_comp:.3e}"
        );
        assert!(
            max_drift_naive > 10.0 * max_drift_comp.max(ulp_final),
            "naive max-snapshot must beat compensated by ≥ 10× of max(comp, ULP); \
             got naive {max_drift_naive:.3e} vs comp {max_drift_comp:.3e}"
        );
    }

    /// Compensated drift stays bounded by 1 ULP of the running total
    /// for every `N ∈ {10⁴, 10⁵, 10⁶}`, confirming the bound is not
    /// an artifact of one specific step count.
    #[test]
    fn compensated_within_one_ulp_across_step_counts() {
        let dt_f64 = 0.1_f64;
        let dt_intent = BigRational::new(BigInt::from(1), BigInt::from(10));

        for &n_steps in &[10_000_usize, 100_000, 1_000_000] {
            let oracle = rational_to_f64(
                &(dt_intent.clone() * BigRational::from_integer(BigInt::from(n_steps))),
            );
            let mut t_comp = CompensatedF64::ZERO;
            for _ in 0..n_steps {
                t_comp += dt_f64;
            }
            let d_c = (t_comp.total() - oracle).abs();
            let ulp_n = oracle.abs() * f64::EPSILON;
            assert!(
                d_c <= ulp_n,
                "compensated drift exceeds 1 ULP at N={n_steps}: {d_c:.3e} > {ulp_n:.3e}"
            );
        }
    }

    /// Order independence under mixed magnitudes. Two compensated
    /// summers consume identical terms in opposite orders; both must
    /// agree with the oracle within ε.
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

        let oracle_rat = f64_to_rational(big)
            + BigRational::from_integer(BigInt::from(n_small as i64)) * f64_to_rational(small);
        let oracle = rational_to_f64(&oracle_rat);

        assert!((a.total() - b.total()).abs() < 1e-9);
        assert!((a.total() - oracle).abs() < 1e-9);
        assert!((b.total() - oracle).abs() < 1e-9);
    }

    /// Sanity for the rational oracle: integer sums round-trip through
    /// `BigRational` to the exact expected value.
    #[test]
    fn rational_oracle_is_exact_for_representable_sums() {
        let xs = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let mut acc = BigRational::zero();
        for &x in &xs {
            acc += f64_to_rational(x);
        }
        assert_eq!(rational_to_f64(&acc), 15.0);

        let xs2 = vec![0.5_f64; 100];
        let mut acc2 = BigRational::zero();
        for &x in &xs2 {
            acc2 += f64_to_rational(x);
        }
        assert_eq!(rational_to_f64(&acc2), 50.0);
    }
}
