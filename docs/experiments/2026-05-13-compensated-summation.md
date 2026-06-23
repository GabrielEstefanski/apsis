# Compensated summation drift in ULPs of total

**Date:** 2026-05-13
**Subject:** `apsis::math::CompensatedF64` drift over 10⁶ additions, measured in ULPs of the running total against an arbitrary-precision intent oracle.

**Status:** Protocol declared a priori; numbers measured on the gate machine.

**Branch:** `feat/compensated-f64`

---

## Setup

Two scenarios sum 10⁶ values into a naive `f64` accumulator and a `CompensatedF64` accumulator. The oracle is the **exact rational sum of the intent values the cadence is meant to represent** — `1/10` for the constant case, `1/20 + (state / u64::MAX) · 1/10` for the adaptive case — projected to the nearest representable `f64`.

Comparing summer outputs to the intent oracle measures both representation error (irreducible in `f64`: each per-term rounding is bounded by `0.5 · ULP(term)` per IEEE-754) and accumulation error (what compensated removes).

| Scenario | `dt_f64` consumed | `dt_intent` summed by oracle |
| --- | --- | --- |
| Constant cadence | `0.1` | `BigRational(1, 10)` |
| Adaptive cadence | `0.05 + 0.10 · r`, `r = state / u64::MAX_f64`, `state` from xorshift64 | `BigRational(1, 20) + BigRational(1, 10) · BigRational(state, u64::MAX)` |

## Results

Drift is reported absolute and as **ULPs of the oracle total** (`ε_machine × |oracle|`). The ULP unit is the natural floor for `f64` accuracy — anything below 1 ULP is indistinguishable from the projected oracle in `f64`.

| Scenario | Oracle | Naive drift | Naive ULPs | Compensated drift | Compensated ULPs |
| ---: | ---: | ---: | ---: | ---: | ---: |
| Constant `dt = 1/10` | 100000.000000 | 1.33 × 10⁻⁶ | 60028 | 0 | 0 (sub-ULP) |
| Adaptive `dt ∈ [1/20, 3/20)` | 100035.273721 | 1.27 × 10⁻⁹ | 57 | 1.46 × 10⁻¹¹ | 0.65 |

### What the constant-cadence result means

`f64(0.1)` differs from `1/10` by `~5.55 × 10⁻¹⁸`, well under `0.5 · ULP(0.1) ≈ 1.4 × 10⁻¹⁷`. Accumulated monotonically over 10⁶ steps the representation error reaches `~5.5 × 10⁻¹²`, still below the 1-ULP envelope of the total (`~2.2 × 10⁻¹¹`). The compensated sum reaches the nearest `f64` to the analytic intent (`100000.0` exactly representable), so the measured drift rounds to `0.0`.

This is the IEEE-754 floor, not a missing test. Naive summation crosses 60000 ULPs at the same `N` because each per-step add additionally drops bits the compensator captures.

### What the adaptive-cadence result means

Random-sign variation in `dt` introduces partial cancellation in both summers. Naive drift settles at 57 ULPs (random-walk scaling `O(√N · ε · |total|)`); compensated drift becomes measurable at `0.65 ULP` — non-zero, well within the per-snapshot bound.

Snapshots every `N / 10` steps confirm the compensated trajectory does not exceed `2.052 × 10⁻¹¹ ≈ 0.92 ULP` at any intermediate point, ruling out a clean cancellation hiding a large mid-run excursion.

## Scaling

`compensated_within_one_ulp_across_step_counts` (constant `dt`):

| N | Compensated drift (ULPs of N · 0.1) |
| ---: | ---: |
| 10⁴ | ≤ 1 |
| 10⁵ | ≤ 1 |
| 10⁶ | ≤ 1 |

## Oracle validation

`rational_oracle_is_exact_for_representable_sums`: `Σ {1, 2, 3, 4, 5}` round-trips through `BigRational` to exactly `15.0`; `Σ 0.5 × 100` round-trips to exactly `50.0`. The oracle uses no compensated-summation primitives, so the three-way comparison (naive ↔ compensated ↔ rational oracle) is non-circular.

## Scope

- Scalar accumulator. Vector accumulation (positions in WHFast) is the integrator's concern; PR B handles it.
- No performance numbers in this notebook. Per-step cost lives in PR B's bench, where `CompensatedF64` runs alongside Kepler / kick / drift work.

## What this experiment does NOT claim

- Compensated drift is universally `0.0`. The `0.0` measurement at constant cadence is the IEEE-754 floor for the specific `(N, dt)` chosen; in adaptive cadence the same algorithm produces a measurable `0.65 ULP` drift.
- Compensated dominates naive in every scenario. For sums where every term is identically `f64(0.0)` or where naive happens to land at the same `f64` by coincidence, the two agree. The generic claim is bounded: compensated drift is `O(ULP(|total|))`; naive grows worst-case `O(N · ε · |total|)`.

## References

Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen. *ZAMM* 54, 39–51.

Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §4.3 ("Compensated summation"). SIAM.

Rein, H., & Tamayo, D. (2015). WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations. *MNRAS* 452, 376–388. §2.7 establishes the unbiased `O(√N · ε)` round-off envelope (Brouwer's law).
