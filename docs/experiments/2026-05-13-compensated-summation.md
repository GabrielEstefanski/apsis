# Compensated summation drift vs arbitrary-precision oracle

**Date:** 2026-05-13
**Subject:** `apsis::math::CompensatedF64` (Neumaier-compensated `f64` accumulator) drift over 10⁶ additions, measured against a `BigRational` oracle that converts each `f64` term losslessly and projects the analytic sum back to the nearest representable `f64`.
**Status:** Protocol declared a priori; numbers measured on the gate machine.
**Branch:** `feat/compensated-f64`

---

## Setup

Two scenarios sum 10⁶ values into a naive `f64` accumulator and a `CompensatedF64` accumulator. The oracle is computed in arbitrary precision (`num_rational::BigRational`), then projected to the nearest `f64`. The oracle does not use any algorithm under test.

| Scenario | `dt` |
|---|---|
| Constant cadence | `0.1` (representation `0.1 + 1.4 × 10⁻¹⁷`) |
| Adaptive cadence | `0.05 + 0.10 · r`, `r` from a 64-bit xorshift, range `[0.05, 0.15)` |

Both summers consume the same `dt` sequence in the same order.

## Results

| Scenario | Oracle | Naive `f64` drift | Compensated drift |
| ---: | ---: | ---: | ---: |
| Constant `dt = 0.1` | 100000.000000 | 1.33 × 10⁻⁶ | 0.0 |
| Adaptive `dt ∈ [0.05, 0.15)` | 100035.273721 | 1.28 × 10⁻⁹ | 0.0 |

Naive `f64` exhibits the expected scaling: worst-case `O(N · ε_machine · |total|)` for the constant-cadence case, random-walk `O(√N · ε_machine · |total|)` for the adaptive case where partial sign cancellation reduces the accumulated error. Both summers achieve the IEEE-754 floor for the analytic sum at 10⁶ steps; compensated drift below `1e-9` is the gate threshold (machine epsilon `× |total|` at this scale ≈ 2.2 × 10⁻¹¹, with two orders of slack for platform variance).

## Oracle validation

`tests::rational_oracle_is_exact_for_representable_sums`:

- `Σ {1, 2, 3, 4, 5}` round-trips through `BigRational` to exactly `15.0`
- `Σ 0.5 × 100` round-trips to exactly `50.0`

The oracle is therefore independent of compensated summation: every term is converted to its exact dyadic-rational representation, summed exactly, and only the final projection back to `f64` rounds (one rounding total, vs `N` for naive).

## Scope

- Scalar accumulator only. Vector accumulation (positions in WHFast) is the integrator's concern; PR B handles it.
- No performance numbers in this notebook. Per-step cost lives in PR B's bench, where `CompensatedF64` runs alongside Kepler / kick / drift work.

## References

Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen. *ZAMM* 54, 39–51.

Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §4.3 ("Compensated summation"). SIAM.

Rein, H., & Tamayo, D. (2015). WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations. *MNRAS* 452, 376–388. § 6 documents REBOUND's compensated summation use.
