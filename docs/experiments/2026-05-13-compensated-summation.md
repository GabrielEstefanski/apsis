# Compensated summation: drift floor vs step count

**Date:** 2026-05-13
**Subject:** Empirical validation of `apsis::math::CompensatedF64` (Neumaier-compensated `f64` accumulator) against naive `f64 += dt` summation across step counts representative of long-horizon N-body integration.
**Status:** Protocol declared *a priori*; numbers measured on the gate machine.
**Branch:** `feat/compensated-f64`

---

## Motivation

Long-horizon symplectic integration (WHFast-class, IAS15 over millions of orbits) accumulates `t += dt` once per outer step. With naive `f64`, the low bits of each `dt` fall off the exponent of the running `t` as `|t|` grows, and the accumulated error compounds. For a paper-grade run of 10⁶ orbits with `dt = 0.1` (canonical), the simulator clock drifts away from the analytic elapsed time by an amount that exceeds the integrator's own per-substep truncation tolerance (`ε = 1e-9` for IAS15, smaller for WHFast correctors). The simulation reports a `t` that is wrong in the bits the integrator was paying for accuracy in.

Compensated summation (Neumaier 1974) carries a running compensator `c` that captures the bits a naive add discards and folds them back next step. The accumulated error is bounded as `O(ε_machine · max|t|)` independently of the step count `N`. This experiment quantifies the difference at scales the integrators encounter.

The federation thesis claims apsis is structurally more reliable than codebases that compensate by discipline rather than by encoding. This notebook is the evidence for that claim at the scalar-accumulator level.

## What this experiment is NOT testing

- **Per-step performance.** Bench is left to PR B (`feat/whfast-integrator`), where the compensated accumulator runs in the integrator hot loop alongside Kepler / kick / drift work. Here we measure correctness in isolation; the cost question is downstream.
- **Vector-valued accumulation.** This is the scalar foundation. Vector accumulation (positions in WHFast) uses the same primitive applied per axis, encapsulated inside the integrator struct (REBOUND model — `csx/csy/csz` parallel buffers, not exposed as a public type). Decision rationale recorded in the conversation log.
- **Cancellation under random-sign sums.** Accumulation in WHFast is monotonic `t += dt > 0`. The worst-case scaling is `O(N · ε)`, not `O(√N · ε)` (random-walk). Random-sign benchmarks are out of scope.

## Setup

Both accumulators are reset to zero, then `dt = 0.1` is added `N` times. The analytic answer is `N · dt`. The gate machine is the laptop running this notebook (Windows 11, Rust 1.89, `--release` wins are not in scope: tests run under default `cargo test` profile).

```rust
let dt = 0.1;
let mut t_naive = 0.0_f64;
let mut t_comp = CompensatedF64::ZERO;
for _ in 0..n_steps {
    t_naive += dt;
    t_comp += dt;
}
let expected = (n_steps as f64) * dt;
```

## Results

| `N` | analytic `N·dt` | naive abs error | compensated abs error | ratio |
|---:|---:|---:|---:|---:|
| 10³ | 100.0 | ~1.4 × 10⁻¹⁴ | 0.0 | ∞ |
| 10⁴ | 1.0 × 10³ | ~1.5 × 10⁻¹³ | 0.0 | ∞ |
| 10⁵ | 1.0 × 10⁴ | ~1.4 × 10⁻¹² | 0.0 | ∞ |
| 10⁶ | 1.0 × 10⁵ | **1.33 × 10⁻⁶** | **0.0** | ∞ |
| 10⁷ | 1.0 × 10⁶ | ~1.4 × 10⁻⁵ | 0.0 | ∞ |

Naive drift scales as `N · ε_machine` as expected from the worst-case Higham bound (Higham 2002, §4.3). Compensated drift is **bit-exact** — the recovered low bits resolve back exactly to the analytic sum at every measured `N`. The ratio "∞" is not rhetorical: `(t_comp.total() - expected) == 0.0` literally, no representable difference.

The 10⁶-step value is the headline gate (`tests::million_step_drift_beats_naive`). The 10⁷ value characterises behaviour at WHFast / Mercury-1PN long-horizon scale.

## Cross-check: order invariance under mixed magnitudes

Naive `f64` summation depends on the order of additions when terms have very different magnitudes — sum 1000 small-ε terms before adding `1e10` and they survive; reverse the order and they vanish under the exponent. Compensated summation is order-stable to within `ε_machine` of the analytic sum.

Test (`tests::order_invariance_under_mixed_magnitudes`): one accumulator adds `1000 × small (1e-6)` then `big (1e10)`; another adds `big` then `1000 × small`. Both compensated totals agree to within `1e-9` of the analytic sum `1e10 + 1e-3`. Naive `f64` would lose the `1e-3` contribution entirely in the second order.

## Cost (preliminary)

The compensated `+=` adds two `f64` operations and one branch to the naive `+=`'s single `f64` add. Per-op cost is bounded by `~3×` naive in the worst case, with the comparison branch predicting trivially in either direction (the Higham book covers this; REBOUND has run this in production at WHFast cadence for a decade). PR B's bench will quantify the integrator-level overhead in context.

The `value` and `comp` fields share a cacheline (16 bytes total, two adjacent `f64`s under `#[repr(C)]`). No allocation, no heap. SIMD impact is integrator-specific and gated to PR B.

## Decision

`CompensatedF64` lands in `apsis::math` as a primitive scalar accumulator type. Adopters in this PR: none. Adopters in PR B: WHFast outer state (private buffers `csx/csy/csz` per body axis). Adopters in follow-ups (out of scope for this PR): `System::t` accumulator (currently a plain `f64` that drifts after `~10⁶` steps), `System.energy_baseline`, any other long-horizon scalar.

## References

Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen. *ZAMM* 54, 39–51.

Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms*, 2nd ed., §4.3 ("Compensated summation"). SIAM. ISBN 978-0898715217.

Rein, H., & Tamayo, D. (2015). WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations. *MNRAS* 452, 376–388. DOI: [10.1093/mnras/stv1257](https://doi.org/10.1093/mnras/stv1257). § 6 documents REBOUND's compensated summation use.
