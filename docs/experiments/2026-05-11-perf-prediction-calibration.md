# Perf prediction calibration — lessons from the MAC and SoA experiments

**Date:** 2026-05-11
**Subject:** Two consecutive perf predictions on the gravity hot path (MAC comparison, SoA layout) overestimated their gains. This doc records the pattern, identifies the misframing root cause, and calibrates the next perf prediction (SIMD) so the same overestimation does not recur.

**Status:** Standing reference. Future perf experiments on the gravity hot path should consult §Calibration rule before declaring a-priori bounds.

---

## Context

The engine ceiling experiment (`2026-05-09-engine-ceiling.md` §Decision, lines 348-354) classified the production engine on the recorded hardware (Ryzen 5 7600X, AVX-512 available) at the v1 target N (10³-10⁵):

> "*Compute-bound* would mean `t_per_interaction` is dominated by FLOPs/cycle and would respond linearly to SIMD throughput. Observed `t_per_interaction = 1.3-2.1 ns` is already low — the per-interaction kernel is well-vectorised by the compiler, and explicit SIMD will yield **engineering-baseline gains, not order-of-magnitude wins**.
>
> *Memory-bound* would show `t_per_interaction` growing materially with N as the working set spills DRAM. Observed variation is 1.24× across two orders of magnitude — mild cache pressure, not the dominant signal.
>
> ***Interaction-bound* is the right framing**: the dominant cost is **how many interactions the walk performs**, not the cost of each."

The roadmap predicted multiplicative ceiling `MAC × SIMD × SoA = 0.7 × 0.5 × 0.9 ≈ 0.32×` (3× speedup) "in the most optimistic scenario where the factors compose without erosion. Realistic central estimate: 2.5-3.5× total."

Two of the three axes have now been measured. Both predictions overestimated.

---

## The two predictions

### MAC

| Quantity | A-priori range | Measured | Status |
| --- | ---: | ---: | --- |
| `t_walk_M1 / t_walk_M0` at N = 10⁴ | ∈ [0.75, 0.95] | **2.73** | FAILED, 2.9× outside upper bound |
| `n_interactions_M1 / n_interactions_M0` at N = 10⁴ | < 1 (some reduction) | **2.71** | wrong direction |

§Decision: defer M1, do not implement M2/M3. Triangle-inequality `δ_max` aggregation accumulates slack at every recursion level; for our regime the geometric MAC tightening was net-negative on interaction count.

### SoA layout

| Quantity | A-priori range | Measured | Status |
| --- | ---: | ---: | --- |
| `t_walk_AoS / t_walk_SoA` at N = 10⁴ | ∈ [1.20, 1.50] | **1.015** | FAILED, below lower bound by 1.18× |
| Tier 1 bit-exact accelerations | == 0 ULP | == 0 ULP | PASS |
| Tier 3 pack overhead | ≤ 1 % of compute | 0.45 % | PASS |

§Decision: ship SoA on structural grounds (Tier 1 + Tier 3 + SIMD pre-requisite), not on Tier 2 measurement. The cache-line waste arithmetic that motivated SoA was correct in isolation but ignored three real effects.

---

## The pattern: three causes for over-prediction

Both predictions assumed effects that the data shows are not load-bearing on the recorded hardware at the v1 target N:

### 1. The compute-bound vs interaction-bound misframing

Both predictions were built on a "compute-bound" or "memory-bound" framing of the regime — load-locality (SoA) or per-interaction efficiency (MAC) would translate directly to wall-time. Engine ceiling §Decision had already classified the system as **interaction-bound** before either prediction was made. In an interaction-bound regime, the dominant lever is **work volume reduction** (fewer interactions), not **per-interaction cost reduction** (faster kernel) or **per-load efficiency** (better layout).

MAC was the volume-reduction axis — but the chosen MAC (Barnes 1990 with triangle-inequality `δ_max`) **increased** interaction count instead of reducing it. SoA was a per-load-efficiency axis — its theoretical 3.2× cache-line waste reduction did not materialise because the actual bottleneck is compute throughput per interaction (already low), not load throughput.

### 2. HW prefetcher coverage

Modern x86 hardware prefetchers (Zen 4 L2 on Ryzen 5 7600X) speculatively load adjacent cache lines on detected sequential access. The SoA prediction's "AoS wastes ~70% of each cache line load" arithmetic assumed the wasted bytes are loaded cold, then discarded. The prefetcher's behaviour changes that accounting: in sequential leaf iteration, the "wasted" bytes are pre-loaded and sit in L1 for the next iteration. The waste is real on paper, prefetch-hidden in practice.

### 3. Plummer kernel is cheap

The per-interaction Plummer kernel is ~30-50 cycles scalar (`sqrt`, `recip`, three multiplies, three adds). Engine ceiling measured `t_per_interaction = 1.3-2.1 ns` ≈ 5-9 cycles on Zen 4 — meaning the compiler is already getting 5-10× speedup over the naive scalar count via auto-vectorisation, FMA, and pipelining. Explicit SIMD has limited remaining headroom because the compiler already covers most of the gain.

---

## Calibration rule for future perf predictions

For any perf experiment on the gravity hot path (BH walk, force kernel, octree build), the a-priori bound construction must:

1. **Treat engine ceiling §Decision bounds as ceilings, not floors.** If engine ceiling predicted `0.5×` per-interaction speedup from SIMD (1.5-2.5× speedup, 3× theoretical ceiling), the next experiment's a-priori range must respect that ceiling. Predicting above the engine-ceiling ceiling without new evidence repeats the over-prediction.

2. **Discount predictions of cache-locality gains.** The two axes that targeted load efficiency (SoA, AoSoA, Morton) all measured below their predicted gains. Future cache-locality predictions must explicitly account for: (a) HW prefetcher coverage of sequential access patterns, (b) compute-bound regime where load improvements don't translate to wall-time gain, (c) per-interaction kernel arithmetic dominating the load.

3. **Discount predictions of per-interaction cost reductions in scalar regimes.** SIMD predictions in particular must respect that the compiler is already getting auto-vectorisation gains; the ceiling for explicit SIMD is the *remaining* headroom, not the lane-count ratio.

4. **Recompute multiplicative ceilings as axes are measured.** Engine ceiling predicted `MAC × SIMD × SoA = 0.32×` (3× speedup). With `MAC = 1.0×` (deferred) and `SoA = 1.0×` (no Tier 2 gain), the realistic ceiling now is just `SIMD ≈ 0.5-0.7×` (1.4-2× speedup), not 3×. Predictions that re-imply the original 3× ceiling without acknowledging the cascade are wrong.

5. **Honest framing in the §Motivation.** Each new perf experiment's notebook should reference this doc by name and include a short paragraph stating which calibration rules apply. Reviewers should be able to verify the prediction was not quietly inherited from an obsolete frame.

---

## Revert discipline (joint criterion)

The accumulated complexity from the perf series is non-trivial:

- `BodyArrays` SoA snapshot type — ~250 LOC including pack helpers and tests
- `Body` field rename `pos_x/y/z, vel_x/y/z` — touches ~36 files mechanically
- (Pending SIMD) AoSoA chunked layout, two-phase walk, scalar/AVX2/AVX-512 kernel paths, runtime dispatch — projected ~1000-1500 LOC

If the SIMD measurement also misses its (recalibrated) bound — Tier 3 walk speedup at N = 10⁴ ≤ 1.0× — the cumulative complexity is unjustified by measurement. The joint revert criterion must fire: revert SoA and SIMD together, document as deferred. The MAC §Decision template applies — negative result with documented mechanism, not engineering failure.

This rule is not symmetric. If SIMD lands gain inside its (recalibrated) range, both ship as planned. The asymmetry is intentional: shipping requires positive evidence, reverting requires accumulated negative evidence + cost-of-complexity argument.

---

## Reference

- Engine ceiling §Decision: `docs/experiments/2026-05-09-engine-ceiling.md` lines 348-394 (interaction-bound classification + 4-axis roadmap with multiplicative ceiling)
- MAC §Decision: `docs/experiments/2026-05-09-octree-mac.md` §Decision (defer M1)
- SoA §Decision: `docs/experiments/2026-05-10-soa-layout.md` §Decision (ship on structural grounds)
- Memory `feedback_no_tuning_to_pass.md`: bound revision is forbidden unless backed by concrete arithmetic
- Memory `project_codebase_neighborhood.md`: apsis is REBOUND-class (small-N high-precision), not GADGET-class — perf predictions calibrated against GADGET-class regime expectations are incompatible
