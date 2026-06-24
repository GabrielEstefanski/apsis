# Engine ceiling re-measurement after SIMD lands

**Date:** 2026-05-12
**Subject:** Re-run the engine-ceiling profiling harness on the production code path with the AVX2 leaf-pair SIMD kernel active. Compare to the pre-SIMD baseline in `docs/experiments/2026-05-09-engine-ceiling.md`. The objective is decomposition, not a ship/revert gate: validate the back-out estimate that leaf-pair was ~67-80 % of pre-SIMD walk time, and quantify how much per-interaction headroom remains for further kernel-arithmetic optimisation (RSQRT-class kernel approximation) vs further vectorisation (accepted-node phase).

**Status:** Diagnostic measurement, not an axis experiment. No production code changes. The result feeds the choice between the next perf optimisation axes — kernel approximation (rsqrt + Newton-Raphson, simultaneously affects leaf and accepted-node phases), rayon scheduling tune, accepted-node SIMD, or AoSoA + Morton.

---

## Motivation

The SIMD lab notebook §Interpretation derived a model where the BH walk is decomposed into

```text
walk = f_dispatch + f_emit + f_kernel_leaf + f_kernel_node
```

with pre-SIMD fractions estimated as `f_dispatch ≈ 0.13`, `f_emit ≈ 0.08`, `f_kernel_leaf ≈ 0.35-0.40`, `f_kernel_node ≈ 0.40-0.45`. That estimate was used to derive the AVX2 walk-speedup envelope `[1.3, 2.0]×` from a kernel-isolated speedup of `S_kernel ∈ [1.8, 2.5]×`.

**The estimate was wrong.** Backing-out from the measured AVX2 walk speedups (1.5× median on Zen 4; 1.67× median on Sapphire Rapids), holding `S_kernel = 2.0×` and assuming SIMD only touches the leaf-pair phase:

```text
walk_post / walk_pre = 1 / S_walk
S_walk = 1.5  →  leaf_pre = 0.67   (Zen 4)
S_walk = 1.67 →  leaf_pre = 0.80   (Sapphire Rapids)
```

Leaf-pair was 67-80 % of pre-SIMD walk time, not 35-40 %. That changes the Amdahl analysis materially:

- Non-leaf fraction (accepted-node + dispatch + emit) is 20-33 %, not 60 %.
- Kernel-only Amdahl ceiling is `1 / 0.20 = 5.0×` (Sapphire Rapids back-out) or `1 / 0.33 = 3.0×` (Zen 4 back-out) walk speedup — substantially more headroom than the original §Tier 3 derivation predicted.
- Current AVX2 sits at 1.5-1.67× — between 33 % and 50 % of the kernel-only ceiling.

This experiment validates the back-out empirically by re-measuring engine-ceiling Cell V with SIMD active and reading off:

1. `t_per_interaction` change vs the May 9 baseline (1.3-2.1 ns).
2. `t_per_body` change at matched N.
3. The relative time inside `t_bh_walk` after SIMD.

The walk-counter ratios (`n_node_visits`, `n_bh_accepted`, `n_leaf_interactions`) are algorithm-only and should be **identical** to the May 9 numbers — they are unaffected by SIMD. Any deviation is a regression in the BH algorithm, not a SIMD effect.

## A-priori predictions

| Quantity | May 9 baseline | Predicted post-SIMD | Rationale |
| --- | ---: | ---: | --- |
| `t_per_interaction` at N = 10⁴ | 1.5 ns | ∈ [0.8, 1.2] ns | AVX2 leaf-pair kernel reduces per-interaction cost on the leaf phase only; weighted across all interactions (leaf + accepted-node), expected drop is ~30-45 %. |
| `t_per_body` at N = 10⁴ | 4.05 µs | ∈ [2.4, 3.1] µs | Same logic propagated to per-body; consistent with measured 1.5× walk speedup. |
| `n_node_visits / n_bh_accepted / n_leaf_interactions` ratios | (see baseline) | unchanged | Algorithm unchanged; counters depend on tree structure + θ + body distribution only. |
| `bh_acceptance_ratio` at N = 10⁴ | 0.29 | unchanged (0.29) | Same reason. |
| Build / walk wall-time fraction at N = 10⁴ | walk dominates ~96 % | walk dominates ~94 % | Build cost flat; walk dropped ~33 %; build's fraction grows slightly. |
| Trail-on / trail-off ratio | (see baseline) | similar | Trail cost is a fixed fraction of step time; should track walk savings only weakly. |

The bound is not a hard gate — this is a diagnostic, not a ship/revert axis. The numbers either confirm the back-out (and inform the next axis choice) or contradict it (in which case the §Interpretation analysis needs revision).

## Methodology

Re-run the existing `engine_ceiling_v` harness on the recorded Zen 4 desktop with the production code at HEAD. Same parameters as May 9:

- Seed: `0x6E63696C`
- Body distribution: sphere log-normal mass
- VV dt: `1e-3`
- θ: `0.5`
- Warmup: 10 steps
- Measured: 100 steps
- N grid: `[100, 1_000, 5_000, 10_000, 50_000, 100_000]`
- Trail variants: `[off, on]`

Output schema unchanged from `crates/apsis/src/physics/engine_ceiling.rs`. Compare side-by-side with the May 9 numbers from `2026-05-09-engine-ceiling.md` §Tier 1 and §Tier 2.

## Run

```text
cargo test --release -p apsis engine_ceiling_v -- --ignored --nocapture
```

CSV output goes to `target/engine-ceiling/profile_v.csv`.

---

## Results

Run completed `2026-05-12` on the recorded Zen 4 desktop, seed `0x6E63696C`, 10 warmup + 100 measured steps per cell. CSV at `target/engine-ceiling/profile_v.csv`.

### Tier 1 — Wall-time + phase decomposition (Cell V, trail-off)

| N | step (ms) | SPS | build % | walk % | integ % | t_per_body (µs) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.088 | 11 302 | 13.3 | 86.1 | 0.4 | 0.88 |
| 1 000 | 1.376 | 727 | 10.6 | 89.0 | 0.3 | 1.38 |
| 5 000 | 19.617 | 51 | 4.4 | 95.3 | 0.3 | 3.92 |
| 10 000 | 35.280 | 28 | 4.2 | 95.4 | 0.4 | 3.53 |
| 50 000 | 243.728 | 4.10 | 4.0 | 95.3 | 0.7 | 4.87 |
| 100 000 | 626.338 | 1.60 | 4.8 | 94.7 | 0.6 | 6.26 |

Walk continues to dominate ≥ 86 % of step time at every N; build fraction stayed within 0.5 percentage points of the May 9 baseline (≈ 3-13 %). Trail-on cells (omitted for brevity, in CSV) show the trail phase at 0.3-1.1 % of step — same as May 9.

### Tier 2 — Cost normalised by work (Cell V) — comparison vs May 9 baseline

| N | t_per_int May 9 | t_per_int May 12 | Δ | t_per_body May 9 | t_per_body May 12 | Δ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 2.8 ns | 4.0 ns | **+43 %** | 0.64 µs | 0.88 µs | +37 % |
| 1 000 | 1.7 ns | 1.1 ns | **−35 %** | 2.06 µs | 1.38 µs | **−33 %** |
| 5 000 | 1.3 ns | 1.0 ns | −23 % | 4.73 µs | 3.92 µs | −17 % |
| 10 000 | 1.5 ns | 1.3 ns | −13 % | 4.05 µs | 3.53 µs | −13 % |
| 50 000 | 1.8 ns | 1.7 ns | −6 % | 4.99 µs | 4.87 µs | −2 % |
| 100 000 | 2.1 ns | 1.9 ns | −10 % | 7.07 µs | 6.26 µs | −11 % |

`t_per_interaction` dropped at every N ≥ 1 000, with the largest relative reduction at N = 1 000 (−35 %) and the smallest at N = 50 000 (−6 %). The N = 100 cell regressed (+43 %) — at that scale absolute step time is sub-millisecond and run-to-run variance dominates the SIMD effect; the trail-on N = 100 row (2.9 ns) is closer to the May 9 baseline than the trail-off row (4.0 ns), confirming small-N noise.

### Tier 3 — Counter ratios (Cell V)

| N | int/body May 9 | int/body May 12 | accept ratio May 9 | accept ratio May 12 |
| ---: | ---: | ---: | ---: | ---: |
| 100 | 189 | 189 | 0.01 | 0.01 |
| 1 000 | 1 110 | 1 110 | 0.05 | 0.05 |
| 5 000 | 3 571 | 3 571 | 0.09 | 0.09 |
| 10 000 | 2 667 | 2 667 | 0.29 | 0.29 |
| 50 000 | 2 724 | 2 724 | 0.41 | 0.41 |
| 100 000 | 3 151 | 3 151 | 0.48 | 0.48 |

All counter ratios match the May 9 baseline exactly. The BH algorithm is unchanged by SIMD (as predicted) — this rules out any "the algorithm shifted under us" confounders for the timing comparison above.

### Tier 4 — Cross-check against `perf_simd_walk` measurements

The previous `perf_simd_walk` harness (since deleted) measured scalar→AVX2 walk speedup directly via A/B on the same seed within a single run. Comparing the two measurement methodologies at matched N:

| N | engine_ceiling May 9→May 12 step ratio | `perf_simd_walk` direct A/B median |
| ---: | ---: | ---: |
| 1 000 | 1.49× | 1.42× |
| 5 000 | 1.21× | 1.29× |
| 10 000 | 1.15× | 1.54× |

At N = 1 000 and N = 5 000 the two methodologies agree within ~5 %. At N = 10 000 they diverge (1.15× vs 1.54×). The seeds differ (`0x6E63696C` vs three perf-canonical seeds), and engine ceiling compares across two runs separated by 3 days while `perf_simd_walk` does in-run A/B. The large-N discrepancy is most plausibly a combination of seed-distribution-dependent leaf-pair fraction and across-day variance in the May 9 baseline. The in-run A/B from `perf_simd_walk` is the cleaner direct measurement of SIMD-as-such.

---

## Interpretation

The §A-priori predictions table verdict:

| Prediction | A-priori | Observed | Verdict |
| --- | --- | --- | --- |
| `t_per_interaction` at N = 10⁴ ∈ [0.8, 1.2] ns | yes | 1.3 ns | **above range** by 0.1 ns |
| `t_per_body` at N = 10⁴ ∈ [2.4, 3.1] µs | yes | 3.53 µs | **above range** by 0.43 µs |
| Counters unchanged | yes | identical | pass |
| Walk dominates ≈ 94 % at N = 10⁴ | yes | 95.4 % | pass |

The two predictions outside their range are both higher than expected — the SIMD walk-speedup at engine-ceiling's seed/distribution at N = 10⁴ is smaller than `perf_simd_walk` measured on the perf-canonical seeds. This is the most important finding from this re-measurement.

### Finding 1 — Leaf-pair fraction is distribution-dependent, not a single number

Backing out `f_leaf` (fraction of pre-SIMD walk in the leaf-pair phase) from each measurement assuming `S_kernel ≈ 2.0×`:

| Source | N | walk speedup | back-out `f_leaf` |
| --- | ---: | ---: | ---: |
| `perf_simd_walk` median (Zen 4 desktop, perf-canonical seeds) | 1 000 | 1.42× | 0.59 |
| `perf_simd_walk` median (Zen 4) | 5 000 | 1.29× | 0.45 |
| `perf_simd_walk` median (Zen 4) | 10 000 | 1.54× | 0.70 |
| `perf_simd_walk` median (Sapphire Rapids, perf-canonical) | 10 000 | 1.67× | 0.80 |
| engine-ceiling delta (Zen 4, seed `0x6E63696C`) | 1 000 | 1.49× | 0.66 |
| engine-ceiling delta (Zen 4, seed `0x6E63696C`) | 5 000 | 1.21× | 0.35 |
| engine-ceiling delta (Zen 4, seed `0x6E63696C`) | 10 000 | 1.15× | 0.26 |

The back-out range across (seed, N) is `f_leaf ∈ [0.26, 0.80]` — a 3× spread. The SIMD lab notebook §Interpretation analysis used a single back-out value (~0.67) and treated it as universal. **It is not.** The dispersion is real and reflects body-distribution-dependent walk topology: distributions where the BH walk traverses more accepted-node interactions (high acceptance ratio at large N, here `0.29-0.48`) have a smaller leaf-pair fraction and therefore see a smaller SIMD walk speedup.

This explains the diminishing-returns curve in the Tier 2 table: at small N (`accept ≈ 0.05-0.09`) the walk is dominated by leaf-pair work and SIMD shows large gains (35 % at N = 1 000); at large N (`accept ≈ 0.41-0.48`) the walk is dominated by accepted-node work which is still scalar, and SIMD's contribution shrinks (2-10 %).

### Finding 2 — Accepted-node phase is the load-bearing residual at scale

The previous SIMD §Interpretation listed three mechanisms eroding the lane-count-ratio prediction. With this finding added: **the accepted-node phase being scalar is the dominant Amdahl floor at large N**, more important than dispatch / emit overhead. At N = 100 000 with `accept = 0.48`, roughly half of the walk visits are accepted-node interactions where the SIMD path is structurally not applied; t_per_interaction at large N is therefore the average of fast leaf (~0.5-0.7 ns SIMD) and slow accepted-node (~3-4 ns scalar quadrupole), weighted by the acceptance ratio.

This changes the relative attractiveness of the deferred next-axis candidates:

- **B (vectorise accepted-node phase)** — leverage **grows** with N. At N = 10⁴ savings would be modest (~17 % of walk); at N = 10⁵ savings would approach ~50 %. For the v0.1 paper target (N ≤ 10³) it stays small; for v0.2 scaling work (N ≤ 10⁵) it becomes the dominant axis.
- **C (AoSoA + Morton)** — leverage **shrunk** post-SIMD because it targets the leaf-pair phase whose fraction is now smaller (especially at large N). At engine-ceiling N = 10⁴ where leaf is only ~26 %, even halving leaf-pair gathers saves <13 % of walk.
- **E (RSQRT + Newton-Raphson)** — leverage **unchanged**. It accelerates per-interaction arithmetic regardless of which phase the interaction is in. This makes E the **only** axis that captures both leaf and accepted-node savings simultaneously.
- **D (rayon chunking)** — leverage **largest at small N**, where rayon scheduling overhead is a fraction of per-task work. The Zen 4 small-N walk speedup variance (1.27-1.81× across seeds at N = 1 000) is consistent with rayon overhead being a meaningful fraction.

### Finding 3 — The N = 100 cell is below the regime SIMD targets

`t_per_interaction` at N = 100 went from 2.8 ns to 4.0 ns (or 2.9 ns trail-on). At this scale the BH walk has overhead from rayon spawn + tree-walk dispatch that exceeds the kernel work; SIMD's per-interaction savings are dwarfed by that fixed overhead. This is below the N regime apsis targets (the smallest validation cell is N ≈ 10² and that path uses exact-eval anyway via the EXACT_THRESHOLD branch — no BH walk for that case in production). The observation is recorded but not actionable — apsis production never enters this regime through the SIMD path.

### Finding 4 — Build phase is now ~5 % of step at large N

Pre-SIMD the build phase was ≤ 4 % of step at N = 10⁴; post-SIMD it is 4.2 %. The walk shrunk so the build proportionally grew, but build's absolute cost (~1.5 ms at N = 10⁴) is unchanged. Build-phase optimisation remains low-priority — would have to halve build cost to gain 2 % of step time.

---

## Decision

This notebook does not have a ship/revert gate. The output informs the next perf axis choice. The four candidates left open by the SIMD §Decision are now ranked by the data above:

| Axis | Pre-experiment estimated leverage | Post-experiment leverage | Recommended priority |
| --- | --- | --- | --- |
| **E** — rsqrt + Newton-Raphson on Plummer kernel arithmetic | "+25-30 % walk speedup, cross-phase" | confirmed cross-phase (only candidate that helps both leaf and accepted-node simultaneously); absolute headroom bounded by per-interaction floor (`~5-7 cycles` at midrange N) | **Next** — small-LOC, low-risk, cross-phase. Bound construction is straightforward (rsqrt+1NR delivers ~24-bit precision; needs Tier 1 against IAS15 conservation). |
| **D** — rayon chunking / scheduler tune | "+5-15 % at small N where overhead pesa" | consistent with the small-N variance observed in `perf_simd_walk` on Zen 4 (1.27-1.81× spread at N = 1 000); not measured directly here | **After E** — separate axis, low risk, modest gain. Worth a 1-day spike with `par_chunks` + thread-pool size tune. |
| **B** — accepted-node phase SIMD (gather Node fields → batched quadrupole) | "+20-30 % at N ≥ 10⁴ if accepted-node is large fraction" | **leverage confirmed and growing with N**: at N = 10⁵ accepted-node is ~half the walk and is currently scalar | **Defer to v0.2 scaling work** — for v0.1 paper (N ≤ 10³) leverage is small; for v0.2 (N ≤ 10⁵) it becomes the dominant axis. |
| **C** — AoSoA + Morton ordering | "+5 % marginal" | **leverage shrunk** because leaf-pair fraction is smaller post-SIMD; engine-ceiling N = 10⁴ has leaf-pair = 26 %, even halving its gathers saves < 13 % of walk | **Drop** — was conditional on leaf-pair being a residual bottleneck; data shows it is not. |

### Recommended next axis: E (rsqrt + Newton-Raphson)

The SIMD §Decision rules require a-priori bound construction before code lands. For E, the bound construction is:

- **Hardware reciprocal-sqrt (`vrsqrtps` / scalar `rsqrtss`)** delivers ~12 bits of precision at single FMA cost. One Newton-Raphson refinement step (3 FMAs) brings it to ~24 bits (single precision). A second Newton step takes it to ~48 bits — close to but short of double-precision (52 bits).
- The Plummer kernel's `(r² + ε²).sqrt().recip()` chain has hardware throughput on Zen 4 of `~12 cyc (sqrt) + 15 cyc (div) = 27 cyc` per 4-lane SIMD chunk. Replacing with `rsqrt + 1 Newton` is `~4 cyc + 3 FMAs ≈ 7 cyc`. Per-interaction kernel critical path shortens by `~20 cyc / 4 lanes = 5 cyc / lane` — that maps to roughly halving the kernel arithmetic cost on the leaf-pair phase, and proportionally on the accepted-node scalar path (which would use the scalar `rsqrtss`).
- Predicted Tier 2a kernel speedup of E over current AVX2: `[1.4, 1.8]×`. Predicted walk speedup over current AVX2: `[1.2, 1.5]×` (smaller because non-vectorisable phases stay constant).
- Tier 1 correctness gate must be tighter than the AVX2 Tier 1's `1e-13` because rsqrt+1NR introduces a systematic ~`5e-8` relative error in `1/r`. The proposed gate: per-body `||a_e[i] - a_avx2[i]|| / ||a_avx2[i]|| ≤ 1e-7` AND IAS15 conservation gate `Δ|L|/|L₀| ≤ 1e-7` over 100 Kepler periods on a 2-body inclined orbit. If the kernel-error is `~5e-8` per interaction and accumulates over `~3000` interactions per body per step, per-step acceleration error is `~1.5e-4` relative — within IAS15's adaptive-dt error budget but at the edge of it. The Kepler conservation gate is the load-bearing one.

If those bounds hold, E ships. If Tier 1 misses, the kernel approximation is dropped and we move to D + reassessment for B/v0.2.

This decision becomes the §Hypothesis of a new lab notebook (`docs/experiments/2026-05-NN-rsqrt-kernel-approx.md`) when the work starts. The notebook `2026-05-11-perf-prediction-calibration.md` standing rules apply: the engine-ceiling envelope (now updated by this experiment to `~3-5×` kernel-only ceiling depending on distribution) is the upper bound on any prediction, not the lane-count ratio.
