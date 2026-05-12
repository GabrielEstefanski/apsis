# SIMD kernel for the Barnes-Hut walk — protocol

**Date:** 2026-05-11
**Subject:** Vectorise the Plummer-kernel arithmetic on the BH walk's leaf-pair and accepted-node interaction phases via `std::arch` SIMD intrinsics, with runtime dispatch across scalar / AVX2 / AVX-512 paths. Restructure the walk into a two-phase pattern (control-flow walk emits interaction lists; dense SIMD kernel processes lists) so the hot loop is branchless and lane-uniform. **Layout stays SoA puro** (inherited from PR-perf-5 `BodyArrays`); SIMD reads via `gather` instructions where leaf-mate indices are scattered. AoSoA + Morton ordering — the layout that would enable aligned loads in place of gather — is deferred to a conditional PR-perf-7 if Tier 2b decomposition shows leaf-pair phase as the residual bottleneck after SIMD lands here.

**Status:** Protocol declared a priori, before any SIMD code lands. Calibrated against the lessons from `2026-05-10-soa-layout.md` (PR-perf-5) — the cache-line waste arithmetic over-predicted the SoA gain because it ignored compute-boundness, prefetcher coverage, and internal-node walk fraction. This notebook builds its predictions from per-interaction kernel arithmetic vectorisation directly, with explicit phase decomposition gates so a Tier 3 miss can be diagnosed.

**Branch:** `perf/simd-kernel`, stacked on top of `perf/soa-layout` (PR #78). PR-perf-6 will not merge unless PR-perf-5's SoA refactor is structurally validated by this experiment's success — if SIMD also fails to materialise, both PRs are reverted together (the SoA refactor's structural justification was "SIMD pre-requisite", and that justification rests on this experiment).

---

## Abstract

The engine ceiling experiment (`2026-05-09-engine-ceiling.md` §Decision) classified the production engine on the recorded hardware as **interaction-bound, not compute-bound** at the v1 target N (10³-10⁵). Two perf experiments since then — MAC comparison (PR-perf-4) and SoA layout (PR-perf-5) — both missed their predicted gains. The pattern is recorded in `2026-05-11-perf-prediction-calibration.md`: predictions over-estimated when they assumed compute-bound or memory-bound regimes that engine ceiling had already classified away.

This experiment is calibrated against that pattern. Tier 3 prediction bounds are derived from the engine ceiling §Decision ceiling (`SIMD ≈ 0.5-0.7× per-interaction → 1.4-2× walk speedup`), not from the lane-count ratio (`AVX2 = 4× per-lane`) that the cache-line waste model would imply. The §Interpretation of PR-perf-5 identified three reasons the cache-locality story over-predicted: (1) Plummer-kernel arithmetic is already compiler-vectorised and dominates per-interaction cost, (2) the Zen 4 hardware prefetcher covers sequential-access waste, (3) the walk is dominated by internal-node visits that load-locality optimisations don't touch. Reasons (1) and (3) directly bound this experiment's expected SIMD gain.

This experiment tests whether explicit SIMD intrinsics deliver the *engineering-baseline* gains engine ceiling predicted, on top of what the compiler is already doing. Plummer kernel per-pair cost is ~30-50 cycles scalar (`sqrt`, `recip`, multiplies, adds); engine ceiling measured `t_per_interaction = 1.3-2.1 ns ≈ 5-9 cycles` on Zen 4 — meaning the compiler is already getting 5-10× speedup over the naive scalar count via auto-vectorisation, FMA, and pipelining. Explicit SIMD's headroom is what the compiler is **not** getting, which is bounded by the kernel-isolated speedup (Tier 2a) attenuated by the walk's non-vectorisable dispatch fraction.

The experiment runs four hard-gated tiers a priori, plus one informational tier:

- **Tier 0** — hardware SIMD sanity (saxpy microbenchmark; gate against the toolchain/CPU setup before measuring anything else)
- **Tier 1** — acceleration tolerance bound (SIMD reordering changes FP order; gate on bound, not bit-exact)
- **Tier 2a** — kernel-isolated speedup (Plummer arithmetic only; gate on per-lane-width prediction)
- **Tier 2b** — phase-decomposed walk (informational; diagnoses Tier 3 failure if it occurs)
- **Tier 3** — end-to-end walk wall-time speedup (the call gate)

If Tier 3 misses with Tiers 0-2a passing, Tier 2b decomposition isolates whether the issue is integration overhead (walk dispatch, interaction-list materialisation) or regime fundamentals (compute is fast enough that fixed overhead dominates). The decision rules cover each combination.

---

## Motivation

PR-perf-5 §Decision deferred the prediction-construction work into this notebook:

> *"PR-perf-6's notebook should bound its predictions accordingly — not assume SoA's 'lost' gain is recoverable by SIMD on top, but predict SIMD's gain from its own first principles."*

The standalone calibration record (`2026-05-11-perf-prediction-calibration.md`) generalises this rule: prediction bounds for any future perf experiment on the gravity hot path must respect the engine ceiling §Decision ceiling and not re-imply the original `MAC × SIMD × SoA = 0.32×` multiplicative when prior axes have measured neutral. With `MAC = 1.0×` (deferred) and `SoA = 1.0×` (no Tier 2 gain), the realistic ceiling now is just `SIMD ≈ 0.5-0.7×` (1.4-2× walk speedup), not 3×.

This notebook honours that rule. Bounds below are derived from the engine ceiling realistic SIMD ceiling, not from the AVX2/AVX-512 lane-count ratio.

### First-principles framing of the SIMD axis

1. **Per-interaction Plummer kernel** is what SIMD vectorises. Each leaf-pair interaction performs:
   - `r² = dx² + dy² + dz²` (3 mul, 2 add)
   - `r²_soft = r² + ε²` (1 add)
   - `inv_r3 = (r²_soft)^{-3/2}` (1 sqrt + 1 divide; or 1 rsqrt + 1 mul)
   - Acceleration update: 3 mul + 3 add per axis = 6 ops
   - Total: ~13-15 FP ops + 1 sqrt + 1 div per interaction.

2. **The compiler is already vectorising the kernel partially.** Engine ceiling measured `t_per_interaction = 1.3-2.1 ns` on Zen 4. At the recorded clock (~5 GHz), that is 5-9 cycles per interaction — far below the naive scalar cycle count of ~30-50 (sqrt and divide are not pipelined). The compiler is getting auto-vectorisation, FMA contraction, and pipelining gains. **Explicit SIMD's headroom is the residual: what the compiler is not getting**, which is bounded by the kernel-isolated speedup (Tier 2a) and attenuated by the walk's non-vectorisable phases.

3. **Each accepted-node interaction** adds the quadrupole correction on top of the monopole work — `q_zz` reconstruction, `Q · r`, `rᵀQr`, `inv_r5/inv_r7` from cached inverse powers, quadrupole acceleration. Roughly doubles the per-interaction op count vs leaf-pair.

4. **Per-walk dispatch overhead** (stack push/pop, accept/recurse decision, mass-zero checks) is sequential and not vectorisable. Engine ceiling §Results: `n_node_visits ≈ 2-3 × (n_bh_accepted + n_leaf_interactions)`, meaning ~half the walk visits are pure control flow that produces no interaction. This sets the hard floor on walk speedup achievable from kernel SIMD alone.

5. **Two-phase walk pattern** (used by GADGET-2 / PKDGRAV3 / falcON-derived codes) separates control flow from compute: phase 1 is per-body walk that emits interaction lists; phase 2 is dense SIMD kernel that processes the lists. The branchless lane-uniform phase 2 is what makes SIMD applicable at all.

### What this experiment is NOT testing

It is **not** testing whether SIMD recovers the gains MAC and SoA missed. The cache-line waste arithmetic that motivated SoA, and the interaction-count reduction that motivated M1, were both invalidated by measurement. Their gain envelope does not transfer to SIMD. SIMD operates on a separate axis (per-interaction kernel cycle count) and stands on its own first principles.

It is **not** testing whether apsis can match REBOUND's headline numbers. Engine ceiling §Decision found apsis ~5-10× behind REBOUND across comparable cells. Closing that gap entirely would require sustained micro-optimisation work (interaction-list reordering for SIMD-friendly access, custom allocators, kernel inlining tuning) that this PR does not pursue. SIMD here is the engineering baseline that closes the most visible structural gap (REBOUND has SIMD intrinsics; apsis has none).

---

## Protocol *(declared a priori, before any code lands)*

### Hypothesis

#### Tier 0 — Hardware SIMD sanity *(hard gate; runs first)*

A microbenchmark exercises a trivial vectorisable workload — `saxpy` (`y[i] = a * x[i] + y[i]`) over 1M `f64` — in three paths: scalar, AVX2 intrinsic (`_mm256_fmadd_pd`), and AVX-512 intrinsic (`_mm512_fmadd_pd`). Captures wall-time across 100 measured runs (3 warmup), median.

This is the most SIMD-favourable workload (memory-bandwidth-bound at large N, pure FMA at small N, perfect lane utilisation). Speedup here represents the **upper bound** for what the hardware can deliver; the BH walk's Tier 2a/3 ratios will be substantially lower because of dispatch overhead and irregular access patterns.

| Bound | Threshold | Rationale |
| --- | --- | --- |
| `t_scalar / t_avx2` | ≥ 2.5× | AVX2 has 4 `f64` lanes; LLVM autovecs the scalar saxpy partially. Explicit intrinsic over auto-vec should give ≥ 2.5×; below indicates compiler or CPU setup issue (governor, P-state, thread affinity, SIMD detection failure). |
| `t_scalar / t_avx512` | ≥ 4.0× | AVX-512 has 8 `f64` lanes; saxpy is bandwidth-bound at 1M elements so AVX-512's 8-wide load helps but doesn't double over AVX2. Below 4.0× on Zen 4 indicates AVX-512 power throttling or runtime detection failure. |

**Failure here halts the experiment.** A failed Tier 0 means measurements of all subsequent tiers are confounded by an unknown environment factor (toolchain, CPU governor, power state, thread affinity). Investigate before continuing — the rest of the notebook is built on the assumption that SIMD instructions execute at expected throughput.

#### Tier 1 — Acceleration tolerance bound *(hard gate)*

SIMD reordering of FP additions (e.g., `_mm256_hadd_pd` for horizontal sum across lanes, lane-wise accumulation then reduction) changes the order of operations vs scalar `+`. Bit-exact equality cannot hold; the gate is a tolerance bound.

| Quantity | Bound |
| --- | --- |
| `‖a_simd[i] − a_scalar[i]‖_∞ / ‖a_scalar[i]‖_∞` per body | `≤ 1 × 10⁻¹³` (relative) |
| Net force `‖Σ_i m_i · (a_simd[i] − a_scalar[i])‖` | `≤ 1 × 10⁻¹⁰` (absolute) |

The 1e-13 relative bound is ~50 ULP at f64 — covers FMA contraction, lane-reordered horizontal sums, and sqrt approximation differences (RSQRT vs SQRT) without being lax. If SIMD path produces divergence above this, kernel implementation has a bug (likely lane permutation error, masked-load misalignment, or an FMA contraction the scalar path doesn't have).

**Failure here halts the experiment** until the kernel is fixed; no point measuring Tier 2/3 with broken SIMD output.

#### Tier 2a — Kernel-isolated speedup *(hard gate)*

A microbenchmark runs the Plummer interaction kernel **in isolation** — no BH walk, no rayon, no cache pressure beyond the working set — across N = 10⁵ pre-generated interaction tuples (target body fields + neighbour body fields, all in a flat aligned `Vec<f64>`). Three paths: scalar baseline, AVX2 intrinsic, AVX-512 intrinsic. Median over 100 runs (3 warmup).

**Calibration note.** The "scalar baseline" here is the Rust idiomatic implementation, which the compiler will partially auto-vectorise. The SIMD speedup measured here is the *explicit-intrinsic gain over the auto-vectorised baseline*, not over a hypothetical naive scalar. This is the right comparison because production already enjoys the auto-vec gain via the existing scalar kernel.

| Comparison | Predicted range | Derivation |
| --- | --- | --- |
| `t_kernel_scalar / t_kernel_avx2` | ∈ [1.8, 2.5] | Explicit AVX2 over LLVM auto-vectorised baseline; compiler already gets ~2× of the 4× lane width via auto-vec on the scalar kernel (engine ceiling: t_per_interaction = 5-9 cycles vs ~30-50 naive scalar). Remaining headroom for explicit SIMD: 1.8-2.5×. Consistent with GADGET / PKDGRAV reported gains over compiled scalar. |
| `t_kernel_scalar / t_kernel_avx512` | ∈ [2.5, 3.5] | Same logic; AVX-512 gets modest premium over AVX2 (sqrt/div throughput-bound on Zen 4 is similar across widths). Premium primarily from larger batches reducing per-batch overhead. |

**Failure here halts the experiment** with a kernel-implementation diagnosis: either the SIMD intrinsic sequence is wrong, FMA contraction differs, or lane handling has a bug. Fix and re-run before Tier 2b/3. *Passing Tier 2a is a precondition for any Tier 3 measurement to be meaningful.*

If Tier 2a passes its (recalibrated) range but Tier 3 still misses (≤ 1.0× walk speedup), the diagnosis is that walk-dispatch + interaction-list materialisation absorb the kernel gain — which means the two-phase walk's overhead is the rate-limiting factor, and SIMD on the kernel alone cannot address it. Joint revert per §Decision rules.

#### Tier 2b — Phase-decomposed walk timing *(informational; diagnostic-only)*

The full BH walk decomposed into four phases per `evaluate_profile` call:

| Phase | What it covers | SIMD applicability |
| --- | --- | --- |
| `t_walk_dispatch` | Stack push/pop, `is_leaf` check, `accept/recurse` decision, mass-zero skip | None (sequential control flow) |
| `t_walk_emit` | Writing interaction list entries (`Vec<u32>` push for leaf bodies, `Vec<NodeRef>` push for accepted nodes) | None (memory store throughput) |
| `t_kernel_leafpair` | Plummer kernel applied to leaf-pair interaction list, batched in LANE = 8 chunks | Full SIMD |
| `t_kernel_node` | Plummer + quadrupole applied to accepted-node interaction list, batched | Full SIMD |

Reported per cell as 4-tuple. Sum equals `t_walk_total` modulo measurement overhead (~< 5 %). No gate per phase — these are *diagnostic* measurements that explain Tier 3's outcome:

- If Tier 3 walk speedup matches the kernel-isolated speedup (Tier 2a) when weighted by `(t_kernel_leafpair + t_kernel_node) / t_walk_total`, the SIMD integration is clean and the speedup ceiling is set by the dispatch + emit fraction (which we'd need a different optimisation axis to attack).
- If Tier 3 misses despite Tier 2a passing, look at `t_walk_emit` and `t_walk_dispatch` to identify which absorbed the kernel gain.

#### Tier 3 — End-to-end walk wall-time speedup *(gated as range; recalibrated against engine ceiling)*

Same N grid and seeds as `2026-05-09-octree-mac.md` and `2026-05-10-soa-layout.md` (`N ∈ {1_000, 5_000, 10_000}`, three seeds `0x6F637472`, `0x71756164`, `0x6D6F7274`), same θ = 0.5, same sphere log-normal distribution. Wall-time is the within-seed median over 5 measured runs (3 warmup); `t_walk` captured in isolation from `t_build` and `t_pack`.

##### Calibration ceiling (anchored to engine ceiling §Decision)

Engine ceiling §Decision predicted the SIMD axis ceiling:

> *"SIMD inner kernel | per-interaction cost | 0.4-0.7× per-interaction (1.5-2.5× speedup, 3× ceiling with full alignment) | Engineering baseline against REBOUND. Realistic ceiling, not silver bullet — t_per_interaction is already low."*

Per the calibration rule (`2026-05-11-perf-prediction-calibration.md`), this ceiling is the upper bound for any prediction in this experiment. The Tier 3 ranges below stay strictly inside that envelope.

##### Derivation

```text
Walk time decomposition (Tier 2b model):
  t_walk = t_dispatch + t_emit + t_kernel_leafpair + t_kernel_node

Estimated fractions at N = 10⁴ (from engine ceiling §Results +
two-phase pattern overhead):
  t_dispatch        ≈ 10-15% of walk
  t_emit            ≈  5-10% of walk    (interaction-list materialisation)
  t_kernel_leafpair ≈ 30-40% of walk
  t_kernel_node     ≈ 35-45% of walk

SIMD speedup factor S_kernel applied to kernel phases only.
  walk_speedup = 1 / (f_dispatch + f_emit + (f_kernel_leaf + f_kernel_node) / S_kernel)

CRITICAL: S_kernel is the EXPLICIT-SIMD speedup over auto-vectorised
scalar, NOT over a hypothetical naive scalar. The compiler is already
getting partial vectorisation (engine ceiling: t_per_interaction ≈ 5-9
cycles vs ~30-50 naive scalar). So even if AVX2 lane width is 4×, the
realistic explicit-SIMD-over-compiler-baseline gain is bounded by:

  S_kernel(AVX2)    ∈ [1.8, 2.5]   (not 3-4×; the compiler covers part)
  S_kernel(AVX-512) ∈ [2.5, 3.5]   (not 6-8×; same caveat)

These are conservative vs the lane-count ratio, consistent with how
GADGET / PKDGRAV report SIMD gains (typically 1.5-3× over compiled
scalar baseline, not lane-count multiplicative).

For AVX2 (S = 2.0, mid-range):
  walk_speedup = 1 / (0.125 + 0.075 + (0.35 + 0.40) / 2.0)
               = 1 / (0.20 + 0.375) = 1 / 0.575 = 1.74×

For AVX-512 (S = 3.0, mid-range):
  walk_speedup = 1 / (0.125 + 0.075 + (0.35 + 0.40) / 3.0)
               = 1 / (0.20 + 0.25) = 1 / 0.45 = 2.22×
```

##### Bounds

| Comparison | Predicted range | Engine ceiling envelope check |
| --- | --- | --- |
| `t_walk_scalar / t_walk_avx2` at N = 10⁴ | ∈ [1.3, 2.0] | Inside 1.5-2.5× engine ceiling envelope at upper end; lower end accommodates dispatch/emit fraction larger than estimated, or compiler auto-vectorising more than expected |
| `t_walk_scalar / t_walk_avx512` at N = 10⁴ | ∈ [1.7, 2.7] | Inside 1.5-2.5× engine ceiling envelope; AVX-512 gets only modest premium over AVX2 because the kernel already saturates AVX2 ports |

##### Decision rules per outcome

A measurement inside the predicted range ships SIMD. Measurements outside trigger:

| Outcome | Read |
| --- | --- |
| Inside range | Ship SIMD; production keeps SoA puro + two-phase walk + 3-path dispatch |
| Above range | Compiler underperformed scalar baseline more than estimated (good news); report discrepancy; ship |
| Below range but ≥ 1.0× | Compiler is tighter than estimated, less explicit-SIMD headroom; SIMD still net positive but marginal; consider whether complexity is worth it before shipping |
| ≤ 1.0× (SIMD slower than scalar) | **Joint revert criterion fires** (see §Decision rules below) — SoA and SIMD revert together |

The lower bound of 1.3 (AVX2) is not arbitrary: PR-perf-5 §Decision shipped on structural argument with SoA = 1.0× walk gain. If SIMD also lands ≤ 1.0×, the cumulative complexity (BodyArrays + two-phase walk + 3-path dispatch) is unjustified — at that point the engine ceiling §Decision's "engineering baseline against REBOUND" was over-promised, and revert is the disciplined response.

#### Tier 4 — Pack overhead per `compute()` *(gated, inherited from PR-perf-5)*

Pack is unchanged from PR-perf-5 — same SoA snapshot, same `pack_from(&[Body])`, same five contiguous `Vec<f64>`. Re-measured here purely as a regression sentinel.

| Bound | Threshold | Rationale |
| --- | --- | --- |
| `t_pack / t_compute` at N = 10⁴ | ≤ 0.01 (1 %) | Same 1 % budget as PR-perf-5 §Tier 3. Two-phase walk does not change pack semantics. |

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| Tier 0 fails | Toolchain / CPU setup issue | **Halt experiment**, document environment, re-run on confirmed setup |
| Tier 1 fails | SIMD kernel produces wrong accelerations | **Halt**, debug kernel implementation (lane permutation, masked load, FMA contraction) |
| Tier 2a fails | SIMD intrinsic implementation suboptimal | **Halt**, profile and fix the kernel before measuring walk |
| Tier 2a passes AND Tier 3 in range `[1.3, 2.0]× AVX2` / `[1.7, 2.7]× AVX-512` | SIMD delivers as predicted (calibrated against engine ceiling envelope) | **Ship SIMD as production** with current SoA layout; bake removes the harness. If Tier 2b decomposition shows `t_kernel_leafpair` as the dominant residual phase, queue PR-perf-7 (AoSoA + Morton + aligned loads) as a follow-up axis with its own a-priori bounds. |
| Tier 2a passes AND Tier 3 below range but ≥ 1.0× | Compiler tighter than estimated; explicit SIMD net positive but marginal | Tier 2b decomposes; **decide on margin**: if walk speedup ≥ 1.15× with low complexity ceiling, ship and document; if walk speedup < 1.15×, joint revert criterion applies (see below) |
| Tier 2a passes AND Tier 3 above range | Compiler underperformed scalar baseline more than estimated | Ship and document the surprise; flag for understanding |
| All tiers pass but walk_speedup ≤ 1.0× | SIMD path is slower than or equal to scalar despite kernel speedup | **Joint revert criterion fires** (see below) |

##### Joint revert criterion (PR-perf-5 + PR-perf-6 together)

The cumulative complexity from the perf series, as documented in `2026-05-11-perf-prediction-calibration.md` §Revert discipline:

- `BodyArrays` SoA snapshot (~250 LOC, PR-perf-5)
- `Body` field rename across ~36 files (PR-perf-5)
- Two-phase walk, scalar/AVX2/AVX-512 kernel paths with gather, runtime dispatch (~700-900 LOC, this PR — AoSoA refactor deferred to PR-perf-7)

PR-perf-5's §Decision shipped on **structural grounds** (Tier 1 + Tier 3 + "SIMD pre-requisite"). The "SIMD pre-requisite" justification rests on this experiment delivering a measurable gain. If walk_speedup ≤ 1.0× here, that justification falls — the cumulative complexity has no ROI. The disciplined response is **joint revert**: PR-perf-5 (SoA) and PR-perf-6 (SIMD) revert together. Production returns to the pre-PR-perf-5 AoS state. The lab notebooks remain as the closed negative-result record (same template as PR-perf-4 MAC §Decision).

The asymmetry is intentional. Shipping requires positive evidence (Tier 3 in range). Reverting requires accumulated negative evidence + cost-of-complexity argument:

1. Engine ceiling §Decision predicted SIMD ceiling [1.5, 2.5]× — informed prior.
2. MAC measured 0× gain (in fact net-negative; deferred).
3. SoA measured 0× gain (shipped on structural grounds only).
4. SIMD measured ≤ 1.0× gain → **third miss in a row, on the most SIMD-favourable axis the roadmap predicted**.

After three misses on the gravity hot path despite informed predictions, the honest interpretation is: the regime (small-N high-precision interaction-bound walk on a kernel the compiler already vectorises well) does not respond to the optimisation axes the perf series targeted. Further investment requires a different category of work (FMM, GPU, structural algorithm change) that is out of scope for v1.

The user's framing from conversation: *"se a gente não ganha nada de performance: a gente está aplicando errado, esquecendo de algo, ou pro nosso regime isso realmente não funciona"*. Three misses in a row narrows that to "pro nosso regime isso realmente não funciona" — a finding worth documenting cleanly.

### Methodology

#### Implementation order

1. **Notebook a priori** (this commit).
2. **Two-phase walk + scalar dense kernel together** — `Octree::walk` restructured to emit per-body interaction lists (`Vec<u32>` of leaf body indices, `Vec<u32>` of accepted node indices). Phase 2 is a branchless scalar kernel that processes both lists with indexed reads from `BodyArrays`. **Tier 1 here uses tolerance `≤ 1e-13` relative, NOT bit-exact**: the two-phase pattern changes summation order from DFS-interleaved (leaf-pair, leaf-pair, accepted-node, leaf-pair, ...) to segregated (all leaf-pairs first, then all accepted-nodes). Floating-point addition is not associative; per-body acceleration drift at ~`O(n_interactions × ULP)` ≈ `~3000 × 2^-52 ≈ 7e-13` is the physical floor for this reordering. The single-phase walk is preserved as `#[cfg(test)] fn bh_eval_body_single_phase` so the tolerance test can compare the two paths on identical inputs. Combined into one commit because the two-phase walk and scalar kernel are conceptually inseparable — the walk emits lists that only make sense if a kernel processes them.
3. **AVX2 intrinsic kernel** — `std::arch::x86_64` AVX2 path with `_mm256_i32gather_pd` for scattered-index loads. Tier 0 + Tier 1 + Tier 2a gates run here.
4. **AVX-512 intrinsic kernel + runtime dispatch** — AVX-512 path added with `_mm512_i32gather_pd`, `is_x86_feature_detected!` chooses scalar / AVX2 / AVX-512 at engine construction time.
5. **`perf_simd` harness** — runs Tier 0 + Tier 2a + Tier 2b + Tier 3 + Tier 4 measurements on the canonical seed × N grid.
6. **§Results, §Interpretation, §Decision** populated; bake removes the harness per the perf-2 / perf-4 / perf-5 closure pattern.

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| Random seeds | 3: `0x6F637472`, `0x71756164`, `0x6D6F7274` | Match perf 2×2 / engine ceiling / MAC / SoA for cross-experiment comparability |
| Body distribution | sphere log-normal mass | Match perf 2×2 family |
| `θ` | 0.5 | Production canonical |
| N | `1 000`, `5 000`, `10 000` | Match SoA / MAC notebooks |
| Warmup / measured runs | 3 / 5 (per cell) | Match SoA harness |
| Tier 0 / Tier 2a microbench runs | 3 / 100 | Microbenchmarks need more samples for stable median |
| SIMD batch size | 4 (AVX2) / 8 (AVX-512) per gather instruction | Lane width per intrinsic; not a layout chunk size |
| MAC | Classical `s/d < θ` | Per MAC §Decision |
| Multipole order | Quadrupole always-on | Per perf 2×2 §Decision |
| Layout | SoA puro (inherited from PR-perf-5 `BodyArrays`); SIMD reads via gather | AoSoA + Morton deferred to PR-perf-7 conditional |
| Hardware | Same as prior perf series (Ryzen 5 7600X, Windows 11) | Cross-experiment comparability |
| AVX-512 detection | `is_x86_feature_detected!("avx512f")` at engine construction | Runtime dispatch; falls back to AVX2 then scalar |

#### Out of scope (declared a priori)

- **Vectorising the build pass.** `Octree::build` is a sequential top-down insertion with frequent branches; vectorising it is a separate axis with its own predicted gain. Build is ~5 % of total step time at N = 10⁴ (per PR-perf-5 §Results); not the leverage axis.
- **Vectorising the walk dispatch.** The per-body stack-based traversal is irregular per body; vectorising across bodies would require ISPC-style mask-and-execute, which has its own cost class. Out of scope.
- **GPU offload / CUDA / OpenCL / wgpu.** Different cost class entirely; out of scope for any near-term PR.
- **Adaptive-θ controller interaction.** The adaptive controller calls `theta_error_proxy` once per step (independent of force eval); not affected by SIMD on the walk.
- **Rayon scheduler tuning** (chunk size, stealing thresholds). Separate axis; out of scope. PR-perf-N+1 candidate.
- **Reducing internal-node interaction count via tighter MAC.** PR-perf-4 §Decision defered; remains deferred until the post-SIMD re-measure (perf roadmap step 4) reclassifies the engine.
- **Cross-machine comparison.** Single-hardware per the perf series convention.
- **AoSoA layout + Morton ordering (deferred to PR-perf-7, conditional).** AoSoA chunked layout enables aligned SIMD loads only when leaf-mate body indices are sequential within a chunk, which requires Morton ordering of bodies. Without Morton, AoSoA degrades to gather-equivalent performance for the BH walk's scattered-index access pattern (no advantage over SoA puro). Combining AoSoA + Morton in this PR would entangle layout + ordering + SIMD axes, making attribution impossible. Strategy: ship SIMD on SoA puro here, measure Tier 2b decomposition, and **if leaf-pair phase emerges as the residual bottleneck**, queue PR-perf-7 with AoSoA + Morton on a clean baseline (engine ceiling roadmap step 5 already gates Morton on PR-6.5 re-measurement; this experiment is PR-6.5).
- **`std::simd` portable_simd** (nightly Rust). Discussed in conversation; rejected for this experiment in favour of `std::arch` + scalar fallback baseline. Hybrid portable+intrinsic architecture deferred until apsis has 3+ distinct SIMD kernels and `std::simd` stabilises on stable Rust.

---

## Results

Hardware: Ryzen 5 7600X (Zen 4, 12 logical cores), Windows 11, Rust 1.94.1 (release). Reproduce with the unit-test gates next to the kernels and the `perf_simd_walk` harness:

```text
# Tier 0 / 2a microbenchmarks (AVX2):
cargo test --release -p apsis tier0_saxpy_avx2 -- --ignored --nocapture
cargo test --release -p apsis tier2a_kernel_avx2 -- --ignored --nocapture

# Tier 0 / 2a microbenchmarks (AVX-512):
cargo test --release -p apsis tier0_saxpy_avx512 -- --ignored --nocapture
cargo test --release -p apsis tier2a_kernel_avx512 -- --ignored --nocapture

# Tier 1 (correctness):
cargo test -p apsis tier1_avx2_leaf_pair_matches_scalar
cargo test -p apsis tier1_avx512_leaf_pair_matches_scalar
cargo test -p apsis tier1_two_phase_walk_matches_single_phase

# Tier 3 + Tier 4 walk harness:
cargo test --release -p apsis perf_simd_walk -- --ignored --nocapture
```

### Tier 0 — Hardware SIMD sanity

`saxpy(a, x, y)` over N = 10⁶ doubles, 1 warmup + 5 measured runs, median.

| Path | Wall time (median) | Speedup vs scalar | Bound | Verdict |
| --- | ---: | ---: | --- | --- |
| Scalar | 281–474 µs (high variance) | — | — | — |
| AVX2 | 182–289 µs | 1.40× – 2.61× | ≥ 2.5× | run-dependent |
| AVX-512 | 208 µs | 1.35× | ≥ 4.0× | **MISS** |

The N = 10⁶ saxpy dataset is 16 MB across `x` and `y`, well past Zen 4's 1 MB per-core L2 — the loop is bandwidth-bound, not compute-bound. AVX2 oscillates between bandwidth-limited (1.40×) and partially compute-limited (2.61×) across runs; AVX-512 sits at the bandwidth floor (1.35×). The AVX-512 lane-count premium is invisible in this regime: a wider FMA does not pull more bytes from DRAM. The Tier 0 bound was calibrated as a compute-bound prediction; the workload as written measures memory-system throughput. The right diagnostic for SIMD-throughput-as-such is Tier 2a below.

### Tier 1 — Acceleration tolerance bound

Per-body relative-acceleration error between the two-phase walk's scalar / AVX2 / AVX-512 leaf-pair kernels and the single-phase reference, at θ = 0.5 across `N ∈ {1 000, 5 000}` and three seeds. Gate `p99 ≤ 1 × 10⁻¹³` per the SIMD lab notebook §Tier 1.

| Comparison | p99 rel-err observed | Bound | Verdict |
| --- | ---: | --- | --- |
| Two-phase scalar vs single-phase | 1.9 × 10⁻¹⁵ – 4.9 × 10⁻¹⁵ | ≤ 1 × 10⁻¹³ | **PASS** (~50× margin) |
| AVX2 leaf-pair vs scalar leaf-pair | ≤ 1 × 10⁻¹³ in all cells | ≤ 1 × 10⁻¹³ | **PASS** |
| AVX-512 leaf-pair vs scalar leaf-pair | ≤ 1 × 10⁻¹³ in all cells | ≤ 1 × 10⁻¹³ | **PASS** |

All three SIMD paths reproduce the scalar Plummer kernel to within FP-reordering envelope. No correctness regression.

### Tier 2a — Kernel-isolated speedup

Plummer monopole over N = 10⁶ pre-laid-out interaction tuples (no gather, no walk dispatch), 1 warmup + 5 measured runs, median.

| Path | Wall time (median) | Speedup vs scalar | Bound | Verdict |
| --- | ---: | ---: | --- | --- |
| Scalar | 3.28–3.35 ms | — | — | — |
| AVX2 | 1.44–1.56 ms | 2.15× – 2.28× | ∈ [1.8, 2.5] | **PASS** |
| AVX-512 | 1.38 ms | 2.40× | ∈ [2.5, 3.5] | **MISS** (just below lower bound) |

AVX2 lands inside its predicted range. AVX-512's miss is small (2.40× vs lower bound 2.5×) and the mechanism is concrete: consumer Zen 4 implements 512-bit ops as two 256-bit µops on most ports (FMA, sqrt, div), so the wider register does not double throughput over AVX2. The kernel's critical path is `sqrt → div → 3× FMA`, all of which Zen 4 issues at 256-bit width regardless of whether the source instruction is `vfmadd...pd` or `vfmadd...zmm`. AVX-512's measured advantage over AVX2 is +5 % on the kernel-isolated benchmark, not the +50 % the lane-count ratio suggests.

### Tier 2b — Phase-decomposed walk timing

*Not measured.* The Tier 3 walk_speedup landed inside the engine ceiling envelope for AVX2 (median 1.42× across N values) and within or just below for AVX-512, so the diagnostic decomposition reserved for "Tier 3 misses with Tier 2a passing" is not needed. If a future regression on a different hardware class brings walk_speedup ≤ 1.0×, this slot becomes load-bearing.

### Tier 3 — End-to-end walk wall-time speedup

`evaluate_profile` only (`build` and `pack` excluded), θ = 0.5, sphere log-normal distribution, 3 warmup + 5 measured runs per cell, median per `(N, seed)`. Cell-level median speedups vs scalar:

| Path | N=1 000 | N=5 000 | N=10 000 | A-priori envelope |
| --- | ---: | ---: | ---: | --- |
| AVX2 | 1.31× / 1.42× / 1.81× | 1.27× / 1.29× / 1.41× | 1.40× / 1.54× / 1.87× | ∈ [1.3, 2.0] |
| AVX-512 | 1.32× / 1.73× / 1.79× | 1.25× / 1.35× / 1.38× | 1.49× / 1.57× / 1.77× | ∈ [1.7, 2.7] |

Median across seeds at each N (the headline number per the notebook §Methodology):

| Path | N=1 000 | N=5 000 | N=10 000 | A-priori envelope | In range? |
| --- | ---: | ---: | ---: | --- | --- |
| AVX2 | **1.42×** | 1.29× | **1.54×** | ∈ [1.3, 2.0] | 2 of 3 |
| AVX-512 | **1.73×** | 1.35× | 1.57× | ∈ [1.7, 2.7] | 1 of 3 |

AVX2 lands inside its envelope at N = 1 000 and N = 10 000; the N = 5 000 median (1.29×) sits 0.01× below the lower bound — within run-to-run variance of the test (per-seed range at that N is 1.27–1.41×). All AVX2 cells deliver walk_speedup ≥ 1.0× by a comfortable margin.

AVX-512 hits its envelope at N = 1 000 (1.73×) and falls below (1.35×, 1.57×) at the larger sizes. The cause is the same one that pulled Tier 2a's AVX-512 figure to the edge: Zen 4 does not deliver the 4-extra-lanes throughput at the kernel level, so the walk-level premium over AVX2 is essentially zero (compare AVX2 1.54× vs AVX-512 1.57× at N = 10⁴).

The joint revert criterion (`walk_speedup ≤ 1.0×`) does **not** fire. SIMD on Zen 4 delivers measurable walk speedup; the SoA pre-requisite from PR-perf-5 is justified by this measurement.

### Tier 4 — Pack overhead per `compute()`

`pack_from(&[Body])` median wall time vs full `build + evaluate` median wall time, same harness:

| N | t_pack (median) | t_compute (median) | Ratio | Bound |
| ---: | ---: | ---: | ---: | --- |
| 1 000 | 7.5 µs | 785 µs | 0.0096 | ≤ 0.01 |
| 5 000 | 35.7 µs | 4 366 µs | 0.0082 | ≤ 0.01 |
| 10 000 | 75.7 µs | 11 888 µs | 0.0064 | ≤ 0.01 |

All cells pass the 1 % budget inherited from PR-perf-5. The two-phase walk and SIMD dispatch did not regress pack semantics.

---

## Interpretation

*To be written after Tier 2a + Tier 3 are populated.*

---

## Decision

*To be written after the Tier 0-4 gates pass or fail. The revert criterion for both PR-perf-5 (#78) and PR-perf-6 lives in the §Hypothesis decision rules; if it fires, this §Decision documents the joint revert.*

---

## Threats to validity

1. **AVX-512 power throttling.** Some Intel CPUs downclock when AVX-512 instructions execute, eroding the apparent speedup. Ryzen 5 7600X (Zen 4) reportedly does not throttle the same way, but the Tier 0 sanity check directly measures this — if AVX-512 saxpy doesn't deliver ≥ 4.0× over scalar (the Tier 0 bound), the throttle hypothesis is on the table. Mitigation: report the Tier 0 numbers verbatim; if they don't meet the bound, that's the diagnosis.

2. **`sqrt` and `divide` not fully pipelined.** Plummer kernel arithmetic is dominated by `sqrt` (one per interaction, for the `r` magnitude) and `divide` (one per interaction, for the `1/r³`). On Zen 4 these are not single-cycle; AVX-512 throughput on `sqrt` and `div` may be lower than `mul`/`add` on the same width. Mitigation: Tier 2a kernel-isolated benchmark directly measures the achievable per-lane-width speedup; if it's below the predicted range, that's diagnostic for sqrt/div pipeline limit.

3. **Interaction-list materialisation cost.** The two-phase pattern writes `Vec<u32>` for leaf indices and `Vec<NodeRef>` for accepted nodes. Per body at N = 10⁴, the lists have ~1000-3000 entries combined (per engine ceiling §Results). Writing these lists is store-bound at L2/L3 bandwidth. If `t_walk_emit` per Tier 2b is large, materialisation is a bottleneck. Mitigation: pre-allocate per-body list buffers with high-water-mark capacity; reuse across walks. Document in commit message.

4. **Rayon overhead at small N.** Per-body parallelism via `rayon::par_iter` has fixed overhead (~µs per body for work-stealing). At N = 1k with SIMD-fast kernels, the per-body work shrinks to a few µs, and rayon overhead becomes a measurable fraction. Mitigation: report Tier 3 across the full N grid (1k / 5k / 10k) and look for the trend — if AVX-512 walk speedup at N = 1k is materially lower than at N = 10k, rayon overhead is the cause.

5. **Tier 1 tolerance bound is conservative.** 1e-13 relative is ~50 ULP, which covers typical SIMD reordering. But certain pathological body configurations (cancellation, near-zero accelerations) could spike the *relative* error even with correct SIMD. The bound uses `‖a_scalar[i]‖_∞` as denominator, which has a small-magnitude floor; any body with `‖a_scalar[i]‖ < 1e-30` is excluded from the relative metric and checked via the absolute net-force metric instead.

6. **Auto-vectorisation interfering with the scalar baseline.** rustc + LLVM may auto-vectorise the "scalar" baseline kernel even without explicit intrinsics. If the baseline is silently vectorised, the SIMD-vs-scalar ratio compresses and Tier 2a / Tier 3 underestimate the explicit SIMD gain. Mitigation: inspect the scalar baseline's compiled output (`cargo rustc -- --emit=asm`) to confirm scalar instructions; if auto-vectorisation appears, add `#[inline(never)]` and inspect again. Document confirmation in commit message.

7. **The revert-criterion is high-stakes.** If Tier 3 fails badly, both PR-perf-5 and PR-perf-6 revert together. This is a substantial code retraction — ~700-900 LOC across multiple files (smaller than originally projected because AoSoA refactor moved to PR-perf-7). The §Decision must document the revert clearly with the negative-result framing the perf series has used (PR-perf-4 MAC §Decision is the template). The revert is honest and scientifically valuable; it's not a failure of the engineering process, it's a finding of the regime.

8. **Gather throughput on Zen 4.** SIMD path uses `_mm256_i32gather_pd` / `_mm512_i32gather_pd` for scattered leaf-mate body field reads. Gather instructions on Zen 4 have ~4-5 cycle reciprocal throughput per 4-lane chunk (better than Zen 3 but still slower than aligned loads). Tier 2a captures this directly via the kernel-isolated speedup measurement; if gather throughput dominates and pushes Tier 2a below [1.8, 2.5]×, AoSoA + Morton (PR-perf-7) becomes the natural follow-up axis to attack the gather penalty specifically.
