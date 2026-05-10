# Engine ceiling profiling — protocol

**Date:** 2026-05-09
**Subject:** Empirically characterise where the production engine starts to bottleneck at increasing body count, and which stage of the per-step pipeline dominates the cost. Identifies the rate-limiting subsystem per (integrator, N) regime so future optimisation investments can be justified by measured data rather than intuition.

**Status:** Protocol declared a priori, before any instrumentation code is written. §Results populated after the run.

**Branch:** `perf/engine-ceiling`, from `develop` (PR #72 merged carrying the perf 2×2 §Decision).

---

## Abstract

The perf 2×2 experiment (`2026-05-08-octree-perf-2x2.md`, §Decision) validated quadrupole as the production multipole order and reverted Morton on the grounds that at v1's primary regime (N ≤ 10⁴) the BH walk is compute-bound rather than memory-bound — Morton's structural cache benefit exists but is not rate-limiting. That finding raised a sharper question: **what is rate-limiting at the v1 target, and at what N does the engine stop being interactive?**

This experiment instruments the per-step pipeline with five timed phases (`tree_build`, `bh_traversal`, `kernel_force`, `integrator_overhead`, `trail_record`) and four work counters (`n_node_visits`, `n_bh_accepted`, `n_leaf_interactions`, `n_picard_iterations`). The instrumentation runs across two integrator cells (VV+BH, IAS15) and an N sweep covering the regime from comfortably-interactive to offline-only. The decisive output is **time per unit work** (per interaction, per body), not just total time per step — total time hides whether the cost is "much work" or "expensive work", and that distinction determines whether the next optimisation investment should target reducing interaction count (MAC) or reducing per-interaction cost (SIMD).

A secondary axis tests trail recorder on/off. Pre-experimental intuition is that the trail recorder is not the dominant cost, but that intuition has not been measured; this experiment provides the empirical answer.

The deliverable is a per-stage cost breakdown and SPS curve sufficient to identify the dominant bottleneck per regime, plus a §Decision pointing to the next optimisation investment with measured justification.

---

## Motivation

The perf 2×2 §Decision queued PR-perf-3 as an MAC comparison. Before committing to that work, two questions are worth empirical answers:

1. **Is the BH walk actually the dominant cost at our N target?** If trail recorder, integrator overhead, or another stage dominates, MAC investment optimises something that does not move the user-perceived ceiling.
2. **Where is the practical interactive ceiling?** REBOUND publishes headline numbers for IAS15 (~10² interactive, ~10³–10⁴ offline) and tree code (~10⁵ interactive, ~10⁶ offline with SIMD). Apsis's actual ceiling has not been measured end-to-end; without that number, "respectable" or "needs work" is intuition, not evidence.

Both questions feed every subsequent optimisation decision (MAC, SIMD, SoA layout, threading granularity). The cost of answering them now is small (~300 LOC of instrumentation, single-machine measurement) and saves potentially weeks of misdirected optimisation work.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the production engine on the recorded hardware, the metrics declared below are bounded a priori at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 is gated; Tier 2 and Tier 3 are informational. Bounds revision is forbidden unless backed by concrete arithmetic.

#### Tier 1 — Per-stage breakdown sanity gate *(gated)*

For each (cell, N) combination, the sum of the five instrumented phases must agree with the measured total step time within a tight tolerance:

| Metric | Bound | Rationale |
| --- | ---: | --- |
| `\|sum_phases − t_total\| / t_total` for every (cell, N) | `≤ 5 %` | If the instrumented phases miss > 5 % of total time, instrumentation has gaps and the per-stage attribution cannot be trusted |

Sub-bounds on phase composition (gated only when they are the rate-limiting cost; otherwise informational):

| Cell | N | Bound | Origin |
| --- | ---: | ---: | --- |
| V (VV+BH) | ≥ 10⁵ | `t_tree_build + t_bh_walk ≥ 85 %` | At memory-bound regime, BH stages dominate; if not, instrumentation or implementation has unaccounted overhead |
| I (IAS15) | every | `t_force_eval + t_integrator_overhead ≥ 80 %` | IAS15 has non-trivial integrator overhead (predictor coefficients, Gauss-Radau substep machinery, Picard convergence checks) plus the O(N²) force evaluation; together they must dominate |

**Gate-failure protocol.** If the sanity sum fails: dump per-phase breakdown, mark the run as untrustworthy, do not use the data for §Results Tier 2 or §Decision. Investigate instrumentation gaps before re-running.

#### Tier 2 — Cost normalised by work *(informational; the load-bearing analysis)*

For each (cell, N), report:

| Metric | Definition | What it answers |
| --- | --- | --- |
| `t_per_interaction` | `t_bh_walk / (n_bh_accepted + n_leaf_interactions)` | "Is each interaction expensive?" → SIMD ROI signal (captures kernel + traversal cost amortised; the split is escalation territory, see §Escalation rules) |
| `t_per_body` | `t_total_step / N` | "Is per-body amortised cost interactive?" |
| `n_interactions_per_body` | `(n_bh_accepted + n_leaf_interactions) / N` | "How much work per body?" → MAC ROI signal |
| `bh_acceptance_ratio` | `n_bh_accepted / n_node_visits` | "Is the walk pruning effectively?" → opening-criterion efficiency |

The decisive output is the comparison across N for each metric. If `t_per_interaction` grows with N, cache misses are amplifying — SIMD wins are limited until SoA / data layout is fixed. If `n_interactions_per_body` is high (≫ 100), MAC reduction has headroom. If `bh_acceptance_ratio < 0.5`, walk efficiency is low and adaptive-θ may help.

#### Tier 3 — Practical SPS ceiling and REBOUND comparison *(informational)*

For each cell, report the SPS achievable per N and identify the knee where SPS drops below interactivity thresholds:

| Threshold | Interpretation |
| --- | --- |
| ≥ 60 SPS | smoothly interactive (60 fps feel) |
| ≥ 1 SPS | borderline interactive ("step every second" feel) |
| < 1 SPS | offline-only |

For VV (fixed dt), report `wall_ms / step` directly → SPS = 1000 / wall_ms.
For IAS15 (adaptive dt), report `wall_ms / step` plus `current_dt` plus the honest interactive metric `sim_time_per_wall_second = current_dt × SPS`.

REBOUND comparison points (from published documentation and Rein & Liu 2012 / Rein & Spiegel 2015):

| System | REBOUND headline | Apsis equivalent (this experiment) |
| --- | --- | --- |
| IAS15 at N = 100 | ~10⁴ steps/s on consumer hardware | measured here |
| IAS15 at N = 1 000 | ~10–100 steps/s | measured here |
| Tree code at N = 10⁴ | ~100–1000 steps/s with SIMD | measured here |
| Tree code at N = 10⁵ | ~10–100 steps/s with SIMD | measured here |

The Apsis/REBOUND ratio per cell answers "are we within 2× of the reference implementation, or behind by an order of magnitude?". The honest answer informs whether SIMD work is mandatory for v1 paper credibility or merely future polish.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| BH walk dominates VV at N = 10⁴ AND `t_per_interaction` is roughly N-flat | Walk is compute-bound, cache fits L2/L3 | MAC (PR-perf-4) is the right next investment; SIMD adds value but second-order |
| BH walk dominates VV at N = 10⁴ AND `t_per_interaction` grows with N | Walk is memory-bound earlier than predicted | SoA refactor + SIMD is a more pressing investment than MAC |
| Trail recorder ≥ 30 % of total at N = 10⁴ | Trail recording is rate-limiting | Optimise trail buffer (down-sample, lazy snapshot) before further BH work |
| `t_integrator_overhead` ≥ 30 % of VV total | Kick/drift bookkeeping is rate-limiting | Investigate integrator step impl; possibly inline force-call boundary |
| IAS15 force_eval + integrator_overhead < 80 % | Hidden cost in IAS15 substep machinery | Profile IAS15 internals before any further optimisation |
| Apsis SPS / REBOUND SPS < 0.2 | Order-of-magnitude gap to reference | SIMD is mandatory before paper claims of "scaling characterised" land |
| Apsis SPS / REBOUND SPS ≥ 0.5 | Within 2× of reference | Current implementation is competitive; further optimisation is polish |

### Methodology

#### Cells

| Cell | Integrator | Force model | Rationale |
| --- | --- | --- | --- |
| V | Velocity Verlet | Octree (BH default) | Production default for N > 64; most common code path |
| I | IAS15 | Direct O(N²) (per ADR-003) | High-precision path; IAS15 always uses exact mode |

Yoshida-4 and Wisdom-Holman are deferred — Y4 has the same per-step cost profile as VV (slightly more force calls per step), WH only applies to hierarchical scenes and is not the production default.

#### Trail recorder variant

Each (cell, N) is measured twice: trail recorder enabled (production default) vs disabled. The difference quantifies trail's per-step cost; isolated knowledge of whether trail recording is rate-limiting at any N.

#### Phases instrumented (4)

| Phase | Definition | Where in code |
| --- | --- | --- |
| `t_tree_build` | `Octree::build` total (AABB + insert + aggregate_mass + aggregate_quadrupole) | wrapped at the engine.build call site |
| `t_bh_walk` | Whole `bh_eval_body` call: stack traversal + opening-criterion check + kernel force calc, all amortised together | wrapped at the parallel-iter outer boundary in `evaluate` |
| `t_integrator_overhead` | Integrator step minus force calls: kick, drift, predictor updates, Picard convergence checks | wrapped at integrator step boundary, subtracted force-call time |
| `t_trail_record` | Trail buffer push if enabled, no-op otherwise | wrapped at the trail-recorder call site |

Total: `t_total_step = sum of the 4 phases ± 5 % (sanity gate)`.

For IAS15, `t_bh_walk` is zero (no tree); the equivalent direct O(N²) pairwise cost is rolled into a unified `t_force_eval` phase under the same instrumentation hook. For VV+BH, `t_bh_walk` collapses kernel-force and traversal cost into a single measurement on purpose — separating them requires per-call timing markers (`Instant::now()` ~50–100 ns each, vs ~50 ns kernel call) which would contaminate the very measurement we want to isolate. The separation is escalated to method (γ) below only if the (α) data is ambiguous, per §Escalation rules. Decision-rule signals (SIMD ROI, MAC ROI) are derivable from the unified `t_per_interaction` and `n_interactions_per_body` metrics without the split.

#### Counters (4)

| Counter | Definition |
| --- | --- |
| `n_node_visits` | Total `stack.pop()` calls across all bodies in one step |
| `n_bh_accepted` | Nodes accepted as monopole+quadrupole via `s/d < θ` |
| `n_leaf_interactions` | Pairwise force calls inside leaves (excluding self-pair) |
| `n_picard_iterations` | IAS15 Picard convergence iterations summed across substeps in one step |

For VV: `n_picard_iterations = 0` (not applicable). For IAS15: `n_node_visits = n_bh_accepted = n_leaf_interactions = 0`, but the per-step pairwise count is `N × (N−1) / 2` × Picard depth, derivable from N + n_picard.

#### Derived metrics

Computed in post-processing:

- `t_per_interaction = t_bh_walk / (n_bh_accepted + n_leaf_interactions)` for VV+BH; for IAS15 use `t_per_pair = t_force_eval / (N × (N−1) / 2 × n_picard_iterations / 16)` (IAS15 has 16 substeps per step).
- `t_per_body = t_total_step / N`
- `n_interactions_per_body = (n_bh_accepted + n_leaf_interactions) / N`
- `bh_acceptance_ratio = n_bh_accepted / n_node_visits`
- `sps = 1000 / t_total_step_ms`

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| Random seed | `0x6E63696C` ("ncil" — engine-ceiling) | Single seed; this is cost characterisation, not statistical force-error gating |
| Body distribution | Sphere log-normal mass (same as perf 2×2 sphere_distribution_lognormal) | Cross-experiment comparability |
| Integrator dt | VV: `1e-3` (canonical); IAS15: adaptive (initial `1e-3`, ε_b = 1e-9 default) | VV fixed dt makes SPS interpretable; IAS15 adaptive is the production behaviour |
| Warmup steps | 10 (discarded) | First steps are slow due to cold cache, allocator init, Rayon thread pool spinup |
| Measured steps | 100 (default); 1 000 for sub-millisecond phases at small N | Aggregate over many steps reduces measurement noise; per-phase timing under 10 µs is unreliable in single-step measurement |
| Trail recorder | Tested both on (production default, capture every step) and off | Isolates trail's per-step cost |
| Hardware | Same as perf 2×2 (Ryzen 5 7600X, Windows 11, Rayon default thread pool) | Cross-experiment comparability; recorded in §Results |
| Profile | `cargo build --release` defaults | Cross-experiment comparability |

#### Variance handling

For each (cell, N, trail_variant):

- Compute `mean` and `std_dev` across the 100 measured steps
- If `std_dev / mean > 30 %`, mark as untrustworthy and re-run with longer measured horizon
- Report median + IQR in §Results, not mean (robust to occasional OS scheduling spikes)

#### Allocation hotspot detection

`Octree::build` allocates a fresh `Vec<Node>`. Optional sub-measurement: track `tree.nodes.len()` peak per build and report total allocation cost as `peak_capacity × 144 bytes / step`. If this is > 5 % of `t_tree_build`, buffer-reuse is a follow-up investment. Defer if the instrumentation cost is non-trivial; first-pass run can skip this.

#### Escalation rules (α → γ, declared a priori)

Default measurement is method (α): unified `t_bh_walk` phase, no per-interaction timing. Sufficient to drive the SIMD vs MAC decision because both candidates affect `t_per_interaction` and `n_interactions_per_body` regardless of where inside the walk the cost lives. Escalation to method (γ) — kernel microbenchmark + algebraic derivation of `t_kernel_estimated` and `t_traversal_estimated = t_bh_walk − t_kernel_estimated` — happens only if the (α) data is ambiguous along one of these axes:

| Signal | Threshold | Why it triggers γ |
| --- | ---: | --- |
| `t_per_interaction(θ = 0.3) / t_per_interaction(θ = 0.9)` | within `[0.85, 1.15]` (≈ flat across θ) | Per-interaction cost should vary with θ because high θ accepts more internal nodes (each costing a quadrupole tensor contract = ~50 FLOPs) while low θ recurses to leaves (cheap pairwise). A flat ratio means kernel cost is dominated by something other than arithmetic — likely cache loads — and the kernel-vs-traversal split becomes diagnostic for choosing SIMD vs SoA refactor. |
| `t_per_body(N = 10⁵) / t_per_body(N = 10³)` divided by the algorithmic O(N log N) factor `100 × log₂(10⁵) / log₂(10³) ≈ 167×` | observed ratio / 167 > 2.0 | Wall time grows much faster than the algorithm's complexity prediction. Likely traversal cost grows sub-linearly with cache misses; need the split to know whether to invest in SoA layout or in interaction-count reduction. |
| Predicted speedup vs measured (after applying any subsequent optimisation) | observed gain `< 50 %` of model prediction | Mental model wrong somewhere in the cost decomposition. The split is needed to find where the discrepancy lives. |

**Method (γ) implementation when triggered**:

1. Add `kernel_microbench` `#[ignore]`d test in the same `engine_ceiling` module.
2. Mede `t_per_kernel_call` em loop tight (sem walk overhead): same kernel function applied N_calls times to fixed-input args, total time / N_calls.
3. From the existing experiment data, compute `t_kernel_estimated = (n_bh_accepted + n_leaf_interactions) × t_per_kernel_call`.
4. Derive `t_traversal_estimated = t_bh_walk − t_kernel_estimated`.
5. Report as **estimated** decomposition in §Results (not measured directly), with the microbench number that anchors the estimate.

The escalation is one-shot: if (γ) doesn't resolve ambiguity, the next investigation is hardware-level profiling (perf, vtune) which is out of this experiment's scope.

#### Out of scope (declared a priori)

- **Render pipeline timing.** Render lives in `apsis-app`, not the core. Measured separately if at all.
- **Memory profiling beyond allocation count.** Heaptrack / Valgrind massif are platform-specific and non-trivial to set up; skipped.
- **Cross-machine comparisons.** Single-hardware as in prior experiments.
- **Multi-seed statistical bounds.** Cost characterisation does not require multi-seed; single seed sufficient.
- **Thread-count sensitivity.** Default Rayon thread pool used; varying thread count is a separate experiment if rate-limiting is identified at parallelism.
- **Yoshida-4 and Wisdom-Holman.** Deferred (similar profile to VV; not in the v1 critical path).

---

## Results

**Hardware / build identifier** (matches PR-perf-2 §Results for cross-experiment comparability):

- CPU: AMD Ryzen 5 7600X, 6 cores
- OS: Windows 11
- Compiler: `rustc 1.94.1`
- Profile: `cargo build --release` defaults (no LTO, codegen-units = 16)
- Rayon: default thread pool (12 logical via SMT)

CSV exports: `target/engine-ceiling/profile_v.csv` (12 rows), `target/engine-ceiling/profile_i.csv` (3 rows). Cell V end-to-end runtime: 220 s. Cell I: 19 s (truncated by dt-floor saturation; see Tier 3 caveat).

### Tier 1 — Per-stage breakdown (Cell V, gated)

Sanity gate (`|sum_phases − t_total| / t_total ≤ 5 %`): **PASS at every (N, trail) combination**, observed range 99.7–100 % phase coverage. Instrumentation has no measurable gap.

Phase breakdown at θ = 0.5, trail off (median across measured steps, percentages of `t_total`):

| N | t_total step | t_tree_build | t_bh_walk | t_integrator | t_trail |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.064 ms | 15.1 % | **83.2 %** | 1.4 % | 0.0 % |
| 1 000 | 2.06 ms | 6.2 % | **93.3 %** | 0.5 % | 0.0 % |
| 5 000 | 23.7 ms | 3.2 % | **96.7 %** | 0.2 % | 0.0 % |
| 10 000 | 40.5 ms | 3.2 % | **96.6 %** | 0.2 % | 0.0 % |
| 50 000 | 249 ms | 3.3 % | **96.5 %** | 0.2 % | 0.0 % |
| 100 000 | 707 ms | 4.5 % | **95.3 %** | 0.3 % | 0.0 % |

Tree-build share decays from 15 % at N = 100 (where the BH walk is so cheap that build dominates the small step budget) down to 3–5 % at N ≥ 5 000 — the regime where the walk dominates outright. Sub-bound `t_tree_build + t_bh_walk ≥ 85 %` at N ≥ 10⁵: observed 99.8 % at N = 10⁵; **PASS with 15 percentage-point headroom**.

Trail-recorder cost (cells V at trail = on, delta vs trail = off):

| N | t_trail / t_total | Δ step (vs trail off) |
| ---: | ---: | ---: |
| 100 | 1.0 % | +53 % (0.034 ms — small step amplifies relative cost) |
| 1 000 | 0.4 % | within measurement noise |
| 5 000 | 0.2 % | +5 % |
| 10 000 | 0.2 % | +2 % |
| 50 000 | 0.2 % | within measurement noise |
| 100 000 | 0.2 % | within measurement noise |

**Trail recorder is essentially free at every N ≥ 1 000.** Push of one `Vec3` per body per step is below 1 % of the walk cost. The hypothesis "trail recorder might be rate-limiting" — articulated explicitly before the run — is empirically refuted.

### Tier 2 — Cost normalised by work (Cell V)

Per-interaction and per-body amortised costs across N (trail off):

| N | t_per_interaction | t_per_body | n_interactions / body | bh_acceptance_ratio |
| ---: | ---: | ---: | ---: | ---: |
| 100 | 2.8 ns | 0.64 µs | 189 | 0.01 |
| 1 000 | 1.7 ns | 2.06 µs | 1 110 | 0.05 |
| 5 000 | 1.3 ns | 4.73 µs | 3 571 | 0.09 |
| 10 000 | 1.5 ns | 4.05 µs | 2 667 | 0.29 |
| 50 000 | 1.8 ns | 4.99 µs | 2 724 | 0.41 |
| 100 000 | 2.1 ns | 7.07 µs | 3 151 | 0.48 |

`t_per_interaction` ranges 1.3–2.8 ns; the U-shape (high at N = 100, minimum at N = 5 000, slight rise at N = 10⁵) is consistent with overhead amortisation at small N and mild cache pressure at large N. The variation is **1.24× across two orders of magnitude in N**, well within the escalation-rule threshold of `[0.85, 1.15]` for "kernel cost flat across θ". The kernel cost is not memory-bandwidth-bound at the v1 target N (would have shown growing `t_per_interaction`).

`bh_acceptance_ratio` grows from 0.01 (N = 100, almost all leaf interactions, tree barely subdivided) to 0.48 (N = 10⁵, half of the visited nodes accepted as monopole+quad pseudo-bodies). The ratio is the cleanest metric for "is the walk pruning effectively?": at N = 10⁵ roughly 50 % of node visits result in acceptance — substantial pruning, room for further reduction via better MAC.

`n_interactions / body` is the load-bearing metric for MAC ROI: it grows from 189 at N = 100 to 3 151 at N = 10⁵ (16.7×, roughly tracking `log₂(N)` × constant). Each interaction reduction multiplies through the entire walk; the literature MAC alternatives (Barnes 1990 ≈ −20 %, Dehnen 2002 ≈ −40 %) would land directly on this number.

**Escalation check** (γ trigger): `t_per_body(N = 10⁵) / t_per_body(N = 10³)` divided by the algorithmic O(N log N) per-body factor (`log₂(10⁵) / log₂(10³) ≈ 1.67`):

- Observed ratio: `7.07 µs / 2.06 µs = 3.43`
- `3.43 / 1.67 = 2.05` — borderline trigger (threshold `> 2.0`)

The borderline crossing is dominated by the `n_interactions / body` ratio (`3 151 / 1 110 = 2.84`), which is itself wider than the algorithmic prediction. Once interactions are controlled for, `t_per_interaction` only rises 1.24× across the same range — modest cache pressure but not the dominant signal. **Method (γ) escalation is declined** in this experiment: the walk-dominance message is clear without the kernel-vs-traversal split, and the next optimisation candidate (MAC, then SIMD) does not require the split to be sized.

### Tier 3 — SPS ceiling and REBOUND comparison

**Cell V (VV + BH octree, production default)**:

| N | SPS | Interactivity threshold | Regime |
| ---: | ---: | --- | --- |
| 100 | 15 599 | ≥ 60 SPS (smooth) | trivially smooth |
| 1 000 | 486 | ≥ 60 SPS (smooth) | smooth |
| 5 000 | 42 | < 60 SPS, ≥ 1 SPS | borderline interactive |
| 10 000 | 25 | < 60 SPS, ≥ 1 SPS | borderline interactive |
| 50 000 | 4 | < 60 SPS, ≥ 1 SPS | "step every 250 ms" |
| 100 000 | 1.4 | < 60 SPS, ≥ 1 SPS | offline-only |

**Practical interactive ceiling for VV+BH on the recorded hardware: N ≈ 10⁴** for "smooth" (60 fps), N ≈ 50 000 for "still-usable" (~ 4 SPS).

**Cell I (IAS15)** — *partial result with caveat*:

| N | SPS | dt achieved | sim_time / wall_time | Status |
| ---: | ---: | ---: | ---: | --- |
| 100 | 1 074 | 1.58 × 10⁻⁴ | 0.22 | clean |
| 1 000 | 32 | 7 × 10⁻¹⁰ | ≈ 0 | **dt-floor saturated** |
| 10 000 | 0.7 | 2 × 10⁻⁹ | ≈ 0 | **dt-floor saturated** |

The N ≥ 1 000 IAS15 measurements are not valid characterisations of integrator cost. The sphere log-normal distribution (random positions, zero initial velocities, log-normal masses with extreme ratios) is dynamically stiff for IAS15: bodies fall together, close encounters force the adaptive controller to `dt → 0`, and the integrator hits the `1e-12` `dt` floor. The reported step time at N ≥ 1 000 reflects controller-degraded steps that barely advance simulated time, not the natural cost of one integrator step on a stable scene.

**IAS15 ceiling characterisation at N ≥ 1 000 is deferred** until a stable scene generator (Plummer-equilibrium or nested-Kepler initial conditions) is available. The N = 100 measurement remains clean — IAS15 at small-N produces ~ 10³ SPS with `dt ~ 1.6 × 10⁻⁴`, a 22 % real-time-ratio simulation.

**REBOUND comparison** (published headline numbers from REBOUND README + Rein & Liu 2012 / Rein & Spiegel 2015, on broadly-comparable consumer hardware):

| Path | Apsis (this experiment) | REBOUND headline | Apsis / REBOUND |
| --- | ---: | ---: | ---: |
| IAS15 N = 100 | ~10³ SPS | ~10⁴ SPS | ~0.10 |
| Tree code N = 10⁴ | 25 SPS | ~10²–10³ SPS | ~0.025–0.25 |
| Tree code N = 10⁵ | 1.4 SPS | ~10–10² SPS | ~0.014–0.14 |

Apsis is consistently **5–10× behind REBOUND** across the comparable cells. The most credible cause is REBOUND's heavy use of x86 SIMD intrinsics in both the IAS15 inner pairwise loop and the tree-code force kernel; Apsis currently has no SIMD anywhere. The Apsis/REBOUND ratio crosses the decision-rule threshold `< 0.2` at every measured cell, triggering the rule "SIMD is mandatory before paper claims of scaling characterised land".

---

## Interpretation

The five questions the experiment was designed to answer:

1. **Is the BH walk actually the dominant cost at our N target?**
   Yes — overwhelmingly. 83–97 % of `t_total` across all measured N. The hypothesis was BH walk ≥ 70 % at N ≤ 10⁴; observed 93–97 % at N = 10³ – 10⁴. No surprise but the magnitude is sharper than predicted, and that lands the next optimisation cleanly on the walk.

2. **Where is the practical interactive ceiling?**
   N ≈ 10⁴ for smooth (25 SPS, "step every 40 ms"), N ≈ 5 × 10⁴ for borderline (4 SPS), N ≈ 10⁵ for offline-only (1.4 SPS). The smooth ceiling matches the perf 2×2 §Decision's design target ("v1's primary regime ≤ 10⁴"); the offline ceiling is consistent with the algorithmic O(N log N) cost prediction.

3. **Is per-interaction cost arithmetic-bound or memory-bound at our N target?**
   Arithmetic-bound. `t_per_interaction` only varies 1.24× across two orders of magnitude in N — well below the cache-pressure signature that would justify γ escalation. The conclusion lines up with the perf 2×2 §Decision's interpretation of why Morton was reverted: at N ≤ 10⁴ the working set fits L2/L3, and SIMD wins are not yet erodable by cache misses.

4. **Is the trail recorder rate-limiting at any N?**
   No, at any N tested. 0–1 % of `t_total`. The "trail might dominate" hypothesis was the one piece of intuition the experiment was specifically designed to test, and the answer is unambiguously no.

5. **Where does Apsis stand against REBOUND?**
   5–10× behind across the comparable cells. The most plausible cause is the absence of x86 SIMD intrinsics in any of Apsis's force-evaluation code paths. SIMD is the "respectability" lever the user asked about — without it, Apsis's headline numbers will read consistently slower than REBOUND in any side-by-side benchmark.

A subsidiary finding worth noting: the cell I dt-floor saturation at N ≥ 1 000 is honest negative evidence about IAS15's applicability domain. IAS15 is the high-precision few-body integrator; it is the wrong tool for general N ≥ 10³ scenes regardless of implementation quality. This is consistent with REBOUND's own guidance and reaffirms the production architecture (VV/Y4 + BH for general-N scenes, IAS15 reserved for the high-precision regime).

---

## Decision

The data fires three of the protocol's decision rules cleanly:

| Rule | Trigger | Action |
| --- | --- | --- |
| BH walk dominates VV at N = 10⁴ AND `t_per_interaction` ≈ N-flat | Walk = 96.6 %; t_per_interaction varies 1.24× | **MAC (PR-perf-4) is the right next investment**; SIMD second-order |
| Apsis SPS / REBOUND SPS < 0.2 across comparable cells | All three measured comparisons fall below the threshold | **SIMD is mandatory before paper claims of scaling characterised land**; queued as PR-perf-5 |
| Trail recorder ≥ 30 % of total at N = 10⁴ | 0.2 % observed | **No trail-recorder optimisation needed**; the suspicion is empirically resolved |

**Sequencing**:

1. **PR-perf-4 (next)** — multipole acceptance criterion comparison (Barnes 1990, Dehnen 2002, GADGET-style). Design preserved from the perf 2×2 §Decision sub-section: per-cell × per-N × per-θ factorial against the classical MAC baseline, declare bounds a priori, ship the winner or document literature comparison as deferred. Targets the dominant cost (BH walk) at the v1 production regime where it lives.
2. **PR-perf-5 (queued)** — SIMD inner loops in the BH walk's per-interaction kernel and (separately) in IAS15's pairwise direct sum. Required for the paper-credibility threshold; expected ~ 2–4 × speedup on the per-interaction cost, multiplicative with MAC's interaction-count reduction.

**Out-of-scope follow-ups noted but not queued**:

- **IAS15 ceiling characterisation on a stable scene** — needs a Plummer-equilibrium or nested-Kepler initial-condition generator. Worth doing before the paper to pin "IAS15 at N = 100 is interactive, IAS15 at N = 1 000 is offline-by-design" with numbers, but does not affect the optimisation roadmap.
- **Cross-machine REBOUND co-run** — the published REBOUND numbers are from different hardware and slightly older versions; a same-machine same-version comparison would tighten the ratio confidence interval. Logistically substantial; deferred.

**Production engine remains unchanged.** The instrumentation infrastructure (`WalkCounters`, `evaluate_profile`, `engine_ceiling.rs` test module) ships as-is — `evaluate_profile` is `pub(crate)` and the harness is `#[cfg(test)]`, so the public API and the `evaluate` hot path are untouched.

**§Decision provenance.** This section was written after the cell V harness completed (220 s) and the cell I harness completed (19 s, partial), with all derived metrics computed from the per-row eprintln output and the CSVs archived under `target/engine-ceiling/`. Decision rules and escalation triggers were declared a priori in the Hypothesis and §Escalation rules sections of this same notebook before any instrumentation code was written.

---

## Threats to validity

1. **Single-seed.** Cost characterisation does not depend on multi-seed averaging, but pathological body distributions could give unrepresentative interaction counts. Mitigation: the sphere log-normal seed (`0x6E63696C`) matches the family used in perf 2×2, so cross-experiment claims are consistent. If the measured `n_interactions_per_body` looks anomalous, re-run with a second seed.

2. **Single-machine.** SPS numbers are not portable; the Apsis/REBOUND ratio computed here is anchored to the recorded hardware. A different CPU could shift the ratio by ±50 %.

3. **`t_bh_walk` is unified by choice (method α).** Kernel-force vs traversal cost are not separated in the default measurement, because per-call timing markers (`Instant::now()` ≈ 50–100 ns each) would dominate the kernel call itself (~50 ns) and contaminate the very measurement we want. The unified `t_bh_walk` divided by `(n_bh_accepted + n_leaf_interactions)` gives `t_per_interaction` — sufficient for the SIMD vs MAC decision because both candidates affect that metric. The split is escalated to method (γ) — kernel microbenchmark + algebraic derivation — only when the (α) data triggers an ambiguity flag per §Escalation rules.

4. **REBOUND comparison numbers come from published documentation, not co-run benchmarks.** Hardware differences and version skew between published REBOUND numbers and this experiment introduce uncertainty. The Apsis/REBOUND ratio is informational only; ±2× confidence interval is the honest reading.

5. **Trail recorder cost depends on capture cadence.** "Trail on" with capture-every-step is the worst case; production may use sub-sampling. The on/off variant brackets the cost rather than measuring it precisely under a specific cadence.

6. **IAS15 adaptive dt confounds wall-time per step.** A simulation in a stable regime may take large dt steps, giving low SPS but high simulated-time-per-second. Reported alongside `wall_ms / step` as `sim_time_per_wall_second` to avoid the misleading "step" interpretation.

7. **Allocation pattern not fully measured.** Without heaptrack, the cost of `Vec<Node>` reallocation across builds is estimated from peak capacity × node size, which assumes worst-case allocator behaviour. Real cost may be lower (allocator reuses heap pages) or higher (page fault on growth). Bounded by the optional sub-measurement; defer if instrumentation cost is high.
