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

*To be populated after the experiment runs.*

### Tier 1 — Per-stage breakdown

*Pending.*

### Tier 2 — Cost normalised by work

*Pending.*

### Tier 3 — SPS ceiling and REBOUND comparison

*Pending.*

---

## Interpretation

*To be written after §Results is populated.*

---

## Decision

*To be written after §Interpretation; identifies the dominant bottleneck per regime and the next optimisation investment with measured justification.*

---

## Threats to validity

1. **Single-seed.** Cost characterisation does not depend on multi-seed averaging, but pathological body distributions could give unrepresentative interaction counts. Mitigation: the sphere log-normal seed (`0x6E63696C`) matches the family used in perf 2×2, so cross-experiment claims are consistent. If the measured `n_interactions_per_body` looks anomalous, re-run with a second seed.

2. **Single-machine.** SPS numbers are not portable; the Apsis/REBOUND ratio computed here is anchored to the recorded hardware. A different CPU could shift the ratio by ±50 %.

3. **`t_bh_walk` is unified by choice (method α).** Kernel-force vs traversal cost are not separated in the default measurement, because per-call timing markers (`Instant::now()` ≈ 50–100 ns each) would dominate the kernel call itself (~50 ns) and contaminate the very measurement we want. The unified `t_bh_walk` divided by `(n_bh_accepted + n_leaf_interactions)` gives `t_per_interaction` — sufficient for the SIMD vs MAC decision because both candidates affect that metric. The split is escalated to method (γ) — kernel microbenchmark + algebraic derivation — only when the (α) data triggers an ambiguity flag per §Escalation rules.

4. **REBOUND comparison numbers come from published documentation, not co-run benchmarks.** Hardware differences and version skew between published REBOUND numbers and this experiment introduce uncertainty. The Apsis/REBOUND ratio is informational only; ±2× confidence interval is the honest reading.

5. **Trail recorder cost depends on capture cadence.** "Trail on" with capture-every-step is the worst case; production may use sub-sampling. The on/off variant brackets the cost rather than measuring it precisely under a specific cadence.

6. **IAS15 adaptive dt confounds wall-time per step.** A simulation in a stable regime may take large dt steps, giving low SPS but high simulated-time-per-second. Reported alongside `wall_ms / step` as `sim_time_per_wall_second` to avoid the misleading "step" interpretation.

7. **Allocation pattern not fully measured.** Without heaptrack, the cost of `Vec<Node>` reallocation across builds is estimated from peak capacity × node size, which assumes worst-case allocator behaviour. Real cost may be lower (allocator reuses heap pages) or higher (page fault on growth). Bounded by the optional sub-measurement; defer if instrumentation cost is high.
