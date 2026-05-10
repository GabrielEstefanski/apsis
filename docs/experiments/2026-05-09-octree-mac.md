# Octree multipole acceptance criterion (MAC) comparison — protocol

**Date:** 2026-05-09
**Subject:** Compare three multipole acceptance criteria from the literature against the production classical `s/d < θ` baseline. Identify whether a tighter or adaptive opening criterion reduces interaction count enough at fixed accuracy to ship as the new production default.

**Status:** Protocol declared a priori, before any MAC implementation lands. §Results populated incrementally; §Decision written after enough cells measure to make the call.

**Branch:** `perf/octree-mac`, from `develop` (PR #73 merged carrying the engine ceiling §Decision that names this experiment).

---

## Abstract

The engine ceiling experiment (`2026-05-09-engine-ceiling.md`, §Decision) classified the engine as **interaction-bound** at the v1 target regime: `t_per_interaction` is already low (1.3–2.1 ns on the recorded hardware) and the dominant cost is the *number* of interactions the BH walk performs (1 110 / body at N = 10³ growing to 3 151 / body at N = 10⁵). The natural compute-bound attack is to reduce that interaction count without sacrificing accuracy, which is exactly what a tighter or adaptive multipole acceptance criterion does.

This experiment runs four MAC cells against the same body distributions, θ values, and N range used in perf 2×2 (`2026-05-08-octree-perf-2x2.md`), holding everything except the opening criterion fixed:

| Cell | MAC | State per node | Opening test |
| --- | --- | --- | --- |
| **M0** | classical (production) | none | `s / d < θ` |
| **M1** | Barnes 1990 | `δ_max` (max body offset from COM) | `s / (d − δ_max) < θ` |
| **M2** | Dehnen 2002 | per-walk error accumulator | predicted multipole error contribution / accumulated error budget |
| **M3** | GADGET-style (Springel 2005 §2.4) | per-walk acceleration accumulator + node `a_internal_scale` | Barnes 1990 + relative-acceleration error budget |

The expected ranking by literature is M2 ≥ M3 > M1 > M0 in interaction-count reduction at fixed accuracy. The expected ranking by implementation cost is the inverse: M1 cheapest (one extra `f64` per node, modified opening test), M3 intermediate (per-walk state), M2 most complex (walk-state with derivable error bound).

The acceptance gates are organised in three tiers: per-cell force-accuracy preserved at the perf 2×2 Tier 1 percentile bounds (gated), wall-time at matched accuracy versus M0 (gated against literature-derived ranges), interaction-count reduction at fixed θ (informational; the structural lever the experiment is testing).

---

## Motivation

The engine ceiling §Decision queued this MAC comparison as the first attack on the interaction-bound regime. The reasoning chain:

1. Engine ceiling measured `n_interactions / body` = 1 110 (N = 10³) → 2 667 (N = 10⁴) → 3 151 (N = 10⁵). Walk cost is dominated by this volume, not by per-interaction expense.
2. Each MAC controls how many of those interactions are *necessary* to maintain the target accuracy. A tighter criterion accepts more cells as monopole+quadrupole pseudo-bodies and recurses less, lowering the count.
3. Literature reports 15–50 % interaction reduction at matched accuracy across the alternative MACs, multiplicative with later SIMD and SoA-layout optimisations. None of those alternatives are in production today.

The classical `s / d < θ` is the cheapest possible MAC — a single division, no per-node state — but it does not account for two known correction factors: (a) a node's body distribution can extend past its COM by a margin `δ_max`, making "distance from body to COM" overstate the actual closest approach; (b) a node's contribution to a body's acceleration can be much smaller than `M / r²` would predict if the body's accumulated force is already large, meaning the relative-error budget for accepting that node is wider than the geometric criterion captures.

Barnes 1990 addresses (a). Dehnen 2002 addresses (b) directly with an error-driven criterion. Springel 2005 (GADGET-2's MAC) combines both with a relative-acceleration tolerance.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the four MAC cells under test, the metrics declared below are bounded a priori at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 is gated; Tier 2 ranges are gated against literature-bound intervals (not point values); Tier 3 is informational. Bound revision is forbidden unless backed by concrete arithmetic. Literature comparison ranges accommodate the empirical-not-derived nature of the cited literature numbers.

#### Tier 1 — Force accuracy preserved at perf 2×2 bounds *(gated; per cell)*

Same body distribution, seeds, and θ as perf 2×2 §Tier 1 (sphere log-normal, N ∈ {1 000, 10 000}, three seeds `0x6F637472`, `0x71756164`, `0x6D6F7274`, θ = 0.5). Per-body acceleration error measured against an independent O(N²) reference (the `exact_pairwise_forces` path established in perf 2×2 commit `d44dfda`).

| Cell | N | Bound p50 | Bound p95 |
| --- | ---: | ---: | ---: |
| M0, M1, M2, M3 | 1 000 | `≤ 1 × 10⁻³` | `≤ 5 × 10⁻³` |
| M0, M1, M2, M3 | 10 000 | `≤ 1 × 10⁻³` | `≤ 5 × 10⁻³` |

Bounds match the perf 2×2 Tier 1 bounds for cell C (quadrupole-on, current production); MAC alternatives must preserve accuracy at least as well as M0 at the same θ. **A MAC variant that lowers wall-time but degrades p95 above 5 × 10⁻³ is rejected** — accuracy floor is not negotiable for the ship decision.

For Dehnen-style MACs (M2, M3), the per-walk error budget tolerance is tuned to land within the same p95 band as M0; that tolerance is reported in §Methodology and locked before the wall-time measurement.

#### Tier 2 — Wall-time at matched accuracy *(gated as ranges; literature-referenced)*

The candidate MACs are **strictly more conservative than M0 at the same θ** (M1's `s / (d − δ_max) < θ` recurses earlier than M0's `s / d < θ`; M2 / M3 use the same δ_max correction plus an error / acceleration budget). Comparing them at fixed θ would show the candidates *slower* and *more accurate* than M0, missing the structural lever — the candidate's strictness is what lets it accept more aggressively at a *higher* θ for the same accuracy budget.

The honest comparison is therefore **at matched p95 accuracy**: pick θ_match for each candidate such that the candidate's p95 at θ_match equals M0's p95 at θ = 0.5 (the production canonical), then compare wall times at those θ values. Same matched-accuracy pattern from `2026-05-08-octree-perf-2x2.md` §Tier 2 (where C / D quadrupole speedup was measured at matched accuracy via θ binary search).

| Comparison | Range bound | Reference |
| --- | --- | --- |
| `t_eval_M1(θ_match) / t_eval_M0(0.5)` at N = 10⁴ | ∈ [0.75, 0.95] | Barnes 1990 §3.2: 5–25 % speedup at matched accuracy |
| `t_eval_M2(θ_match) / t_eval_M0(0.5)` at N = 10⁴ | ∈ [0.50, 0.80] | Dehnen 2002 §5: 20–50 % speedup at matched accuracy |
| `t_eval_M3(θ_match) / t_eval_M0(0.5)` at N = 10⁴ | ∈ [0.60, 0.85] | Springel 2005 §2.4 / GADGET-2: 15–40 % speedup |

θ_match is found by binary search over θ ∈ [0.5, 1.0], target tolerance ±0.01, accepted when `|p95_candidate(θ_match) − p95_M0(0.5)| / p95_M0(0.5) ≤ 0.05` (5 % p95 agreement). The search runs once per (candidate, N, seed) before the timed measurement; the resulting θ_match is locked and reported alongside the wall-time ratio.

A measurement outside its range is investigated; a measurement at the edge is reported with the discrepancy. Decision-rule "ship" requires the cell to land inside its range AND for the ratio to be lower than at least one cheaper cell (otherwise the cheaper cell wins on parsimony grounds).

#### Tier 3 — Interaction-count reduction at fixed θ *(informational)*

For each cell at θ = 0.5, N ∈ {1 000, 10 000, 100 000}, report:

- `n_interactions_per_body` (median across seeds)
- `n_bh_accepted / n_node_visits` (acceptance ratio)
- Reduction ratio vs M0: `n_int_M_k / n_int_M_0` per cell

This is the structural lever the experiment is testing. A MAC that reduces wall-time without reducing interactions is suspect (likely an instrumentation defect or a different bottleneck shifting). A MAC that reduces interactions but not wall-time has hit a sub-bottleneck — investigate before the ship decision.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| M1 lands inside [0.75, 0.95] AND preserves Tier 1 | Barnes 1990 captures meaningful speedup at minimal complexity | **Ship M1 as the new production default**; consider M2/M3 only if the residual gap to REBOUND demands it later |
| M1 lands at upper edge (0.92–0.95) AND Tier 1 preserved | Barnes 1990 marginal at our distribution; leave on table | Implement and measure M3 (intermediate complexity); revisit decision |
| M2 or M3 lands inside its range AND ≥ 5 % better than M1 at matched accuracy | The structural Dehnen/GADGET refinement pays off beyond Barnes 1990 | **Ship the cheaper-by-complexity winner** between M2 and M3; document the other for posterity |
| Any cell exceeds Tier 1 p95 bound | Implementation defect or MAC over-relaxes accuracy | Investigate per-body error distribution; tune tolerance for Dehnen-style MACs; never relax the Tier 1 bound |
| Tier 3 shows interaction reduction but Tier 2 wall-time gain < 50 % of expected ratio range | Walk has a sub-bottleneck eating the savings (likely traversal cost growing with depth) | Escalate to per-walk node-visit instrumentation before deciding; possibly defer MAC pending SoA / SIMD work |
| All cells measured and none lands inside its range | The literature ranges are wrong for our regime — possible distribution-specific effect | Document honestly, defer ship until investigated; do not relax bounds to fit |

### Methodology

#### Implementation order: cheapest first, escalate by ROI

1. **M0 baseline measurement** — current production code, no changes. One measurement pass to establish the ratios.
2. **M1 (Barnes 1990)** — implement, measure, decide. If M1 captures the upper end of its expected range and Tier 1 holds, the experiment may stop here pending §Decision logic.
3. **M2 / M3** — implemented only if M1 does not capture enough of the literature-predicted gain. Order between M2 and M3 to be decided based on which one's tolerance-tuning task looks more tractable when the work is reached.

The incremental escalation prevents over-investing in complex MACs when a cheaper one already lands inside the ship-decision range.

#### Cell M0 (classical, production)

No code change. Measurement uses the existing `BarnesHutEngine::evaluate_profile` path.

#### Cell M1 (Barnes 1990)

`Node` gains one `f64` field `delta_max`: the maximum distance from any body in the subtree to the node's COM. Computed during `aggregate_mass` (single bottom-up pass, cheap):

- Leaf: `delta_max = max over body in leaf of |body.pos - leaf.com|`.
- Internal: `delta_max = max over children c of (|c.com - parent.com| + c.delta_max)`. (Triangle inequality bound; tight enough for practical use.)

Opening criterion in `bh_eval_body`:

```text
M0:  s / d < θ
M1:  s / (d - delta_max) < θ
```

with a guard `d - delta_max > 0` (recurse if not — the body is inside the node's bounding sphere, no acceptance possible).

Implementation cost: ~60 LOC, no walk-state changes, no API surface. Per-node memory: +8 bytes (`f64`). Production engine field-bumped from ~144 bytes to ~152 bytes per node — small, well within the perf 2×2 layout budget.

#### Cell M2 (Dehnen 2002)

Multipole acceptance based on a derivable bound for the contribution of an unaccepted node to the body's accumulated force error. For a node with multipole moments at distance `d`, the error contribution scales as approximately `|Q| / d⁴` (quadrupole next-order) plus higher-order terms; Dehnen gives a closed-form bound.

Per-walk state: an accumulator `eps_acc` tracking the body's current force-error budget. Opening criterion: accept node if its predicted error contribution `e_node ≤ tolerance × |a_accumulated|` (or similar — exact form to be locked in implementation, with the citation reference).

Tolerance parameter: tuned per Tier 1 to land within the p95 ≤ 5 × 10⁻³ band, then locked for Tier 2 measurement.

Implementation cost: substantive — walk-state accumulator, tolerance tuning loop, possibly a separate test that the tolerance landed in the right p95 band. ~200 LOC + tuning experiment. Defer until M1 / M3 don't suffice.

#### Cell M3 (GADGET-style)

Combines Barnes 1990's geometric correction with an acceleration-relative tolerance:

```text
M3:  (s / (d - delta_max) < θ) AND (a_internal_scale / |a_accumulated| < eta)
```

where `a_internal_scale` is a per-node estimate of "the force a body just outside this node would feel from the node's mass" and `eta` is the relative-acceleration tolerance (per Springel 2005 ≈ 0.005 for production GADGET-2 runs).

Per-walk state: same `a_accumulated` accumulator as M2.

Implementation cost: intermediate — `a_internal_scale` is per-node (one additional `f64`, computed in aggregate_mass); the `a_accumulated` walk-state mirrors M2's.

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| Random seeds | 3: `0x6F637472`, `0x71756164`, `0x6D6F7274` | Match perf 2×2 for cross-experiment comparability |
| Body distribution | sphere log-normal mass | Match perf 2×2 family |
| `θ` (M0, M1) | 0.5 | Production canonical |
| Tolerance (M2, M3) | tuned per Tier 1 to match M0 p95 | Locked before Tier 2 measurement |
| N (Tier 1 + Tier 2) | 1 000, 10 000 | Reference O(N²) feasible at both |
| N (Tier 3 informational) | 1 000, 10 000, 100 000 | Sampled reference at N = 10⁵ (K = 512, same as perf 2×2) |
| Warmup / measured | 1 / 10 (per perf 2×2 convention) | Robust to OS scheduling jitter |
| Multipole order | Quadrupole always-on | Per perf 2×2 §Decision |
| Hardware | Same as perf 2×2 / engine ceiling (Ryzen 5 7600X, Windows 11) | Cross-experiment comparability |

#### Out of scope (declared a priori)

- **Adaptive θ controller** — the existing `ThetaController` infrastructure is not part of this experiment; θ stays fixed at 0.5 across all cells. Adaptive opening per-body is a different axis from MAC selection.
- **Higher-order multipoles (p ≥ 3)** — quadrupole-only per perf 2×2 §Decision; octupole and Dehnen-FMM proper are post-MAC investigations.
- **SIMD / SoA layout** — separately queued as PR-perf-5/6.
- **Cross-machine comparison** — single-hardware as in prior experiments.
- **MAC interaction with Morton sortation** — Morton was reverted in perf 2×2; engine ceiling §Decision queues a re-evaluation after SoA + SIMD land. MAC × Morton cross-product is post-PR-perf-7 if both still relevant.

---

## Results

*To be populated incrementally as cells are implemented and measured.*

### Tier 1 — Force accuracy preserved

*Pending.*

### Tier 2 — Wall-time at matched accuracy

*Pending.*

### Tier 3 — Interaction-count reduction at fixed θ

*Pending.*

---

## Interpretation

*To be written after enough cells land to make a structural reading.*

---

## Decision

*To be written after Tier 2 ranges are populated and the decision rules can fire.*

---

## Threats to validity

1. **Single distribution family.** Sphere log-normal is the perf 2×2 canonical and provides cross-experiment comparability, but a different distribution (Plummer cluster, hot disk, hierarchical binary) could produce different `δ_max` values, different `a_internal_scale` distributions, and different cell rankings. Mitigation: report the seed and distribution explicitly; flag any cell whose ranking depends on this choice.

2. **Tolerance tuning loop for M2 / M3 risks circular validation.** If we tune the tolerance to land at p95 = 5 × 10⁻³ exactly, the Tier 1 bound passes by construction. Mitigation: tolerance tuning targets `p95 ≤ 0.5 × bound = 2.5 × 10⁻³` (half-bound headroom), so the Tier 1 bound passes with margin and the wall-time measurement is taken under genuinely accuracy-preserving conditions.

3. **Walk-state cost (M2, M3) not measured directly.** The per-walk `eps_acc` / `a_accumulated` accumulator adds a few f64 ops per walk step. If walks are very short (small N), this overhead can erode the interaction-count gain. Mitigation: report `t_per_interaction` per cell — if it grows materially in M2/M3 vs M0/M1, the walk-state cost is showing up and needs separation.

4. **Literature comparison ranges are empirical, not derived.** Cited speedups come from REBOUND / GADGET / falcON measurements on different distributions, hardware, and code-base maturity. The ranges in Tier 2 carry ±50 % confidence; "outside the range" is a flag for investigation, not automatic rejection.

5. **MAC ranking can flip with N.** Barnes 1990's `δ_max` correction matters most when nodes have wide-spread bodies (deep trees, large N). Dehnen / GADGET-style accumulator-based criteria matter most when accumulated force is large (also deep trees, large N). Small-N cells (N = 10³) may show all four MACs converging; the interesting differentiation is at N ≥ 10⁴.

6. **Production engine continues to use M0 throughout the experiment.** The MAC variants live behind `pub(crate)` runtime knobs (similar pattern to the perf 2×2 multipole and Morton toggles). The final commit removes all knobs and bakes the chosen winner; if no winner emerges, M0 stays as production and the alternatives are documented as deferred.
