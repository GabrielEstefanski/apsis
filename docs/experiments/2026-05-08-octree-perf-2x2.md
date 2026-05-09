# Octree perf — quadrupole + Morton 2×2

**Date:** 2026-05-08
**Subject:** Establish what makes the basic Barnes-Hut octree fast at the project's target scale (N ∈ [10³, 10⁵]) by adding (a) symmetric traceless quadrupole expansion at internal nodes and (b) Morton (Z-order) body sort. Measure both factors in a 2×2 factorial (mono/quad × no-Morton/Morton) at matched accuracy across N and θ. Decide which factors ship as always-on baked-in based on observed cost-vs-error frontier, not theory.

**Status:** Protocol declared a priori, before any code is written. §Results populated incrementally — PR-perf-1 fills cells A and C; PR-perf-2 fills cells B and D and writes §Decision.

**Branch base:** PR-perf-1 stacked on `feat/octree-port` (PR #69, open at protocol time). Both PRs rebase onto `master` once #69 lands.

---

## Abstract

Barnes-Hut force evaluation cost depends on three first-order knobs:

1. the multipole order to which each node's gravitational signature is expanded,
2. the spatial layout of body insertion / walk traversal,
3. the opening criterion θ that decides when a cell is monopole-acceptable.

This experiment treats (1) and (2) as two binary factors and (3) as the swept axis of the cost-vs-error frontier. The 2×2 factorial design separates the contribution of each factor from confounds:

|  | Monopole only | Monopole + quadrupole |
| --- | --- | --- |
| **Body-array order (no Morton)** | A — current production (PR #69 baseline) | C — quadrupole isolated |
| **Morton-sorted insert + walk** | B — Morton isolated | D — combined production candidate |

Decision rule: at any matched per-body force accuracy on the (cost, error) Pareto frontier, the configuration with strictly lower wall time at N = 10⁵ ships as always-on. If multiple configurations dominate at different N regimes, the crossover N is documented and the engine baked-in to the dominant configuration of the project's primary regime (10⁴–10⁵).

The acceptance gates are organised in three tiers: per-body force agreement against an independent O(N²) reference (Tier 1, gated against Salmon-Warren and Hernquist-Katz error bounds), wall-time at matched accuracy (Tier 2, gated against literature-bound ranges from Hernquist & Katz 1989, Dehnen 2002, Springel 2005), and cache-effect characterisation across N (Tier 3, informational).

---

## Motivation

The octree-port (PR #69) §Results Tier 3 measured wall time scaling with the basic monopole octree. Two limitations were left as follow-up:

1. **Per-body error at θ = 0.5 sits at ~5%** (Salmon-Warren bound). Acceptable for production at N ≤ 10⁴ but pinches downstream perturbation work — a 5% force noise floor masks first-order PN corrections of similar magnitude in close encounters.

2. **Worst-case N-doubling wall-time ratio of 3.43×** at the 500 → 1000 step, against ≈ 2.3× expected from O(N log N). The signature suggests cache-locality loss when the working set crosses L2 / L3 thresholds — testable but not tested.

The quadrupole factor addresses (1) directly: literature reports ≈ 10× per-body error reduction at fixed θ, equivalent to using a relaxed θ ∈ [0.7, 0.9] for the same accuracy as monopole at θ = 0.5. The relaxed θ accepts more cells as monopole-acceptable and reduces interaction count by a factor that more than compensates the per-interaction cost increase (~2.5–3× monopole due to tensor contraction).

The Morton factor addresses (2) by ordering body insertion and walk traversal so spatially adjacent bodies process consecutively. The walk for body i + 1 finds tree nodes warm in cache from body i's walk.

The 2×2 design is critical for two reasons:

- **Cache pressure interaction.** Quadrupole increases per-node memory by ~28 % (5 independent traceless tensor components added to the existing ~144-byte node). At N where the node array straddles L2, this can amplify cache penalty — exactly the regime Morton mitigates. Whether the combination is super-additive (Morton recovers quadrupole-amplified cache loss + delivers its own gain) or sub-additive (quadrupole's memory pressure overflows what Morton can hide) is an empirical question, not derivable from individual-factor analysis.

- **Honest attribution.** Without isolating each factor, a measured speedup of the full quad+Morton configuration cannot be apportioned. Shipping a single-PR combination that wins overall could mask a regression in one factor that the other compensates — bad for the codebase's long-term maintainability.

PR-perf-1 lands quadrupole isolated (cells A, C with toggle for Morton off), letting C-vs-A be reviewed and gated against literature before adding the Morton confound. PR-perf-2 lands Morton (cells B, D), populates the full table, and writes §Decision. This staging is what makes the 2×2 attribution honest at PR-review granularity.

Both factors are scoped as research-validated optimisations; neither is a refactor pretending to be perf, neither carries hidden semantic changes. Bodies arrive in their original order, accelerations are written back in the original order; no public API change.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the four configurations under test (cells A, B, C, D), the metrics declared below are bounded a priori at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are gated; failure of any gated metric reproves the affected cell and implementation is fixed before merge — bound revision is forbidden unless backed by concrete arithmetic (FP-reorder floor) or by re-derivation of the literature estimate from interaction-count arithmetic. Literature comparison bounds use **ranges**, not point values, to acknowledge that the references are empirical reports, not derived bounds. Tier 3 is informational and never reproves.

#### Tier 1 — Force accuracy *(gated; per cell + cross-cell consistency)*

Random-sphere distribution: N ∈ {1000, 10000} bodies sampled uniformly inside a unit sphere with masses drawn from `LogNormal(μ = 0, σ = 1)`. Three seeds: `0x6F637472`, `0x71756164`, `0x6D6F7274` ("octr", "quad", "mort"). Per-body acceleration error is measured against an independent O(N²) reference computed at the same body distribution.

**Per-cell error bounds** at θ = 0.5:

| Cell | Configuration | Bound | Origin |
| --- | --- | ---: | --- |
| A | mono, no Morton | `max_i \|Δa_i\| / \|a_i^exact\| ≤ 5 × 10⁻²` | Salmon-Warren 1994 monopole-at-θ=0.5 bound; baseline from octree-port §Results |
| B | mono, Morton | `≤ 5 × 10⁻²` | Same as A; Morton is permutation-invariant for the FINAL tree structure, so accuracy must not change |
| C | quad, no Morton | `≤ 5 × 10⁻³` | Hernquist & Katz 1989 quadrupole-at-θ=0.5; ≈ 10× monopole improvement |
| D | quad, Morton | `≤ 5 × 10⁻³` | Same as C; Morton is permutation-invariant |

**Cross-cell consistency** (Morton-on cells must produce identical forces to their no-Morton counterparts up to FP-reorder floor):

| Comparison | Bound | Origin |
| --- | ---: | --- |
| `max_i \|a_B(i) − a_A(i)\| / \|a_A(i)\|` at θ = 0.5, N = 1000 | `≤ 1 × 10⁻¹²` | Per-body BH walk visits ≈ 8 · log₂(N) ≈ 80 nodes; relative drift bound `80 · ε ≈ 1.8 × 10⁻¹⁴`; with N = 1000 leaf accumulation extending it to ≈ 2 × 10⁻¹³; 5× headroom |
| `max_i \|a_D(i) − a_C(i)\| / \|a_C(i)\|` at θ = 0.5, N = 1000 | `≤ 1 × 10⁻¹²` | Same FP-reorder floor; Morton must change the order of force computation, never the computed forces |

#### Tier 2 — Wall-time at matched accuracy *(gated as ranges; literature-referenced)*

The honest comparison is **at matched per-body accuracy**: pick θ_quad such that quadrupole at θ_quad has the same `max_i error_per_body` as monopole at θ = 0.5; then compare wall times at that θ_quad.

θ_match for each cell is determined by binary search over θ ∈ [0.5, 1.0], target tolerance ±0.01 on θ, accepted when `|error_quad(θ_match) − error_mono(0.5)| / error_mono(0.5) ≤ 0.05`.

Wall-time bounds are stated as **ranges**. Literature references (Hernquist & Katz 1989, Dehnen 2002 §5, Springel 2005 §2.4) are empirical reports; a measurement outside the range is investigated, a measurement at the edge is reported with the discrepancy.

| Comparison | Range bound | Reference |
| --- | --- | --- |
| `t_eval_C(θ_match_C) / t_eval_A(0.5)` at N = 10⁴ | ∈ [0.30, 0.70] | Quadrupole-at-matched-accuracy speedup vs monopole; Dehnen 2002 §5 Table 1 (≈ 2–3× faster); Springel 2005 §2.4 (≈ 2× in GADGET-2) |
| `t_eval_B(0.5) / t_eval_A(0.5)` at N = 10⁴ | ∈ [0.40, 0.95] | Morton spatial-sort speedup; literature reports 1.5–3× (Springel 2005 §5.2; Wang et al. 2018); 0.95 upper accommodates work-stealing fragmentation; 0.40 lower is the optimistic 2.5× regime |
| `t_eval_D(θ_match_D) / t_eval_C(θ_match_C)` at N = 10⁴ | ≤ 0.90 | The decisive D-vs-C bound: Morton's contribution **on top of quadrupole**. Floor of 10 % gain at matched accuracy required to justify Morton always-on |
| `t_eval_D(θ_match_D) / t_eval_A(0.5)` at N = 10⁵ | ≤ 0.30 | Combined production-target gain at the size class that motivates this work |

**Build vs walk decomposition** (informational, not gated, but required reporting):

| Comparison | Expected sign | Why the decomposition matters |
| --- | --- | --- |
| `t_build_B(N) / t_build_A(N)` | ≤ 1.0 (Morton helps build) | Morton orders the `Vec<Node>` growth pattern; build cache benefit |
| `t_build_C(N) / t_build_A(N)` | ≥ 1.0 (quadrupole adds tensor cost) | `aggregate_mass` computes traceless tensor at every node |
| `t_walk_C(θ_match_C) / t_walk_A(0.5)` | < 1.0 (the actual quadrupole win) | Walk visits fewer nodes at relaxed θ; net reduction even with ~2.5–3× per-interaction cost |

The decomposition allows separating "build slower but walk faster" from "both slower" — the latter would indicate the implementation has a defect, not a tradeoff.

Reading the user's framing: **Morton attacks build + walk; quadrupole only alters walk** (it shrinks the interaction count if θ rises but makes each interaction more expensive). Without the decomposition the reader cannot tell which knob did the work.

#### Tier 3 — Cache-effect characterisation *(informational, NOT gated)*

Re-run the octree-port Tier 3 wall-time table (N ∈ {100, 250, 500, 1000, 2500}) for all four cells at θ = 0.5. Report:

- Per-N median evaluate wall time (4 cells × 5 N values = 20 measurements per seed).
- Worst N-doubling time ratio per cell.
- Empirical exponent of `t(N) = c · N^k` from log-log regression, per cell.

Expected: the worst N-doubling ratio drops below the baseline 3.43× for cells with Morton on (B, D); empirical exponent moves toward the theoretical 1.0–1.2 range. A measurement showing Morton-on cells still ≥ 3.0 indicates the cache effect was misdiagnosed (different bottleneck — likely thread contention on `accels` writes), warranting separate investigation.

#### Pareto frontier *(reported, not gated)*

For each of the 4 cells, sweep θ ∈ {0.3, 0.5, 0.7, 0.9} at N ∈ {1000, 10000, 100000}. Plot (max-per-body-error, t_eval) per (cell, θ, N). The shipped configuration is the cell whose frontier dominates at the project's target accuracy budget (5 × 10⁻³ per-body error, the "PN-friendly" floor implied by the perturbation framework's 1PN test scenarios).

CSV export: `target/perf-2x2/octree_pareto_<seed>.csv` with columns `cell,N,theta,seed,max_rel_err,t_build_ms,t_walk_ms,t_eval_ms,std_err_t_eval`. Plotting deferred to a separate `docs/experiments/2026-05-08-octree-perf-frontier.py` script alongside the CSV.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| All Tier 1 + Tier 2 ranges pass; D-vs-C ≤ 0.90 at N = 10⁴ AND Morton gain grows with N | Quadrupole hits literature; Morton contributes ≥ 10 % under quad memory pressure; gain scales | Both ship always-on; toggle removed in PR-perf-2 final commit; quadrupole becomes default mode; cells A and B documented as historical reference |
| All Tier 1 pass; quadrupole hits range; Morton gain D-vs-C < 0.90 at N = 10⁴ but unchanged at N = 10⁵ | Quadrupole alone delivers; Morton marginal at our scales | Quadrupole ships; Morton reverted; documented as deferred (revisit if N target moves to 10⁶) |
| Quadrupole below range (`t_eval_C / t_eval_A > 0.70`) | Implementation defect or cache pressure eroding theoretical gain | Investigate: (a) profile the tensor evaluation; (b) check Morton's separate measurement (if D-vs-C shows large gain, cache is the cause); fix or document; never relax the bound to merge |
| Morton causes Tier 1 force discrepancy above FP-reorder floor | Permutation changed final tree structure (tie-breaking bug at exact octant boundaries) or Morton encoding has a sign / interleave bug | Investigate `child_octant` against bodies at quantisation boundaries; verify Morton encoding round-trips on a synthetic 8-corner test; fix and re-run; never relax the bound |
| D-vs-C **sub-additive** (`t_walk_D > t_walk_C` despite both Morton enabled and quadrupole enabled) | Quadrupole memory pressure overflows what Morton can hide; the AoS layout assumption is wrong | This is a structural signal, not a feature decision. Stop the perf branch; investigate node layout (consider splitting hot fields and cold fields, or SoA); revisit before merging anything |

### Methodology

#### Quadrupole expansion (p = 2)

**Tensor representation.** Symmetric traceless quadrupole tensor `Q` at every node, stored as 5 independent components `Q_xx, Q_xy, Q_xz, Q_yy, Q_yz` with `Q_zz = −(Q_xx + Q_yy)` enforced by traceless constraint. Adds 40 bytes per node (5 × 8). Node footprint grows from ~144 to ~184 bytes (~28 % increase).

**Construction during `aggregate_mass`.** For each node, after children's masses and COMs are known:

```text
For internal node with children {c}:
  M_node, COM_node ← already aggregated from children
  Q_node ← Σ_c [Q_c + M_c · (3 · D_c ⊗ D_c − I · |D_c|²)]
  where D_c = COM_c − COM_node, ⊗ is outer product, I is identity tensor
For leaf node with bodies {b}:
  Q_node ← Σ_b m_b · (3 · d_b ⊗ d_b − I · |d_b|²)
  where d_b = pos_b − COM_node
```

The internal-node formula is the parallel-axis theorem for second moments (translation of a child's tensor from its own COM to the parent's COM). Reference: Goldstein, Poole & Safko §11.3.

**Force evaluation in `bh_eval_body`.** When a node passes the BH opening criterion, the acceleration on the target at vector displacement `r` from the node's COM (magnitude `r`, unit vector `r̂`):

```text
a_mono = −G · M / r³ · r
a_quad = −G / r⁵ · [Q · r̂ − (5/2) · (r̂ᵀ · Q · r̂) · r̂]
a_total = a_mono + a_quad
```

Standard expansion derived from the gradient of `Φ = −G · [M/r + (1 / (2 r³)) · Σ_ij Q_ij · n_i · n_j]`. References: Hernquist & Katz 1989 eq. (2.11); Dehnen 2002 §3.

Per-interaction cost ratio quad/mono ≈ 2.5–3× (consistent with Dehnen 2002 §5 and Springel 2005 §2.4). The Tier 2 bound range [0.30, 0.70] for `t_eval_C / t_eval_A` accommodates this; values outside that range trigger investigation per Decision Rules.

#### Morton encoding (Z-order)

- Per body: normalise position to `[0, 1]³` against the build-time AABB, quantise each axis to 21 bits, interleave bits via standard "magic-number" sequence to a 63-bit `u64` Morton code.
- Sort `Vec<(u64, u32)>` (code, body index) via `sort_unstable_by_key`.
- Insert in Morton order during `Octree::build`.
- Walk in Morton order during `BarnesHutEngine::evaluate`: `perm.par_iter().map(|&i| ...).for_each(|(i, acc)| accels[i] = acc)`. Output ordering is preserved (writes go to original-index slots).

If `core/system/tests::replay::*` fails post-Morton due to non-deterministic equal-Morton-code ordering, switch to `sort_by_key` (stable). 21-bit quantisation gives 2²¹ ≈ 2.1 million cells per axis, so equal-code collisions are negligible in realistic distributions.

#### Matched-accuracy θ search

Binary search on θ ∈ [0.5, 1.0] for each Morton×multipole cell:

1. Measure `error_mono(0.5)` once per (cell, N, seed) — the reference accuracy.
2. For quadrupole cells, bisect θ to find θ_match such that `|error_quad(θ_match) − error_mono(0.5)| / error_mono(0.5) ≤ 0.05`.
3. Convergence in ≤ 6 iterations (θ tolerance ±0.01 across [0.5, 1.0]).
4. **Pre-flight monotonicity check.** Sample θ ∈ {0.3, 0.5, 0.7, 0.9} at the cell; verify strict monotonic error growth before invoking the binary search. If monotonicity fails, the body distribution has degeneracies the search assumes away — investigate before reporting θ_match.

The search lives in the bench harness (`crates/apsis/benches/perf_2x2.rs`), not the production engine. Production engine ships with a single fixed θ baked-in once §Decision is written.

#### Toggle: `pub(crate)` runtime knob, removed in final commit

`BarnesHutEngine` gains two `pub(crate)` setters:

- `set_multipole_order(order: MultipoleOrder)` accepting `Monopole | Quadrupole`
- `set_morton_enabled(enabled: bool)`

Tests and benches in PR-perf-1 / PR-perf-2 exercise the 2×2 matrix via these setters.

After §Decision is written, the final commit of PR-perf-2 removes both setters and bakes the chosen configuration in. The toggle exists only during the experiment; the shipped engine has no dynamic switch. Reproducibility for future research: the experiment commits remain in `git log` linked from §Results CSV.

#### Frozen variables

| Variable | Pinned value | Why pinned |
| --- | --- | --- |
| Compiler | rustc 1.85+ (workspace `rust-version`) | Per workspace `Cargo.toml`; exact `rustc --version` recorded in §Results |
| Profile | `cargo build --release` defaults: `opt-level = 3`, `lto = false`, `codegen-units = 16`, `incremental = false` | Confirmed via inspection: no `[profile.release]` overrides at workspace or crate level. Pinned so future LTO-enabled measurements are clearly cross-experiment, not cross-config |
| Allocator | Rust stdlib default (Windows 11 default heap; `mimalloc` / `jemallocator` confirmed absent from `Cargo.lock` and source tree) | Morton sort interacts with allocator behaviour; pinning ensures B-vs-A and D-vs-C are not contaminated by allocator-side effects |
| OS / hardware | Windows 11, same machine as octree-port Tier 3 | Wall-time numbers are not portable; recorded in §Results |
| Rayon thread pool | Default (`rayon::current_num_threads()`) | A control on this would defeat the cache-locality story; recorded in §Results |
| Warm-up evaluations | 1 per cell per N per θ before timed runs | Excludes first-touch effects (CPU frequency scaling, allocator warm-up, page faults) |
| Timed evaluations | 10 per cell per N per θ; report median + 1σ | Robust to outliers from OS scheduling jitter; σ feeds the variance-stability decision |
| Seeds | 3: `0x6F637472`, `0x71756164`, `0x6D6F7274` | Multi-seed addresses the single-seed threat-to-validity from octree-port |

#### Out of scope (declared a priori)

- **Higher-order multipoles (p ≥ 3).** Quadrupole is the canonical first improvement; octupole and Dehnen FMM proper are post-Morton work, not part of this experiment.
- **SIMD inner loops in `bh_eval_body`.** Marginal-future bucket per the perf-categorisation plan; gated on Morton landing first to provide spatial coherence.
- **Adaptive θ controllers.** `ThetaController` consumes the θ-error proxy unchanged; tuning it for quadrupole is a separate experiment.
- **GPU offload.** Niche bucket; out of scope.
- **Radix sort for Morton codes.** `std::sort_unstable_by_key` is the chosen sort; radix is a follow-up only if profiling shows the standard sort dominating (predicted not — for N = 10⁵ the sort is estimated < 5 % of build cost).
- **Cargo features for the toggles.** Cargo features must be additive (Rust API guideline C-FEATURE); a `morton` or `quadrupole` feature would be exclusive/negative, breaking composability for downstream crates. The runtime knobs are `pub(crate)` only and removed in the final commit. No Cargo features, no public setters.
- **Hot/cold field splitting in `Node`.** Only triggered if D-vs-C sub-additivity is observed (Decision Rule "structural signal"); otherwise the AoS layout matches the quadtree's and the diff stays focused.

---

## Results

*To be populated incrementally — PR-perf-1 fills cells A and C; PR-perf-2 fills cells B and D and writes §Decision.*

### Tier 1 — Force accuracy

*Pending.*

### Tier 2 — Wall-time at matched accuracy

*Pending.*

### Tier 3 — Cache-effect characterisation

*Pending.*

### Pareto frontier

*Pending.*

---

## Interpretation

*To be written after both PRs land and §Results is fully populated.*

---

## Threats to validity

1. **Multi-seed but single-machine.** Cache effects are hardware-sensitive; the gain measured on the development machine (Windows 11, Rayon-default thread pool, identified in §Results) may not reproduce on machines with different L2/L3 sizes, prefetcher behaviour, or core counts. Tier 2 ranges are conservative against literature spread to accommodate ±50 % variance from this; cross-machine reproducibility would require re-baselining all four cells.

2. **Rayon work-stealing fragmentation.** Morton-ordered iteration into Rayon's parallel iterator does not guarantee consecutive iterations process on the same thread. Work-stealing can fragment Morton blocks. Mitigation: if the Tier 2 walk-time bound for Morton-on cells fails, investigate `with_min_len` chunk-size tuning before declaring a regression.

3. **Quadrupole tensor algebra correctness.** The parallel-axis-theorem combination of children's tensors is error-prone (sign conventions, traceless enforcement). Mitigation: a synthetic 2-body test (one child, one body) exercises the leaf path; a 4-body symmetric test (3-child internal node with bodies symmetric about COM) exercises the parallel-axis path; both gated in Tier 1.

4. **Matched-accuracy θ search convergence.** Binary search on θ assumes monotonic error growth with θ. This holds for any well-behaved BH implementation but can break if the body distribution has degeneracies (e.g. all bodies on a plane, which the random-sphere distribution rules out). Pre-flight monotonicity check at θ ∈ {0.3, 0.5, 0.7, 0.9} per cell is required before invoking the binary search.

5. **Sub-additivity false positive.** D being slower than C does not always mean the AoS layout is wrong; it could mean the build-time Morton encoding cost is large relative to the walk gain at the measured N. Mitigation: Tier 2's separated build/walk decomposition isolates this — sub-additivity is concluded only when `t_walk_D > t_walk_C`, not when `t_eval_D > t_eval_C`.

6. **Pareto frontier visualisation gap.** The frontier is reported as CSV; without a plotting step in the notebook itself, future readers must produce the figures externally. Mitigation: a separate `docs/experiments/2026-05-08-octree-perf-frontier.py` script alongside the CSV in PR-perf-2 produces matplotlib figures from the CSV; the §Decision references the figure.
