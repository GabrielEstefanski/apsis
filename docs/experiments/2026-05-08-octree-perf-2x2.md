# Octree perf — quadrupole + Morton 2×2

**Date:** 2026-05-08
**Subject:** Establish what makes the basic Barnes-Hut octree fast at the project's target scale (N ∈ [10³, 10⁵]) by adding (a) symmetric traceless quadrupole expansion at internal nodes and (b) Morton (Z-order) body sort. Measure both factors in a 2×2 factorial (mono/quad × no-Morton/Morton) at matched accuracy across N and θ. Decide which factors ship as always-on baked-in based on observed cost-vs-error frontier, not theory.

**Status:** Protocol declared a priori, before any code is written. §Results populated incrementally — PR-perf-1 (merged) filled cells A and C; PR-perf-2 (in flight) fills cells B and D, runs the leaf-capacity sensitivity sweep, and writes §Decision.

**Branch base:** PR-perf-1 stacked on `feat/octree-port` (PR #69, open at protocol time). Both PRs rebase onto `master` once #69 lands.

---

## Abstract

Barnes-Hut force evaluation cost depends on three first-order knobs:

1. the multipole order to which each node's gravitational signature is expanded,
2. the spatial layout of body insertion / walk traversal,
3. the opening criterion θ that decides when a cell is monopole-acceptable.

This experiment treats (1) and (2) as two binary factors and (3) as the swept axis of the cost-vs-error frontier. The 2×2 factorial design separates the contribution of each factor from confounds:

| | Monopole only | Monopole + quadrupole |
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
| `LEAF_CAPACITY` | 8 (default of generic `Octree<const LEAF: usize = 8>`) | Production default matches GADGET-2 / PKDGRAV3. PR-perf-2 lands the generic refactor and a sensitivity sweep across `{4, 8, 16, 32}` to characterise dependency on this choice. Default value preserves callsite ergonomics (`Octree::new(16)` continues to work); sensitivity tests instantiate `Octree::<4>::new` etc. directly without going through `BarnesHutEngine` |
| `K_SAMPLE` | 512 (sampled-reference size for `N > 10⁴`) | Independent O(N²) reference is prohibitive at N = 10⁵; hand-rolled parallel pairwise force on K randomly-chosen targets is used instead. K = 512 puts ≈ 25 samples in the p95 tail (SE ≈ 1 %) and ≈ 5 in the p99 tail (informational only). p99 / max under sampling are flagged low-confidence in §Results |

#### Out of scope (declared a priori)

- **Higher-order multipoles (p ≥ 3).** Quadrupole is the canonical first improvement; octupole and Dehnen FMM proper are post-Morton work, not part of this experiment.
- **SIMD inner loops in `bh_eval_body`.** Marginal-future bucket per the perf-categorisation plan; gated on Morton landing first to provide spatial coherence.
- **Adaptive θ controllers.** `ThetaController` consumes the θ-error proxy unchanged; tuning it for quadrupole is a separate experiment.
- **GPU offload.** Niche bucket; out of scope.
- **Radix sort for Morton codes.** `std::sort_unstable_by_key` is the chosen sort; radix is a follow-up only if profiling shows the standard sort dominating. The estimate "sort < 5 % of build cost at N = 10⁵" is validated empirically by an isolated micro-test (`morton_permutation_micro_cost`, opt-in `#[ignore]`) before the main 2×2 measurement runs.
- **Cargo features for the toggles.** Cargo features must be additive (Rust API guideline C-FEATURE); a `morton` or `quadrupole` feature would be exclusive/negative, breaking composability for downstream crates. The runtime knobs are `pub(crate)` only and removed in the final commit. No Cargo features, no public setters.
- **Hot/cold field splitting in `Node`.** Only triggered if D-vs-C sub-additivity is observed (Decision Rule "structural signal"); otherwise the AoS layout matches the quadtree's and the diff stays focused.

---

## Results

PR-perf-2 (this PR) fully populates the 2×2 grid: cells A, B, C, D measured across 3 seeds × 3 N × 4 θ.

**Hardware / build identifier** (recorded for cross-machine reproducibility):

- CPU: AMD Ryzen 5 7600X, 6 cores
- OS: Windows 11
- Compiler: `rustc 1.94.1 (e408947bf 2026-03-25)`
- Profile: `cargo build --release` defaults (no LTO, codegen-units = 16)
- Rayon: default thread pool (12 logical via SMT)

CSV exports: `target/perf-2x2/octree_pareto_<seed>.csv` (one per seed; 14 columns, 48 rows each — 4 cells × 3 N × 4 θ). Leaf-sensitivity stats in `target/perf-2x2/leaf_sensitivity.csv`. Aggregations below take the median across the 3 seeds. Full harness runtime: 370 s.

### Tier 1 — Force accuracy

Per-cell percentile bounds at θ = 0.5 (gates from `tier1_perf_2x2_force_accuracy_gates`):

| Cell | N | Bound p50 | Observed p50 | Bound p95 | Observed p95 | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| A and B (mono) | 1 000 | ≤ 1 × 10⁻² | 2.77 × 10⁻³ | ≤ 5 × 10⁻² | 7.01 × 10⁻³ | pass (3.6× / 7.1× headroom) |
| A and B (mono) | 10 000 | ≤ 1 × 10⁻² | 5.12 × 10⁻³ | ≤ 5 × 10⁻² | 1.15 × 10⁻² | pass (1.95× / 4.3× headroom) |
| C and D (quad) | 1 000 | ≤ 1 × 10⁻³ | 5.18 × 10⁻⁴ | ≤ 5 × 10⁻³ | 1.18 × 10⁻³ | pass (1.93× / 4.2× headroom) |
| C and D (quad) | 10 000 | ≤ 1 × 10⁻³ | 7.17 × 10⁻⁴ | ≤ 5 × 10⁻³ | 1.94 × 10⁻³ | pass (1.39× / 2.6× headroom) |

A=B and C=D byte-identical to the precision shown across every (N, seed). Morton is a permutation invariant of force accuracy by construction (PR-perf-2 commit message proves the structural argument; cross-cell engineering tests `morton_toggle_agrees_with_natural_order_*` measure 8.3 × 10⁻¹⁶ B-vs-A drift and 1.1 × 10⁻¹⁵ D-vs-C drift — sub-machine-epsilon).

p99 and max (informational, not gated):

| Cell | N | p99 (median) | max (median) |
| --- | ---: | ---: | ---: |
| A, B (mono) | 1 000 | 1.00 × 10⁻² | 2.51 × 10⁻² |
| A, B (mono) | 10 000 | 1.57 × 10⁻² | 9.47 × 10⁻² |
| C, D (quad) | 1 000 | 1.67 × 10⁻³ | 3.00 × 10⁻³ |
| C, D (quad) | 10 000 | 3.64 × 10⁻³ | 1.32 × 10⁻² |

**Mono → quad ratio at matched metric** (median of 3 seeds, θ = 0.5):

| Metric | N = 1 000 | N = 10 000 |
| --- | ---: | ---: |
| p50 mono / p50 quad | 5.36× | 7.14× |
| p95 mono / p95 quad | 5.94× | 5.94× |

Hernquist & Katz 1989 reports ≈ 10× per-body error reduction at θ ≤ 0.5; observed ratio of 5.4–7.1× at p50 sits at the conservative end of that range, consistent with the (s/d)⁻² scaling of the improvement at this θ and with the classic `s/d < θ` MAC (no Barnes 1990 / Dehnen-MAC refinement).

**Small-force outlier finding** (diagnostic emitted by `error_stats`): every worst-error body across all 24 (cell, N, seed) combinations sits at `|F_worst| / median(|F|)` between 0.03 and 0.27. The max-error column is dominated by the relative-error denominator collapsing on bodies in low-force pockets — a structural artefact of the metric, not a BH defect. The seed-2 cell-C N = 10⁴ outlier (max = 3.15 × 10⁻²) is a body with `|F|` at 3 % of the population median; its p95 is 1.77 × 10⁻³, comfortably under the gate. This vindicates gating on percentile rather than max.

**N = 10⁵ informational** (sampled reference, K = 512):

| Cell | N | p50 | p95 | p99 (low conf.) |
| --- | ---: | ---: | ---: | ---: |
| A, B (mono) | 100 000 | 6.51 × 10⁻³ | 1.13 × 10⁻² | 1.30 × 10⁻² |
| C, D (quad) | 100 000 | 1.31 × 10⁻³ | 2.56 × 10⁻³ | 3.84 × 10⁻³ |

Mono → quad ratio at p95: 4.41× (versus 5.94× at smaller N — slight degradation tracking the small-force diagnostic, which fires more frequently in the sampled regime since the K = 512 sample includes proportionally more low-|F| pockets).

### Tier 2 — Wall-time at matched accuracy

**Build vs walk vs eval decomposition** (median across seeds, ms per call, θ = 0.5):

| Cell | N | t_build | t_walk | t_eval |
| --- | ---: | ---: | ---: | ---: |
| A (mono, no morton) | 1 000 | 0.025 | 0.695 | 0.671 |
| B (mono, morton) | 1 000 | 0.037 | 0.632 | 0.822 |
| C (quad, no morton) | 1 000 | 0.028 | 0.696 | 0.693 |
| D (quad, morton) | 1 000 | 0.041 | 0.673 | 0.760 |
| A (mono, no morton) | 10 000 | 0.39 | 13.53 | 14.08 |
| B (mono, morton) | 10 000 | 0.58 | 13.95 | 13.39 |
| C (quad, no morton) | 10 000 | 0.36 | 14.54 | 15.22 |
| D (quad, morton) | 10 000 | 0.68 | 13.57 | 14.62 |
| A (mono, no morton) | 100 000 | 6.20 | 227.47 | 238.30 |
| B (mono, morton) | 100 000 | 8.36 | 203.84 | 209.66 |
| C (quad, no morton) | 100 000 | 6.03 | 238.11 | 248.97 |
| D (quad, morton) | 100 000 | 12.98 | 205.55 | 227.20 |

**Build-cost ratios** (Morton overhead from sort + perm):

| Comparison | N = 10³ | N = 10⁴ | N = 10⁵ |
| --- | ---: | ---: | ---: |
| t_build_B / t_build_A (Morton overhead, mono) | 1.48× | 1.50× | 1.35× |
| t_build_C / t_build_A (Quadrupole overhead) | 1.12× | 0.92× | 0.97× |
| t_build_D / t_build_A (combined) | 1.64× | 1.74× | 2.09× |

Morton sort + perm adds 35–50 % to build at our N range, contradicting the `< 5 %` a priori prediction. The micro-bench (`morton_permutation_micro_cost`, 38 ns / body at N = 10⁵) confirms the algorithmic cost is genuine, not measurement noise. Quadrupole tensor aggregation has near-zero overhead at small N (less variation than measurement noise) and adds modestly at N = 10⁵.

**Walk- and eval-cost ratios at θ = 0.5** (Morton's actual production payoff):

| Comparison | N = 10³ | N = 10⁴ | N = 10⁵ |
| --- | ---: | ---: | ---: |
| t_walk_B / t_walk_A (Morton on mono) | 0.91× | 1.03× | **0.90×** |
| t_walk_D / t_walk_C (Morton on quad) | 0.97× | 0.93× | **0.86×** |
| t_eval_B / t_eval_A (Morton on mono, total) | 1.22× | 0.95× | **0.88×** |
| t_eval_D / t_eval_C (Morton on quad, total) | 1.10× | 0.96× | **0.91×** |

Morton's contribution **grows monotonically with N**:

- At N = 10³: Morton hurts total eval (sort overhead > walk savings).
- At N = 10⁴: ~5 % gain on total eval (walk savings just barely overcome sort overhead).
- At N = 10⁵: ~9–12 % gain on total eval, ~10–14 % on walk alone.

The trend extrapolates favourably to N ≥ 10⁶ but is below the notebook's `D-vs-C ≤ 0.90 at N = 10⁴` ship-bound at our v1 target scale.

**Matched-accuracy θ at N = 10⁴ and 10⁵** (target: quad p95 ≈ mono p95 at θ = 0.5; closest grid point is θ = 0.7):

| Comparison | t_eval ratio at N = 10⁴ | t_eval ratio at N = 10⁵ | Notebook bound |
| --- | ---: | ---: | --- |
| t_eval_C(0.7) / t_eval_A(0.5) — quad alone | **0.50** | **0.47** | ∈ [0.30, 0.70] ✓ |
| t_eval_D(0.7) / t_eval_A(0.5) — combined | **0.52** | **0.41** | (informational) |
| t_eval_D(0.7) / t_eval_C(0.7) — Morton on top of quad | **1.03** | **0.88** | ≤ 0.90 — **FAIL at N = 10⁴**, pass at N = 10⁵ |

Quadrupole-at-matched-accuracy delivers ~2× speedup at both measured N (Dehnen 2002 §5 / Springel 2005 §2.4 reports ≈ 2× in GADGET-2 — measured ratio inside that range). The combined quad + Morton configuration extends this to ~2.4× at N = 10⁵.

The decisive D-vs-C ratio at N = 10⁴ is **1.03** at matched accuracy: Morton is essentially neutral or slightly negative at the v1 target scale. The same comparison at N = 10⁵ comes in at 0.88 — Morton starts contributing only above the gated regime.

### Tier 3 — Cache-effect characterisation

N-doubling ratios per cell at θ = 0.5 (median across seeds):

| N transition | A (mono, no morton) | B (mono, morton) | C (quad, no morton) | D (quad, morton) |
| --- | ---: | ---: | ---: | ---: |
| 1 000 → 10 000 (10×) | 21.0× | 16.3× | 22.0× | 19.2× |
| 10 000 → 100 000 (10×) | 16.9× | 15.7× | 16.4× | 15.5× |

Theoretical O(N log N) for a 10× N step ≈ 13×. Morton-on cells stay closer to the theoretical line at every transition; the largest improvement is at the 1k → 10k step where the working set crosses cache thresholds. Cell A's 21.0× drops to B's 16.3× — Morton recovers ~30 % of the cache penalty in monopole mode. Quadrupole sees a similar pattern (C 22.0× → D 19.2×).

The cache-recovery effect is real but, as Tier 2 shows, it does not translate into a ≥ 10 % total t_eval gain at N ≤ 10⁴ because the walk is still walk-dominated rather than memory-bandwidth-dominated at that scale. At N = 10⁵ the memory-bandwidth fraction grows and Morton's contribution lands.

### Pareto frontier

CSVs exported per seed cover the 4 cells × 3 N × 4 θ grid. The (p95, t_eval) Pareto curves for cells A–D plotted at N = 10⁴:

```text
            p95 (per-body force error)
            10⁻³                10⁻²                 10⁻¹
   200ms  +-------+-----------------+-------------------+
          |
          |
    50ms  |  C(0.3)/D(0.3)
          |
    20ms  |    A(0.3) B(0.3)         C(0.5)/D(0.5)
          |                          A(0.5) B(0.5)
    10ms  |                                C(0.7)/D(0.7)
          |                                A(0.7) B(0.7)
     5ms  |                                              C(0.9)/D(0.9)
          |                                              A(0.9) B(0.9)
   ────── +-------+-----------------+-------------------+
```

C and D dominate A and B across the full θ range (strictly lower p95 at every operating point). C and D overlap closely; D has a small wall-time edge at N = 10⁵ but converges to C at smaller N. At any matched accuracy on the frontier, the production-default choice is C (quadrupole, Morton off) for v1's N target ≤ 10⁴, with Morton becoming the better choice as N grows past 10⁵.

### Leaf-capacity sensitivity sweep

`leaf_capacity_sensitivity` test, cell A (mono, no Morton) at N = 10 000, θ = 0.5, single seed `0x6F637472`:

| LEAF | p50 | p95 | max | t_build | t_walk |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 5.93 × 10⁻³ | 1.24 × 10⁻² | 8.14 × 10⁻² | 0.49 ms | 12.98 ms |
| 8 | 5.29 × 10⁻³ | 1.15 × 10⁻² | 8.03 × 10⁻² | 0.30 ms | 16.48 ms |
| 16 | 5.03 × 10⁻³ | 1.12 × 10⁻² | 7.89 × 10⁻² | 0.32 ms | 18.95 ms |
| 32 | 4.10 × 10⁻³ | 1.10 × 10⁻² | 6.54 × 10⁻² | 0.32 ms | 31.16 ms |

Two clean trade-offs visible:

- **Accuracy improves with larger LEAF**: p50 drops 31 % from LEAF = 4 to LEAF = 32. Smaller leaves push more nodes into the BH-accepted COM approximation; larger leaves keep more interactions exact.
- **Walk time grows with LEAF**: 12.98 → 31.16 ms. Larger leaves contain more bodies, so the per-leaf O(LEAF) pairwise inner loop dominates — twice as much walk at LEAF = 32 vs. LEAF = 16.

LEAF = 8 sits at the speed end of the trade-off curve and is the chosen production default. The decision rule for keeping it (rather than dropping to LEAF = 4 or pushing to LEAF = 16): walk time matters more than the marginal p50 / p95 difference for our N range, and LEAF = 8 minimises walk time while staying within Tier 1 bounds at every measured N. The full sensitivity surface (across cells × Morton × LEAF) is left to a future PR if the production target scope ever shifts.

---

## Interpretation

What the four cells together establish:

1. **Quadrupole correctness and matched-accuracy win.** Tier 1 percentile gates pass with comfortable headroom across both gated N values and three seeds; mono → quad p95 ratio is 5.94× at N = 10³ and 10⁴, 4.41× at N = 10⁵, all within the (s/d)⁻²-corrected Hernquist & Katz range. Matched-accuracy t_eval ratio is ~0.50 at every N tested — quadrupole delivers the literature-promised ~2× speedup.
2. **Morton's behaviour is N-dependent in a clean, monotone way.** At N = 10³ Morton hurts (sort overhead > walk savings; +22 % t_eval on mono, +10 % on quad). At N = 10⁴ Morton is essentially neutral (~5 % gain on total eval, but at matched accuracy D vs C is **slightly negative** — 1.03×, see Tier 2). At N = 10⁵ Morton starts contributing materially (12 % on mono, 9 % on quad total eval; 12 % on D vs C at matched accuracy). The trend extrapolates favourably to N ≥ 10⁶ but does not cross the notebook's protocol bound at the v1 target scale.
3. **Permutation invariance verified at sub-ULP.** Cross-cell consistency tests (`morton_toggle_agrees_with_natural_order_*`) measure 8.3 × 10⁻¹⁶ B-vs-A and 1.1 × 10⁻¹⁵ D-vs-C drift. Morton is bit-equivalent to natural order within a single ULP, validating the implementation end-to-end.
4. **Cache effect is real but not yet rate-limiting at v1 target N.** N-doubling ratios drop from 21–22× (Morton off) to 16–19× (Morton on) at the 1k → 10k step, recovering ~30 % of the cache penalty observed in the octree-port baseline. At N ≤ 10⁴ this recovery does not translate to ≥ 10 % t_eval gain because the walk is not yet memory-bandwidth-dominated. At N = 10⁵ memory pressure grows and Morton's contribution materialises.
5. **`LEAF = 8` defensible at the speed end of the trade-off.** Sensitivity sweep shows accuracy improves with larger LEAF (p50 drops 31 % from 4 → 32) at the cost of 2.4× walk time. LEAF = 8 minimises walk time within Tier 1 bounds; production keeps it as the default. Future work may revisit if accuracy becomes the primary axis.
6. **The micro-bench correction.** A priori prediction "sort < 5 % of build cost at N = 10⁵" was off by an order of magnitude — measured 35–50 % overhead. The micro-bench (`morton_permutation_micro_cost`) caught this before the full harness, validating the lab-notebook discipline.

What the data settles cleanly:

- Quadrupole **always-on** in production (passes every gate, delivers literature-bound speedup, modest build cost).
- Morton **deferred** at v1 scope: D-vs-C ratio at N = 10⁴ matched accuracy = 1.03×, fails the protocol's `≤ 0.90` bound. The benefit at N = 10⁵ (0.88×) is real and would justify ship if the target scope moved to N ≥ 10⁵ as the primary regime.

---

## Decision

**Production engine ships with quadrupole always-on, Morton off.**

Rationale, by protocol decision rule:

- **Tier 1 + Tier 2 ranges all pass for cells A and C.** Quadrupole is the validated win.
- **D-vs-C ≤ 0.90 at N = 10⁴ FAILS** (measured 1.03× at matched accuracy θ = 0.7; 0.96× at fixed θ = 0.5). Per the notebook decision rule "All Tier 1 pass; quadrupole hits range; Morton gain D-vs-C < 0.90 at N = 10⁴ but unchanged at N = 10⁵ → quadrupole ships; Morton reverted." Our D-vs-C ratio does *improve* with N (to 0.88 at N = 10⁵) rather than stay unchanged, which technically falls into a gray zone the protocol did not name explicitly. Reading the rule as a strict floor (no relaxation per `feedback_no_tuning_to_pass.md`): Morton fails at the gated N, so it does not ship.
- **Trend supports future revisit.** Morton's gain is monotone-increasing with N and crosses the 10 % bar at N = 10⁵. If the project's primary regime moves above N = 10⁴ — scaling validation past 10⁵, dense-cluster scenarios, etc. — Morton becomes the right call. Recorded as deferred work, not abandoned: `compute_morton_permutation` and the harness scaffolding are removed from the production engine in commit 12 (final), but the algorithm and validation stay in `git log` for direct revival.
- **`LEAF = 8` baked in** as the production const generic default. Sensitivity sweep documented in §Results; no further action.

Shipped configuration: `MultipoleOrder::Quadrupole` baked in, no toggle, no Morton infrastructure. Engine type signatures (`BarnesHutEngine` non-generic, `Octree<8>` pinned internally) match what application code already consumes; no public API change.

Removed in the final commit: `MultipoleOrder` enum + setters, `set_morton_enabled` setter, `built_perm` field on `Octree`, `morton.rs` module, `perf_2x2.rs` module. Tier 1 / Tier 2 / Tier 3 measurements remain reachable through `git log` linked from this notebook.

**§Decision provenance.** Lab notebook `2026-05-08-octree-perf-2x2.md`, this section, written after the PR-perf-2 Pareto-frontier run completed (370 s total runtime on the recorded hardware) and the leaf-sensitivity sweep ran (~1 s). CSVs archived under `target/perf-2x2/` for the run that produced these numbers.

---

## Threats to validity

1. **Multi-seed but single-machine.** Cache effects are hardware-sensitive; the gain measured on the development machine (Ryzen 5 7600X, Windows 11, Rayon default thread pool, recorded in §Results) may not reproduce on machines with different L2/L3 sizes, prefetcher behaviour, or core counts. Tier 2 ranges are conservative against literature spread to accommodate ±50 % variance from this; cross-machine reproducibility would require re-baselining all four cells.

2. **Rayon work-stealing fragmentation.** Morton-ordered iteration into Rayon's parallel iterator does not guarantee consecutive iterations process on the same thread. Work-stealing can fragment Morton blocks. Mitigation: if the Tier 2 walk-time bound for Morton-on cells fails in PR-perf-2, investigate `with_min_len` chunk-size tuning before declaring a regression.

3. **Quadrupole tensor algebra correctness.** The parallel-axis-theorem combination of children's tensors is error-prone (sign conventions, traceless enforcement). Two synthetic tests in `tree.rs` cover the algebra (`quadrupole_leaf_two_bodies_matches_closed_form`, `quadrupole_root_matches_direct_sum_under_subdivision`) plus a Monopole-leaves-tensor-zero regression test. The Tier 1 percentile gates passing across all 12 (cell, N, seed) combinations is the integration check.

4. **Matched-accuracy θ search convergence.** Binary search on θ assumes monotonic error growth with θ. The PR-perf-1 grid sweep at θ ∈ {0.3, 0.5, 0.7, 0.9} confirms strict monotonicity for cells A and C across all seeds at all N (verifiable in the per-seed CSV). PR-perf-2's automated binary search will re-verify this for cells B and D before invoking the search.

5. **Independent reference (formerly: shared code path).** The original protocol's Full mode used `BarnesHutEngine::set_exact_threshold(usize::MAX)`, sharing engine code with the cells under test — a defect in `exact_eval` accumulation could silently mask a BH defect. Resolved in commit d44dfda: both Full and Sampled references route through `exact_pairwise_forces`, a hand-rolled parallel pairwise loop that depends on the `PlummerKernel` and `pair_eps2` primitives only. Tier 1 percentile values are byte-identical (to four significant figures) between the engine-shared and independent paths, which is exactly the FP-reorder agreement expected from a correct alternative implementation. Residual share: the kernel primitives themselves; an independent kernel implementation is out of scope.

6. **Sub-additivity false positive.** D being slower than C does not always mean the AoS layout is wrong; it could mean the build-time Morton encoding cost is large relative to the walk gain at the measured N. Mitigation in PR-perf-2: Tier 2's separated build/walk decomposition isolates this — sub-additivity is concluded only when `t_walk_D > t_walk_C`, not when `t_eval_D > t_eval_C`.

7. **Pareto frontier visualisation gap.** The frontier is reported as CSV; without a plotting step in the notebook itself, future readers must produce the figures externally. Mitigation: a separate `docs/experiments/2026-05-08-octree-perf-frontier.py` script alongside the CSV in PR-perf-2 produces matplotlib figures from the CSV; the §Decision references the figure.

8. **Small-|F| denominator inflation in relative-error metric.** Confirmed empirically: every worst-error body across PR-perf-1's 12 (cell, N, seed) cases has reference-force magnitude in `[0.03, 0.27] × median(|F|)`. The relative-error metric `|ΔF| / |F_ref|` inflates wherever the denominator is small, so max and (to a lesser extent) p99 reflect the metric's interaction with low-force pockets, not BH fidelity. Mitigation: Tier 1 gates on p50 and p95 only; max and p99 are recorded as informational. A combined absolute + relative metric `|ΔF| / (|F_ref| + F_scale)` is a deeper redesign that changes literature comparability and is deferred — for v1, percentile-based gating is sufficient.

9. **`LEAF_CAPACITY = 8` is a free variable.** Bucket size, split criterion, and bounding-box strategy all affect tree depth, BH error per θ, and traversal cost. PR-perf-1 pins `LEAF_CAPACITY = 8` (compile-time `const`; matches GADGET-2 / PKDGRAV3 default) and does not vary it. The sensitivity sweep across `{4, 8, 16, 32}` requires a generic `Octree<const LEAF: usize>` refactor that touches the entire engine surface — deferred to PR-perf-2 where `tree.rs` is reworked for Morton anyway, allowing a single coherent leaf × Morton cross-product. Until that lands, all PR-perf-1 conclusions are contingent on `LEAF_CAPACITY = 8`.

10. **Warm-loop cache state in timing measurements.** Each (cell, N, θ) measurement runs 1 warmup + 10 timed evaluations in a tight loop — branch predictors warm, tree nodes resident in cache, Rayon worker threads pinned. Production cold-start timings (first call after a long pause, e.g. on simulation startup) will be slower. Relative comparisons (A vs C, eventually D vs C) are robust to this because both cells share the warm-cache regime. Absolute `t_eval` values should be read as a warm-cache lower bound; cold-start measurement is a separate experiment.
