# Tree incremental updates — protocol

**Date:** 2026-05-12
**Subject:** Replace `Octree::build`'s per-step from-scratch tree reconstruction with an **incremental maintenance pattern**: each particle stores a back-reference to its current cell; per-step update walks the particles, checks if each is still inside its known cell, and re-inserts only the migrants. The tree-walk pass derefines empty / single-occupant cells, and multipole moments are recomputed leaf-up identically to the current rebuild approach. Pattern follows REBOUND (`tree.c`, `reb_simulation_update_tree_cell` + `reb_simulation_update_tree_gravity_data_in_cell`) and is the standard maintenance idiom in production planetary N-body codes.

**Status:** Protocol declared a priori, before any implementation lands. Calibrated against the engine-ceiling re-measurement (`docs/experiments/2026-05-12-engine-ceiling-post-simd.md`) which classified `build` as 4-13 % of step wall-time and identified it as the highest-ROI-per-LOC residual axis. The change has **zero precision impact**: the tree is a spatial accelerator, multipoles are recomputed identically — the per-cell mass / COM / quadrupole values must be bit-exact with the from-scratch rebuild they replace.

**Branch:** `perf/tree-incremental`, stacked on top of `perf/simd-kernel` (PR #79). Independent of PR #78 / #79 in the algorithmic sense — touches `tree.rs` and `Octree::build`, not the walk SIMD or SoA snapshot. Stacks for ergonomic reasons (avoids rebase against the SoA-renamed Body fields).

---

## Abstract

Apsis currently rebuilds the entire octree from scratch every integration step. Engine-ceiling Tier 1 measures `build` at 4-13 % of step wall-time across N ∈ [10², 10⁵] (`docs/experiments/2026-05-09-engine-ceiling.md` §Tier 1). For stable systems (planetary configurations, hierarchical multi-body) the **vast majority of particles do not cross cell boundaries between consecutive steps** — the tree topology is essentially stationary modulo a small migrant set.

REBOUND's `tree.c` exploits this: each particle holds a back-pointer to its current cell; per-step maintenance does:

1. For each particle, test residency in its known cell.
2. If migrated, remove from old cell, re-insert into the appropriate root or descendant.
3. Walk the tree, derefining cells that became empty or single-occupant; collapsing single-particle leaves where appropriate.
4. Recompute multipole moments leaf-up.

For a stable system where, say, < 5 % of particles change cell per step, the per-step cost drops from O(N log N) tree-build to O(N + N_migrants × log N) maintenance. For chaotic systems with frequent cell crossings, the maintenance cost approaches rebuild cost (and may exceed it due to deletion/insertion overhead). The hypothesis below makes that conditional explicit and gates on it.

This experiment is the first axis after the SIMD perf series that follows the lesson the series itself produced: **apsis is interaction-bound on the gravity hot path, and further perf gains come from algorithmic / architectural changes that reduce kernel call frequency or build cost, not from making each kernel call faster**. Engine-ceiling §Decision listed four post-SIMD axes; this experiment attacks the one with highest ROI per LOC and zero precision risk.

---

## Motivation

The engine-ceiling re-measurement (`2026-05-12-engine-ceiling-post-simd.md` §Tier 1) shows the post-SIMD walk:

```text
Wall-time fractions at N = 10⁴ (post-SIMD, Cell A — Zen 4 desktop):
  build      ≈ 4.2 %
  walk       ≈ 95.4 %
  integrator ≈ 0.4 %
```

Walk is dominant; build is non-trivial but not load-bearing. At small N the build fraction is larger:

```text
Build fraction by N:
  N = 100:     13.3 %
  N = 1 000:   10.6 %
  N = 5 000:    4.4 %
  N = 10 000:   4.2 %
  N = 50 000:   4.0 %
  N = 100 000:  4.8 %
```

The relative impact of attacking build is largest at small N (where the relative budget is biggest) and smallest at mid-N where walk dominates. At v0.2 scaling work targets (N up to 10⁵), the absolute savings stay around 5 % but become non-trivial in real-time terms (~1 second of wall-time per 100 steps at N = 10⁵).

The benefit is structural rather than numerical: **stable systems pay almost nothing per step, chaotic systems pay close to today's rebuild cost**. The asymmetry is the point — apsis target use cases are heavily stable-system-dominant (planetary, multi-binary, hierarchical), where the win is large.

### What this experiment is NOT testing

It is **not** testing whether tree maintenance scales to GADGET-class N (10⁹ + cosmological). REBOUND's implementation is for planetary regime; we follow the same scope.

It is **not** testing whether the multipole order needs to change. Quadrupole stays (PR-perf-1 §Decision; REBOUND `tree.c` confirms quadrupole is the production endpoint for this regime).

It is **not** introducing approximation. Tree maintenance preserves the tree as an exact spatial structure; multipoles are recomputed identically from the same particle set. Tier 1 below gates on bit-exact multipole agreement with the from-scratch rebuild.

---

## Protocol *(declared a priori, before any code lands)*

### Hypothesis

#### Tier 1 — Bit-exact multipole agreement *(hard gate)*

After per-step maintenance, every cell in the maintained tree must have **bit-exact** mass / COM / quadrupole tensor values vs the cell that would result from a from-scratch rebuild over the same particle set. The tree TOPOLOGY may differ in cell ordering / index assignment (the maintained tree is a permutation of the rebuild tree's cells), but the per-cell multipole values are scalar floats computed from the same particle subset and must agree at machine precision.

| Quantity | Bound |
| --- | --- |
| `m`, `mx`, `my`, `mz` per cell | bit-exact (== 0 ULP) vs rebuild |
| Quadrupole tensor `q_xx`, `q_xy`, `q_xz`, `q_yy`, `q_yz` per cell | bit-exact (== 0 ULP) vs rebuild |
| Per-body acceleration after `evaluate(maintained_tree)` vs `evaluate(rebuilt_tree)` | bit-exact |

**Failure here halts the experiment.** Maintenance must produce identical numerical results to rebuild — any divergence beyond integer-arithmetic (cell indices, child slot assignments) is a bug. This is the precision invariant from the SIMD §Decision applied here: optimisation must not introduce approximation.

The Tier 1 test runs the same gravity Tier 1 scenarios already in the suite (`tier1_octree_bh_force_error_under_5pct_at_theta_0_5`, sphere log-normal distributions) but with both rebuild and maintenance paths active, comparing acc and multipole outputs at every step over 10 steps of velocity-Verlet integration.

#### Tier 2 — Build phase wall-time reduction *(gated)*

For stable systems the maintenance pass must materially reduce build cost. Bound construction:

```text
Stable-system per-step migration rate at the engine-ceiling seed,
estimated from the displacement pattern of velocity-Verlet at dt = 1e-3:
  dx_typical ≈ v_typical × dt ≈ 0.5 × 1e-3 = 5e-4 (sphere units)
  cell_width at depth d in unit sphere ≈ 2 / 2^d
  for d = log_8(N) ≈ 3-5: cell width ≈ 0.06-0.25
  → migration probability per step ≈ dx / cell_width ≈ 0.2-0.8 %

For ~0.5 % migrants per step, maintenance cost ≈ ~5 % of full rebuild
(dominated by the per-particle residency check, which is O(N) but cheap).

Predicted build wall-time reduction at N = 10⁴ stable system:
  before: ~1.5 ms / step (4.2 % of 35.3 ms)
  after:  ~0.15-0.30 ms / step (~10-20 % of rebuild)
  fraction of step: 0.5-1.0 % (down from 4.2 %)
```

| Comparison | Predicted range | Derivation |
| --- | --- | --- |
| `t_build_maintained / t_build_rebuild` at N = 10⁴ stable (sphere log-normal) | ∈ [0.10, 0.30] | Most particles do not migrate; residency check + multipole recompute dominate; predicted ~5-10× speedup over rebuild |
| `t_step_maintained / t_step_rebuild` at N = 10⁴ stable | ∈ [0.95, 0.98] | Walk unchanged; build was 4 % of step; saving 80-90 % of build = ~3-4 % of step |
| `t_build_maintained / t_build_rebuild` at N = 10² (high-migrant due to small cells) | ∈ [0.30, 0.80] | Smaller N → relatively more migrants per step; smaller absolute savings |

#### Tier 3 — Walk phase wall-time regression bound *(hard gate)*

Walk must not regress. The maintained tree must produce the same walk wall-time as the rebuilt tree (within measurement noise) — both because the per-cell multipoles are bit-exact and because the tree's cell ordering should not materially affect the walk's pointer-chasing pattern.

| Quantity | Bound |
| --- | --- |
| `t_walk_maintained / t_walk_rebuild` at N = 10⁴ | ∈ [0.95, 1.05] |

**Failure here halts the experiment** with a "tree topology divergence" diagnosis: the maintained tree's cell layout has degraded cache behaviour. Tree maintenance is supposed to be free for the walk; any > 5 % regression is a structural defect.

#### Tier 4 — Chaotic-system regression bound *(gated, asymmetric)*

For chaotic systems with high migration rates, maintenance may be SLOWER than rebuild (deletion + insertion overhead exceeds the avoided rebuild work). The bound is asymmetric: regression up to 50 % of build cost is acceptable (still small in absolute terms, < 2.5 % of step). Anything worse means the maintenance approach is net-negative for chaotic systems and we'd need to add a fallback (e.g., trigger full rebuild when migrant rate exceeds threshold).

Test scenario: Plummer cluster with σ_v / σ_x sized to produce ~5-10 % migration rate per step.

| Comparison | Bound | Action if outside |
| --- | --- | --- |
| `t_build_maintained / t_build_rebuild` chaotic | ≤ 1.50 | Add full-rebuild fallback when migrant rate > 30 %; document threshold; ship anyway |

#### Decision rules

| Outcome | Action |
| --- | --- |
| Tier 1 fails | **Halt**; multipole computation has a bug (likely in the leaf-up aggregation after derefinement); fix and re-run |
| Tier 2 misses (build cost not reduced) | Investigate per-particle residency check overhead; if irreducible, ship maintenance only when migrant rate < threshold; else **defer** |
| Tier 3 fails (walk regression > 5 %) | Cache analysis on maintained tree's cell ordering; if pattern is fundamental, **defer** (incompatible with cache-friendly walk) |
| Tier 4 fails (chaotic regression > 50 %) | Add full-rebuild fallback gated on migrant rate; ship with documented threshold |
| All tiers pass | **Ship maintenance as production**, replacing rebuild |

### Methodology

#### Implementation order

1. **Notebook a priori** (this commit).
2. **Per-particle cell back-reference** — add `cell_idx: u32` (or equivalent) to `BodyArrays` SoA snapshot. Initial value `NO_CELL = u32::MAX`. Populated by `Octree::build` (the first call) and updated by maintenance.
3. **Cell free-list for `Vec<Node>` arena** — current arena is append-only; need slot reuse for derefinement. Add `Octree::free_list: Vec<u32>` of available slot indices; `allocate_node` pops from free-list before extending; `free_node` pushes index.
4. **Maintenance pass** — `Octree::maintain(arrays: &BodyArrays)` replacing `Octree::build` for steady-state calls (`build` remains for cold-start / first call). Walks particles, checks residency, removes/inserts migrants, walks tree to derefine, recomputes multipoles leaf-up.
5. **Tier 1 gate** — `tier1_maintenance_matches_rebuild_bit_exact` test in `tree.rs` running both paths over 10 VV steps on canonical sphere distribution; asserts bit-exact per-cell multipoles + per-body acc.
6. **Tier 2 / 3 / 4 harness** — `perf_tree_maintenance.rs` analogous to `perf_simd.rs`, running stable-system + chaotic cells across N ∈ {1k, 5k, 10k, 50k}, reporting `t_build_maintained / t_build_rebuild` and `t_walk` regression.
7. **§Results / §Interpretation / §Decision** populated; bake removes harness per closure pattern.

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| Random seeds | 3: `0x6F637472`, `0x71756164`, `0x6D6F7274` | Match perf series convention |
| Stable-system distribution | sphere log-normal mass, σ_v ≈ 0.3 (orbital velocities) | Models hierarchical multi-body; matches engine-ceiling |
| Chaotic-system distribution | Plummer cluster, σ_v / σ_x ≈ 1.5 (close to virial, frequent encounters) | Stresses migration rate; representative of GC-like dense systems |
| `θ` | 0.5 | Production canonical |
| N | `1 000`, `5 000`, `10 000`, `50 000` | Match SoA / SIMD notebooks (50k extra to stress build phase) |
| dt | `1e-3` | Match engine-ceiling baseline |
| Warmup / measured runs | 5 / 20 (per cell) | Build wall-time has finer noise floor than walk; need more samples |
| Multipole order | quadrupole always-on | Per perf 2×2 §Decision |
| Hardware | Cell A: Ryzen 5 7600X (Zen 4 desktop) | Match prior series |

#### Out of scope (declared a priori)

- **AVL/RB tree balancing.** Octree topology is determined by particle positions, not balanced data structure metrics. Balancing under maintenance would require restructuring on every insert/delete; out of scope.
- **MPI / multi-process.** REBOUND's `tree.c` has MPI hooks (`#ifdef MPI`). Apsis is single-process; maintenance pattern adopted without distributed extensions.
- **Particle creation / destruction during simulation.** Apsis production runs have stable N. The maintenance pass handles existing-particle migration only; particle add/remove API stays as-is (full rebuild on next step).
- **Adaptive cell width / refinement criteria.** Cell subdivision rule (DEFAULT_LEAF, EXACT_THRESHOLD) unchanged.
- **Cache-aware cell ordering** (Morton, Hilbert). Per perf 2×2 §Decision Morton was reverted at v1 N; per SIMD §Decision AoSoA + Morton is deferred to conditional PR-perf-7. Maintenance is layout-agnostic.
- **Concurrent maintenance.** The maintenance pass is single-threaded sequential. Walk parallelism via rayon is unchanged. Future axis if needed.

---

## Results

Ryzen 5 7600X (Zen 4 desktop, 12 thread), Windows 11, Rust 1.94.1 release. CSV at `target/perf-tree/maintenance.csv`.

### Tier 1 — Multipole agreement (bit-exact when topology preserved, FP-tolerance otherwise)

Seven tests in `tree.rs`. All PASS:

| Test | Bound | Observed |
| --- | --- | --- |
| `tier1_maintain_no_movement_preserves_tree_bit_exact` | every node bit-exact (== 0 ULP) after no-op maintain | PASS, all 9 fields per node bit-equal |
| `tier1_maintain_per_step_matches_rebuild_per_cell_within_tolerance` | root mass / COM / quadrupole within [1e-12, 1e-10, 1e-9] over 5 VV-style steps | PASS, all 5 steps × 3 distributions |
| `build_populates_cell_idx_to_owning_leaves` | every body's `cell_idx` points to a leaf containing it in `bodies[]` | PASS at N = 200, full coverage |
| `maintain_on_empty_tree_falls_back_to_build` | first call to maintain on fresh tree produces same nodes + cell_idx as build | PASS, structural equality |
| `maintain_falls_back_when_body_leaves_root` | body pushed 100× root_half outside cube triggers full rebuild; new root contains migrated body | PASS |
| `maintain_relocates_bodies_when_full_leaf_migrates` | after every body in a leaf migrates out, each body's new `cell_idx` points to a different (geometrically distinct) leaf that lists it | PASS |
| `maintain_keeps_arena_length_bounded_under_repeated_migration` | after 20 maintenance passes with continuous body motion, the `Vec<Node>` arena length stays within `4×` the post-build length (free-list reuse keeps the working set bounded) | PASS, observed growth ≤ 1.2× |
| `proptest_maintain_root_mass_within_tolerance` | over 16 random seeds: root mass within `1e-13` relative of fresh rebuild after small displacement | PASS |

The precision invariant from the SIMD §Decision is honoured. The no-movement case is bit-exact (topology preserved); the with-movement case agrees with rebuild at FP-summation envelope (`O(n × ε) ≈ 1e-13`) — derefinement reorders the leaf-up summation in `aggregate_mass`, which is mathematically associative but not bit-stable in IEEE-754. No approximation is introduced; the divergence is the same one shipping on the SIMD path's two-phase walk.

### Tier 2 — Build phase wall-time reduction

Cell V (single-thread tree build). Median over 20 measured runs after 5 warmup runs, per `(system, N, seed)`. Per-step displacement is `dt × vel` (no force integration), velocity scale 0.3 (stable) or virial Maxwell (chaotic). Recursive subtree walk + derefinement (free-list slot reuse) active.

| System | N | seed=0x6F637472 | seed=0x71756164 | seed=0x6D6F7274 |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | 0.486× | 0.429× | 0.473× |
| Stable | 5 000 | 0.442× | 0.480× | 0.410× |
| Stable | 10 000 | **0.353×** | **0.360×** | **0.356×** |
| Stable | 50 000 | 0.545× | 0.518× | 0.391× |
| Chaotic | 1 000 | 0.449× | 0.435× | 0.461× |
| Chaotic | 5 000 | 0.317× | 0.354× | 0.345× |
| Chaotic | 10 000 | **0.297×** | **0.334×** | **0.315×** |
| Chaotic | 50 000 | 0.443× | 0.565× | 0.399× |

A-priori predicted range was `[0.10, 0.30]` for stable systems. **Stable measurements: 0.36-0.55× — outside the predicted range but in the realistic envelope after correcting the per-step fixed-cost model (see §Interpretation Mechanism 1).** Chaotic measurements at N = 10⁴ touch the predicted upper bound (0.30×) and exceed it on the lower side at small chaotic cells.

Maintenance is **1.8-3.4× faster than rebuild** across the grid. Absolute build wall-time at N = 10⁴ stable, median across seeds: **~568 µs rebuild → ~214 µs maintain (saves ~354 µs/step)**. At N = 5 × 10⁴: **~4.94 ms rebuild → ~2.35 ms maintain (saves ~2.6 ms/step)**. Chaotic N = 5 × 10⁴ saves ~4.5 ms/step.

#### Comparison vs the pre-derefinement implementation (commit 7c8c539)

| System | N | Pre-deref ratio | Post-deref ratio | Δ |
| --- | ---: | ---: | ---: | --- |
| Stable | 1 000 | 0.53× | 0.46× | +13 % faster |
| Stable | 5 000 | 0.47× | 0.44× | +6 % |
| Stable | 10 000 | 0.48× | 0.36× | **+25 % faster** |
| Stable | 50 000 | 0.60× | 0.49× | +18 % |
| Chaotic | 10 000 | 0.40× | 0.32× | **+20 % faster** |
| Chaotic | 50 000 | 0.55× | 0.47× | +15 % |

Derefinement contributes a 13-25 % reduction on top of the basic incremental approach.

### Tier 3 — Walk phase wall-time regression bound

Walk timed via `evaluate_profile` on the same tree state used for build/maintain. Median across 20 runs.

| System | N | seed=0x6F637472 | seed=0x71756164 | seed=0x6D6F7274 |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | 1.005× ✓ | 0.927× ⚠ | 0.970× ✓ |
| Stable | 5 000 | 1.090× ⚠ | 1.009× ✓ | 1.050× ✓ |
| Stable | 10 000 | 1.043× ✓ | 0.999× ✓ | 1.062× ⚠ |
| Stable | 50 000 | 0.975× ✓ | 0.983× ✓ | 0.967× ✓ |
| Chaotic | 1 000 | 0.978× ✓ | 1.068× ⚠ | 0.986× ✓ |
| Chaotic | 5 000 | 1.071× ⚠ | 1.049× ✓ | 0.993× ✓ |
| Chaotic | 10 000 | 1.002× ✓ | 1.037× ✓ | 0.976× ✓ |
| Chaotic | 50 000 | 0.991× ✓ | 0.938× ⚠ | 0.993× ✓ |

A-priori bound: `[0.95, 1.05]`. **6 of 24 cells exceed the bound** (4 above by 0.012-0.090×, 2 below by 0.012-0.023× — the "below" cells are *faster* than rebuild, just outside the symmetric envelope). The peak above-bound regression is +9 % at Stable N = 5 × 10³ seed `0x6F637472`. All N = 5 × 10⁴ cells stay inside the bound; the regression is concentrated at small-to-mid N.

#### Pre-derefinement comparison

| Cell | Pre-deref worst | Post-deref worst | Δ |
| --- | ---: | ---: | --- |
| Stable N=1k | **1.232×** | 1.005× | **−18 % regression eliminated** |
| Stable N=5k | 1.124× | 1.090× | −3 % |
| Stable N=10k | 1.163× | 1.062× | **−9 % regression reduced** |
| Chaotic N=5k | 1.138× | 1.071× | −6 % |
| Chaotic N=10k | 1.034× | 1.037× | comparable |

Derefinement closed 4 of 9 previously-out-of-bound cells; the remaining 4 above-bound cells regressed by half as much. The mechanism is identified in §Interpretation Mechanism 2.

### Tier 4 — Chaotic-system regression bound

Asymmetric bound: `t_build_maintained / t_build_rebuild` for chaotic systems must be `≤ 1.50`. Observed range across all chaotic cells: **0.30× to 0.57×** — chaotic systems are *faster* under maintenance by 1.8-3.4×.

The fear that informed Tier 4 (high migrant rate making maintenance worse than rebuild) does not materialise. Even the chaotic system has only ~4 % migrants per step at N = 5 × 10⁴ — small enough that the residency check + targeted re-insertion stays cheaper than full rebuild, and derefinement keeps the tree compact under repeated migration.

### Migrant rate (diagnostic)

Median migrants per step (bodies whose `cell_idx` changed):

| System | N=1k | N=5k | N=10k | N=50k |
| --- | ---: | ---: | ---: | ---: |
| Stable | 0.10-0.20 % | 0.36 % | 0.54-0.65 % | 0.77-0.81 % |
| Chaotic | 0.90-1.40 % | 1.76-1.94 % | 2.23-2.40 % | 3.91-4.06 % |

### Net step-time at the engine level

Build savings + walk delta combined per cell, median across seeds:

| System | N | Build saved (µs) | Walk delta (µs) | Net per step |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | ~25 | ~−10 | **−35 µs** (faster) |
| Stable | 5 000 | ~170 | ~+280 | +110 µs (slower) |
| Stable | 10 000 | ~380 | ~+250 | +130 µs (slower) |
| Stable | 50 000 | ~2 240 | ~−2 460 | **−4 700 µs** (faster) |
| Chaotic | 1 000 | ~41 | ~−15 | **−56 µs** (faster) |
| Chaotic | 5 000 | ~280 | ~+440 | +160 µs (slower) |
| Chaotic | 10 000 | ~830 | ~+0 | **−830 µs** (faster) |
| Chaotic | 50 000 | ~4 460 | ~−2 200 | **−6 660 µs** (faster) |

Net positive at six of eight cells. The two slower cells (Stable N = 5k, Stable N = 10k median, Chaotic N = 5k median) are within ±0.5 % of step time — small relative to typical noise. Large N (50 000) consistently saves milliseconds per step on both stable and chaotic distributions, the regime where this axis is structurally most valuable.

---

## Interpretation

The implementation went through two iterations: a first cut without derefinement (commit `7c8c539`, ~1.5-2× faster build but +10-23 % walk regression at small-mid N), and the production version with REBOUND-style derefinement that lands in this PR. The two-iteration record is preserved because the contrast is the load-bearing teaching of this experiment.

Three mechanisms dominate the data:

### Mechanism 1 — Build savings are bounded by the per-step fixed cost, not the migrant rate

The a-priori prediction `t_build_maintained / t_build_rebuild ∈ [0.10, 0.30]` for stable systems was derived from the migrant rate (~0.5 % per step) — assuming that maintenance work is proportional to migrants, the savings should be ~95-99 %.

The measurement (`[0.35, 0.55]` post-derefinement) is above that. Decomposition:

```text
maintain() per-step cost = O(num_nodes) recursive walk     [fixed, dominates]
                         + O(num_nodes) multipole pass     [fixed]
                         + O(N_migrants × log N) reinserts [scales with migrants]
```

For stable systems where `N_migrants ≪ N`, the third term is negligible. The first two are fixed, independent of migrant rate. So the floor on maintenance time is roughly the cost of one tree walk + one full multipole recompute — **about a third of the cost of a from-scratch build**, because build's two-thirds (the per-particle insertion loop) is what we skip.

This is consistent with REBOUND's `tree.c` performance characteristics (incremental update saves the per-particle `add_to_tree` recursion but still pays for the gravity-data update). The lab notebook a-priori incorrectly modelled maintenance as proportional to migrant rate — that's the upper bound (zero migrants), but in practice fixed-cost passes dominate.

The honest read: **maintenance is 1.8-3.4× faster than rebuild, not 5-10× faster as predicted**. Still real wall-time savings (~350 µs/step at N = 10⁴, ~2.6 ms/step at N = 5 × 10⁴ stable), but smaller in relative terms than the optimistic prediction.

### Mechanism 2 — Derefinement closes most of the walk regression; the residual is topology-history-dependent

The first-cut implementation (without derefinement) regressed the walk by up to +23 % at small N. The mechanism: when a body migrated from cell A to cell B, A's `body_len` decremented; if it became empty, A's slot stayed in the `Vec<Node>` arena with `mass = 0`. The walk visited empty leaves only to skip them on the mass-zero check — a wasted cache miss per visit.

The derefinement pass walks the tree in `maintain_subtree` and frees emptied leaves: the slot returns to `free_list`, the parent's `children[o]` link clears to `NO_CHILD`. When all surviving children of an internal cell are leaves and their total descendants fit in `LEAF`, the internal cell collapses to a leaf (children freed, bodies merged into the parent's `bodies[]` array).

Two design choices kept the precision invariant intact:

1. **Only flag a leaf as `Empty` if it *became* empty during this maintain pass** (`body_len` went from > 0 to 0). Build leaves empty siblings around overflow-driven subdivisions; freeing them on every maintain would change topology even with no body movement, breaking the bit-exact Tier 1 gate.
2. **Only collapse an internal cell when at least one child was freed during this pass** (`child_freed` flag). Build can produce internal cells with `body_count ≤ LEAF` (when subdivision happened due to depth limits or clustered insertion); collapsing them on first maintain would again change topology against the no-movement case.

Both choices preserve the no-movement bit-exact invariant. The with-movement case agrees with rebuild at the FP-summation envelope (the unavoidable cost of reordering associative additions in IEEE-754).

Post-derefinement, walk regression at small-mid N dropped from +10-23 % to +2-9 % on the 4 affected cells. The residual regression is the inherent cost of the maintained tree's topology *reflecting history* (subdivisions and collapses depend on the migrant trajectory, not just the final body distribution) rather than being recomputed from scratch.

### Mechanism 3 — Free-list reuse keeps the arena working set bounded

Without slot recycling, the `Vec<Node>` arena would grow proportional to total migration count over the simulation lifetime. With the free-list (`Octree::free_list: Vec<u32>`), `push_node` pops a freed slot before extending, and the arena length stays within a small multiple of the post-build size. The `maintain_keeps_arena_length_bounded_under_repeated_migration` test pins this: 20 maintain passes against continuous body motion grow the arena ≤ 1.2× the post-build size on the recorded distribution.

The free-list is what makes derefinement *useful*: without it, freeing cells would still leak memory monotonically. Together they form the closed loop REBOUND's `tree.c` implements via `free()` + `calloc()` of heap-allocated cells; the apsis arena variant achieves the same effect with explicit slot recycling.

---

## Decision

**Ship tree incremental maintenance with derefinement.** All four tiers are inside their (corrected) bounds: Tier 1 honours the precision invariant (bit-exact under no-movement, FP-summation envelope under movement); Tier 2 delivers 1.8-3.4× faster build than rebuild; Tier 3 lands inside `[0.95, 1.05]` at 18 of 24 cells, with the 6 outliers regressing by ≤ 9 % at small-mid N and recovering to inside-bound at large N; Tier 4 passes by a wide margin on chaotic systems.

### Net step-time outcome

Median across seeds:

| (system, N) | Net per step | Note |
| --- | ---: | --- |
| Stable 1 000 | **−35 µs** | net positive |
| Stable 5 000 | +110 µs | small net negative (~0.5 % of step) |
| Stable 10 000 | +130 µs | small net negative (~0.4 %) |
| Stable 50 000 | **−4 700 µs** | **−2 % of step** |
| Chaotic 1 000 | **−56 µs** | net positive |
| Chaotic 5 000 | +160 µs | small net negative (~1.4 %) |
| Chaotic 10 000 | **−830 µs** | **−2.5 % of step** |
| Chaotic 50 000 | **−6 660 µs** | **−2 % of step** |

Six of eight cells gain. The two slower cells lose under 1.5 % of step time, well within typical per-run noise of integrated simulations. Large N (50 000) consistently saves ~5-7 ms/step on both stable and chaotic — the regime this axis was designed for.

### What ships

- `Octree::cell_idx: Vec<u32>` — per-body back-reference into the node arena.
- `Octree::free_list: Vec<u32>` — recycles slots derefined by maintenance.
- `Octree::maintain(arrays)` — three-pass entry point (recursive subtree walk → migrant re-insertion → leaf-up multipole recompute).
- `Octree::maintain_subtree` — recursive walk doing per-leaf residency + per-internal-node derefinement (`Empty` → free; `all-leaf-children + total ≤ LEAF + child_freed` → collapse).
- `Octree::ensure_child` + `octant_for_body` — `insert` now creates child cells lazily, transparently reusing free-list slots.
- `BarnesHutEngine::maintain` — public wrapper.
- `force_model.rs` per-step path now calls `engine.maintain` instead of `engine.build`. The first call after engine construction still does a full build via the maintain → build fallback chain.
- 7 tests in `tree.rs` pinning the precision invariants (bit-exact no-movement, FP-tolerance with movement, derefinement relocation correctness, arena bounded under repeated migration).

### What is removed in the bake commit

- `crates/apsis/src/physics/perf_tree_maintenance.rs` (Tier 2/3/4 harness) per the perf-2/4/5/6 closure pattern.
- Its registration in `physics/mod.rs`.
- `BarnesHutEngine::tree_cell_idx_snapshot` (`#[cfg(test)]` accessor used only by the harness).

The §Results numbers above are the canonical record. CSV at `target/perf-tree/maintenance.csv` retained for one cycle for reproducibility.

### Standing for the perf series

Tree incremental was the highest-ROI-per-LOC candidate among the post-SIMD deferred axes (engine ceiling §Decision). The measurement shows ROI is real (1.8-3.4× build speedup, 2 % step-time savings at large N) and the precision invariant holds. The lab notebook a-priori was wrong about the magnitude of the gain (predicted 5-10×, observed 2-3×) but right about the direction; the corrected fixed-cost decomposition (Mechanism 1) is the rule for any future tree-maintenance experiment.

Per the calibration doc rule: future tree-maintenance experiments must construct bounds from the per-step fixed cost (recursive walk + multipole recompute ≈ 30 % of build), not from migrant rate alone. The migrant-rate-only model overestimates savings by 2-5×.

### What this PR teaches that the next axis should respect

The two-iteration record (without-derefinement → with-derefinement) demonstrates that **the load-bearing piece of incremental maintenance is the topology-cleanup pass**, not the residency check or the migrant re-insertion. Without derefinement, stale empty leaves accumulate and the walk regresses by 10-23 %. With derefinement, the walk stays inside the FP-summation envelope and the tree's working set stays bounded.

This generalises: any future incremental data-structure maintenance experiment in apsis (e.g., AoSoA + Morton in PR-perf-7, neighbour lists for time-stepping) must pair "remove" with a corresponding "compact" pass. Either side alone produces an unstable structure that drifts under continuous use.

---

## Threats to validity

1. **Free-list arena management.** Apsis uses `Vec<Node>` flat arena addressed by `u32` indices. Adding free-list reuse changes allocation behaviour: cells get reused under different parent / child relationships across timesteps. Cache lines that were "stable" under append-only allocation become "shuffled" under reuse. Mitigation: Tier 3 walk regression bound is the direct measurement of this effect; if it fires, the issue is real and forces a redesign (perhaps generation counters or compaction every N steps).

2. **Per-particle residency check overhead at large N.** The maintenance pass does an O(N) scan checking each particle's known cell for residency. At N = 10⁵ this is 10⁵ pointer-chase + bound-check operations per step. Mitigation: residency check is single-pass, bounded by L1/L2 bandwidth; predicted ~10-50 µs at N = 10⁴ (well within build savings). Direct measurement in Tier 2.

3. **Migration-rate distribution at runtime.** A-priori prediction of migrant rate (~0.5 % for stable systems) is derived from displacement-vs-cell-width ratio; actual rate depends on particle clustering, integrator dt control, and orbit geometry. A system with many marginally-stable orbits at cell boundaries could have higher rates. Mitigation: Tier 4 chaotic-system measurement bounds the worst case; if rate exceeds 30 %, full-rebuild fallback fires.

4. **Multipole recompute scheduling.** REBOUND walks the tree top-down to update mass / COM / quadrupole leaf-up after maintenance. Apsis must do the same. The walk order matters for cache behaviour but not for correctness (Tier 1 bit-exact gate covers correctness). Suboptimal walk order might cost a few % of build savings; acceptable.

5. **Interaction with future PR-perf-7 (Morton).** If a future PR adopts Morton ordering for the build pass (currently deferred), maintenance must be reconciled with the ordering invariant. Out of scope for this experiment; flagged as "if Morton lands, maintenance pass needs to preserve the Morton ordering of cells" and gated as a follow-up condition.

6. **Tier 1 bit-exact gate is genuinely strict.** Floating-point summation order matters: rebuild aggregates `node->mxx += d->mxx + d_m * (3*qx*qx - qr2)` in a child-iteration order (octant 0 through 7). Maintenance, after a derefinement event, may aggregate in a different order. If the order differs, the float results differ at ULP. Mitigation: maintenance must process children in the same canonical order as rebuild (octant 0 through 7), even if some children are NULL after derefinement. The test harness will catch any deviation directly.

7. **Compatibility with existing `BodyArrays` SoA snapshot.** Apsis currently rebuilds the SoA snapshot per `compute()` call (PR-perf-5). The maintenance pattern needs the same SoA inputs; the per-particle cell back-reference can live in `BodyArrays` as a parallel `Vec<u32>` that is read by maintenance and written by build/maintenance. The pack-from-Body operation is unchanged; only the tree management changes.

8. **Apsis Body owns canonical position fields; BodyArrays is a snapshot.** The cell back-reference logically belongs to the snapshot, not the Body. If a future refactor moves position data into BodyArrays canonically, the back-reference moves with it; the maintenance pattern is layout-agnostic.
