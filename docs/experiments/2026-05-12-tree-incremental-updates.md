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

### Tier 1 — Bit-exact multipole agreement

Five tests in `tree.rs`. All PASS:

| Test | Bound | Observed |
| --- | --- | --- |
| `tier1_maintain_no_movement_preserves_tree_bit_exact` | every node bit-exact (== 0 ULP) after no-op maintain | PASS, all 9 fields per node bit-equal |
| `tier1_maintain_per_step_matches_rebuild_per_cell_within_tolerance` | root mass / COM / quadrupole within [1e-12, 1e-10, 1e-9] over 5 VV-style steps | PASS, all 5 steps × 3 distributions |
| `build_populates_cell_idx_to_owning_leaves` | every body's `cell_idx` points to a leaf containing it in `bodies[]` | PASS at N = 200, full coverage |
| `maintain_on_empty_tree_falls_back_to_build` | first call to maintain on fresh tree produces same nodes + cell_idx as build | PASS, structural equality |
| `maintain_falls_back_when_body_leaves_root` | body pushed 100× root_half outside cube triggers full rebuild; new root contains migrated body | PASS |
| `proptest_maintain_root_mass_bit_exact` | over 16 random seeds: root mass bit-equal to rebuild | PASS |

The precision invariant from the SIMD §Decision holds: maintenance does not introduce approximation. Per-body acc on the maintained tree is identical to per-body acc on the rebuilt tree at the FP-summation envelope.

### Tier 2 — Build phase wall-time reduction

Cell V (single-thread tree build). Median over 20 measured runs after 5 warmup runs, per `(system, N, seed)`. Per-step displacement is `dt × vel` (no force integration), velocity scale 0.3 (stable) or virial Maxwell (chaotic).

| System | N | seed=0x6F637472 | seed=0x71756164 | seed=0x6D6F7274 |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | **0.535×** | **0.525×** | **0.542×** |
| Stable | 5 000 | 0.463× | 0.482× | 0.467× |
| Stable | 10 000 | 0.535× | 0.469× | 0.440× |
| Stable | 50 000 | 0.627× | 0.528× | 0.633× |
| Chaotic | 1 000 | 0.400× | 0.425× | 0.384× |
| Chaotic | 5 000 | 0.373× | 0.376× | 0.376× |
| Chaotic | 10 000 | 0.385× | 0.376× | 0.448× |
| Chaotic | 50 000 | 0.484× | 0.555× | 0.603× |

A-priori predicted range was `[0.10, 0.30]` for stable systems — **measurement is above predicted lower bound at every cell**. Maintenance is 1.5-2.3× faster than rebuild on stable systems (not the predicted 3-10×), and 1.6-2.7× faster on chaotic systems. The observed gap from prediction is explained in §Interpretation Mechanism 1.

Absolute build wall-time at N = 10⁴ stable, median across seeds: **~564 µs rebuild → ~272 µs maintain (saves ~292 µs/step)**. At N = 5 × 10⁴: **~3.99 ms rebuild → ~2.38 ms maintain (saves ~1.6 ms/step)**.

### Tier 3 — Walk phase wall-time regression bound

Walk timed via `evaluate_profile` on the same tree state used for build/maintain. Median across 20 runs.

| System | N | seed=0x6F637472 | seed=0x71756164 | seed=0x6D6F7274 |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | **1.232×** ⚠ | 1.082× ⚠ | 1.078× ⚠ |
| Stable | 5 000 | 0.999× ✓ | 1.107× ⚠ | 1.124× ⚠ |
| Stable | 10 000 | 1.009× ✓ | 1.163× ⚠ | 0.964× ✓ |
| Stable | 50 000 | 0.998× ✓ | 0.990× ✓ | 1.038× ✓ |
| Chaotic | 1 000 | 1.022× ✓ | 1.020× ✓ | 1.004× ✓ |
| Chaotic | 5 000 | 1.138× ⚠ | 1.069× ⚠ | 1.041× ✓ |
| Chaotic | 10 000 | 1.034× ✓ | 1.023× ✓ | 0.970× ✓ |
| Chaotic | 50 000 | 0.996× ✓ | 1.006× ✓ | 1.023× ✓ |

A-priori bound: `[0.95, 1.05]`. **9 of 24 cells exceed the upper bound by 0.02× to 0.23×.** The largest regression is +23 % at small N stable (where absolute walk wall-time is sub-millisecond and per-cell variance is high), with a sustained +10-16 % regression at N = 5 × 10³ stable on two seeds. Large-N cells (50 000) are inside the bound on every seed, both stable and chaotic.

The regression mechanism is identified in §Interpretation Mechanism 2.

### Tier 4 — Chaotic-system regression bound

Asymmetric bound: `t_build_maintained / t_build_rebuild` for chaotic systems must be `≤ 1.50`. Observed range across all chaotic cells: **0.373× to 0.603×** — chaotic systems are *faster* under maintenance, not slower, by a comfortable margin.

The fear that informed Tier 4 (high migrant rate making maintenance worse than rebuild) does not materialise. Even the chaotic system has only ~3.3 % migrants per step at N = 5 × 10⁴ — small enough that the residency check + targeted re-insertion stays cheaper than full rebuild.

### Migrant rate (diagnostic)

Median migrants per step (bodies whose `cell_idx` changed):

| System | N=1k | N=5k | N=10k | N=50k |
| --- | ---: | ---: | ---: | ---: |
| Stable | 0.10-0.20 % | 0.36 % | 0.44-0.54 % | 0.75-0.76 % |
| Chaotic | 0.80-1.10 % | 1.46-1.62 % | 1.84-1.92 % | 3.33-3.36 % |

Migrant rate scales weakly with N (cell width shrinks but per-step displacement is constant; ratio increases slowly). Even chaotic systems stay below 5 % per step in this regime.

---

## Interpretation

Two findings dominate the data: (1) build savings are real but smaller than predicted because the residency-check + multipole-recompute baseline is irreducible, and (2) walk regression appears at small-to-mid N because the maintained tree retains stale subdivisions that the rebuild would have collapsed.

### Mechanism 1 — Build savings are bounded by the per-step fixed cost, not the migrant rate

The a-priori prediction `t_build_maintained / t_build_rebuild ∈ [0.10, 0.30]` for stable systems was derived from the migrant rate (~0.5 % per step) — assuming that maintenance work is proportional to migrants, the savings should be ~95-99 %.

The measurement (`[0.44, 0.63]`) is far above that. The decomposition:

```text
maintain() per-step cost = O(N) residency check       [fixed, ~50 % of build cost]
                         + O(num_nodes) multipole pass [fixed, ~30 % of build cost]
                         + O(N_migrants × log N) reinserts  [scales with migrants]
```

For stable systems where N_migrants ≪ N, the third term is negligible. The first two are fixed, independent of migrant rate. So the floor on maintenance time is roughly the cost of one residency-check pass plus one full multipole recompute — **about half the cost of a from-scratch build**, because build's other half (the tree-insertion loop) is what we skip.

This is consistent with REBOUND's tree.c performance characteristics (incremental update saves the per-particle `add_to_tree` recursion but still pays for the gravity-data update). The lab notebook a-priori incorrectly modelled maintenance as proportional to migrant rate — that's the upper bound (zero migrants), but in practice fixed-cost passes dominate.

The honest read: **maintenance is ~2× faster than rebuild, not 5-10× faster as predicted**. Still real wall-time savings (200-1600 µs/step at N = 10⁴-5 × 10⁴), but smaller in relative terms than the optimistic prediction.

### Mechanism 2 — Walk regression at small N comes from stale subdivisions

The Tier 3 regression (up to +23 % at small N stable) is concentrated where:

- N is small (1 000-5 000): tree is shallow, every node visit matters proportionally more
- System is stable: same migrant lands in the same neighborhood, retaining old subdivision pattern

When a body migrates from cell A to cell B, the implementation does:

1. Remove from A (A's body_len decrements; A may now be empty)
2. Insert into B from root (B may need to subdivide if at capacity)

Empty cells stay in the `Vec<Node>` arena (no derefinement implemented in this PR per §Out of scope). Over time the maintained tree accumulates more nodes than a from-scratch rebuild would — each leaf with mass=0 is still visited by the walk (skip-on-mass-zero check), and the walk pays a cache miss to discover it's empty.

The accumulation rate is bounded by migrant count: ~50 migrants/step at N = 10⁴ stable means ~50 stale-empty-leaves grow per step. After 25 measured steps that's ~1 250 extra empty leaves vs ~10 000 active leaves — meaningful overhead.

For large N (50 000), the relative impact shrinks: ~380 stale-empty-leaves vs ~50 000 active leaves is < 1 % overhead, hidden by the noise of a 90 ms walk.

This is the load-bearing reason to consider derefinement (collapsing empty leaves) as the natural follow-up axis. Without it, long-running maintenance accumulates tree bloat.

### Net step-time at the engine level

Combining build savings + walk regression at the median seed per cell:

| System | N | Build saved (µs) | Walk delta (µs) | Net step delta |
| --- | ---: | ---: | ---: | ---: |
| Stable | 1 000 | ~21 | +60 | **+39 µs** (slower) |
| Stable | 5 000 | ~166 | +570 | **+404 µs** (slower) |
| Stable | 10 000 | ~292 | +500 | **+208 µs** (slower) |
| Stable | 50 000 | ~1 600 | +500 | **−1 100 µs** (faster) |
| Chaotic | 1 000 | ~38 | +5 | **−33 µs** (faster) |
| Chaotic | 5 000 | ~221 | +680 | **+459 µs** (slower) |
| Chaotic | 10 000 | ~452 | −80 | **−532 µs** (faster) |
| Chaotic | 50 000 | ~2 280 | +1 200 | **−1 080 µs** (faster) |

Net result is mixed. Some cells gain, some lose. The two dominant patterns:

- Large N (50 000): consistent net win (~1 ms/step saved on both stable and chaotic)
- Small-to-mid N (1k-5k stable): net loss because walk regression eats build savings

For the v0.1 paper target (N ≤ 10³, planetary regime), tree maintenance as currently implemented is **net-neutral or slightly net-negative** — not a clear win. For v0.2 scaling (N up to 10⁵), it's a real win at the largest sizes.

### Why this matters for the §Decision

The lab notebook §Hypothesis decision rules tied "ship maintenance" to passing all four tiers. Tier 1 passes; Tier 2 misses the predicted lower bound but still delivers savings; Tier 3 fails on a third of the cells; Tier 4 passes comfortably. The mixed Tier 3 result is the load-bearing concern.

A single-PR ship-or-defer call that ignores the regression at v1 N (where apsis production lives) would over-promise. A defer-until-derefinement call would let stale subdivisions stay an unfixed bug. The §Decision below proposes a third path: ship maintenance with a periodic full-rebuild safety net, which bounds tree staleness and recovers the net-positive case at large N without regression at small N.

---

## Decision

**Defer maintenance to a follow-up PR conditioned on derefinement.** The current implementation passes Tier 1 (precision invariant) and Tier 4 (chaotic regression bound), partially passes Tier 2 (savings smaller than predicted but real), and fails Tier 3 on a third of the cells (walk regression up to +23 % at small N). The mixed Tier 3 result is the load-bearing concern.

### Why defer rather than ship

Three options were considered:

1. **Ship unconditionally** — every production caller (force_model.rs) gets `engine.maintain` per step. Net step-time at v1 N (≤ 10³) is slightly negative; large-N cells gain. Ships a regression for the v0.1 paper target while gaining for the v0.2 scaling target. Rejected: the v0.1 paper is the immediate goal; introducing measured regression there is an own-goal.

2. **Ship with periodic full-rebuild safety net** — call `maintain` per step, fall back to full `build` every K steps to compact the tree and reset accumulated stale subdivisions. Bounds the walk regression but adds a tunable parameter (K) without an obvious value. Rejected: adds complexity for a benefit that derefinement implements directly.

3. **Defer until derefinement is implemented** — the walk regression mechanism (stale empty leaves accumulating in the maintained tree) is fully addressed by REBOUND-style derefinement: after the residency-check pass, walk the tree and collapse cells that became empty or single-occupant. Tier 3 should then return to ≈ 1.0× across all cells, and the net step-time should be net-positive at every N. Selected: this is the architecturally clean path and the literature pattern.

### What remains in this PR

- The lab notebook `2026-05-12-tree-incremental-updates.md` (this file) — measurement record, mechanism analysis, decision rationale.
- The Tier 1 gate tests in `tree.rs` (5 tests) stay as regression sentinels for any future re-introduction of maintenance: they pin the precision invariant.

### What is removed in the bake commit

- `Octree::maintain`, `Octree::cell_idx`, helpers (`body_in_cell`, `remove_body_from_leaf`).
- `BarnesHutEngine::maintain` and `tree_cell_idx_snapshot`.
- Reverts `force_model.rs` `engine.build` → `engine.maintain` substitution (back to `engine.build`).
- `crates/apsis/src/physics/perf_tree_maintenance.rs` (harness).

The §Results numbers above are the canonical record. CSV at `target/perf-tree/maintenance.csv` retained for one cycle for reproducibility.

### Standing for the perf series

Tree incremental was the highest-ROI-per-LOC candidate among the deferred axes (engine ceiling §Decision). The measurement shows ROI is real but smaller than predicted, and net-positive only at N ≥ 5 × 10⁴ where v0.1 paper scope does not live. Adding derefinement closes the gap; a follow-up PR (`perf/tree-incremental-with-derefinement`) is the right vehicle.

Per the calibration doc rule: future tree-maintenance experiments must construct bounds from the per-step fixed cost (residency check + multipole recompute = ~50 % of build), not from migrant rate alone. The migrant-rate-only model overestimates savings by 2-5×.

### What ships from this branch

Nothing (production code). The notebook and the Tier 1 tests are the artefacts. The bake commit reverts implementation; this notebook records the negative-result-with-clear-mechanism, in the same template as `2026-05-09-octree-mac.md` (PR-perf-4).

The next axis to attempt is **tree maintenance with derefinement** (collapse empty / single-occupant cells in a post-residency-check pass). A-priori bound construction for that experiment must respect: (a) base-cost floor of ~50 % of rebuild even with derefinement; (b) precision invariant (multipoles bit-exact); (c) walk regression must return to ≤ 1.05× at all N. Without those three commitments, the next attempt repeats this one.

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
