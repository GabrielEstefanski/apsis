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

*To be populated incrementally as commits land.*

### Tier 1 — Bit-exact multipole agreement

*Pending.*

### Tier 2 — Build phase wall-time reduction (stable system)

*Pending.*

### Tier 3 — Walk phase wall-time regression bound

*Pending.*

### Tier 4 — Chaotic-system regression bound

*Pending.*

---

## Interpretation

*To be written after Tier 1-4 are populated.*

---

## Decision

*To be written after the Tier 1-4 gates pass or fail. The precision invariant (Tier 1 bit-exact gate) is the load-bearing one for ship/revert; the wall-time gates inform whether to ship unconditionally, gate on migrant-rate threshold, or defer.*

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
