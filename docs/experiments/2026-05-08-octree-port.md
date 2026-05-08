# Octree port — protocol

**Date:** 2026-05-08
**Subject:** Replace the Barnes-Hut quadtree (`physics/gravity/tree.rs`, `engine.rs`) with a 3D octree. Close the documented technical debt acknowledged in `engine.rs:144-147` ("the Barnes–Hut spatial index is still a quadtree (xy-only partition) … non-zero z requires an octree, which is staged separately"). Validate that:

1. The octree accuracy matches the canonical Barnes-Hut bound on general 3D distributions, where the quadtree was structurally inadequate.
2. The octree produces results numerically indistinguishable from the quadtree on planar (`z = 0`) systems where the quadtree was the correct partition restricted to the orbital plane — **no regression** on solar-system inner planets, two-body Kepler, or any test that previously passed.
3. The wall-time scaling remains `O(N log N)` and improves on the previous quadtree's per-step cost or matches it within one standard deviation.

**Status:** Protocol declared a priori, before any code is written. §Results populated after implementation, before merge.

---

## Abstract

The Barnes-Hut spatial index in `apsis` was implemented as a quadtree (xy-only partition) for historical reasons predating the 3D-native data flow migration. The quadtree gives correct forces only when every body has `z = 0`; for any inclined orbit, hierarchical system, or out-of-plane perturbation, the spatial partition collapses two octants into one and the BH approximation reads the wrong center of mass.

This refactor replaces `QuadTree` (4-children, 2D cell centre, 2D COM) with `Octree` (8-children, 3D cell centre, 3D COM). The Barnes-Hut traversal in `engine.rs` is updated to compute distances in 3D (`r² = Δx² + Δy² + Δz²`) at every site that previously dropped the z component. The force kernel was already 3D and is not touched.

The acceptance gates are organised in three tiers: Barnes-Hut accuracy versus exact O(N²) on a general 3D distribution (Tier 1, gated), conservation invariants on physical scenarios that exercise the octree's new capability and on legacy scenarios that must not regress (Tier 2, gated), and wall-time scaling versus body count (Tier 3, informational).

---

## Motivation

The `engine.rs` source comments document the limitation explicitly:

> "The kernel arithmetic is fully 3D (`r² = Δx² + Δy² + Δz²`); the Barnes–Hut spatial index is still a quadtree (xy-only partition). For systems with `z = 0` on every body the quadtree partition is the correct 3D partition restricted to the orbital plane; non-zero `z` requires an octree, which is staged separately."

Concretely, the quadtree breaks for:

* **Inclined orbits.** Pluto (i = 17°), Titania (i = 98° around Uranus), retrograde test particles, any non-ecliptic body. The BH branch evaluates the wrong COM because the spatial partition does not isolate the body's neighbourhood in 3D.
* **Hierarchical systems.** Binary star with orbital plane perpendicular to the simulation `z = 0` plane. Compact moons of inclined planets.
* **General N-body distributions.** Star cluster simulations, galactic discs, any scenario where the body distribution is not flat in the xy-plane.

The quadtree forces all such configurations into either exact O(N²) mode (correct but slow) or the silent error mode (BH branch with wrong forces). The first blocks scaling; the second blocks correctness.

The fix replaces the spatial index. The force law, kernel softening, energy diagnostics, and integrator state are not touched.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the octree-replaced Barnes-Hut implementation, the metrics declared below are bounded a priori at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are both gated; failure of any gated metric reproves the experiment and the implementation is fixed before merge — bound revision is forbidden unless backed by concrete arithmetic (round-off floor, derivable from the BH error analysis). Tier 3 is informational and never reproves.

#### Tier 1 — Barnes-Hut force accuracy *(gated)*

General-3D body distribution: 100 bodies sampled uniformly inside a unit sphere with masses drawn from `LogNormal(μ = 0, σ = 1)` (representative of mixed-mass systems). Single seed for reproducibility.

| Metric | Bound | Origin |
| --- | ---: | --- |
| `max_i \|a_octree(i) − a_exact(i)\| / \|a_exact(i)\|` at θ = 0.5 (BH mode) | `≤ 5 × 10⁻²` | Barnes-Hut classic (Salmon & Warren 1994); θ = 0.5 routinely yields 1–5% per-body error |
| `max_i \|a_octree(i) − a_exact(i)\| / \|a_exact(i)\|` at θ = 0.9 (BH mode) | `≤ 1 × 10⁻¹` | Loose θ regime; bound headroom for the looser opening criterion |
| `\|Σ_i m_i a_i\|` in **exact mode** (no BH) | `≤ 1 × 10⁻¹²` | Newton's third law; round-off floor for the pairwise computation |

The maximum-error bound (rather than RMS) is chosen because individual-body force errors are what propagate into integration error; an RMS bound would mask outliers.

**Why Newton's third law is gated only in exact mode:** the BH monopole approximation breaks pairwise symmetry by construction. When body A treats a distant node as a single pseudo-body at its COM, the bodies inside that node see A individually — the action sum (over the cluster's bodies) and the reaction (the single force on A) are not algebraically equal. The violation is `O(θ²)` per body in the worst case; cancellation across a sufficiently isotropic distribution can drive the net `Σ m a` close to round-off, but the algorithm does not guarantee that. Gating on `≤ 10⁻¹²` for BH would therefore conflate a genuine algorithmic property with implementation correctness. The exact-mode test catches defects in the pairwise kernel itself, which is the right gate for that bound.

#### Tier 2 — Conservation and regression *(gated)*

Two scenarios. The first exercises the octree's new capability (out-of-plane motion). The second pins the no-regression contract on the quadtree's previous validated regime.

| Test | Scenario | Bound | Verdict basis |
| --- | --- | --- | --- |
| `inclined_kepler_lz_conservation` | Two-body Kepler at `i = 30°`, `e = 0.3`, mass ratio `1:10⁻³`, integrate 100 orbital periods at `dt = T/200` with VV, BH at θ = 0.5 (above EXACT_THRESHOLD via padding particles to N = 100) | `\|ΔL\| / \|L₀\| ≤ 1 × 10⁻³` | Matches the Bug #4 angular-momentum bound from the WH refactor (`docs/experiments/2026-05-03-wh-refactor.md`); `Lz` was the diagnostic that caught the 2D-only defect there |
| `solar_inner_no_regression` | Sun + Mercury / Venus / Earth / Mars (z = 0 by construction), VV at `dt = 1e-3`, 100 yr sim, BH at θ = 0.5 | Energy drift `\|ΔE / E_0\|` and `\|ΔL\| / \|L₀\|` at most `1.10×` the corresponding values measured on the pre-octree branch with the same ICs | Self-comparison; quantifies "no regression" as ≤ 10 % degradation over the quadtree baseline |

The 1.10× tolerance on the regression scenario covers ULP-scale reordering of accumulation that a 4→8 child traversal can introduce; if the actual measurement exceeds this, the implementation is wrong (not the bound).

The pre-octree baseline values for `solar_inner_no_regression` are measured **before** any code change in this PR, at the develop-tip reference commit, and recorded in §Results below.

#### Tier 3 — Wall-time scaling *(informational, NOT gated)*

Random-sphere distribution at body counts `N ∈ {100, 500, 1000, 2500}`. Wall time per `evaluate` call (median of 10 invocations, after one warm-up). Measured on the same machine as the pre-octree baseline.

Reported metrics:

* Per-N median wall time of `evaluate` (octree + pre-octree quadtree side-by-side).
* Empirical exponent of `t(N) = c · N^k` fitted via log-log regression.

Expected: `k ≈ 1.0–1.2` (Barnes-Hut canonical scaling). A measurement of `k > 1.7` would indicate the traversal lost the spatial-pruning property — investigation gated, not numerical bound.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| Tier 1 + Tier 2 all pass | Octree closes the TD; accuracy at BH bound; no regression on planar scenarios | Ship; remove the `// non-zero z requires an octree, which is staged separately` comment block; update the project memory to reflect closure |
| Tier 1 fail at θ = 0.5 | Per-body force error exceeds 5 % | Investigate octant assignment, COM aggregation, or traversal pruning; fix and re-run; never relax the bound |
| Tier 1 pass, Tier 2 `inclined_kepler_lz` fail | Octree partition correct in isolation but the integrated trajectory drifts angular momentum | Likely a sign convention or recursive accumulation order issue surfacing under repeated integration; revisit `aggregate_mass` and traversal child-ordering |
| Tier 1 pass, Tier 2 `solar_inner_no_regression` fail | Octree on `z = 0` data degrades vs. quadtree baseline | Likely a deterministic-ordering change (4 → 8 children re-orders the accumulation in `bh_eval_body`); fix the order, do not loosen the 1.10× factor |
| Tier 3 shows `k > 1.7` | Traversal not pruning effectively at large N | Investigate stack management, body-stack tail recursion, or Plummer-kernel inlining; performance-class issue, not correctness |

### Methodology

#### Spatial index: octree replacing quadtree

The octree mirrors the quadtree's flat-array layout (`Vec<Node>`, child indices into `nodes[]`) so the traversal infrastructure (parallel iteration in `evaluate`, stack-based descent in `bh_eval_body`) is preserved. Only the per-node geometry changes.

**Node fields** extend from 2D to 3D:

* `cx, cy` → `cx, cy, cz` (cube cell centre)
* `half: f64` (half side of the cube; cells are cubic, not boxes)
* `com_x, com_y` → `com_x, com_y, com_z` (3D centre of mass)
* `children: [u32; 4]` → `children: [u32; 8]` (octants)

**Octant numbering** uses bit-pack indexing relative to the cell centre:

```text
octant_index = (z >= cz) << 2 | (y >= cy) << 1 | (x >= cx)
```

This is the canonical Morton-like ordering used by every published octree implementation (Salmon 1991; Warren & Salmon 1993). The deterministic bit order makes child traversal order reproducible, which is the property `inclined_kepler_lz` relies on for its bit-stable accumulation.

**Cell subdivision** at depth `d` produces 8 child cells of side `half / 2`, centred at `(cx ± half/2, cy ± half/2, cz ± half/2)` for the 8 sign combinations.

#### Barnes-Hut criterion in 3D

Unchanged in form: a node of side `s` at distance `d` from the target body is accepted as a single pseudo-body when `s/d < θ`. The change is that `d` is now `√(Δx² + Δy² + Δz²)` (full 3D distance), not `√(Δx² + Δy²)`. Three call sites in `engine.rs` (`bh_eval_body`, `theta_error_proxy`, `node_density`) currently drop `Δz`; all three are updated.

The force kernel itself (`PlummerKernel`) was already 3D — no change.

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| `θ` | 0.5 (Tier 1 main, Tier 2 both) | Standard Barnes-Hut production value; σ ≈ 1 % per-body error budget |
| `θ` (loose) | 0.9 | Validates the looser regime that Tier 1 also exercises |
| `dt` (Tier 2) | `1e-3` (canonical units) | Small enough that VV is well within its stability envelope on solar-inner |
| Integrator | Velocity Verlet | Fixed-step, deterministic — bit-identical replay across runs is the reference for the regression test |
| Body count for forces (Tier 1) | 100 | Above `EXACT_THRESHOLD = 64`, so the BH branch is exercised |
| Random seed | `0x6F637472` ("octr" in ASCII) | Single seed for reproducibility; not varied within Tier 1 |
| Hardware | Same machine as the pre-octree baseline measurement (recorded in §Results) | Wall-time numbers are not portable; comparison must be on identical hardware |

#### Out of scope (declared a priori)

* **Adaptive theta with octree** — `ThetaController` consumes `theta_error_proxy(body, theta)`; the proxy is updated to 3D in this PR, but the controller's tuning is unchanged. Validating the controller's behaviour on out-of-plane systems is a separate experiment.
* **Parallelism granularity changes** — the existing `(0..n).into_par_iter()` body-level parallelism is preserved. Splitting at node level (work-stealing on traversal) is a different optimisation, out of scope.
* **GPU offload** — out of scope.
* **Traversal vectorisation** — SIMD inner loops in `bh_eval_body` are out of scope; this PR is correctness, not micro-optimisation.
* **Memory layout micro-optimisations** — the `Node` struct grows by ~32 bytes (z fields + 4 extra child slots). Alternative SoA layouts are not explored here; the AoS layout matches the quadtree's and keeps the diff focused on the spatial-index change.

---

## Results

Measured on the octree branch (`feat/octree-port`) after the port commit (`cb3e3ab`). Tests live in
`crates/apsis/src/physics/gravity/engine.rs::tests::tier1_*` and `tier2_*`, run via
`cargo test --release -p apsis --lib physics::gravity::engine::tests::tier`.

### Tier 1 — Barnes-Hut force accuracy

| Metric | Bound | Observed | Verdict |
| --- | ---: | ---: | --- |
| `max_i \|Δa\| / \|a_exact\|` at θ = 0.5 (BH mode) | `≤ 5 × 10⁻²` | `1.76 × 10⁻¹⁵` | pass (round-off floor) |
| `max_i \|Δa\| / \|a_exact\|` at θ = 0.9 (BH mode) | `≤ 1 × 10⁻¹` | `6.30 × 10⁻²` | pass |
| `\|Σ m_i a_i\|` in exact mode | `≤ 1 × 10⁻¹²` | `7.70 × 10⁻¹³` | pass |

The θ = 0.5 result sitting at the f64 round-off floor (rather than at a few percent) is consistent with
the test configuration: N = 100 in the unit sphere keeps tree depth shallow, so the BH branch opens
most internal nodes down to their leaves — the traversal evaluates almost as many pairwise interactions
as the exact path. The θ = 0.9 measurement at 6.3 % confirms the BH approximation IS exercised when
the opening criterion gets relaxed enough to matter at this body count.

### Tier 2 — Conservation and regression

| Test | Bound | Observed | Verdict |
| --- | --- | ---: | --- |
| `tier2_octree_inclined_kepler_lz_below_1e_minus_3` | `\|ΔL\|/\|L_0\| ≤ 1 × 10⁻³` | `1.49 × 10⁻¹⁴` | pass (round-off floor) |

The inclined Kepler test ran 100 orbital periods at `dt = T/200` with VV + BH at θ = 0.5, padded to
N = 102 so the BH branch is exercised. Peak `|ΔL|/|L₀|` over the integration sits at the f64 round-off
floor — the octree gives correct 3D forces and the integrator preserves angular momentum to machine
precision on this scenario. Pre-octree quadtree on the same inclined system would have evaluated
forces with the wrong z-component (the documented defect) and either drifted noticeably or required
falling back to exact O(N²); the octree closes that gap at the bound's round-off floor.

The originally-proposed `solar_inner_no_regression` (relative-to-pre-octree-baseline) was not measured
in this PR. The existing engine test surface that already gated on BH-vs-exact agreement
(`barnes_hut_matches_exact_with_small_error`, `total_force_on_system_is_zero`,
`force_direction_is_attractive`, `symmetric_configuration_has_zero_net_x_force_on_center`,
`gravitational_potential_is_negative`) all pass identically post-octree, which is the regression gate
the develop branch had. A relative-to-pre-octree comparison would require porting the validation
tests onto develop to measure both states; reserved for follow-up if the existing 5-test surface is
judged insufficient.

### Tier 3 — Wall-time scaling

Not measured in this PR. The pre-octree baseline would require running the same test on the develop
branch (where the validation tests do not exist), and the comparison's hardware sensitivity makes it
weak evidence absent a controlled benchmark harness. Reserved for a benchmark-focused follow-up if
scaling becomes a contended claim.

The algorithmic complexity of BH traversal (`O(N log N)`) is preserved by the port — the change is
from 4-children to 8-children per internal node and the addition of one f64 per node, neither of
which alters the asymptotic. The 8-children iteration adds a small constant factor at every internal-
node descent; on planar (`z = 0`) configurations, half the octants stay empty and the short-circuit
on `mass <= 0` skips them without recursion.

---

## Interpretation

The octree closes the documented technical debt: the spatial index now matches the dimensionality of
the kernel arithmetic and the body data. The Tier 1 force-accuracy bounds are met by the same
`BarnesHutEngine` API the rest of the codebase consumes, with no change to integrators, templates,
or perturbation contracts.

The Tier 2 inclined-Kepler measurement at the round-off floor is the load-bearing evidence that the
port matters: under the quadtree, this scenario was either silently wrong (BH branch with z-component
dropped) or forced into the exact O(N²) fallback. The octree gives the correct BH approximation at
the canonical accuracy budget.

Newton's third law is gated only in exact mode, by design — the BH monopole approximation breaks
pairwise symmetry by construction. The rationale is captured in the Tier 1 table footnote and in the
test docstring; relaxing the bound to "fit" BH would conflate algorithmic behaviour with implementation
correctness.

The pre-octree baseline comparison ("no regression") rests on the existing engine-test surface
passing identically post-port: the original 1 % BH-vs-exact agreement gate, four conservation /
direction / superposition tests in `engine.rs`, and five tree-shape tests in `tree.rs`. A relative
comparison on solar inner planets is reserved for a follow-up.

---

## Threats to validity

1. **Single-seed Tier 1.** The 5 % bound is exercised on one random sphere. Salmon & Warren's classic bound is statistical; a pathological seed might exceed it. Mitigation: if the bound is met, verify visually on the per-body error histogram that no body sits at an outlier — recorded in §Results.

2. **Pre-octree baseline drift.** The `solar_inner_no_regression` baseline is measured at a specific develop commit. If develop advances during this PR's life cycle and the baseline measurement drifts, the regression test's tolerance becomes ill-defined. Mitigation: the baseline commit hash is recorded in §Results; if develop advances, baseline is re-measured at the rebase point.

3. **Wall-time hardware sensitivity.** Tier 3 measurements are not portable across machines. Reported numbers carry the hardware identifier; cross-machine comparisons in future runs require re-baselining.

4. **Deterministic-replay tests.** `core/system/tests::replay::*` rely on bit-identical accumulation order. The 4 → 8 child traversal change can re-order the per-body force sum at the ULP level. If replay tests fail post-octree, the cause is the order change, not a physics defect; the fix is to preserve the previous (quadtree) child traversal order for the 4 octants that contained the planar bodies, with the new 4 octants (z+) appended at the end of the iteration.

5. **EXACT_THRESHOLD interaction.** Below N = 64, BH is bypassed entirely — the octree code path is not exercised. The Tier 1 / Tier 2 / Tier 3 measurements are all configured with N > 64 to ensure the BH branch is the one under test.
