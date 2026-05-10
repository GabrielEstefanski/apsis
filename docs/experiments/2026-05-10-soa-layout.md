# Struct-of-Arrays body layout — protocol

**Date:** 2026-05-10
**Subject:** Introduce a per-step immutable Struct-of-Arrays snapshot of the body state for the gravity hot path, leaving the integrator and the public API on the existing `Vec<Body>` Array-of-Structs. Measure the cache-locality wall-time gain on the BH walk at the v1 target N range.

**Status:** Protocol declared a priori, before the type lands. §Results populated incrementally; §Decision written after Tier 1 + Tier 2 are populated.

**Branch:** `perf/soa-layout`, from `origin/develop` (parallel to `perf/octree-mac` / PR #74; no shared hot files).

---

## Abstract

The engine ceiling experiment (`2026-05-09-engine-ceiling.md` §Decision) classified the engine as interaction-bound at the v1 target regime. The MAC comparison experiment (`2026-05-09-octree-mac.md` §Decision) defered the MAC axis and shifted SoA to the lead, with SIMD as the immediate follow-up.

The current `Body` struct is approximately 104 bytes. The BH walk's hot path reads only five fields per leaf-pair interaction (`pos_x, pos_y, pos_z, mass, softening` = 40 bytes per body), so each `Body` load brings ~88 bytes of cold data through the cache (Body straddles two 64-byte cache lines). This is structural waste that no algorithmic change can address — the fix is to lay the five hot fields out contiguously so the cache line carries useful payload only.

This experiment introduces `BodyArrays`, a five-field SoA snapshot (`pos_x, pos_y, pos_z, mass, softening`), packed once per step from `Vec<Body>` immediately before the force evaluation, consumed by the tree build and the BH walk in indexed `for i in 0..n` form, and discarded after the force eval returns. The integrator continues to operate on `Vec<Body>` (it is compute-bound in its own internal arrays, not Body-bound). The public API surface — templates, save format, render layer, perturbation API, inspector — continues to read `Body`.

The acceptance gates are organised in three tiers: bit-exact accelerations (Tier 1, hard gate), BH walk wall-time speedup at fixed θ (Tier 2, gated against arithmetic-derived ranges), pack overhead per step (Tier 3, gated at ≤ 1 % of total step wall time).

---

## Motivation

The MAC §Decision baked an updated four-axis roadmap with SoA at axis 2. The reasoning chain for SoA leading:

1. The current `Body` layout wastes ~70 % of every cache line load in the gravity hot path (40 useful bytes of 128 loaded).
2. SIMD optimisation (axis 3) is only worthwhile when the loads it vectorises are dense — vectorising scattered AoS loads wastes most of the SIMD lane bandwidth.
3. SoA captures a measurable gain *before* SIMD lands and is the structural pre-requisite for SIMD to pay its full speedup.
4. REBOUND (the closest-regime production codebase to apsis — small-N, error-control rigorous, IAS15-class; cf. `project_codebase_neighborhood`) keeps `reb_particle` AoS in its public API but has SIMD-friendly internal buffers; the layered approach taken here mirrors that design.

The arithmetic above gives a theoretical upper bound of `128 / 40 ≈ 3.2×` for the leaf-pair phase of the walk. Internal-node loads (`Node.com_x/y/z, mass, q_xx..q_yz` = 72 bytes of a ~144-byte node) are not refactored by this PR and improve only when SIMD lands. Engine ceiling §Results reports `n_leaf_interactions ≈ 0.6 × n_bh_accepted` at N = 10⁴ with quadrupole, which puts the leaf phase at roughly 35–50 % of total walk time depending on per-interaction cost ratio. The blended walk speedup is bounded a priori in §Hypothesis below.

### Design constraint: SoA as a write-once snapshot, not a domain type

The SoA buffer is the *execution* state of one phase (force eval); it is not a domain object. To keep ownership unambiguous and rule out an entire class of staleness bugs by construction:

- **Lifecycle.** `pack_from(&[Body])` at the start of each step → tree build reads SoA → BH walk reads SoA → force eval returns. Conceptually discarded; in practice the buffer is reused by the next step's `pack_from`. **Never mutated during the step.**
- **Sole consumer is the gravity hot path.** Integrator reads/writes `Vec<Body>`. API consumers read `Vec<Body>`. SoA is never observed outside the force eval window.
- **No invalidation flags, no dirty bits.** Because SoA is rebuilt every step from authoritative `Vec<Body>`, there is no concurrent writer to track. The pack itself is the synchronisation point.
- **Integration stays AoS.** Not because AoS is intrinsically better for integrators, but because IAS15 is compute-bound in its `b/g/e/csb` coefficient arrays — the Body load is not the dominant cost. If a future integrator profile shows Body-bound integrator hot paths, the analysis re-opens; the rule today is empirical, not aspirational.
- **Indexed kernel, no iterator abstraction.** The walk and tree build read `arrays.pos_x[i]` directly in `for i in 0..n` loops. No `Iterator<Item = Body>` wrapper, no trait-object. This keeps the inner loop in a shape that PR-perf-6 (SIMD) can unroll and vectorise without structural refactor.

Operational rule per step:

```text
pack AoS → SoA           (System.body_arrays.pack_from(&self.bodies))
build tree using SoA     (Octree::build(&self.body_arrays))
BH walk using SoA        (engine.evaluate(&self.body_arrays, ...))
integrate in AoS         (integrator.step(&mut self.bodies, &acc, dt))
SoA stale until next pack
```

Simple. Deterministic. No shared mutable state across the step boundary.

---

## Protocol *(declared a priori, before any code lands)*

### Hypothesis

The metrics declared below are bounded a priori at the values stated. Bound revision is forbidden unless backed by concrete arithmetic (cf. the `feedback_no_tuning_to_pass` rule applied throughout the perf series).

#### Tier 1 — Accelerations bit-exact *(hard gate)*

For the same body distribution, the same θ, and the same kernel, the accelerations returned by `BarnesHutEngine::evaluate` reading `BodyArrays` must equal those returned by the current `&[Body]` path **bit-for-bit per component**. The walk does no cross-body reduction in the per-body acceleration sum (each body's acceleration is computed independently in the parallel iter), so floating-point reordering cannot be invoked as an excuse for any divergence.

| Quantity | Bound |
| --- | --- |
| `‖a_soa[i] − a_aos[i]‖_∞` per body | `== 0.0` (bit-exact) |
| `Σ_i m_i · a_soa[i] − Σ_i m_i · a_aos[i]` (system net force) | `< 1 × 10⁻¹²` |

Net-force tolerance is non-zero only because the `phi` reduction in `evaluate_profile` sums in non-deterministic rayon order and the SoA path may interleave differently — this is the *only* legitimate source of divergence; per-body accelerations must remain bit-exact regardless of scheduling.

**Failure here is non-negotiable**: any per-body acceleration divergence indicates an indexing bug in the SoA path and the experiment halts until the bug is found.

#### Tier 2 — BH walk wall-time speedup *(gated as range)*

Same N grid and seeds as `2026-05-09-octree-mac.md` (`N ∈ {1_000, 5_000, 10_000}`, three seeds `0x6F637472`, `0x71756164`, `0x6D6F7274`), same θ = 0.5, same sphere log-normal distribution. Wall-time is the within-seed median over 5 measured runs (3 warmup); `t_walk` captured in isolation from `t_build` and `t_pack` via `evaluate_profile` instrumentation already in place.

Predicted ranges derived from the leaf-phase fraction of total walk time:

| Comparison | Range | Derivation |
| --- | --- | --- |
| `t_walk_AoS / t_walk_SoA` at N = 10⁴, leaf phase = 35 % of walk | ∈ [1.18, 1.32] | `1 / (0.65 + 0.35 / 3.2) = 1.32×` upper, `1 / (0.65 + 0.35 / 2.0) = 1.18×` lower |
| `t_walk_AoS / t_walk_SoA` at N = 10⁴, leaf phase = 50 % of walk | ∈ [1.33, 1.52] | `1 / (0.50 + 0.50 / 3.2) = 1.52×` upper, `1 / (0.50 + 0.50 / 2.0) = 1.33×` lower |
| **Combined Tier 2 range** | **∈ [1.20, 1.50]** | union of the two scenarios; covers leaf-phase uncertainty |

Lower bound `2.0×` for SoA leaf-phase gain (vs the theoretical `3.2×`) accounts for the prefetcher already amortising part of the AoS waste under sequential leaf iteration. Upper bound assumes the realistic inner loop hits the L1-bound regime the arithmetic predicts.

A measurement inside the range ships SoA. Below the range triggers root-cause investigation (likely: prefetcher already covers more than estimated, or another bottleneck dominates) before deciding. Above the range is reported with the discrepancy and shipped — over-performance is welcome but flagged for understanding.

#### Tier 3 — Pack overhead per step *(gated)*

`pack_from(&[Body]) → BodyArrays` is a linear-in-N memory copy: `5N × 8 bytes` written sequentially across five `Vec<f64>`. At ~10 GB/s memory bandwidth, expected cost is ~4 µs per 1 000 bodies, ~40 µs per 10 000 bodies.

| Bound | Threshold | Rationale |
| --- | --- | --- |
| `t_pack / t_step` at N = 10⁴ | `≤ 0.01` (1 %) | Pack overhead must not erase the walk gain. With expected `t_step ≈ 4 ms` post-SoA at N = 10⁴, 1 % budget is ~40 µs — exactly the bandwidth-limited estimate. |

Pack frequency is **once per step** by construction (rebuilt every step from `Vec<Body>`; no caching). The harness logs `t_pack` per N + per seed so future regressions are caught — if `t_pack` ever crosses the 1 % threshold in a later experiment, that PR has to root-cause before merge.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| Tier 1 fails | Indexing bug in SoA path | **Halt experiment, fix bug, restart Tier 1** |
| Tier 1 passes AND Tier 2 inside [1.20, 1.50] AND Tier 3 ≤ 1 % | SoA captures the predicted cache-locality gain at acceptable overhead | **Ship SoA as production**; the refactor itself is the bake — no toggle to remove |
| Tier 1 passes AND Tier 2 below 1.20 | Prefetcher covers more than estimated, or another bottleneck dominates | Root-cause via cache-miss profiling before shipping; if root cause is benign (e.g. compiler already reordered loads), still ship — Tier 1 is the gate |
| Tier 1 passes AND Tier 2 above 1.50 | Either leaf phase larger than estimated OR a second-order win (fewer branch mispredicts, etc.) | Ship and document; flag for follow-up understanding |
| Tier 3 above 1 % | Pack is slower than memory-bandwidth-bound estimate | Investigate before shipping (likely: copy is going through a non-contiguous path or hitting allocator) |

### Methodology

#### Implementation order — Tier 1 gate first, then measurement

1. **Notebook a priori** (this commit).
2. **Rename `Body` fields** to `pos_x/y/z, vel_x/y/z` for nomenclatural alignment with `BodyArrays`. Mechanical refactor across ~36 files; behaviour-identical.
3. **Introduce `BodyArrays`** in `crate::domain::body_arrays`. Type + `pack_from(&[Body])` + roundtrip / pack-correctness tests. No consumers. **No `unpack_into`** — the snapshot is write-once.
4. **`Octree::build` and `BarnesHutEngine::evaluate` consume `&BodyArrays`** with indexed `for i in 0..n` inner loops. **Tier 1 gate runs here** — bit-exact accelerations vs the AoS path on a fixed seed before the step orchestration changes.
5. **`System::step` packs the snapshot once per step** before calling the force eval. Buffer reused across steps (allocated once at `System::new`, `clear()` + extend on each `pack_from`); integrator continues to read/write `Vec<Body>` unchanged.
6. **`perf_soa` harness + run + populate §Results**.
7. **§Decision** — ship if gates clear; halt-and-fix if Tier 1 fails.

#### Run parameters

| Parameter | Value | Justification |
| --- | --- | --- |
| Random seeds | 3: `0x6F637472`, `0x71756164`, `0x6D6F7274` | Match perf 2×2 / engine ceiling / MAC for cross-experiment comparability |
| Body distribution | sphere log-normal mass | Match perf 2×2 family |
| `θ` | 0.5 | Production canonical |
| N | `1 000`, `5 000`, `10 000` | Match MAC experiment |
| Warmup / measured runs | 3 / 5 (per cell) | Match MAC harness |
| Multipole order | Quadrupole always-on | Per perf 2×2 §Decision |
| MAC | Classical `s/d < θ` | Per MAC §Decision |
| Hardware | Same as prior perf series (Ryzen 5 7600X, Windows 11) | Cross-experiment comparability |

#### Out of scope (declared a priori)

- **Velocity in the SoA snapshot.** BH walk does not read velocity; including `vel_x/y/z` would add 24 bytes per body to the SoA payload (40 → 64) and reduce the per-cache-line density that motivates the refactor. Integrator stays on AoS for velocity reads.
- **AoSoA / chunked layout (4 or 8 bodies per chunk)**. Pure SoA first; AoSoA is a SIMD-era optimisation deferred to PR-perf-6 if the SIMD profile asks for it.
- **Single-allocation `Vec<f64>` with offset slices**. Separate `Vec<f64>` per field is simpler and sufficient for cache locality. Migration to single-allocation aligned buffers is a small follow-up if SIMD aligned-loads need it.
- **AVX-strict alignment (32 / 64 byte).** `Vec<f64>` provides natural 8-byte alignment. AVX2 / AVX-512 aligned loads are a PR-perf-6 concern.
- **Universal SoA (deleting `Body`)**. Layered approach keeps `Body` as the API surface; perturbations, templates, save format, render, and inspector continue to read `Body`. SoA is execution state, not a domain type.
- **Refactoring `Octree::Node` layout**. Internal nodes stay AoS in this PR. SIMD-era refactor candidate.
- **Integrator SoA refactor**. IAS15 is compute-bound in its internal coefficient arrays; Body load is not the dominant cost. Re-evaluate if a future integrator profile shows otherwise.
- **MAC re-evaluation**. Per the MAC §Decision, MAC re-enters scope only after SoA + SIMD land and the re-measure step (axis 4) reclassifies the engine.
- **Cross-machine comparison**. Single-hardware as in prior experiments.

---

## Results

*To be populated incrementally as commits land.*

### Tier 1 — Accelerations bit-exact

*Pending (commit 4 gate).*

### Tier 2 — BH walk wall-time speedup

*Pending (commit 6).*

### Tier 3 — Pack overhead per step

*Pending (commit 6).*

---

## Interpretation

*To be written after Tier 2 is populated.*

---

## Decision

*To be written after Tier 1 + Tier 2 + Tier 3 gates pass or fail.*

---

## Threats to validity

1. **Single distribution family.** Sphere log-normal matches the perf series canonical and provides cross-experiment comparability; a clustered distribution would have heavier leaf phases (deeper trees, more pairs per leaf), pushing both AoS and SoA walk costs higher and potentially shifting the leaf-phase fraction. Mitigation: report the seed and distribution explicitly; flag any cell whose ranking depends on this choice.

2. **Bit-exactness assumption.** The walk is per-body and computes accelerations independently; no cross-body reduction touches the per-body result. If a future kernel introduces a reduction in the per-body path (shouldn't happen, but worth noting), the bit-exact gate would need to relax to `< 1 × 10⁻¹⁵` per component. The current Plummer kernel and the truncated-Plummer counter-test kernel both have purely per-body computations.

3. **Cache effects are machine-dependent.** Ryzen 5 7600X has 32 KB L1d per core, 1 MB L2 per core, 32 MB L3 shared. The leaf phase at N = 10⁴ with `LEAF = 8` involves ~40 × 8 = 320 bytes per leaf pair-set — fits in L1 with margin even before SoA. The cache-locality gain is therefore L2/L3-driven (whole-distribution loads), not L1-driven. Different cache hierarchies (Apple M-series, server EPYCs, older Intel) may show different gain magnitudes. Single-machine measurement is honest about this.

4. **Compiler may already do part of the optimisation.** rustc + LLVM can sometimes hoist field loads and reorder to improve locality even with AoS — this is the lower bound of the Tier 2 range. If Tier 2 lands at the lower edge or below, that is the most likely explanation; the gain remains real because the compiler cannot reorder *across* heap allocations.

5. **Layered approach leaves AoS waste in non-hot paths.** Templates, save format, render, perturbations, and inspector still read `Body`. None of these are hot at N ≤ 10⁵ (per the engine ceiling Cell V profile, render is < 5 % of step time and templates run once per session). If a future perf experiment shows one of these hitting cache pressure, Universal SoA re-enters scope — but the threshold for that is empirical, not aspirational.

6. **Pack overhead measurement assumes warm cache.** Pack frequency is once per step, but the harness measures pack alongside the step and the cache may already be warm from the previous step's walk. Cold-cache pack would be slower. Mitigation: report `t_pack` per cell and flag if it shows N-superlinear scaling (which would indicate allocator or layout pathology).
