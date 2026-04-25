# Operational Domain — Speed, Conservation, and Pareto Benchmarks

**Date:** 2026-04-24
**Subject:** Quantifying the `(integrator, scenario, N) → regime` surface for `apsis`
**Baseline commit:** `354f82f` (modular `scaling_benchmark` example landed)
**Tooling:** `examples/scaling_benchmark.rs`, `examples/conservation_scan.rs`,
planned `examples/integrator_pareto.rs` and `examples/template_matrix.rs`
**Status:** in progress — this note accumulates findings across the four
benchmarks and becomes the source material for the paper's "Operational
domain" section.

---

## Abstract

The paper currently declares `N ≤ 10³` as deliberate scope and
describes the integrator zoo (VV, Y4, WH, IAS15) without quantitative
justification for the choice of four. Both statements underestimate
what the library can demonstrably do. This note establishes the
empirical surface that the paper's Operational-domain section will
draw from:

1. **Speed** — wall-time per step and simulation-time per wall-second,
   measured across four physically-distinct scenarios and three
   integrators (§1).
2. **Conservation** — energy and angular-momentum drift over
   scenario-specific dynamical windows, measured across the same
   matrix (§2).
3. **Pareto** — single scenario, `(integrator × dt)` scan revealing
   the cost/precision frontier. Plans: §3.
4. **Template coverage** — every `TemplateKind` in the library
   catalogue runs for N periods under a recommended integrator,
   reporting drift and wall-time. Plans: §4.

The paper is then reformulated with `(integrator, scenario, N)` as
the axes of the regime description, not `N` alone.

---

## 1. Speed scaling — findings

**Source:** `examples/scaling_benchmark.rs` at `354f82f`.
**Run on:** 12-core workstation, release build, warm-up 10 steps,
measured 50 steps per cell.

### Method

Four scenarios — `friendly_cluster`, `hierarchical_kepler`,
`clustered_substructure`, `multiple_binaries` — each declare a natural
`dt_hint(N)` derived from the scenario's shortest dynamical
timescale. VV, Y4, and IAS15 are swept across `N ∈ {128, …, 65536}`
(IAS15 capped at 4096 because its per-step cost makes larger N
exceed the 10-s-per-step abort threshold).

### Headline numbers

Interactive ceiling (≤ 33 ms/step) per (scenario, integrator):

| scenario                | VV ceiling | Y4 ceiling | IAS15 ceiling |
|-------------------------|-----------:|-----------:|--------------:|
| friendly_cluster        |    32 768  |    16 384  |           512 |
| hierarchical_kepler     |    16 384  |     8 192  |           512 |
| clustered_substructure  |    32 768  |    16 384  |           512 |
| multiple_binaries       |    32 768  |    16 384  |           512 |

Batch-realtime ceiling (≤ 1 s/step) is typically `2×` the interactive
ceiling for VV and Y4, limited by tree traversal cost scaling as
O(N log N).

### Interpretation

The paper's `N ≤ 10³` framing is off by at least one order of
magnitude for VV/Y4 on three of four scenarios. The correct claim
form is a table, not a scalar: the ceiling depends on
`(integrator, scenario)`. `hierarchical_kepler` is the tightest
scenario by ~2× because its `dt_hint` is set by the innermost Kepler
period (at `a = 0.3 AU`), which is smaller than the cluster
dynamical times at comparable N.

IAS15's ceiling is fundamentally lower — roughly `N = 512` for
interactive — driven by two structural factors, both present
regardless of scenario:

1. The integrator requires a deterministic force, which in `apsis`
   forces the automatic switch to direct O(N²) evaluation.
2. Each IAS15 step performs ~7 Picard force evaluations instead of
   the 1 (VV) or 6 (Y4) of the fixed-step schemes.

These multiply: direct O(N²) grows as `N²`, and each of those
evaluations costs 7× the VV single evaluation. At `N = 1024`, the
combined multiplier vs VV is roughly `N/log(N) × 7 ≈ 700×`. For the
dense-random-cluster scenarios we also see `dt`-floor events at
`N ≥ 1024` — the adaptive controller shrinks `dt` to its floor
(1e-12) trying to resolve close encounters that are frequent by
construction.

### Consequence for the paper

The Operational-domain section will carry:

- A scenario-indexed ceiling table for VV, Y4.
- A standalone note on IAS15's regime (small-N, high-precision: the
  Mercury/1PN test case is exactly the regime it excels in, since
  `N = 2` reduces both above multipliers to ~7×).

The paper's former `N ≤ 10³` claim is replaced by
`N ≤ [scenario-specific ceiling] for VV/Y4; smaller for IAS15, see
§Operational domain`.

---

## 2. Conservation — findings (pending)

**Source:** `examples/conservation_scan.rs`.
**Status:** bench run in flight as of this draft; numbers will
populate below once it completes.

### Method

For each cell, integrate `10 × t_characteristic` simulation units at
the scenario's natural `dt_hint` and record `|dE/E|` and `|ΔL_z|`
between the post-warm-up state and the final state. `t_characteristic`
is the scenario-declared dynamical time: disk crossing time for
cluster scenarios, inner-orbit period for `hierarchical_kepler`,
binary period for `multiple_binaries`.

### Expected numbers

Placeholder rows pending the bench run. What the table should reveal:

- VV: leading-order `dE/E ~ (dt)² · N_orbits`, growing with N
  because steps accumulate. Absolute numbers: 1e-4 to 1e-6 range
  expected.
- Y4: leading-order `dE/E ~ (dt)⁴ · N_orbits`, orders of magnitude
  smaller than VV at the same `dt`.
- IAS15: machine precision (1e-14 to 1e-12) expected when the
  adaptive controller is not in dt-floor distress.

### Interpretation (pending)

Once populated, the table answers "is my simulator accurate at N = X
for integrator Y over 10 orbits?" — the quality half of the
`N`-ceiling question.

---

## 3. Integrator Pareto — planned

**Target source:** `examples/integrator_pareto.rs`.

A single scenario (planned: `hierarchical_kepler` at N = 1024) with
each integrator swept across a range of `dt` values. Plot
`|dE/E|` vs wall-time on a log-log axis; the resulting Pareto
frontier demonstrates which integrator dominates at which cost
budget.

Expected result: VV wins at low precision (fast, lossy), Y4
dominates the middle, IAS15 is only cheap-per-precision at very
high accuracy targets. No single integrator dominates across the
frontier — which is why the library ships all four.

This supplants the paper's current qualitative "integrator zoo"
language with an explicit design-space map.

---

## 4. Template matrix — planned

**Target source:** `examples/template_matrix.rs`.

Run every `TemplateKind` in the library's built-in catalogue
(~19 entries) under a scenario-appropriate integrator for N
periods, and record drift and wall-time. Evidence that the
framework is general — the Operational-domain argument does not
rest on ad-hoc scenarios invented for the paper.

---

## 5. Open questions and next actions

- **Dense random clusters stress IAS15.** `dt`-floor events at
  N = 1024 on `friendly_cluster` and `clustered_substructure`.
  Expected (IAS15 is not designed for that regime), but worth a
  footnote in the paper rather than a silent anomaly.
- **Should WH be included in Pareto?** Wisdom–Holman is absent
  from the current benches because it is only applicable to
  hierarchical configurations with a dominant primary. The
  Pareto experiment planned in §3 is the natural place to bring
  it back in — `hierarchical_kepler` is the scenario that fits.
- **Softening asymmetry (N3L violation in BH).** Deferred to v0.2
  with its own fix. No impact on v0.1 Operational-domain claims;
  relative force error at our operational regime is bounded by
  the existing `ε/r` diagnostic in `docs/softening.md`.
