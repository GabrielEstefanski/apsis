---
date: 2026-05-01
status: current
---

# `recommended_dt` heuristic — validation note

`System::recommended_dt()` is the utility that returns a physics-justified timestep for fixed-step integrators (Velocity Verlet, Yoshida-4, Wisdom–Holman). This note records the formula, the validation harness, and the current verdict.

## Heuristic

```text
dt = min(dt_dynamic, dt_pair, dt_softening)
     clamped to [1e-9, 1e6]

  dt_dynamic   = 0.05 · min(√(r_min / a_max), a_max / |jerk|)
  dt_pair      = 0.01 · min_ij 2π · √(r_ij³ / μ_ij)
  dt_softening = 0.05 · √(ε / a_max)    when kernel ε > 0, else +∞
```

`r_min` is the closest-pair separation, `μ_ij = G · (m_i + m_j)`. References: Aarseth (2003) §2 for the dynamical timescale family; Power et al. (2003) for the softening-based criterion.

## Validation harness

Lives at [`validation/recommended-dt/`](../../validation/recommended-dt/). Two Cargo examples:

- `recommended_dt_validation` — drives 13 templates × 3 fixed-step integrators at `recommended_dt()` for 100 substeps; writes per-step trajectory to `out/runs.csv`.
- `recommended_dt_compare` — loads `runs.csv`, computes peak `|ΔE/E_0|` and `|ΔLz|` per cell, applies the gates, emits `out/comparison.json`, exits non-zero if any gated cell fails.

### Gates

- VV: `|ΔE/E_0| ≤ 1 × 10⁻³`
- Y4: `|ΔE/E_0| ≤ 1 × 10⁻⁶`
- VV + Y4: `|ΔLz| ≤ max(10⁻¹⁰ · |Lz_0|, 10⁻¹⁰)` (isclose-style two-sided bound)
- WH: informational; its symplectic structure depends on dt–period commensurability, which the heuristic does not encode.

### Out-of-regime scenarios (informational)

Four templates fall outside the heuristic's quasi-regular envelope and report informationally rather than gated:

- `alpha_centauri_ab` — `e = 0.52` binary (Kervella et al. 2017)
- `hd_80606_b_system` — `e = 0.93` eccentric hot Jupiter (Naef et al. 2001)
- `three_body_figure_eight` — periodic 3-body with tight pair passes
- `three_body_pythagorean` — chaotic; template suggests Mercurius

## Current verdict

`18/18 gated cells pass + 21 informational cells`. Tightest gated cell: `hot_jupiter` Y4 at `8.7 × 10⁻³` of bound; all gated cells within an order of magnitude of their bounds, none at the edge.

## Contract

`recommended_dt` targets quasi-regular orbital regimes — hierarchical, low-to-moderate eccentricity, no strong close encounters. Outside that envelope the function returns a finite value but neither the dynamical nor the pair-Kepler term guarantees that fixed-step VV / Y4 hold the conservation bounds above. The harness records `gated: false` for the four scenarios outside the envelope rather than tuning the formula to a regime it was not derived for.

## History

Prior phases of the experiment live in the git log:

- Phase A — original Tier 2 bound `|ΔLz/Lz_0| ≤ 10⁻¹⁰` (sub-round-off in scenarios with small `|Lz_0|`); commit `0d71547`.
- Phase B — isclose-style two-sided revision; commit `b30d278`. Verdict `26/26 gated cells pass` against the heuristic available at that commit (pre-PR-#121).
- Phase C — kernel-softening removal (PR #121) silently swapped the active branch; per-template IC audit surfaced the `binary` `v_body` √2 mismatch; heuristic refactored to `min(dt_dynamic, dt_pair, dt_softening)`; four scenarios reclassified as out-of-regime. Commits `3f7b809`, `86e5802`, `f064346`.
