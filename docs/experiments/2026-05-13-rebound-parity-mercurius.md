# REBOUND parity — Mercurius

**Date:** 2026-05-13
**Subject:** Numerical agreement between apsis Mercurius (`crates/apsis/src/physics/integrator/mercurius.rs`) and REBOUND `MERCURIUS` (`integrator_mercurius.c`) on a Solar-System scenario with a Jupiter-crossing test particle. The paired Mercurius lab notebook (`docs/experiments/2026-05-13-mercurius-hybrid.md` §Tier 2) deferred this gate to a separate validation PR; this is that PR's notebook.

**Status:** Protocol declared *a priori*, before any code lands. Bounds and scenario locked from the Mercurius lab notebook §Tier 2 specification, refined to the practical horizon and dt that fit the existing `validation/rebound-parity/` infrastructure.

**Branch:** `validation/mercurius-parity`, branched from `develop` post-merge of PR #83. Independent of any other in-flight work.

---

## Abstract

Mercurius shipped in PR #83 with structural (Tier 1) and cost (Tier 4) validation but no cross-implementation parity. The federated FPM thesis treats integrators as first-class scientific artifacts that age well — the strongest substantiation for that claim is bit-level numerical agreement with the canonical REBOUND implementation on a scenario that actively engages the close-encounter sub-integration. This experiment is that gate.

The scenario is a Sun + 4 outer planets + 1 Jupiter-crossing test particle, integrated for 10⁴ years (~840 Jupiter orbits, ~1170 test-particle orbits) under Mercurius with REBOUND's default `r_crit_hill = 3` on both sides. Output is sampled at 1-year cadence; conservation diagnostics (ΔE / E₀, Δ Lz / Lz₀) and the test particle's osculating elements (a, e, i) are compared between implementations.

---

## Motivation

The Mercurius lab notebook §Decision flagged Tier 2 cross-implementation parity as the prerequisite for MERCURIUS appearing in the v0.1 paper §Validation table. Without it, the implementation is "passes its own tests" — REBOUND parity is what closes the gap to "matches the canonical reference within numerical precision." The other rebound-parity scenarios (Kepler, Pythagorean, figure-8, retrograde) validate IAS15; this one is the first to validate a federated multi-integrator construction.

The framing is **validation**, not competition: tolerances are set by the precision the physics admits given two independent implementations of the same algorithm. Any divergence beyond those tolerances signals an algorithmic bug in the apsis port, not a performance ranking.

---

## Protocol *(declared a priori, before any code lands)*

### Initial conditions

Units: solar AU-year (`UnitSystem::solar()`), G ≈ 4π². All bodies start at t = 0 with the centre of mass at the origin and zero net momentum (apsis applies the COM shift; REBOUND uses `sim.move_to_com()`).

| Body | Mass (M_sun) | Heliocentric a (AU) | e | i (rad) | True anomaly at t=0 (rad) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sun | 1.0 | — | — | — | — |
| Jupiter | 9.55 × 10⁻⁴ | 5.20 | 0.0 | 0.0 | 0.0 |
| Saturn | 2.86 × 10⁻⁴ | 9.58 | 0.0 | 0.0 | π/2 |
| Uranus | 4.37 × 10⁻⁵ | 19.18 | 0.0 | 0.0 | π |
| Neptune | 5.15 × 10⁻⁵ | 30.07 | 0.0 | 0.0 | 3π/2 |
| Test particle | 1.0 × 10⁻⁹ | 4.20 | 0.40 | 0.05 | 0.0 |

The test particle's `a = 4.20 AU`, `e = 0.40`, `i = 0.05 rad` give an apoapsis at `5.88 AU` (above Jupiter's circular orbit) and periapsis at `2.52 AU` (below). Inclination of 0.05 rad keeps the encounter geometry 3D-non-trivial (out-of-plane). This guarantees the test particle's orbit crosses Jupiter's during the 1000-year integration; the exact encounter times are dictated by the resonance structure and are not specified — the parity gates on conservation and final orbital elements, both of which are insensitive to encounter phase.

The four outer planets are placed on circular coplanar orbits at named heliocentric distances (Murray & Dermott §3 canonical values). Eccentricity and inclination set to zero so the IC are bit-reproducible across both implementations from the (a, M, ν) → (q, v) Keplerian conversion.

### Integrator settings

| Parameter | apsis Mercurius | REBOUND `MERCURIUS` |
| --- | --- | --- |
| Outer dt | 0.01 yr (≈ T_J / 1186) | `sim.dt = 0.01` |
| Hill multiplier α (`r_crit_hill`) | 3.0 (`with_alpha(3.0)` or default) | `sim.ri_mercurius.r_crit_hill = 3` |
| Changeover function | `L_mercury` (built-in) | default `L_mercury` |
| Inner integrator (encounter step) | `Ias15` with default tolerance | REBOUND IAS15 default |
| Coordinate convention | democratic-heliocentric (built-in) | built-in DH |
| Force model | direct O(N²) (Mercurius bypasses ctx.force) | direct |
| Exact finish time | not enforced (consumed_dt = dt every call) | `sim.exact_finish_time = 1` for the per-sample query |

### Run parameters

- Total integration: 1000 years (≈ 84 Jupiter orbits, ≈ 117 test-particle orbits at a=4.20).
- Output cadence: 1 sample per year, plus initial state — 1001 samples per body per run.
- Output format: per-body (pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, mass) + total energy + total angular momentum z-component, at every sample time.

### Hypothesis

#### Tier 1 — Conservation parity *(hard gate)*

| Metric | Bound | Rationale |
| --- | --- | --- |
| Cross-implementation ΔE/E₀ peak (Mercurius − REBOUND, sample-wise absolute) | ≤ 5 × 10⁻⁹ | Mercurius lab notebook §Tier 2 spec. Two independent implementations of a 2nd-order symplectic-class scheme should agree on energy at the IAS15-saturated floor times the encounter-step engagement rate. |
| Cross-implementation ΔLz/Lz₀ peak (absolute) | ≤ 10⁻¹⁰ | Angular momentum is exactly conserved by the analytical Kepler drift; the encounter step's IAS15 conserves Lz at machine precision. Cross-impl floor sits below energy. |
| Per-side ΔE/E₀ peak (each implementation independently, absolute) | ≤ 10⁻⁸ | Mercurius is a 2nd-order method; energy oscillates around the initial value at scale set by the K-weighted kick truncation. |

#### Tier 2 — Test-particle orbital element parity *(hard gate)*

After 1000 years, the test particle has experienced at least one Jupiter encounter (the orbital geometry guarantees it). Compare its osculating elements at t = 1000 yr:

| Element | Bound (relative) |
| --- | --- |
| Semi-major axis a | ≤ 10⁻⁵ |
| Eccentricity e | ≤ 10⁻⁵ |
| Inclination i | ≤ 10⁻⁵ |

The bound is loosened from the Mercurius lab notebook's 10⁻⁶ to 10⁻⁵ to absorb cross-implementation drift in the IAS15 controller's adaptive step sequence during the encounter — both implementations make IAS15-precision-floor decisions about when to refine, but those decisions are not bit-deterministic across independent codebases. The element values themselves should agree; the path through phase space may differ.

#### Tier 3 — Reference-side conservation sanity *(soft gate)*

REBOUND-side ΔE/E₀ and ΔLz/Lz₀ peak should match the published REBOUND MERCURIUS conservation behaviour on Solar-System scenarios (Rein et al. 2019 §3 reports ~10⁻⁹ to 10⁻¹⁰ on similar runs). If REBOUND-side conservation fails, the parity test is uninformative — either the scenario is mis-specified or one of the integrator settings is wrong.

### Methodology

Three-side test infrastructure following the existing `validation/rebound-parity/{kepler,figure8,pythagorean,retrograde}/` pattern:

1. **apsis side** (`crates/apsis/examples/rebound_parity_mercurius.rs`): instantiates the scenario, runs Mercurius for 1000 years, writes `out/apsis.csv` with one row per sample.
2. **REBOUND side** (`validation/rebound-parity/mercurius/rebound_side.py`): reads apsis's actual sample times from `out/apsis.csv`, sets up the same scenario in REBOUND, runs MERCURIUS landing at apsis's sample times via `exact_finish_time = 1`, writes `out/rebound.csv`.
3. **Comparator** (`validation/rebound-parity/mercurius/compare.py`): loads both CSVs, computes the Tier 1/2/3 metrics, exits 0 iff every gated metric is within tolerance, writes `out/comparison.json` with the structured report.

The orchestrator `validation/rebound-parity/mercurius/run.py` chains the three steps and exits with the comparator's code.

#### Why orbital invariants for the test particle, not |Δr|

Same argument as the Kepler parity notebook (`docs/experiments/2026-04-25-rebound-parity-kepler.md` §Pilot Analysis): IAS15 inside the encounter step is not bit-deterministic across independent implementations, and any cross-impl phase drift accumulated through the encounter would dominate the point-wise position drift Δr at the comparison time. Orbital elements (a, e, i) are the physically meaningful signals — they are constants of pure Kepler motion in the no-encounter regime and change deterministically (per the encounter geometry) when the test particle crosses Jupiter's orbit. Reporting Δr informationally is acceptable for transparency; gating on it is not.

---

## Results

*Populated after the run completes on the Ubuntu environment with REBOUND 4.x installed.*

### Tier 1 — Conservation parity

| Metric | apsis | REBOUND | Δ | Bound | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| ΔE/E₀ peak (per side, absolute) | TBD | TBD | — | ≤ 10⁻⁸ each | TBD |
| Cross-impl ΔE/E₀ peak (absolute) | — | — | TBD | ≤ 5 × 10⁻⁹ | TBD |
| Cross-impl ΔLz/Lz₀ peak (absolute) | — | — | TBD | ≤ 10⁻¹⁰ | TBD |

### Tier 2 — Test-particle orbital element parity

| Element | apsis | REBOUND | Δ relative | Bound | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| a | TBD | TBD | TBD | ≤ 10⁻⁵ | TBD |
| e | TBD | TBD | TBD | ≤ 10⁻⁵ | TBD |
| i | TBD | TBD | TBD | ≤ 10⁻⁵ | TBD |

### Tier 3 — Reference-side sanity

| Diagnostic | REBOUND | Expected | Status |
| --- | ---: | --- | --- |
| ΔE/E₀ peak (absolute) | TBD | ~10⁻⁹ to 10⁻¹⁰ (Rein et al. 2019 §3) | TBD |
| ΔLz/Lz₀ peak (absolute) | TBD | ~10⁻¹² | TBD |

---

## Interpretation

*Populated after Results land.*

---

## Decision

*Populated after Interpretation lands. Possible outcomes:*

- *Pass all gates → Mercurius enters the v0.1 paper §Validation table; close PR.*
- *Pass Tier 1, fail Tier 2 → IAS15 controller drift inside encounter step is dominating; loosen Tier 2 bound with explicit drift accounting, or re-derive the bound from first principles.*
- *Fail Tier 1 → algorithmic bug in apsis Mercurius port; bisect against the REBOUND source until the divergence is traced.*

---

## References

- Rein, H., Hernandez, D. M., Tamayo, D., & Brown, G. (2019). *Hybrid symplectic integrators for planetary dynamics.* MNRAS, 489(4), 4632–4640.
- Mercurius implementation lab notebook: `docs/experiments/2026-05-13-mercurius-hybrid.md`.
- Kepler parity protocol template: `docs/experiments/2026-04-25-rebound-parity-kepler.md`.
- REBOUND `integrator_mercurius.c` (canonical reference; copy at `~/Desktop/rebound/mercurius.c`).
