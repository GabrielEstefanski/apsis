# Long-horizon Mercury 1PN — millennium-scale precession reproducibility

**Date:** 2026-05-13
**Subject:** Sun + Mercury under apsis IAS15 with the apsis-1pn perturbation registered, integrated for 1000 years (~4150 orbits) and compared against the closed-form Schwarzschild test-particle perihelion advance `Δω = 6πGM/(c²a(1−e²))` per orbit. Cross-validated against `Mercurius + apsis-1pn` once the Mercurius perturbation contract fix (PR #86) merges — this is the headline FPM federation demonstration: integrator-of-integrators composing with a first-class perturbation produces the same physics as a homogeneous integrator + same perturbation, both matching the analytical GR prediction at the f64 precision floor.

**Status:** Protocol declared *a priori*, before the long-horizon run. Builds on the existing 500-orbit gate (`crates/apsis-1pn/tests/mercury_precession_gate.rs`, ~1 ppm developer-hardware) and Issue #29 (long-horizon Mercury 1PN).

**Branch:** `feat/mercury-1pn-long-horizon`, branched from `develop`. The Mercurius half of Tier 2 stacks on PR #86 (Mercurius perturbation hole fix) and is added once that merges.

---

## Abstract

The 500-orbit Mercury 1PN gate establishes that apsis-1pn reproduces the GR perihelion advance at ~1 ppm on developer hardware over ~120 years. This experiment extends the horizon by a factor of ~8.3× (to 1000 years, ~4150 orbits) to demonstrate that the precision is **stable** over horizons relevant to GR observation campaigns, and adds a cross-integrator parity check (IAS15 vs Mercurius, both with the same apsis-1pn perturbation registered) to validate that the federation contract — any apsis integrator composes with any apsis perturbation — holds at the precision floor of the lead empirical result.

For the v0.1 paper this is the headline FPM demonstration: the federation thesis is not "two integrators happen to give similar answers", it is "the federation **contract** delivers the same physics regardless of which first-class artifact occupies each slot, to the precision the field admits." The simplest perturbation (1PN) composing with the simplest hybrid integrator (Mercurius + 1PN, where Sun-Mercury is non-encountering and Mercurius reduces to its WH-like outer step) is the floor of the federation; if this fails, the thesis fails.

---

## Motivation

Three claims chain into this experiment:

1. **apsis-1pn reproduces GR.** Established at 500 orbits / ~1 ppm by the existing gate.
2. **The precision is stable.** Issue #29's open question — does the agreement hold at 4× the horizon, or does cumulative round-off swamp the signal? Answer here.
3. **The federation contract works at the precision floor.** Not just "WH + 1PN works" or "IAS15 + 1PN works" — *any* integrator + *any* perturbation. Mercurius is the most architecturally complex integrator in the apsis zoo (rewind hybrid with embedded IAS15 sub-integrator); if 1PN composes cleanly through Mercurius's K-weighted half-kicks (post PR #86), the contract is demonstrated at the level the paper needs.

The paper positioning ([[project_paper_positioning]]) makes the federation an explicit point of differentiation against monolithic codebases. A clean Mercurius + 1PN parity result is the demonstration that the "federated FPM" framing is operationally true, not just rhetorical.

### What this experiment is NOT testing

- Not a new physics derivation. apsis-1pn shipped under PR #?? and is validated by the existing gate; this extends the horizon, not the formula.
- Not a Mercurius implementation test. Mercurius shipped under PR #83 with structural / cost validation and PR #85 with REBOUND parity at the conservation level.
- Not a WHFast integrator scope. Adding WHFast to the integrator zoo is queued as a separate axis (Rein & Tamayo 2015 — needed for 10⁹-orbit horizons; the 4150-orbit horizon here is comfortably inside both IAS15 and apsis WH-class envelopes).

---

## Locked design decisions

| Question | Decision | Rationale |
| --- | --- | --- |
| Units | Canonical Hénon (`UnitSystem::canonical`, G=1) | Matches the existing 500-orbit gate convention; `PostNewtonian1PN::solar_units()` carries `c = 10065.13` in time-units of year/(2π), which is exactly the canonical Hénon time unit. Re-using these conventions removes "did the units rescale correctly" as a source of error. |
| Bodies | Sun (m=1, unsoftened) + Mercury (m=1.66 × 10⁻⁷, unsoftened) | Sun + 1 planet is the cleanest Schwarzschild test-particle case; matches the analytical formula's assumptions. Mercury mass set to physical value so the orbital frequency is exactly Mercury's. Both unsoftened so the Newtonian baseline is bit-Keplerian (Plummer softening would swamp the 1PN signal by ~2000×, per `apsis-1pn` `kernel_requirements`). |
| Mercury orbital elements | a = 0.387098, e = 0.20563, i = 0 (2D), starting at periapsis | a, e are physical Mercury values (Murray & Dermott Table A.2). i = 0 keeps the precession measurement in 2D where the periapsis-orientation extraction is least ambiguous. Starting at periapsis sets `ω_initial` cleanly to 0 by symmetry. |
| Horizon | 1000 years = 4153 Mercury orbits = 6283.19 canonical time units | 8.3× the existing 500-orbit gate. Long enough to demonstrate stability; short enough that one IAS15 run takes ~1 minute on Cell A. Issue #29 framed as "century-scale", the choice here is "millennium-scale" to make the cumulative drift signal big enough to dominate any per-orbit noise. |
| Outer dt (IAS15 first-call seed) | 1 × 10⁻⁴ canonical | Matches the existing 500-orbit gate. ~15 134 IAS15 sub-steps per Mercury orbit, ~6.3 × 10⁷ over the full integration. |
| Outer dt (Mercurius) | 1 × 10⁻³ canonical (10 sub-steps per Mercury orbit) | Mercurius is fixed-step. Smaller dt is wasted work; larger dt loses K-weighted-kick resolution. 10 sub-steps per Mercury orbit gives ~41 530 outer steps over 1000 yr — fast even at the per-step cost of the analytical Kepler solver. |
| Sample cadence | 1 sample per Mercury orbit (4154 samples per side) | Captures the cumulative ω trajectory at orbit-resolution, plenty for the linear-in-time precession signal. |
| Perturbation | `PostNewtonian1PN::solar_units()` registered via `System::add_perturbation` | Same instance the existing 500-orbit gate uses. Single point of physics truth — any change here would invalidate both the gate and this experiment together. |

---

## Algorithm

### Reference: closed-form GR perihelion advance

Schwarzschild test-particle 1PN expansion (Will 1993 §6, Murray & Dermott §8.1):

$$
\Delta\omega_\text{GR per orbit} = \frac{6\pi G M_\odot}{c^2 \, a \, (1 - e^2)}
$$

In the canonical-Hénon unit system: `G = 1`, `M_sun = 1`, `c = C_SOLAR_UNITS = 10065.13`. Cumulative precession over `N` orbits: `Δω_GR(N) = N · Δω_per_orbit`.

For the locked Mercury IC (`a = 0.387098`, `e = 0.20563`):

$$
\Delta\omega_\text{per orbit} = \frac{6\pi}{(10065.13)^2 \cdot 0.387098 \cdot (1 - 0.20563^2)} = 5.02 \times 10^{-7} \text{ rad/orbit}
$$

Over 4153 orbits: `Δω_GR(end) = 2.084 × 10⁻³ rad ≈ 7.16 arcmin`. Easily resolvable above the f64 noise floor.

### Measured: osculating ω from the eccentricity vector

Per sample, extract the osculating periapsis orientation from the Sun-Mercury relative state via `apsis::physics::orbital::compute_elements`:

$$
\vec{e} = \frac{\vec{v} \times \vec{h}}{\mu} - \hat{r}, \quad \omega = \mathrm{atan2}(e_y, e_x)
$$

Subtract `ω_initial = 0` (Mercury starts at periapsis along +x by construction); unwrap the resulting trajectory to remove `2π` jumps. The cumulative measured precession at sample `k` is `ω_unwrapped(k)`.

---

## Protocol

### Hypothesis

#### Tier 1 — IAS15 + apsis-1pn precision *(hard gate)*

| Metric | Bound | Rationale |
| --- | --- | --- |
| `\|Δω_measured(end) − Δω_GR(end)\| / \|Δω_GR(end)\|` | ≤ 10⁻⁵ (10 ppm) | The existing 500-orbit gate hits ~1 ppm developer-side. Cumulative round-off scales as `√N_substeps · ε ~ 10⁻¹²` absolute, i.e. ~10⁻⁹ relative on the 2 × 10⁻³ rad signal. The 10 ppm bound has ~4 orders of margin against the round-off floor and ~1 order against the 1-ppm developer-hardware result. |
| Per-orbit precession trajectory linearity | R² ≥ 0.99999 against `Δω_GR_per_orbit · k` | Independent secondary check: GR predicts a strictly linear Δω(t); any non-linear residual signals a non-1PN perturbation (numerical drift, controller chatter, etc.). |

#### Tier 2 — Mercurius + apsis-1pn precision *(hard gate, post PR #86)*

Same bounds as Tier 1, applied to a separate Mercurius+1PN run on the same scenario. For Sun + Mercury (N = 2) the encounter step never fires; Mercurius reduces to its WH-like outer step (K-weighted half-kicks with 1PN folded in via the perturbation contract + analytical Kepler drift around the Sun). The fact that *no encounter ever fires* is itself part of the test — Mercurius's structural overhead must not introduce drift that IAS15 doesn't have.

#### Tier 3 — Cross-integrator parity *(hard gate, post PR #86)*

| Metric | Bound | Rationale |
| --- | --- | --- |
| `\|Δω_IAS15(end) − Δω_Mercurius(end)\| / \|Δω_GR(end)\|` | ≤ 5 × 10⁻⁵ (50 ppm) | Cumulative cross-integrator drift bound: each integrator agrees with GR at ≤ 10⁻⁵; cross-integrator drift is bounded above by the sum of independent drifts (~2 × 10⁻⁵), with one decade slack for second-order coupling between adaptive (IAS15) and fixed-step (Mercurius) controllers. |

### Methodology

1. **apsis-side IAS15 run** (`crates/apsis/examples/mercury_1pn_long_horizon_ias15.rs`): builds the scenario, registers `PostNewtonian1PN::solar_units()`, integrates 1000 yr at the locked dt, writes per-orbit (t, x, y, vx, vy, a_osc, e_osc, ω_osc) to `validation/mercury-1pn-long-horizon/out/ias15.csv`.

2. **apsis-side Mercurius run** (`crates/apsis/examples/mercury_1pn_long_horizon_mercurius.rs`, added post PR #86): same scenario, same perturbation, Mercurius integrator, writes to `out/mercurius.csv`.

3. **Comparator** (`validation/mercury-1pn-long-horizon/compare.py`): loads both CSVs, unwraps ω, computes the three Tier metrics against the analytical GR prediction, exits 0 iff all gates pass.

4. **Orchestrator** (`validation/mercury-1pn-long-horizon/run.py`): runs the cargo examples then the comparator. Same pattern as the existing `validation/rebound-parity/{kepler,figure8,...}` scenarios.

---

## Results

*Populated post-run.*

### Tier 1 — IAS15 + apsis-1pn

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| Δω relative error vs GR (end) | TBD | ≤ 10⁻⁵ | TBD |
| Per-orbit linearity R² | TBD | ≥ 0.99999 | TBD |

### Tier 2 — Mercurius + apsis-1pn

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| Δω relative error vs GR (end) | TBD | ≤ 10⁻⁵ | TBD |
| Per-orbit linearity R² | TBD | ≥ 0.99999 | TBD |

### Tier 3 — Cross-integrator parity

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| `\|Δω_IAS15 − Δω_Mercurius\| / \|Δω_GR\|` (end) | TBD | ≤ 5 × 10⁻⁵ | TBD |

---

## Interpretation

*Populated post-results.*

---

## Decision

*Populated post-interpretation. Possible outcomes:*

- **All three tiers pass** → headline result for `paper.md` §Validation: "apsis Mercury 1PN reproduces the GR perihelion advance at X ppm over 1000 years under both IAS15 and Mercurius, with cross-integrator parity at Y ppm. The federation contract operates at the precision the field admits."
- **Tier 1 passes, Tier 2 fails** → bug in PR #86's Mercurius perturbation wiring; bisect against the regression test `ctx_perturbations_are_honored_by_interaction_step`.
- **Tier 1 fails** → either the IAS15 long-horizon precision claim was over-stated by the 500-orbit gate (cumulative drift exceeds expectation), or a regression in `apsis-1pn` since the gate landed.
- **Tier 3 fails despite Tiers 1 and 2 passing** → IAS15 vs Mercurius drift in `Δω` larger than independent-integrator-sum bound predicts; investigate adaptive-vs-fixed step coupling with the perturbation.

---

## References

- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics.* Cambridge University Press, §6 (Schwarzschild perihelion advance).
- Murray, C. D., & Dermott, S. F. (1999). *Solar System Dynamics.* Cambridge University Press, §8.1 + Table A.2 (Mercury orbital elements).
- Rein, H., & Spiegel, D. S. (2015). *MNRAS* 446, 1424 (IAS15).
- Rein, H., Hernandez, D. M., Tamayo, D., & Brown, G. (2019). *MNRAS* 489, 4632 (MERCURIUS).
- Existing 500-orbit gate: `crates/apsis-1pn/tests/mercury_precession_gate.rs`.
- Mercurius implementation lab notebook: `docs/experiments/2026-05-13-mercurius-hybrid.md`.
- Mercurius REBOUND parity: `docs/experiments/2026-05-13-rebound-parity-mercurius.md`.
- Issue #29 (long-horizon Mercury 1PN, century-scale).
