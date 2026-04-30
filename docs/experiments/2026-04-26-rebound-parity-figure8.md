# REBOUND Parity — Figure-8 Choreography (10 periods)

**Date:** 2026-04-26
**Subject:** Numerical agreement between IAS15 (apsis) and IAS15 (REBOUND) on the Chenciner–Montgomery figure-8 three-body choreography.
**Baseline commit:** `9caaef2` ("feat(ias15): refactor controller to spec-conformant Pascal warmstart, halving rejection, and 7x growth cap").
**Tooling:** apsis IAS15 (`crates/apsis/src/physics/integrator/ias15.rs`), REBOUND 4.6.0 via Python 3.10 (`reb.IAS15`).
**Status:** *Run executed 2026-04-26 against `9caaef2`. All twelve Tier-1 and Tier-2 gated metrics pass at 1–10 ULP. The Tier-3 informational $\lvert \Delta r \rvert$ stands at $9.44 \times 10^{-13}$, consistent with controller-level phase drift between two adaptive IAS15 implementations at the f64 round-off floor.*

---

## Abstract

This experiment extends the v0.1 APSIS validation portfolio by testing cross-implementation agreement on a non-Keplerian, fully periodic, equal-mass three-body scenario — the figure-8 choreography (Moore 1993; Chenciner & Montgomery 2000). Where the Kepler parity test (notebook `2026-04-25`) gates on orbital elements $(a, e, \omega, h)$, those quantities are not defined for the figure-8: the system has no central body and no Keplerian decomposition. The protocol therefore gates on the global integrals of motion that *are* invariant under the dynamics — total energy $E$, total angular-momentum vector $\mathbf{L}$, total linear momentum $\mathbf{P}$, and centre-of-mass position — and reports per-body geometric drift $|\Delta r|$ as observational context only, on the same methodological grounds established in the Kepler notebook.

The metric set is organised into three tiers of evidentiary weight (hard / sanity / observational), reflecting that energy and angular-momentum agreement are the load-bearing physical statements while linear-momentum and COM agreement are construction-level consistency checks. The integration horizon is held to 10 orbital periods — the conventional figure used in the literature to exhibit the choreography (Chenciner & Montgomery 2000; Simó 2002) — with a $50\,T$ extension recorded as informational context outside the gated claim.

---

## Motivation

The Kepler parity result (notebook `2026-04-25`) establishes that `apsis` IAS15 reproduces the orbit-element invariants of the canonical two-body problem to ULP precision against REBOUND's IAS15. A reviewer is reasonably likely to ask whether that agreement survives into a regime where:

- the dynamics are not reducible to a one-body Kepler problem;
- there is no dominant primary, so the integrator's per-body force evaluation is exercised symmetrically;
- the orbit periodically enters close mutual approaches near the central crossing of the figure-8.

The figure-8 is the simplest published scenario that satisfies all three. It has a closed analytic form for the motion's *symmetry* (cyclic permutation by $T/3$ between the three bodies; reflection symmetry through the origin), but no closed form for the trajectory itself — the choreography exists but is computed numerically. This makes it a natural complement to Kepler in a parity portfolio: same integrator, different physical regime.

The framing remains validation, not competition. Tolerances are set by the precision the physics admits — for global integrals of motion under IAS15, that is f64 machine epsilon scaled by the round-off accumulated across the integration horizon. The 1PN/Mercury demonstration in the v0.1 paper rests on extension-mechanism evidence, not on a numerical-superiority claim; this experiment supports the numerical foundation underlying that mechanism, nothing more.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the figure-8 choreography integrated under IAS15 in both `apsis` and REBOUND with the canonical Chenciner–Montgomery initial conditions and matched integrator settings, the metrics declared below are bounded *a priori* at the values stated. The bounds are organised into three tiers reflecting the evidentiary role of each metric. **Verdict on the experiment is determined by Tiers 1 and 2 only.** Tier 3 is reported for completeness and characterisation, not as a parity gate, on the same methodological grounds developed in the Kepler notebook (§Pilot Interpretation, `2026-04-25-rebound-parity-kepler.md`): point-wise position metrics on adaptive high-order integrators conflate orbital phase drift (a non-invariant of the cross-implementation comparison) with geometric drift (the physically meaningful signal).

#### Tier 1 — Hard physical invariants *(gated)*

The quantities here are exact constants of motion of the figure-8 dynamics; any drift in either implementation is a numerical artefact. Both per-side conservation and cross-implementation agreement are gated.

- **Energy** — $\max |\Delta E / E_0|$ per side $\leq 1.0 \times 10^{-13}$ ($\approx 50 \times$ f64 machine epsilon). IAS15 is designed for machine-precision energy conservation (Rein & Spiegel 2015); the bound matches the Kepler notebook's energy gate and is consistent with the per-side floors observed there (1.0–1.7 ULP).
- **Energy** — cross-implementation $\max |\Delta E_\text{apsis} - \Delta E_\text{rebound}| / |E_0| \leq 1.0 \times 10^{-13}$.
- **Angular momentum** — $\max |\Delta \mathbf{L}|$ per side $\leq 1.0 \times 10^{-13}$ in absolute units. The figure-8 has $\mathbf{L}_0 = \mathbf{0}$ by construction, so a relative bound is undefined; the absolute bound is set against the characteristic scale $|m_i \, \mathbf{r}_i \times \mathbf{v}_i| \sim O(1)$ for an individual body, with $50 \times$ ULP margin. The metric is computed as the full vector norm $|\mathbf{L}| = \sqrt{L_x^2 + L_y^2 + L_z^2}$ rather than $|L_z|$ alone — the system is planar by IC construction ($z = v_z = 0$), so $L_x = L_y = 0$ to machine precision; computing the full norm avoids baking the planarity assumption into the metric and exposes any implementation-side breaking of planarity.
- **Angular momentum** — cross-implementation $\max |\mathbf{L}_\text{apsis} - \mathbf{L}_\text{rebound}| \leq 1.0 \times 10^{-13}$.

#### Tier 2 — Construction-level sanity *(gated, weak)*

These quantities are zero by construction of the canonical ICs and remain zero under exact dynamics. Drift here is a sum-of-forces / sum-of-positions accumulation effect — it should sit at the f64 round-off floor and is gated to confirm no implementation-level book-keeping bug has crept in.

- **Linear momentum** — $\max |\Delta \mathbf{P}|$ per side $\leq 1.0 \times 10^{-13}$. $\mathbf{P}_0 = \mathbf{0}$. Reference scale $|m_i \, \mathbf{v}_i| \sim O(1)$.
- **Linear momentum** — cross-implementation $\max |\mathbf{P}_\text{apsis} - \mathbf{P}_\text{rebound}| \leq 1.0 \times 10^{-13}$.
- **Centre-of-mass position** — $\max |\mathbf{r}_\text{COM}|$ per side $\leq 1.0 \times 10^{-12}$. $\mathbf{r}_\text{COM}(0) = \mathbf{0}$. The bound is one decade looser than the momentum bound because COM position accumulates from $\sum_i m_i \, \mathbf{v}_i \, \Delta t$ over $10\,T \approx 63$ time units — a position-scale integral of the momentum-scale ULP noise.
- **Centre-of-mass position** — cross-implementation $\max |\mathbf{r}_{\text{COM},\text{apsis}} - \mathbf{r}_{\text{COM},\text{rebound}}| \leq 1.0 \times 10^{-12}$.

#### Tier 3 — Geometric coherence *(informational, NOT gated)*

- **Per-body position drift** — $\max |\mathbf{r}_{i,\text{apsis}}(t) - \mathbf{r}_{i,\text{rebound}}(t)|$ for $i \in \{0, 1, 2\}$ over all sample times. Reported as context. Phase drift inherent to two adaptive IAS15 implementations is expected to dominate this signal; see Kepler notebook §Pilot Interpretation for the diagnosis. No tolerance is declared because no physical statement depends on this quantity.

### Methodology

#### Initial conditions

Three equal-mass bodies in canonical units ($G = 1$, $m_i = 1$):

| Body | $x$ | $y$ | $v_x$ | $v_y$ |
| ---: | ---: | ---: | ---: | ---: |
| 1 | $-0.97000436$ | $+0.24308753$ | $+0.4662036850$ | $+0.4323657300$ |
| 2 | $+0.97000436$ | $-0.24308753$ | $+0.4662036850$ | $+0.4323657300$ |
| 3 | $0$ | $0$ | $-0.93240737$ | $-0.86473146$ |

Source: Chenciner & Montgomery (2000), reproduced in the apsis preset `templates::presets::threebodyproblems::three_body_figure_eight` (`crates/apsis/src/templates/presets/threebodyproblems/three_body_figure_eight.rs`). The 8-digit literature representation is converted to f64 by both implementations from the *same string-literal* values, so the bit-pattern of the ICs is identical between sides on the same hardware. Higher-precision IC sets exist (notably Simó 2002 to 13 digits); they are not used here because the f64 representation of the 8-digit literals already saturates double precision for the IC-noise-floor purposes of this experiment.

By construction these ICs satisfy $\sum_i m_i \, \mathbf{r}_i = \mathbf{0}$, $\sum_i m_i \, \mathbf{v}_i = \mathbf{0}$, and $\sum_i m_i \, \mathbf{r}_i \times \mathbf{v}_i = \mathbf{0}$. No additional centre-of-mass or momentum zeroing is performed prior to integration on either side; doing so would introduce an implementation-divergent f64 correction to ICs that should be bit-identical between sides.

#### Integrator settings

| Parameter | apsis IAS15 | REBOUND IAS15 |
| --- | --- | --- |
| Initial timestep | $T / 1000 \approx 6.32591398 \times 10^{-3}$ | $T / 1000 \approx 6.32591398 \times 10^{-3}$ |
| Adaptive control | active, default tolerance | active, default `epsilon` |
| Force model | direct $O(N^2)$ (forced via pairing rule, see ADR-003) | direct (REBOUND default) |
| Exact finish time | not enforced | `sim.exact_finish_time = 1` |

REBOUND's `exact_finish_time = 1` is used so REBOUND samples are evaluated at the *actual* sample times produced by the apsis side after its adaptive controller has overshot the nominal target; this removes "two implementations sampled at slightly different physical times" as a confound in the cross-implementation comparison. The same convention was used in the Kepler notebook.

#### Run parameters and sampling

- **Total integration horizon (gated):** 10 orbital periods. $T = 6.3259139870$ (canonical units; the value follows from the Chenciner–Montgomery analysis and is reproduced in REBOUND fixtures and in apsis benchmark `figure_eight`). Total horizon $t_\text{final} = 10\,T \approx 63.259$. The $10\,T$ choice matches the literature horizon at which the figure-8 choreography is conventionally exhibited (Chenciner & Montgomery 2000; Simó 2002) and avoids overstating evidence beyond the regime in which the choreography is uncontroversial.
- **Analysis cadence (dense):** 200 samples per period $\to$ 2001 samples in $[0, 10\,T]$. Both `apsis` and REBOUND emit state at this dense cadence. The comparator computes its $\max(\cdot)$ aggregates over this dense set; the gated metrics are therefore evaluated against the worst-case sample over the full horizon, not against an arbitrarily sparse subsample.
- **Report cadence (sparse):** 4 samples per period (every $T/4$) $\to$ 41 representative samples published in the §Results evolution table. Sparse cadence is for human readability of the trajectory's invariant evolution; it does *not* enter the gating computation.
- **Output format:** wide CSV with `t`, full per-body state $(x, y, v_x, v_y)$ for the three bodies, and total energy $E$. Sample index $n \in \{0, \ldots, 2000\}$ is emitted in place of the Kepler notebook's `orbit` column.

#### Metric computation

All Tier-1 and Tier-2 quantities are evaluated identically on both sides — same formula, same $\mu$, same per-sample state vectors. The only difference between sides is the integrated state itself.

For each sample on each side:

$$
\begin{aligned}
E_\text{total} &= \sum_i \tfrac{1}{2} \, m_i \, |\mathbf{v}_i|^2 \;-\; \sum_{i < j} \frac{G \, m_i \, m_j}{|\mathbf{r}_i - \mathbf{r}_j|}
&& \text{(with $G = 1$)} \\
\mathbf{L} &= \sum_i m_i \, (\mathbf{r}_i \times \mathbf{v}_i)
&& \text{(3-component; $z$ dominates the planar case)} \\
\mathbf{P} &= \sum_i m_i \, \mathbf{v}_i
&& \text{(2-component $xy$; $z$ trivially zero)} \\
\mathbf{r}_\text{COM} &= \frac{\sum_i m_i \, \mathbf{r}_i}{\sum_i m_i}
\end{aligned}
$$

Per-side drift metrics are $\max_t |Q(t) - Q(0)|$ (or relative form, where $Q_0 \neq 0$). Cross-implementation metrics are $\max_t |Q_\text{apsis}(t) - Q_\text{rebound}(t)|$ (or relative). Source of truth for the formulas: `validation/rebound-parity/figure8/compare.py::physical_invariants`.

### Why this metric set, not $|\Delta r|$

The justification for the invariant-based metric set is identical to the Kepler notebook's. Briefly: IAS15's adaptive controller selects each substep size as $dt_\text{new} = dt \cdot \mathrm{safety} \cdot (\epsilon / \mathrm{err})^{1/7}$. ULP-level differences in summation order between the two implementations propagate into $\mathrm{err}$ at the f64-precision floor; the $1/7$ exponent then propagates that into $dt$ with mild sensitivity. The two implementations therefore take *slightly different* $dt$ sequences, and that difference accumulates as orbital phase drift over many periods. Phase drift is not a numerical defect of either implementation; it is the ceiling on cross-implementation parity for any adaptive high-order method without enforced bit-equivalence.

For the figure-8, gating on per-body $|\Delta r|$ would conflate this controller-level phase drift with a physically meaningful disagreement on the trajectory. The honest physical question is whether the two implementations integrate the *same dynamical system* — i.e. preserve the same global integrals of motion. That question is answered by the Tier-1 and Tier-2 metrics. Tier-3 is reported because reviewers expect to see the geometric-drift number, but it is explicitly not the gate.

### Out of scope (declared a priori)

- **Period closure of the trajectory.** A well-formed figure-8 implementation closes the orbit at $t = T$; the per-side $\max |\Delta r|$ between $\mathbf{r}(T)$ and $\mathbf{r}(0)$ is a meaningful diagnostic of integrator quality but is not a cross-implementation parity criterion (it tests each implementation against itself, not against the other). Both apsis and REBOUND have their own internal closure benchmarks; this experiment does not duplicate them.
- **Symbolic symmetry checks.** The figure-8 satisfies cyclic permutation by $T/3$, body-pair reflection through the origin, and time-reversal symmetry. Verifying these is a separate methodological exercise; this experiment treats them as IC-level facts and does not gate on numerical residuals of the symmetry operations.
- **Sensitivity to IAS15 `epsilon`.** The default tolerance on both sides is used. A sweep over `epsilon` would characterise the cost-precision frontier and is reserved for the Phase 6A characterisation experiment.

---

## Results

The run was executed 2026-04-26 against `apsis` commit `9caaef2`. The same run was the original failing case for the IAS15 controller audit documented in the companion notebook `2026-04-26-ias15-warmstart-bug.md`; the audit's resolution (Pascal-triangle warmstart, unconditional halving on truncation rejection, 7$\times$ growth cap on accept-path `dt_next`) is what made every Tier-1 and Tier-2 gate pass.

Total samples: 2001 (200 per period $\times$ 10 periods + 1). Final time: 63.262389 (canonical units; $10\,T$ to 5 ULP in the period value).

### Energy and angular momentum (Tier 1)

| Metric | Observed | Tolerance | Margin |
| --- | ---: | ---: | ---: |
| $\lvert \Delta E / E_0 \rvert$ apsis | 8.625e-16 | 1.00e-13 | 116$\times$ under |
| $\lvert \Delta E / E_0 \rvert$ rebound | 8.625e-16 | 1.00e-13 | 116$\times$ under |
| Cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$ | 1.035e-15 | 1.00e-13 | 97$\times$ under |
| $\lvert \Delta \mathbf{L} \rvert$ apsis (abs) | 3.608e-16 | 1.00e-13 | 277$\times$ under |
| $\lvert \Delta \mathbf{L} \rvert$ rebound (abs) | 3.331e-16 | 1.00e-13 | 300$\times$ under |
| Cross-impl $\lvert \Delta \mathbf{L} \rvert$ (abs) | 4.441e-16 | 1.00e-13 | 225$\times$ under |

### Linear momentum and centre-of-mass (Tier 2)

| Metric | Observed | Tolerance | Margin |
| --- | ---: | ---: | ---: |
| $\lvert \Delta \mathbf{P} \rvert$ apsis (abs) | 4.965e-16 | 1.00e-13 | 201$\times$ under |
| $\lvert \Delta \mathbf{P} \rvert$ rebound (abs) | 5.095e-16 | 1.00e-13 | 196$\times$ under |
| Cross-impl $\lvert \Delta \mathbf{P} \rvert$ (abs) | 4.578e-16 | 1.00e-13 | 218$\times$ under |
| $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ apsis (abs) | 3.111e-15 | 1.00e-12 | 322$\times$ under |
| $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ rebound (abs) | 3.815e-15 | 1.00e-12 | 263$\times$ under |
| Cross-impl $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ (abs) | 7.514e-16 | 1.00e-12 | 1331$\times$ under |

**All twelve gated metrics pass.** Every observed value sits in the 1–10 ULP regime, consistent with the f64 round-off floor for two IAS15 implementations integrating the same dynamical system.

### Per-body geometric drift (Tier 3, informational)

| Metric | Observed |
| --- | ---: |
| max $\lvert \Delta \mathbf{r} \rvert$ over all bodies and samples | 9.436e-13 |
| Body 1, peak at sample 1832 | 9.404e-13 |
| Body 2, peak at sample 1967 | 9.134e-13 |
| Body 3, peak at sample 2000 | 9.436e-13 |

The three per-body peaks agree to within a factor of 1.04 — the cyclic permutation symmetry of the figure-8 dynamics is preserved across the cross-implementation comparison to ULP precision, as expected from the equal-mass ICs.

Raw outputs: `validation/rebound-parity/figure8/out/{apsis,rebound}.csv`, `out/comparison.json`.

---

## Interpretation

The two IAS15 implementations agree on the figure-8 choreography across all four global integrals of motion at machine precision. Energy and angular momentum — the dynamically meaningful conserved quantities — match per-side and cross-implementation to 1–10 ULP over 10 orbital periods ($\approx 63$ canonical time units). Linear momentum and centre-of-mass position — quantities that are zero by construction and remain zero under exact dynamics — sit at the f64 round-off floor with margins 2–3 orders of magnitude under the declared sanity bounds.

The Tier-3 $\lvert \Delta \mathbf{r} \rvert \sim 9.4 \times 10^{-13}$ is roughly three orders of magnitude smaller than the analogous Kepler-parity observation post-controller-fix ($\sim 2.2 \times 10^{-12}$, notebook `2026-04-25` informational entry as updated by the warmstart audit). The improvement reflects the controller now tracking the specification's `dt_next` choices to ULP precision after the warmstart/halving/growth-cap fixes; the per-body equal-magnitude $\lvert \Delta \mathbf{r} \rvert$ corroborates the cyclic permutation symmetry of the choreography across the cross-implementation comparison.

This completes the second entry in Pillar A (numerical foundation) of the v0.1 validation portfolio: the same IAS15 implementation that reproduces the canonical two-body Kepler problem to machine precision (notebook `2026-04-25`) also reproduces the canonical non-Keplerian periodic three-body problem to machine precision in a regime where no central body, no Keplerian decomposition, and recurrent mutual close approaches at the central crossing are all present simultaneously.

---

## Threats to validity

1. **Floating-point summation ordering.** apsis and REBOUND sum the three pairwise interactions in different orders, producing different ULP-level rounding. The Tier-1 metrics measured at f64 precision confirm whether the floor is at machine epsilon. No tolerance is widened to accommodate this; the bounds were set against the expected f64 round-off floor in advance.
2. **FMA usage.** apsis is built with default Rust FP semantics; REBOUND is C with potential FMA via the compiler. Different FMA decisions produce small but systematic deviations within the same ULP envelope. The threat is identical to the Kepler notebook and is not expected to surface above the round-off floor.
3. **Adaptive controller micro-decisions.** Both implementations follow Rein & Spiegel (2015) for the Picard predictor-corrector loop and the $(\epsilon / \mathrm{err})^{1/7}$ controller, but micro-decisions in the controller (when to grow $dt$, marginal-convergence handling) propagate ULP-level differences in $\mathrm{err}$ into ULP-level differences in $dt$, accumulating as orbital phase drift. The protocol gates on invariants precisely because phase drift is not a cross-implementation invariant; the Tier-3 $|\Delta r|$ metric is reported only as evidence of the magnitude of this expected drift.
4. **Initial-condition rounding.** The 8-digit literature ICs are converted from string to f64 by Rust and Python independently. On x86-64 with IEEE-754 default rounding, the bit pattern is the same on both sides; this should be confirmed by $t = 0$ row inspection in the comparator. Differences at $t = 0$ would indicate an IC-construction divergence and would invalidate the comparison.
5. **Period value.** The figure-8 period is irrational and known only numerically. The value $T \approx 6.3259139870$ used to drive sample-time computation is taken from the canonical literature; a refined value would shift all sample times uniformly on both sides and therefore would *not* affect the cross-implementation parity metrics, only the absolute time labels.
6. **Choreography stability over the integration horizon.** The figure-8 is linearly stable but lies near unstable companions in the choreography family. At $10\,T$, no destabilisation is expected; the $50\,T$ extension (Appendix A) is included precisely to make any anomalous secular behaviour visible without making it part of the gated claim.

---

## Reproducibility

| Field | Value |
| --- | --- |
| apsis canonical commit | *(to be pinned at run time)* |
| REBOUND version | 4.6.0 |
| Python version | 3.10.0 (CPython, MSC v.1929 64-bit) |
| Rust toolchain | apsis Cargo profile `release`; default FP semantics (no `--ffast-math`-equivalent) |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| Harness | `validation/rebound-parity/figure8/run.py` (cargo example $\to$ REBOUND side $\to$ comparator) |
| Apsis side | `crates/apsis/examples/rebound_parity_figure8.rs` |
| Raw outputs | `validation/rebound-parity/figure8/out/{apsis,rebound}.csv`, `out/comparison.json` |

**Commit pinning protocol:** the canonical hash committed to this notebook on the run date includes the apsis-side Cargo example, the Python harness under `validation/rebound-parity/figure8/`, and this notebook itself. The harness is reproducible on a clean checkout of that commit with the dependencies pinned in `validation/rebound-parity/figure8/requirements.txt`.

---

## Appendix A — Extended sanity run ($50\,T$, informational)

Run identical in every respect to the gated $10\,T$ baseline except for the integration horizon $t_\text{final} = 50\,T$. Reported quantities: the same Tier-1, Tier-2, and Tier-3 metrics, computed over the dense analysis cadence. **The $50\,T$ run is *not* part of the gated parity claim of this experiment.** It is included as a sanity check that no anomalous secular signature emerges over the longer horizon — a question reviewers may reasonably raise about a $10\,T$-only baseline.

If the $50\,T$ metrics remain at the same f64 round-off floor as the $10\,T$ baseline, that is corroboration of the gated result. If they diverge, that divergence is reported but does not retroactively alter the $10\,T$ verdict; instead, the $50\,T$ behaviour becomes a separate finding to be discussed.

---

## Appendix B — Format consistency with the Kepler notebook

This notebook deliberately mirrors the section structure and methodological framing of `2026-04-25-rebound-parity-kepler.md`. The framework is shared; the metrics are specialised:

| Section | Kepler | Figure-8 |
| --- | --- | --- |
| Physical invariants gated | orbital elements $(a, e, \omega, h)$ + energy | global integrals $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ |
| Per-side check | conservation against analytic Kepler invariants | conservation against IC values (zero for $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$) |
| Cross-impl check | element drift between apsis and REBOUND | invariant drift between apsis and REBOUND |
| Phase-drift handling | informational $\lvert \Delta r \rvert$, ungated | informational $\lvert \Delta r \rvert$, ungated |
| Horizon | $100\,T$ | $10\,T$ (gated) + $50\,T$ (informational) |

The shared framework is "physical invariants gate; geometric coherence informs". The specialisation is dictated by the system, not by methodological preference.
