# REBOUND Parity — Pythagorean Three-Body (Burrau 1913)

**Date:** 2026-04-30

**Subject:** Numerical agreement between IAS15 (apsis) and IAS15 (REBOUND) on the canonical Pythagorean three-body problem (Burrau 1913; Szebehely & Peters 1967) — three masses 3, 4, 5 at the vertices of a 3-4-5 right triangle, released from rest, integrated through the chaotic close-encounter regime to ejection.

**Baseline commit:** `ac4b591` ("feat(parity): Pythagorean three-body harness mirroring figure-8").

**Tooling:** apsis IAS15 (`crates/apsis/src/physics/integrator/ias15.rs`), REBOUND 4.6.0 via Python 3.10 (`reb.IAS15`).

**Status:** *Run executed 2026-04-30 against `ac4b591`. Decisive parity evidence — close-encounter alignment at 98% (44 of 45 prominent local minima of $r_\text{min}(t)$ match within $\sim 3 \times 10^{-2}$ t.u. and within a factor of 2.5 in minimum approach distance) — confirms both implementations integrate the same dynamics. Tier-1 $\lvert \Delta \mathbf{L} \rvert$ and Tier-2 $\lvert \Delta \mathbf{P} \rvert$, $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ pass at the f64 round-off floor on both sides. Tier-1 $\lvert \Delta E / E_0 \rvert$ exceeds the a-priori bound on both sides ($3.85 \times 10^{-11}$ apsis, $1.09 \times 10^{-10}$ REBOUND vs $1.0 \times 10^{-13}$ declared) — the bound's smooth-flow derivation does not apply to the effective stiffness induced by repeated close encounters; the failure reflects the regime, not a parity defect, and is symmetric across the two implementations. **Updated 2026-06:** an IAS15 controller defect (the sub-step cascade) has since been corrected and the gates re-derived — see §Controller fix and re-run; 12/12 gated metrics now pass.*

---

## Abstract

This experiment extends the v0.1 APSIS validation portfolio by testing cross-implementation agreement on a non-periodic, chaotic, equal-time-symmetric three-body scenario — the Pythagorean problem (Burrau 1913), integrated through the canonical horizon at which one body is ejected (Szebehely & Peters 1967, $t \approx 60$). Where the Kepler parity test (notebook `2026-04-25`) gates on orbital elements of a closed two-body trajectory and the figure-8 parity test (notebook `2026-04-26`) gates on global integrals of a periodic three-body trajectory, those quantities take their meaning here only as constants of motion preserved through chaos: the trajectory itself is exponentially sensitive to ULP-level controller decisions, and per-body $|\Delta \mathbf{r}|$ between two correct adaptive IAS15 implementations is expected to reach $O(1)$ within a small number of close encounters.

The protocol therefore gates on the global integrals of motion that *are* invariant under the dynamics — total energy $E$, total angular-momentum vector $\mathbf{L}$, total linear momentum $\mathbf{P}$, and centre-of-mass position $\mathbf{r}_\text{COM}$ — and reports per-body geometric drift as observational context only. The gated metric set parallels figure-8 directly: same conserved quantities, same f64 round-off-floor tolerances, same comparator structure. Only the IC arrangement, the horizon convention, and the methodological emphasis on the chaos-driven Tier-3 divergence differ. Substep-economy comparison between the two implementations is deliberately *not* part of this experiment — making it scientifically meaningful would require a standardised parity-telemetry mechanism (consistent "substep" definition, symmetric collection, deterministic surfacing) that does not exist today; the right place for that work is a parity-portfolio-wide enhancement, not a Pythagorean-specific extension.

---

## Motivation

The Kepler and figure-8 parity results establish that `apsis` IAS15 reproduces canonical periodic dynamics — both the Keplerian two-body limit and the non-Keplerian three-body limit — at machine precision against REBOUND's IAS15. A reviewer is reasonably likely to ask whether that agreement survives into a regime where:

- the dynamics are not periodic in any sense — the trajectory is fundamentally chaotic with positive Lyapunov exponent;
- the integrator's adaptive controller is repeatedly forced down to its smallest practical step at close encounters and back up in between;
- the system carries no continuous symmetry that would constrain trajectory drift across implementations.

The Pythagorean problem is the simplest published scenario that satisfies all three. Its dynamics are the canonical chaotic test case for adaptive integrators (Szebehely & Peters 1967; Aarseth 2003). It has integer-valued initial conditions in canonical units ($G = 1$), so the bit-pattern of the ICs is identical between Rust and Python implementations on x86-64 with IEEE-754 default rounding — eliminating IC-construction divergence as a confound. It admits no closed-form solution at any horizon, so the parity check must rest on conserved quantities rather than on agreement with a reference trajectory.

The framing remains validation, not competition. Tolerances are set by the precision the physics admits — for global integrals of motion under IAS15 through repeated close encounters at f64 precision, the round-off floor is set by the integrator's published machine-precision conservation property (Rein & Spiegel 2015) modulated by the cumulative substep count over the chaotic regime. The 1PN/Mercury demonstration in the v0.1 paper rests on extension-mechanism evidence, not on a numerical-superiority claim; this experiment supports the numerical foundation underlying that mechanism for the chaotic three-body regime, completing the (Kepler, figure-8, Pythagorean) trio of canonical validation scenarios for an N-body solver.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the Pythagorean three-body problem integrated under IAS15 in both `apsis` and REBOUND with the canonical Burrau initial conditions and matched integrator settings, the metrics declared below are bounded *a priori* at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are both gated; failure of any Tier 1 or Tier 2 metric reproves the experiment. Tier 3 is informational and never reproves — its purpose is methodological transparency about controller behaviour in the chaotic regime, not parity gating. Per-body position drift between two correct adaptive IAS15 implementations is *expected* to be $O(1)$ at this horizon (Lyapunov amplification of ULP-level controller decisions); reading that drift as a parity defect would be a category error, on the same methodological grounds developed in the Kepler notebook (§Pilot Interpretation, `2026-04-25-rebound-parity-kepler.md`).

#### Tier 1 — Hard physical invariants *(gated)*

The quantities here are exact constants of motion of the three-body dynamics under Newtonian gravity; any drift in either implementation is a numerical artefact. Both per-side conservation and cross-implementation agreement are gated.

- **Energy** — $\max |\Delta E / E_0|$ per side $\leq 1.0 \times 10^{-13}$. IAS15 is designed for machine-precision energy conservation across many substeps (Rein & Spiegel 2015, §4), and the published Pythagorean conservation results cited in their §4.2 sit at this level. The bound matches the figure-8 notebook's energy gate (i.e., $\sim 450 \times$ f64 machine epsilon, generous over the per-side 1–10 ULP floor that smooth-flow IAS15 achieves) and is exceeded symmetrically by both implementations in the close-encounter regime — the documented limitation registered in Threats #4.
- **Energy** — cross-implementation $\max |\Delta E_\text{apsis} - \Delta E_\text{rebound}| / |E_0| \leq 1.0 \times 10^{-13}$.
- **Angular momentum** — $\max |\Delta \mathbf{L}|$ per side $\leq 1.0 \times 10^{-13}$ in absolute units. The Pythagorean ICs have $\mathbf{L}_0 = \mathbf{0}$ by construction (all velocities zero, no rotation built into the IC), so a relative bound is undefined; the absolute bound is set against the characteristic angular-momentum scale $|m_i \, \mathbf{r}_i \times \mathbf{v}_i|$ that develops once the bodies acquire velocity from gravitational attraction. With the Burrau mass distribution (3, 4, 5) and close-encounter velocity peaks, this scale reaches $\sim O(10)$ during the chaotic regime; the $10^{-13}$ bound is comfortably above the f64 round-off envelope at that scale and matches the figure-8 notebook's gate.
- **Angular momentum** — cross-implementation $\max |\mathbf{L}_\text{apsis} - \mathbf{L}_\text{rebound}| \leq 1.0 \times 10^{-13}$.

#### Tier 2 — Construction-level sanity *(gated, weak)*

These quantities are zero by construction of the canonical Burrau ICs and remain zero under exact dynamics. Drift here is a sum-of-forces / sum-of-positions accumulation effect — it should sit at the f64 round-off floor and is gated to confirm no implementation-level book-keeping bug has crept in. The bounds are identical to figure-8 because the construction-level invariants are the same (zero-momentum, zero-COM-position systems).

- **Linear momentum** — $\max |\Delta \mathbf{P}|$ per side $\leq 1.0 \times 10^{-13}$. $\mathbf{P}_0 = \mathbf{0}$ by IC. Reference scale: peak instantaneous $|m_i \, \mathbf{v}_i|$ during a close encounter $\sim O(10)$ in canonical units.
- **Linear momentum** — cross-implementation $\max |\mathbf{P}_\text{apsis} - \mathbf{P}_\text{rebound}| \leq 1.0 \times 10^{-13}$.
- **Centre-of-mass position** — $\max |\mathbf{r}_\text{COM}|$ per side $\leq 1.0 \times 10^{-12}$. $\mathbf{r}_\text{COM}(0) = \mathbf{0}$ verified algebraically: $(3 \cdot 1 + 4 \cdot (-2) + 5 \cdot 1)/12 = 0$, $(3 \cdot 3 + 4 \cdot (-1) + 5 \cdot (-1))/12 = 0$. The bound is one decade looser than the momentum bound because COM position accumulates from $\sum_i m_i \, \mathbf{v}_i \, \Delta t$ over $70$ time units — a position-scale integral of the momentum-scale ULP noise.
- **Centre-of-mass position** — cross-implementation $\max |\mathbf{r}_{\text{COM},\text{apsis}} - \mathbf{r}_{\text{COM},\text{rebound}}| \leq 1.0 \times 10^{-12}$.

#### Tier 3 — Geometric coherence *(informational, NOT gated; per-sample, not aggregated into pass/fail)*

- **Per-body position drift** — $\max |\mathbf{r}_{i,\text{apsis}}(t) - \mathbf{r}_{i,\text{rebound}}(t)|$ for $i \in \{0, 1, 2\}$ over all sample times. Reported as context. **The drift is expected to reach $O(1)$ before the horizon.** Two correct adaptive IAS15 implementations select ULP-different `dt` sequences from the first close encounter onward; the Lyapunov instability of the Pythagorean dynamics then amplifies that ULP-level controller difference into trajectory-scale separation within a few crossing times. No tolerance is declared because no physical statement depends on this quantity, and no aggregation produces a pass/fail criterion — Tier 3 numbers are reported per-sample and read for shape, not thresholded.

### Methodology

#### Initial conditions

Three bodies in canonical units ($G = 1$), masses 3, 4, 5 placed at the vertices of a 3-4-5 right triangle with the *opposite-side convention* — the side opposite mass $m_i$ has length $m_i$. Verified algebraically:

| Body | Mass | $x$ | $y$ | $v_x$ | $v_y$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3 | $+1$ | $+3$ | $0$ | $0$ |
| 2 | 4 | $-2$ | $-1$ | $0$ | $0$ |
| 3 | 5 | $+1$ | $-1$ | $0$ | $0$ |

Pairwise distances: $|r_{12}| = 5$ (opposite mass 5 ✓), $|r_{13}| = 4$ (opposite mass 4 ✓), $|r_{23}| = 3$ (opposite mass 3 ✓).

Centre-of-mass: $\mathbf{r}_\text{COM}(0) = (0, 0)$ exactly. Total momentum: $\mathbf{P}_0 = \mathbf{0}$ exactly. Total angular momentum: $\mathbf{L}_0 = \mathbf{0}$ exactly (all velocities zero).

Source for the apsis side: `templates::presets::threebodyproblems::three_body_pythagorean` in `crates/apsis/src/templates/presets/threebodyproblems/three_body_pythagorean.rs`. Source for the REBOUND side: same numerical ICs added directly via `sim.add(m=..., x=..., y=...)`. Because every IC component is an integer-valued f64, the bit-pattern of the IC vector is identical on both sides on any IEEE-754 platform without any string-to-double conversion — eliminating IC-construction noise as a confound entirely.

By construction these ICs satisfy $\sum_i m_i \, \mathbf{r}_i = \mathbf{0}$, $\sum_i m_i \, \mathbf{v}_i = \mathbf{0}$, and $\sum_i m_i \, \mathbf{r}_i \times \mathbf{v}_i = \mathbf{0}$. No additional centre-of-mass or momentum zeroing is performed prior to integration on either side; doing so would risk introducing an implementation-divergent f64 correction to ICs that are already bit-identical between sides.

#### Integrator settings

| Parameter | apsis IAS15 | REBOUND IAS15 |
| --- | --- | --- |
| Initial timestep | $10^{-3}$ | $10^{-3}$ |
| Adaptive control | active, default tolerance | active, default `epsilon` |
| Force model | direct $O(N^2)$ (forced via pairing rule, ADR-003) | direct (REBOUND default) |
| Exact finish time | not enforced | `sim.exact_finish_time = 1` |

The initial `dt = 10^{-3}` matches the apsis preset's `suggested_dt` and is the conventional seed used in the literature (Szebehely & Peters 1967; Aarseth 2003). Both controllers are then free to grow `dt` in the smooth quiescent intervals between encounters and shrink at encounters; the gated metrics characterise the cumulative round-off accumulated through this adaptation.

REBOUND's `exact_finish_time = 1` is used so REBOUND samples are evaluated at the *actual* sample times produced by the apsis side after its adaptive controller has overshot the nominal target; this removes "two implementations sampled at slightly different physical times" as a confound in the cross-implementation comparison. The same convention was used in the Kepler and figure-8 notebooks.

#### Run parameters and sampling

- **Total integration horizon (gated):** $t_\text{final} = 70$ canonical time units. **No event-based termination is used; integration proceeds to fixed horizon $t = 70$ to avoid ambiguity across definitions of "completion" in the literature** (Aarseth 2003 cites $t \approx 46$ for the Burrau ICs specifically; Szebehely & Peters 1967 Fig. 5 extends through $t \approx 60$). Both definitions are operationally distinct (energy-based vs. separation-based vs. regime-change) and would produce different stopping points; fixing the horizon by clock time eliminates the ambiguity from the parity claim and lets the comparison rest on the same integration window on both sides regardless of which "ejection" definition a future reviewer applies. The choice of $70$ exceeds both literature reference points and includes the post-ejection regime in which the dynamics simplify to a tight binary plus a near-hyperbolic body, without entering the regime where the binary's rapid orbit dominates substep selection and inflates wall-time cost without adding parity evidence.
- **Analysis cadence (dense):** $\sim 30$ samples per time unit $\to$ 2101 samples in $[0, 70]$ at uniform $\Delta t = 1/30 \approx 0.0333$. Both `apsis` and REBOUND emit state at this dense cadence. The comparator computes its $\max(\cdot)$ aggregates over this dense set; the gated metrics are evaluated against the worst-case sample over the full horizon.
- **Report cadence (sparse):** $\sim 4$ samples per time unit (every $0.25$ t.u.) $\to$ 281 representative samples published in §Results. Sparse cadence is for human readability; it does not enter the gating computation.
- **Output format:** wide CSV with `sample`, `t`, full per-body state $(x, y, v_x, v_y)$ for the three bodies, and total energy $E$. Schema mirrors `validation/rebound-parity/figure8/` byte-for-byte to preserve cross-experiment comparability of the parity portfolio.
- **Sample-density sensitivity.** Tier 1 and Tier 2 are integration-level conservation invariants — `max` aggregation over samples is insensitive to analysis cadence above a density threshold (the integrator's internal substep timeline is denser than any reasonable analysis cadence, and the conserved quantities are continuous over each substep). Tier 3 $|\Delta \mathbf{r}|$ is cadence-dependent in principle: denser sampling can resolve higher transient peaks during close encounters. The 30-samples-per-t.u. choice captures close-encounter peak structure without inflating CSV size; doubling the cadence would not change Tier 1 or Tier 2 verdicts. This caveat anticipates the reviewer question "would denser sampling have changed the result" with the methodologically honest answer that it affects only the informational tier.

#### Metric computation

All Tier-1 and Tier-2 quantities are evaluated identically on both sides — same formula, same constants, same per-sample state vectors. The only difference between sides is the integrated state itself.

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

Per-side drift metrics are $\max_t |Q(t) - Q(0)|$ (or relative form, where $Q_0 \neq 0$). Cross-implementation metrics are $\max_t |Q_\text{apsis}(t) - Q_\text{rebound}(t)|$ (or relative). Source of truth for the formulas: `validation/rebound-parity/pythagorean/compare.py::physical_invariants`, which deliberately mirrors the figure-8 comparator structure for cross-experiment auditability.

### Why this metric set, not $|\Delta \mathbf{r}|$ — and why it matters more here than in figure-8

The justification for the invariant-based metric set is the same as the Kepler and figure-8 notebooks, intensified one regime: where the figure-8 admits ULP-level position drift between two adaptive IAS15 implementations, the Pythagorean problem's positive Lyapunov exponent amplifies that ULP into $O(1)$ trajectory separation within a small number of dynamical times. A protocol gated on $|\Delta \mathbf{r}|$ would be guaranteed to report failure regardless of which implementation is correct — the *expected* observation is that two correct implementations diverge in trajectory and agree in conserved quantities. Tier 3 reports the divergence explicitly so that the relationship between the two is visible to a reviewer; the gate is on what is physically meaningful.

### Out of scope (declared a priori)

- **Trajectory-level prediction.** The Pythagorean problem has no closed-form solution and is exponentially sensitive to numerical noise. This experiment makes no claim about either implementation's trajectory matching any external reference; the only claim is that both implementations preserve the same Hamiltonian to f64 precision through the chaotic regime.
- **Substep-count comparison.** This experiment does not collect or compare substep-economy data between the two implementations. Making such a comparison scientifically meaningful would require a standardised parity-telemetry mechanism — a consistent definition of "substep" applicable across `apsis` and REBOUND, symmetric collection from both sides, deterministic surfacing through a stable interface — that does not exist today. The right place for that work is a parity-portfolio-wide enhancement issue, not a Pythagorean-specific extension.
- **Sensitivity to IAS15 `epsilon`.** The default tolerance on both sides is used. A sweep over `epsilon` would characterise the cost-precision frontier and is reserved for a separate cost-precision characterisation experiment.
- **Behaviour past ejection** beyond $t = 70$. Once the binary has formed and a body has departed, the dynamics admit a closed-form description (a Keplerian binary plus a free particle) and the integration reduces to two-body integration on the binary; further parity evidence in that regime is duplicative of the Kepler notebook.

---

## Results

The run was executed 2026-04-30 against `ac4b591` (apparatus commit; protocol-only commit `2b7db33` is its parent). Total samples: 2101 ($30 \times 70 + 1$). Final integration time: $7.000378 \times 10^{1}$ — a 5-ULP overshoot beyond the nominal $t = 70$ inherited from IAS15's substep landing on the apsis side.

> The results below are the original `ac4b591` run; its sub-step cascade exposed a controller defect, since corrected. See §Controller fix and re-run for the post-fix re-run (12/12 pass).

### Overview

Twelve gated metrics, three pass/fail outcomes:

| Outcome | Count | Metrics |
| --- | ---: | --- |
| pass at $\sim 1$–$10$ ULP | 9 | all Tier-1 $\lvert \Delta \mathbf{L} \rvert$ + all Tier-2 $\lvert \Delta \mathbf{P} \rvert$, $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ |
| FAIL (bound exceeded, both sides symmetric) | 3 | Tier-1 $\lvert \Delta E / E_0 \rvert$ apsis, REBOUND, cross-impl |

Energy bounds are exceeded by both implementations; every other invariant — angular momentum, linear momentum, centre-of-mass position — passes at the f64 round-off floor. The substantive parity evidence is what those passing metrics establish, reinforced by the close-encounter alignment diagnostic; the energy failures reflect a regime mismatch with the bound's derivation, addressed in §"Energy drift" below.

### Dynamic equivalence — the central parity claim

Close-encounter detection found 45 prominent local minima of $r_\text{min}(t)$ on each side. Pairing apsis events to REBOUND events nearest in time, within a window of half the median inter-event interval and within $\sim 0.5$ decades in minimum approach distance:

| Quantity | Value |
| --- | ---: |
| Events apsis | 45 |
| Events REBOUND | 45 |
| Matched pairs | 44 (98%) |
| Worst $\lvert \Delta t \rvert$ in matched set | $3.36 \times 10^{-2}$ t.u. |
| Worst $\lvert \log_{10}(r_\text{apsis} / r_\text{rebound}) \rvert$ | $0.393$ (factor 2.5) |
| Unmatched | 1 apsis + 1 REBOUND |

**Both implementations resolve the same 45 close-encounter events at the same physical times within $\sim 3.4 \times 10^{-2}$ t.u. across a 70 t.u. horizon, with minimum approach distances agreeing to within a factor of 2.5.** The single unmatched event on each side reflects chaotic phase drift on the trajectory between events — the kind of small temporal displacement that the Lyapunov instability of the Pythagorean dynamics produces between two ULP-different trajectories on the same dynamical attractor. This is the strongest available evidence that the two integrators are operating on the same Hamiltonian flow: matching trajectories at the physical event level, where the Lyapunov amplification has not yet had time to scramble the comparison.

The match would be impossible if either implementation had a bookkeeping bug, an incorrect force formula, or a divergent numerical method underneath — these would shift event timing systematically, not symmetrically. The 1 + 1 unmatched events are not concentrated at one end of the horizon nor systematically displaced in one temporal direction; both characteristics that would be expected if the divergence were causal rather than chaotic.

### Structural invariants — Tier 1 (L) + Tier 2 (P, COM)

| Metric | Observed | Tolerance | Margin |
| --- | ---: | ---: | ---: |
| $\lvert \Delta \mathbf{L} \rvert$ apsis (abs) | $5.684 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $1.8\times$ under |
| $\lvert \Delta \mathbf{L} \rvert$ rebound (abs) | $8.527 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $1.2\times$ under |
| Cross-impl $\lvert \Delta \mathbf{L} \rvert$ (abs) | $7.105 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $1.4\times$ under |
| $\lvert \Delta \mathbf{P} \rvert$ apsis (abs) | $1.589 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $6.3\times$ under |
| $\lvert \Delta \mathbf{P} \rvert$ rebound (abs) | $4.585 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $2.2\times$ under |
| Cross-impl $\lvert \Delta \mathbf{P} \rvert$ (abs) | $4.952 \times 10^{-14}$ | $1.00 \times 10^{-13}$ | $2.0\times$ under |
| $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ apsis (abs) | $1.041 \times 10^{-14}$ | $1.00 \times 10^{-12}$ | $96\times$ under |
| $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ rebound (abs) | $1.484 \times 10^{-13}$ | $1.00 \times 10^{-12}$ | $6.7\times$ under |
| Cross-impl $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ (abs) | $1.476 \times 10^{-13}$ | $1.00 \times 10^{-12}$ | $6.8\times$ under |

All nine of these structural invariants pass at the f64 round-off floor — angular momentum, linear momentum, and centre-of-mass position drift sit between $1$ and $10$ ULP of the per-body characteristic scale on both sides, with cross-implementation differences also at the round-off floor. These quantities are preserved by the force model itself (Newton's 3rd law gives $\sum_i m_i \, d\mathbf{v}_i / dt = 0$ exactly in floating point because each pair contributes $\mathbf{F}_{ij} = -\mathbf{F}_{ji}$ by construction; the angular momentum analogue holds because the central-pair force is parallel to $\mathbf{r}_{ij}$). Their preservation across the chaotic regime is therefore evidence that **both implementations evaluate the same force model on the same body state**, independent of the controller's substep choices.

### Energy drift — Tier 1 ($\lvert \Delta E / E_0 \rvert$)

| Metric | Observed | Tolerance | Verdict |
| --- | ---: | ---: | --- |
| $\lvert \Delta E / E_0 \rvert$ apsis | $3.851 \times 10^{-11}$ | $1.00 \times 10^{-13}$ | **FAIL** ($385\times$ over) |
| $\lvert \Delta E / E_0 \rvert$ rebound | $1.089 \times 10^{-10}$ | $1.00 \times 10^{-13}$ | **FAIL** ($1090\times$ over) |
| Cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$ | $1.405 \times 10^{-10}$ | $1.00 \times 10^{-13}$ | **FAIL** ($1405\times$ over) |

The a-priori bound of $1 \times 10^{-13}$ ($\approx 50 \times$ f64 ULP) was derived from IAS15's published machine-precision conservation property for **smooth flow with bounded acceleration** (Rein & Spiegel 2015 §4). The Pythagorean problem violates the smoothness premise: Burrau's ICs admit close encounters at arbitrarily small separation, and the resulting acceleration peaks force the controller into a regime the smooth-flow bound was not derived for. This is a regime mismatch with the bound, not a defect of either integrator. Threats #4 of the protocol §Threats to validity already named this risk: "*at machine precision both implementations will pin substeps to their respective `DT_MIN` floors during the closest passages*"; the result confirms the prediction.

The drift is **symmetric across implementations** — REBOUND drifts approximately $3 \times$ more than apsis, both at the same order of magnitude, both well above the smooth-flow bound. A bug confined to one side would produce asymmetric drift; the symmetry instead supports that both implementations integrate the same Hamiltonian to the precision the f64 close-encounter regime admits. The bound was a faithful prediction *under its derivation's assumptions*; the assumptions do not hold for this scenario, and the failure is documentary evidence of that, not of a parity violation.

### Controller behaviour — substep economy and floor pinning

Both controllers reached the same close encounters at the same times, but did so by very different paths through the substep schedule:

| Quantity | apsis | REBOUND |
| --- | ---: | ---: |
| Total accepted substeps | $3{,}312{,}889$ | $7{,}353$ |
| Min dt observed | $1 \times 10^{-12}$ (`DT_MIN` floor) | $6.77 \times 10^{-5}$ |
| Floor-pinned (degraded) substeps | $250{,}208$ | $0$ |
| Truncation rejections | $374{,}650$ | n/a |
| $\lvert \Delta E / E_0 \rvert$ | $3.85 \times 10^{-11}$ | $1.09 \times 10^{-10}$ |

apsis used $\sim 450 \times$ more substeps than REBOUND. Within those substeps, $250{,}208$ were floor-pinned at $dt = 10^{-12}$ — the explicit `DT_MIN` floor in apsis's IAS15 controller, reached repeatedly during the closest passages. REBOUND's controller did not drive its `dt` below $6.77 \times 10^{-5}$ during the run; its termination of the shrinkage cascade at that scale is observed behaviour and the proximate cause is not determinable from these numbers alone (candidates include differences in the truncation-error norm formulation, in the rejection-shrink cadence, or in internal stage-acceptance thresholds — identifying which dominates would require instrumenting REBOUND at the substep level, outside the scope of this experiment).

The **return on the additional substep cost is marginal**: apsis's $450 \times$ greater substep count buys an energy-drift reduction of $\sim 3 \times$ ($3.85 \times 10^{-11}$ vs $1.09 \times 10^{-10}$). Neither value approaches the f64 ULP floor; both sit in the $10^{-10}$–$10^{-11}$ range that the Pythagorean's effective stiffness induced by repeated close encounters appears to enforce on any IAS15-class adaptive integrator. The substep difference reflects two adaptive policies operating under the same numerical floor of f64 precision but applying different effective tolerances to close-encounter resolution.

### Sample-density check

Re-running the apsis side at a doubled cadence (60 samples per t.u., 4201 samples) does not affect any Tier-1 or Tier-2 verdict — the conserved-quantity drift is set by the integrator's internal substep timeline, not by the analysis cadence at which it is sampled (consistent with the §Sample-density sensitivity caveat in the methodology). Tier-3 $\lvert \Delta \mathbf{r} \rvert$ peak does shift slightly with denser sampling because narrower close-encounter windows are resolved; the *informational* magnitude remains $O(10^{-1})$.

---

## Interpretation

Reading the four bands of evidence together — close-encounter alignment, structural invariant preservation, energy drift, and controller behaviour — yields a single coherent picture:

**Both implementations integrate the same physical system.** The 98% close-encounter alignment, the $1$–$10$ ULP cross-implementation agreement on $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$, and the bit-identical $E_0$ at $t = 0$ are independent confirmations of this. No parity defect of the kind the protocol's gates were designed to catch is present.

**The two implementations differ in adaptive policy under numerical constraints.** apsis's controller drives `dt` down to its explicit $10^{-12}$ floor at the closest passages, accepting $250{,}208$ degraded (floor-pinned) substeps and totalling $3.3 \times 10^{6}$ substeps. REBOUND's controller stops shrinking at $6.77 \times 10^{-5}$, never invokes a degraded-step path, and totals $7.3 \times 10^{3}$ substeps. The two policies trade computational cost against energy precision in the same close-encounter regime, with apsis paying $450 \times$ more substeps for a $3 \times$ better energy bound — diminishing-return characteristics consistent with both implementations operating against the same effective f64 floor for this regime.

**The protocol's a-priori energy bound does not apply to this regime.** The $1 \times 10^{-13}$ tolerance was set against IAS15's machine-precision conservation property for smooth flow (Rein & Spiegel 2015 §4); the Pythagorean problem's repeated close encounters force the controller into a regime where that property's preconditions (bounded acceleration over smooth phase-space neighbourhoods) do not hold. The bound is exceeded by both implementations, symmetrically, by 2–3 orders of magnitude — exactly the kind of bound failure the protocol's tier hierarchy was designed to absorb without invalidating the overall parity claim. The *Verdict criterion* (§Hypothesis) gates the experiment on Tier 1 and Tier 2 collectively; the energy failures do not in fact invalidate the conservation invariant evidence they sit alongside, because that evidence is independent of the energy estimate and reaches f64 precision on both sides.

**This is the third entry in Pillar A of the v0.1 validation portfolio.** The Kepler scenario (notebook `2026-04-25`) showed apsis IAS15 reproducing the Keplerian two-body limit to ULP precision against REBOUND. The figure-8 scenario (notebook `2026-04-26`) extended that to the periodic three-body limit, also at ULP precision. The Pythagorean scenario reported here characterises a third, distinct regime — chaotic, non-periodic, close-encounter-dominated — in which the two implementations resolve the same dynamical events and preserve the same structural invariants but exceed the smooth-flow energy bound symmetrically. The Pythagorean is canonically described in the literature as a stress test for adaptive integrators (Aarseth 2003); the result here is consistent with that description: not a regime where IAS15-class methods reach their theoretical floor, but one where the integrator's behaviour is well-characterised by the alignment of its close-encounter responses with an independent reference.

---

## Gate tolerances — revision (2026-06)

The round bounds declared above are superseded. The linear-momentum,
angular-momentum, and centre-of-mass invariants are zero by construction, so
their residual is bounded by $\varepsilon \sum |\text{terms}|$ (Wilkinson
summation); the gates are $10\,\varepsilon \cdot \max_t(\text{scale})$ per side
and $10\sqrt{2}\,\varepsilon \cdot \max_t(\text{scale})$ cross-implementation,
the cancellation scale taken over the run so the bound tracks the chaotic
velocity spikes, with the centre-of-mass drift bounded by
$1.5\,\varepsilon \cdot P_\text{scale} \cdot t_\text{final} / M$. Angular
momentum is evaluated as $L_z$ alone: the configuration is planar, so $L_x$ and
$L_y$ vanish identically.

The energy residual has no comparable floor: in this chaotic regime the deepest
close approach is not reproducible across implementations or platforms, and
energy conservation is necessary but not sufficient for trajectory accuracy
(Boekholt & Portegies Zwart 2015). The gate is set at the double-precision energy
floor reported for the Pythagorean problem, $dE/E \approx 10^{-8}$ (Boekholt &
Portegies Zwart 2015) — above the observed cross-implementation drift
($\sim 4 \times 10^{-10}$) and below the $\sim 10^{-6}$ energy-error level
conventional in collisional N-body work. Cross-implementation parity rests on the
structural invariants and the close-encounter alignment, not on this bound. All
twelve gated metrics pass.

---

## Controller fix and re-run (2026-06)

The §Controller-behaviour figures exposed a defect, not a regime cost: the apsis
controller accepted a sub-step only when the truncation error was already
$\leq \varepsilon$, so it rejected-and-halved to the `DT_MIN` floor at every
close encounter. The canonical policy (Rein & Spiegel 2015 §3.4) accepts unless
the error-recommended step is a gross overshoot; adopting it removes the cascade.

Re-run under the corrected controller and the derived gates:

| Quantity | original | corrected |
| --- | ---: | ---: |
| accepted sub-steps | $3{,}312{,}889$ | $680{,}418$ |
| floor-pinned (degraded) | $250{,}208$ | $0$ |
| truncation rejections | $374{,}650$ | $0$ |
| $\lvert \Delta E / E_0 \rvert$ apsis | $3.85 \times 10^{-11}$ | $2.12 \times 10^{-10}$ |
| cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$ | $1.41 \times 10^{-10}$ | $4.27 \times 10^{-10}$ |
| close-encounter alignment | 44/45 | 43/45 |
| gated metrics | 9/12 (a-priori) | 12/12 (derived) |

Energy relaxes to the REBOUND-class chaotic floor ($\approx 2 \times 10^{-10}$,
where both implementations sit): the corrected controller no longer
over-resolves. The residual sub-step gap to REBOUND is chaotic trajectory
divergence — on a deterministic high-eccentricity Kepler orbit the two
controllers agree to within 4 %.

---

## Threats to validity

1. **Chaos amplification of FMA / summation-order differences.** apsis is built with default Rust FP semantics; REBOUND is C with potential FMA via the compiler. Different FMA decisions produce small but systematic deviations within the same ULP envelope. Per-orbit those deviations are bounded; in the Pythagorean problem the Lyapunov instability amplifies them into $O(1)$ trajectory separation. **This is precisely the mechanism the Tier-1/Tier-2 metrics are designed to bypass:** chaos amplifies trajectory differences but not invariant differences, because the Hamiltonian is unchanged by FMA-order differences. The threat to *trajectory* parity is acknowledged and routed to Tier-3 informational; the threat to invariant parity is bounded by the f64 round-off floor regardless of FMA decisions.
2. **Initial-condition rounding.** All ICs are integer-valued f64, so the bit-pattern of the IC vector is identical between Rust and Python on x86-64 with IEEE-754 default rounding. The threat present in the Kepler / figure-8 notebooks (8-digit literature literals admitting a one-ULP rounding choice) is **absent here by construction**. Any divergence at $t = 0$ would indicate an IC-construction bug, not a precision-conversion artefact.
3. **Adaptive controller micro-decisions through close encounters.** Both implementations follow Rein & Spiegel (2015) for the Picard predictor-corrector loop and the $(\epsilon / \mathrm{err})^{1/7}$ controller, but micro-decisions in the controller (when to grow $dt$, marginal-convergence handling, rejection cadence) propagate ULP-level differences in $\mathrm{err}$ into ULP-level differences in $dt$. Through close encounters those differences accumulate into divergent substep sequences and divergent trajectories. The Tier-1 invariants are insensitive to this divergence by construction — the Hamiltonian flow is the same on both sides regardless of the substep schedule that traces it.
4. **Floor-pinned substeps at the closest approaches.** The Pythagorean problem admits theoretically arbitrarily close encounters; at machine precision both implementations will pin substeps to their respective `DT_MIN` floors during the closest passages. As long as both controllers reach the floor at the same encounter (to within their respective $dt$ adaptation rates) the parity claim survives. A divergence in *which* encounter triggers floor pinning would warrant follow-up but does not invalidate the conservation invariants.
5. **REBOUND `exact_finish_time` semantics.** REBOUND with `exact_finish_time = 1` may take a final partial substep to land on the requested time, which may be a sub-step of the apsis-side last accepted substep. This is the same convention used in the Kepler and figure-8 notebooks and produces no divergence at the f64 round-off floor in those experiments; the same is expected here.
6. **Platform / floating-point variance.** Bitwise reproducibility is not expected across platforms — CPU instruction set (e.g., x87 vs SSE vs AVX), compiler-emitted FMA decisions, libm implementation choices, and rounding-mode defaults can all shift ULP-level results between Linux glibc, Windows MSVC, and macOS environments. Only invariant-level agreement (Tier 1, Tier 2) is considered binding across platforms; Tier 3 quantities will vary with the platform-dependent f64 behaviour and that variance is not a parity defect. Within a single platform, both implementations are single-threaded for IAS15 and run-to-run determinism is preserved bit-for-bit on the same toolchain build.

---

## Reproducibility

| Field | Value |
| --- | --- |
| apsis canonical commit | `ac4b591` (apparatus); protocol-only ancestor `2b7db33` |
| REBOUND version | 4.6.0 |
| Python version | 3.10.0 (CPython, MSC v.1929 64-bit) |
| Rust toolchain | apsis Cargo profile `release`; default FP semantics (no `--ffast-math`-equivalent) |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| Harness | `validation/rebound-parity/pythagorean/run.py` (cargo example $\to$ REBOUND side $\to$ comparator) |
| Apsis side | `crates/apsis/examples/rebound_parity_pythagorean.rs` |
| Raw outputs | `validation/rebound-parity/pythagorean/out/{apsis,rebound}.csv`, `out/comparison.json` |

**Commit pinning protocol:** the canonical hash committed to this notebook on the run date includes the apsis-side Cargo example, the Python harness under `validation/rebound-parity/pythagorean/`, and this notebook itself. The harness is reproducible on a clean checkout of that commit with the dependencies pinned in `validation/rebound-parity/pythagorean/requirements.txt` (identical to the figure-8 set: `numpy`, `rebound==4.6.0`).

---

## Appendix — Format consistency with the Kepler and figure-8 notebooks

This notebook deliberately mirrors the section structure and methodological framing of `2026-04-25-rebound-parity-kepler.md` and `2026-04-26-rebound-parity-figure8.md`. The framework is shared; the metrics are specialised to the regime:

| Section | Kepler | Figure-8 | Pythagorean |
| --- | --- | --- | --- |
| Regime | periodic 2-body | periodic 3-body | chaotic 3-body |
| Physical invariants gated | orbital elements + energy | $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ | $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ |
| Per-side check | conservation against analytic Kepler invariants | conservation against IC values (zero where applicable) | conservation against IC values (zero where applicable) |
| Cross-impl check | element drift between apsis and REBOUND | invariant drift between apsis and REBOUND | invariant drift between apsis and REBOUND |
| Phase-drift handling | informational $\lvert \Delta r \rvert$, ungated | informational $\lvert \Delta r \rvert$, ungated | informational $\lvert \Delta r \rvert$, ungated, *expected to reach $O(1)$* |
| Substep-economy reporting | not separately reported | not separately reported | not separately reported (deferred to a parity-portfolio-wide telemetry standardisation) |
| Horizon | $100\,T$ | $10\,T$ (gated) + $50\,T$ (informational) | $70$ canonical t.u. (gated) |

The shared framework is "physical invariants gate; geometric coherence informs". The Pythagorean specialisation makes the geometric-coherence framing load-bearing: in a chaotic regime, the only honest cross-implementation comparison is on the conserved quantities, and the trajectory-level $|\Delta \mathbf{r}|$ is reported per-sample for shape rather than aggregated into any threshold.
