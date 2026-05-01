# REBOUND Parity — Pythagorean Three-Body (Burrau 1913)

**Date:** 2026-04-30
**Subject:** Numerical agreement between IAS15 (apsis) and IAS15 (REBOUND) on the canonical Pythagorean three-body problem (Burrau 1913; Szebehely & Peters 1967) — three masses 3, 4, 5 at the vertices of a 3-4-5 right triangle, released from rest, integrated through the chaotic close-encounter regime to ejection.
**Baseline commit:** *(to be pinned at run time)*
**Tooling:** apsis IAS15 (`crates/apsis/src/physics/integrator/ias15.rs`), REBOUND 4.6.0 via Python 3.10 (`reb.IAS15`).
**Status:** *Protocol declared a priori; no run executed at the time of writing. Results section to be populated once the harness is invoked on the pinned commit.*

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

The Pythagorean problem is the simplest published scenario that satisfies all three. Its dynamics are the canonical chaotic test case for adaptive integrators (Szebehely & Peters 1967; Aarseth & Lecar 1971; Hairer, Nørsett & Wanner 1993 §II.10). It has integer-valued initial conditions in canonical units ($G = 1$), so the bit-pattern of the ICs is identical between Rust and Python implementations on x86-64 with IEEE-754 default rounding — eliminating IC-construction divergence as a confound. It admits no closed-form solution at any horizon, so the parity check must rest on conserved quantities rather than on agreement with a reference trajectory.

The framing remains validation, not competition. Tolerances are set by the precision the physics admits — for global integrals of motion under IAS15 through repeated close encounters at f64 precision, the round-off floor is set by the integrator's published machine-precision conservation property (Rein & Spiegel 2015) modulated by the cumulative substep count over the chaotic regime. The 1PN/Mercury demonstration in the v0.1 paper rests on extension-mechanism evidence, not on a numerical-superiority claim; this experiment supports the numerical foundation underlying that mechanism for the chaotic three-body regime, completing the (Kepler, figure-8, Pythagorean) trio of canonical validation scenarios for an N-body solver.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the Pythagorean three-body problem integrated under IAS15 in both `apsis` and REBOUND with the canonical Burrau initial conditions and matched integrator settings, the metrics declared below are bounded *a priori* at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are both gated; failure of any Tier 1 or Tier 2 metric reproves the experiment. Tier 3 is informational and never reproves — its purpose is methodological transparency about controller behaviour in the chaotic regime, not parity gating. Per-body position drift between two correct adaptive IAS15 implementations is *expected* to be $O(1)$ at this horizon (Lyapunov amplification of ULP-level controller decisions); reading that drift as a parity defect would be a category error, on the same methodological grounds developed in the Kepler notebook (§Pilot Interpretation, `2026-04-25-rebound-parity-kepler.md`).

#### Tier 1 — Hard physical invariants *(gated)*

The quantities here are exact constants of motion of the three-body dynamics under Newtonian gravity; any drift in either implementation is a numerical artefact. Both per-side conservation and cross-implementation agreement are gated.

- **Energy** — $\max |\Delta E / E_0|$ per side $\leq 1.0 \times 10^{-13}$ ($\approx 50 \times$ f64 machine epsilon). IAS15 is designed for machine-precision energy conservation across many substeps (Rein & Spiegel 2015, §4), and the published Pythagorean conservation results cited in their §4.2 sit at this level. The bound matches the figure-8 notebook's energy gate.
- **Energy** — cross-implementation $\max |\Delta E_\text{apsis} - \Delta E_\text{rebound}| / |E_0| \leq 1.0 \times 10^{-13}$.
- **Angular momentum** — $\max |\Delta \mathbf{L}|$ per side $\leq 1.0 \times 10^{-13}$ in absolute units. The Pythagorean ICs have $\mathbf{L}_0 = \mathbf{0}$ by construction (all velocities zero, no rotation built into the IC), so a relative bound is undefined; the absolute bound is set against the characteristic angular-momentum scale $|m_i \, \mathbf{r}_i \times \mathbf{v}_i|$ that develops once the bodies acquire velocity from gravitational attraction (peak $\sim O(10)$ during the chaotic regime), with $50 \times$ ULP margin against that peak.
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

The initial `dt = 10^{-3}` matches the apsis preset's `suggested_dt` and is the conventional seed used in the literature (Szebehely & Peters 1967; Aarseth 2003 §3). Both controllers are then free to grow `dt` in the smooth quiescent intervals between encounters and shrink at encounters; the gated metrics characterise the cumulative round-off accumulated through this adaptation.

REBOUND's `exact_finish_time = 1` is used so REBOUND samples are evaluated at the *actual* sample times produced by the apsis side after its adaptive controller has overshot the nominal target; this removes "two implementations sampled at slightly different physical times" as a confound in the cross-implementation comparison. The same convention was used in the Kepler and figure-8 notebooks.

#### Run parameters and sampling

- **Total integration horizon (gated):** $t_\text{final} = 70$ canonical time units. **No event-based termination is used; integration proceeds to fixed horizon $t = 70$ to avoid ambiguity across definitions of "completion" in the literature** (Aarseth 2003 §3 cites $t \approx 46$ for the Burrau ICs specifically; Szebehely & Peters 1967 Fig. 5 extends through $t \approx 60$). Both definitions are operationally distinct (energy-based vs. separation-based vs. regime-change) and would produce different stopping points; fixing the horizon by clock time eliminates the ambiguity from the parity claim and lets the comparison rest on the same integration window on both sides regardless of which "ejection" definition a future reviewer applies. The choice of $70$ exceeds both literature reference points and includes the post-ejection regime in which the dynamics simplify to a tight binary plus a near-hyperbolic body, without entering the regime where the binary's rapid orbit dominates substep selection and inflates wall-time cost without adding parity evidence.
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
- **Sensitivity to IAS15 `epsilon`.** The default tolerance on both sides is used. A sweep over `epsilon` would characterise the cost-precision frontier and is reserved for the Phase 6A characterisation experiment.
- **Behaviour past ejection** beyond $t = 70$. Once the binary has formed and a body has departed, the dynamics admit a closed-form description (a Keplerian binary plus a free particle) and the integration reduces to two-body integration on the binary; further parity evidence in that regime is duplicative of the Kepler notebook.

---

## Results

*To be populated post-run. The §Results section of the figure-8 notebook is the format template: tables for Tier-1 and Tier-2 gated metrics with observed value, declared a priori tolerance, and margin; a separate table for Tier-3 informational per-body $|\Delta \mathbf{r}|$ peaks; a final interpretation paragraph framing the result as the third entry in Pillar A of the v0.1 validation portfolio.*

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
| apsis canonical commit | *(to be pinned at run time)* |
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
