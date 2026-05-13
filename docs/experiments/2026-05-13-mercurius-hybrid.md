# Mercurius — close-encounter hybrid integrator — protocol

**Date:** 2026-05-13
**Subject:** First-class implementation of the MERCURIUS hybrid integrator (Rein, Hernandez, Tamayo & Brown 2019, MNRAS 489, 4632–4640): WH symplectic outer step + IAS15 close-field sub-integration with a smooth `K(r/r_cross)` changeover. Federates the existing `WisdomHolman` and `Ias15` integrators into a single algorithm whose far-field cost stays O(N) per step while close encounters get IAS15-grade precision *only on the encountering pairs*.

**Status:** Protocol declared a priori, before any implementation lands. Builds on the design proposal `docs/proposals/close-encounter-hybrid.md` (PR #32) and locks the six §6 open design questions with explicit rationale below.

**Branch:** `feat/mercurius`, branched from `develop`. Independent of the perf series (PRs #78 / #81 / #82 already merged); touches the integrator stack, not the gravity hot path.

---

## Abstract

Apsis today handles close encounters reactively through IAS15's adaptive controller: when the local truncation error grows during a close approach, `dt` shrinks, and the next sub-step runs at higher cost across **every body in the system**. For sparse hierarchical scenarios — Sun + N planets where encounters are rare and pairwise — the per-step cost gets amplified by N when only one pair is encountering. A symplectic far-field + non-symplectic high-precision near-field hybrid (MERCURIUS, Rein et al. 2019) localises the encounter cost to the encountering pair while preserving secular stability over arbitrary horizons elsewhere.

This is the first non-IAS15 integrator that ships as a paper-defining federated artifact: the algorithm itself is composed from existing first-class apsis integrators (WH for the far half, IAS15 for the close half), demonstrating the federated-FPM thesis at the integrator layer rather than the perturbation layer.

---

## Motivation

Three regimes drive the case for a hybrid:

1. **Sparse-encounter planetary** (Sun + ≤ 10 planets, occasional close approaches). IAS15 alone forces every body to pay the encounter cost. Mercurius keeps the far-field at WH cost (one analytical Kepler step per planet) and pays IAS15 only on the encountering pair.

2. **Long-horizon hierarchical with rare encounters** (10⁴ – 10⁶ orbits, planet–planetesimal). IAS15 saturates at the Brouwer-law floor for energy drift over very long horizons; WH preserves secular stability arbitrarily long. Today there is no apsis integrator that handles both regimes within one run.

3. **Federated-FPM demonstration.** The thesis ([[project_thesis_anchor]], [[project_paper_positioning]]) frames apsis as a federation where physics components are first-class. Mercurius federates two existing integrators into a single algorithm — the federation idea applied to the temporal axis instead of the perturbation axis. Worth the implementation cost on positioning grounds alone.

### What this experiment is NOT testing

- **Not a WHFast refactor.** TD-008 is closed (commit `ed12fe1`'s notebook references it); the existing `WisdomHolman` is fit for purpose as the symplectic outer half. A future WHFast upgrade will substitute transparently.
- **Not introducing a new force-evaluation regime.** Direct O(N²) evaluation of the residual close-field forces is appropriate for the planetary N (≤ 10² typically). The Barnes-Hut tree is irrelevant inside Mercurius and is bypassed.
- **Not changing the precision invariant.** Mercurius is a Hamiltonian-split symplectic-class scheme: at K = 1 (no encounters anywhere) it reduces algebraically to the existing WH + Kepler drift; at K = 0 (deep encounter) the close-field IAS15 carries the full force. Both limits are testable algebraic identities.

---

## Algorithm

Following Rein et al. 2019 §2. Hamiltonian decomposition:

$$
H = H_K + H_\text{indirect} + H_\text{far} + H_\text{close}
$$

with

| Term | Form | Operator |
| --- | --- | --- |
| $H_K$ | $\sum_{i \ge 1} (p_i^2 / 2 m_i - G m_0 m_i / |q_i|)$ | Analytical Kepler drift |
| $H_\text{indirect}$ | $(\sum_i p_i)^2 / (2 m_0)$ | Uniform position drift on planets |
| $H_\text{far}$ | $-\sum_{i < j} K(r_{ij}/r_\text{cross}^{(ij)}) \, G m_i m_j / r_{ij}$ | K-weighted planet-planet kick |
| $H_\text{close}$ | $-\sum_{i < j} (1 - K(r_{ij}/r_\text{cross}^{(ij)})) \, G m_i m_j / r_{ij}$ | IAS15 sub-integration |

with $q_i = r_i - r_0$ heliocentric and $K$ the changeover function. Per-pair changeover scale $r_\text{cross}^{(ij)} = \alpha \, R_H^{(ij)}$ from the mutual Hill radius.

Symplectic 2nd-order split (DKD-form):

$$
\Psi_\tau = \Phi_{H_K}^{\tau/2} \; \Phi_{H_\text{indirect}}^{\tau/2} \; \Phi_{H_\text{far}}^{\tau/2} \; \Phi_{H_\text{close}}^{\tau} \; \Phi_{H_\text{far}}^{\tau/2} \; \Phi_{H_\text{indirect}}^{\tau/2} \; \Phi_{H_K}^{\tau/2}
$$

The inner $\Phi_{H_\text{close}}^{\tau}$ is integrated by the embedded IAS15 instance, which adaptively sub-steps within $[0, \tau]$ until the residual force is exhausted. The IAS15 controller's own state persists across Mercurius outer steps so it can learn the close-field truncation scale.

Total composition: 7-stage symplectic split. Outer truncation: $O(\tau^3)$ per step, $O(\tau^2)$ over the integration. Identical to standard WH order when $K \equiv 1$ (no encounters).

---

## Locked design decisions *(resolves design proposal §6 a priori)*

| Question | Decision | Rationale |
| --- | --- | --- |
| §6.1 Hill radius $M_\star$ for non-hierarchical systems | Refuse with `MercuriusError::NonHierarchical` at construction time | Mercurius assumes a dominant central body for the analytical Kepler drift. Non-hierarchical fallback would be a different algorithm (different far-field operator). Honest failure → user routes to IAS15 directly. Matches REBOUND. |
| §6.2 Changeover function shape | $K(y) = y^2 (3 - 2y)$ on $y \in [0, 1]$, clipped at endpoints | $C^1$ polynomial. Standard cubic Hermite smoothstep, matches Rein et al. 2019 §2.2 default and REBOUND `mercurius` switching. $C^2$ alternatives (Rein et al. mention but do not adopt) cost more in evaluation without measurable order improvement on the planetary regime. |
| §6.3 Default $\alpha$ (Hill multiplier for $r_\text{cross}$) | $\alpha = 3$ | REBOUND default; validated by Rein et al. 2019 §3 against several planetary scattering scenarios. |
| §6.4 `fast` integrator selection | Wisdom-Holman | Already implemented (3D-native, dense output, hierarchy classification, 619 tests passing). Yoshida-4 would lose the analytical Kepler step that makes Mercurius cheap on the far field. |
| §6.5 Per-pair vs global changeover | Per-pair $r_\text{cross}^{(ij)}$ | Canonical in Rein et al. 2019 and physically necessary — different planet pairs have different $R_H^{(ij)}$ scales. Bookkeeping cost is $O(N^2)$ memory for the table, recomputed per outer step; trivial at planetary $N$. |
| §6.6 Outer step shrink near encounters | Outer $\tau$ stays fixed; inner IAS15 absorbs cost | Canonical Rein et al. 2019 design. The whole point of Mercurius is to *avoid* shrinking the outer step globally; encounter cost stays local to the IAS15 sub-integration. |

Additional decisions (not in §6, locked here):

| Decision | Choice | Rationale |
| --- | --- | --- |
| Public name | `IntegratorKind::Mercurius`, slug `"mercurius"` | Matches Rein et al. 2019 algorithm name and REBOUND. Self-documenting. |
| Force model interaction | Mercurius bypasses `IntegratorContext::force` for its internal pair forces; computes direct O(N²) K-weighted kicks and (1-K)-weighted residuals internally | The Barnes-Hut tree provides no value at planetary $N$ and the per-pair changeover requires explicit pair iteration. The internal close-field IAS15 gets a wrapping `ChangeoverForceModel` so its FSAL / Picard contracts are honoured without changes. |
| `compute_closeness` 2D defect | Fix in this PR (extends to 3D) | Phase 1 EncounterFlag reads `r_min`; latent 2D bug affects Phase 1 correctness on out-of-plane configurations. Single-line fix; in scope. |
| `requires_deterministic_force` | `true` | Inherits from the embedded IAS15. Mercurius itself is deterministic given a deterministic close-field force model (which the internal `ChangeoverForceModel` is, by construction). |
| `controls_own_step_size` | `false` | Outer $\tau$ is consumed exactly per call (`consumed_dt == dt_hint`). The internal IAS15 *does* control its own step within $[0, \tau]$, but that is not visible at the Mercurius trait boundary. |
| `execution_profile` | `Realtime` | Outer per-step wall time bounded by IAS15-on-encountering-pair cost. Worst-case adversarial (deep encounter at high N) would push into Precision territory; flagged for empirical re-classification after Tier 4 measurement. |

---

## Phasing

Implementation lands in one PR with the following commit history:

1. **Lab notebook a-priori** *(this commit)*. Locks decisions before code.
2. **Phase 1 — Encounter diagnostic surfacing.** `PhysicsConfig::close_encounter_threshold: Option<f64>`, `EncounterFlag` enum (`Far` / `Approaching` / `Close`), `compute_closeness` extended to 3D, `warn_diag!` on `Close` transitions. No behavioural change to existing integrators; observability only.
3. **Phase 3 — Mercurius integrator.** `Mercurius` struct + `Integrator` impl. Internal `ChangeoverForceModel` (close + far variants). Mutual Hill radius computation. Per-pair $r_\text{cross}$ table. 7-stage symplectic step. Tier 1 algebraic-identity tests (K=1 → matches WH; K=0 → matches IAS15-only).
4. **Phase 4 — Public surface.** `IntegratorKind::Mercurius` variant + `make_integrator` factory entry + slug + label + description. `apsis-py` exposure (constructor arg). `apsis-app` integrator-pick UI entry.
5. **Cross-implementation parity + bake.** REBOUND MERCURIUS parity entry on a planetary-scattering scenario. §Results / §Interpretation / §Decision populated post-measurement.

Phase 2 of the design proposal (hard-switch hybrid as a pedagogical stepping stone) is documented in §Results below as a measurement entry — the swap-boundary energy drift signature — but not shipped as a public integrator option, because Phase 3 supersedes it entirely. Recording the hard-switch numerical signature validates the §Decision rationale for smooth changeover.

---

## Protocol *(declared a priori, before any code lands)*

### Hypothesis

#### Tier 1 — Algebraic identity *(hard gate)*

Two limits of the changeover function reduce Mercurius to existing first-class integrators. Both must hold to within f64 round-off floor.

| Limit | Reduction | Bound |
| --- | --- | --- |
| $K \equiv 1$ everywhere ($r_\text{cross} \to 0$) | Mercurius = WH (close-field IAS15 sees zero force; far-field K-kick = full planet-planet kick) | Per-step trajectory difference vs `WisdomHolman` over 100 steps on `solar_system_inner`: $|\Delta r| \le 10 \, \varepsilon \cdot |r|$, $|\Delta E / E_0| \le 10^{-14}$ |
| $K \equiv 0$ everywhere ($r_\text{cross} \to \infty$) | Mercurius = WH outer drift + IAS15 close-field on full force (effectively Kepler-aware IAS15) | Per-step trajectory difference vs reference IAS15-direct on `kepler_circular_e0`: $|\Delta a / a| \le 10^{-12}$ over 100 orbits |

**Failure here halts the experiment.** Either limit failing means the K-weighting or sub-integration plumbing has a sign / scaling error that no amount of cross-implementation parity will catch.

#### Tier 2 — Cross-implementation parity vs REBOUND MERCURIUS *(hard gate)*

Reference implementation: REBOUND `MERCURIUS` integrator (`integrator_mercurius.c`), same algorithm and parameters. Run apsis Mercurius and REBOUND MERCURIUS on identical initial conditions, compare conservation diagnostics.

Scenario: `solar_system_outer_with_test_particle`, a Sun + 4 outer planets + 1 close-passing test particle constructed to encounter Jupiter at $r \sim 0.5 \, R_H$ at $t \approx 100$ years.

| Diagnostic | Bound | Horizon |
| --- | --- | --- |
| $|\Delta E / E_0|$ vs REBOUND | $\le 5 \times 10^{-9}$ peak, $\le 10^{-9}$ pre-encounter | $10^4$ orbits of Jupiter |
| $|\Delta L / L_0|$ vs REBOUND | $\le 10^{-10}$ | same |
| Test-particle orbital elements ($a$, $e$, $i$) parity | $\le 10^{-6}$ relative | post-encounter |

The bound is two orders of magnitude looser than the algebraic-identity test because cross-implementation drift accumulates (different rounding pathways, different IAS15 internal state).

#### Tier 3 — Smooth-changeover validation *(gated, documents §Decision)*

Compare two Mercurius modes on the same encounter:

| Mode | Description | Expected energy signature at encounter boundary |
| --- | --- | --- |
| Hard switch (control) | Step-function $K = 0$ inside $r_\text{cross}$, $K = 1$ outside. Implemented as a feature-gated test variant only, not public. | Visible energy jump $\sim 10^{-7}$ at encounter-boundary crossing |
| Smooth $K(y) = y^2 (3 - 2y)$ (production) | Production Mercurius. | Energy drift continuous, $\le 10^{-12}$ at boundary crossing |

The two-mode comparison documents the §Decision rationale for choosing smooth changeover. A static measurement, not a parity gate.

#### Tier 4 — Cost characterisation *(reported, no gate)*

Wall-time per step on three scenarios at fixed outer $\tau$:

| Scenario | Description | Expected vs IAS15-only |
| --- | --- | --- |
| `quiet_solar_system` | Sun + 8 planets, no encounters | Mercurius cheaper: WH-cost outer + zero close-field |
| `single_pair_encounter` | Sun + 4 planets + 1 scattering body in close approach | Mercurius cheaper: close cost local to encountering pair |
| `cluster_chaos` | Sun + 20 mutually-perturbing low-mass bodies, multiple simultaneous encounters | Mercurius probably **not** cheaper — close-field set saturates, IAS15-on-everything competitive |

The third scenario is reported as a calibration row rather than a gate: it characterises the regime where Mercurius stops paying off. Important for the integrator-selection guidance section in `paper.md`.

### Methodology

- Cross-implementation parity (Tier 2) generates REBOUND reference trajectories via `validation/rebound-parity/` infrastructure. New scenario file added under that directory.
- Tier 1 limits are tested via two new feature-gated `Mercurius` constructor knobs (`with_changeover_threshold(0.0)` and `with_changeover_threshold(f64::INFINITY)`) that are public for testing but not surfaced through `IntegratorKind` config.
- All wall-time measurements (Tier 4) run on Cell A (Zen 4 desktop) only; cross-vendor not relevant for an algorithmic / structural experiment.
- The §Results table below is empty until measurements are run.

---

## Results

*Populated after implementation lands, in the same PR's bake commit.*

### Tier 1 — Algebraic identity

| Limit | $|\Delta r| / |r|$ peak | $|\Delta E / E_0|$ peak | Status |
| --- | ---: | ---: | --- |
| $K \equiv 1$ vs WH | TBD | TBD | TBD |
| $K \equiv 0$ vs IAS15-direct | TBD | TBD | TBD |

### Tier 2 — REBOUND MERCURIUS parity

| Diagnostic | apsis | REBOUND | Δ | Bound | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| $|\Delta E / E_0|$ peak | TBD | TBD | TBD | $\le 5 \times 10^{-9}$ | TBD |
| $|\Delta L / L_0|$ peak | TBD | TBD | TBD | $\le 10^{-10}$ | TBD |
| $\Delta a$ test-particle | TBD | TBD | TBD | $\le 10^{-6}$ | TBD |

### Tier 3 — Smooth-changeover documentation

| Mode | Energy signature at boundary |
| --- | --- |
| Hard switch | TBD |
| Smooth $C^1$ | TBD |

### Tier 4 — Cost characterisation

| Scenario | IAS15-only step time | Mercurius step time | Ratio | Note |
| --- | ---: | ---: | ---: | --- |
| `quiet_solar_system` | TBD | TBD | TBD | TBD |
| `single_pair_encounter` | TBD | TBD | TBD | TBD |
| `cluster_chaos` | TBD | TBD | TBD | TBD |

---

## Interpretation

*Populated after Results land.*

---

## Decision

*Populated after Interpretation lands. The decision narrates: ship as-is with the locked decisions, or revise specific decisions based on measured behaviour.*

---

## References

- Rein, H., Hernandez, D. M., Tamayo, D., & Brown, G. (2019). *Hybrid symplectic integrators for planetary dynamics.* MNRAS, 489(4), 4632–4640.
- Rein, H., & Spiegel, D. S. (2015). *IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits.* MNRAS, 446, 1424–1437.
- Wisdom, J., & Holman, M. (1991). *Symplectic maps for the n-body problem.* AJ, 102, 1528–1538.
- Chambers, J. E. (1999). *A hybrid symplectic integrator that permits close encounters between massive bodies.* MNRAS, 304(4), 793–799.
- Design proposal: `docs/proposals/close-encounter-hybrid.md` (PR #32).
