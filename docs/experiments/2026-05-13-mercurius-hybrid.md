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

Faithful port of REBOUND's `integrator_mercurius.c` (Rein et al. 2019). The
implementation mirrors REBOUND's structure rather than the simplified
"clean Hamiltonian split" sketch the design proposal carried — that
sketch is mathematically incomplete: a pure $(1-K)\,V$ flow generates a
kick (no position evolution), so feeding it to IAS15 (which integrates
$\ddot x = a$) produces a $v\cdot \tau$ free drift on top of the
analytical Kepler step and double-counts position evolution. REBOUND's
actual structure is a **rewind hybrid**: the WH branch advances all
particles, an encounter detector identifies close pairs, those pairs
rewind to pre-Kepler state, and IAS15 re-integrates them over the same
$\tau$ window with the full Sun-pull plus the $(1-K)$-weighted
planet-planet residual. Non-encountering particles keep their
analytical-Kepler positions.

### Per-step structure

In democratic-heliocentric coordinates (planets relative to Sun position;
inertial momenta in COM rest frame):

1. **Interaction kick** $\mathrm{int}(\tau/2)$: $v_i \mathrel{+}= \tau/2 \cdot a_i^K$ where $a_i^K$ is the K-weighted planet-planet acceleration (Sun pull excluded — handled analytically by stage 4).
2. **Jump** $\mathrm{jmp}(\tau/2)$: $q_i \mathrel{+}= \tau/2 \cdot (\sum_j m_j v_j) / m_0$ on every planet.
3. **COM drift** $\mathrm{com}(\tau)$: track inertial COM motion over the step.
4. **Backup + Kepler drift** $\mathrm{kep}(\tau)$: snapshot $(q_i, v_i)$ pre-Kepler; analytically advance every planet around the Sun for $\tau$.
5. **Encounter predict**: scan all pairs, compute $r_\text{min}$ over the step via cubic-Hermite interpolation between pre- and post-Kepler positions; flag pairs with $r_\text{min} < \mathrm{dcrit}_{ij}$ where $\mathrm{dcrit}_{ij} = \max(\mathrm{dcrit}_i, \mathrm{dcrit}_j)$.
6. **Encounter step** $\mathrm{enc}(\tau)$: for every flagged particle, restore the pre-Kepler snapshot. Run IAS15 over $\tau$ on the encountering subset only. IAS15's force model returns full Sun-pull + $(1-K)$-weighted planet-planet (gravity in REBOUND's `mode = 1`); the K-weighted far-field contribution from non-encountering planets onto the encountering pair is implicit because $K \to 1$ at large separation makes the residual zero anyway.
7. **Jump** $\mathrm{jmp}(\tau/2)$: same as stage 2.
8. **Interaction kick** $\mathrm{int}(\tau/2)$.

Outer truncation: $O(\tau^3)$ per step, $O(\tau^2)$ over the integration.
The K-weighted interaction kicks are the only second-order operator in
the split; the rest is either analytical (Kepler), exact (jump, COM),
or high-precision (IAS15 on encountering pairs).

### Changeover function (REBOUND `L_mercury`)

$C^2$ quintic Hermite polynomial with a `0.1 · dcrit` deadband:

$$
y = \frac{d - 0.1 \, \mathrm{dcrit}}{0.9 \, \mathrm{dcrit}}, \qquad
L(y) = \begin{cases}
0 & y \le 0 \\
10 y^3 - 15 y^4 + 6 y^5 & 0 < y < 1 \\
1 & y \ge 1
\end{cases}
$$

The deadband ensures $L \equiv 0$ for $d < 0.1 \, \mathrm{dcrit}$ — IAS15
gets full responsibility deep inside the encounter, with no leakage of
K-weighted force into the close regime.

### Critical radius (REBOUND `dcrit_for_particle`)

Per-particle, max of four criteria:

$$
\mathrm{dcrit}_i = \max \begin{cases}
v_c \cdot 0.4 \, \tau & \text{average velocity} \\
|v_i| \cdot 0.4 \, \tau & \text{current velocity} \\
\alpha \cdot a_i \cdot \sqrt[3]{m_i / (3 m_0)} & \text{Hill radius} \\
2 r_i^\text{phys} & \text{physical radius}
\end{cases}
$$

with $a_i$ the osculating semi-major axis, $v_c = \sqrt{G m_0 / |a_i|}$ the
circular velocity. Pair-wise scale: $\mathrm{dcrit}_{ij} = \max(\mathrm{dcrit}_i, \mathrm{dcrit}_j)$.

---

## Locked design decisions *(resolves design proposal §6 a priori)*

| Question | Decision | Rationale |
| --- | --- | --- |
| §6.1 Hill radius $M_\star$ for non-hierarchical systems | Refuse with a `warn_diag!` event + `used_fallback = true`; step does not advance | Mercurius assumes a dominant central body for the analytical Kepler drift. Non-hierarchical fallback would be a different algorithm. Honest failure → user routes to IAS15 directly. Matches REBOUND. |
| §6.2 Changeover function shape | REBOUND `L_mercury`: $C^2$ quintic Hermite $L(y) = 10 y^3 - 15 y^4 + 6 y^5$ with $y = (d - 0.1\,\mathrm{dcrit})/(0.9\,\mathrm{dcrit})$ | The 0.1·dcrit deadband ensures $L \equiv 0$ deep in the encounter — IAS15 carries the full force without leakage from the K-weighted kick. Earlier draft used a $C^1$ cubic; updated after reading `integrator_mercurius.c` (mid-experiment revision documented in [[feedback_research_commit_discipline]]). |
| §6.3 Default $\alpha$ (Hill multiplier for `dcrit`) | $\alpha = 3$ | REBOUND default; validated by Rein et al. 2019 §3 against several planetary scattering scenarios. |
| §6.4 `fast` integrator selection | Wisdom-Holman analytical Kepler drift, K-weighted planet-planet half-kicks | The K-weighted kick is the only second-order operator; Kepler is analytical, jump and COM are exact, encounter step is high-precision. Yoshida-4 would lose the analytical Kepler advantage. |
| §6.5 Per-pair vs global changeover | Per-pair, via $\mathrm{dcrit}_{ij} = \max(\mathrm{dcrit}_i, \mathrm{dcrit}_j)$ | REBOUND structure — per-particle critical radius derived from 4 criteria (avg velocity, current velocity, Hill radius, physical radius), pair-wise reduced by max. Not the "mutual Hill radius per pair" formulation the design proposal sketched; the 4-criterion form is what REBOUND ships and what Rein et al. 2019 measure against. |
| §6.6 Outer step shrink near encounters | Outer $\tau$ stays fixed; the inner IAS15 (rewind + restart per encounter step) absorbs cost | Canonical Rein et al. 2019 design. Encounter cost stays local to the encountering subset. |

Additional decisions (not in §6, locked here):

| Decision | Choice | Rationale |
| --- | --- | --- |
| Public name | `IntegratorKind::Mercurius`, slug `"mercurius"` | Matches Rein et al. 2019 algorithm name and REBOUND. Self-documenting. |
| Algorithmic structure | Rewind hybrid (REBOUND-faithful), not clean Hamiltonian split | A pure $(1-K)\,V$ Hamiltonian generates a kick (no position evolution); feeding it to IAS15 (which integrates $\ddot x = a$) introduces a $v\cdot\tau$ free drift that double-counts Kepler. REBOUND structure is correct: WH advances all → encounter detect → encountering pairs rewind → IAS15 re-integrates them with full Sun-pull + $(1-K)$ planet-planet. |
| Force model interaction | Mercurius bypasses `IntegratorContext::force` for its own pair forces (direct O(N²) K-weighted kicks + close-field IAS15 with its own internal force model) | The Barnes-Hut tree adds no value at planetary $N$, and the per-pair K weighting requires explicit pair iteration. |
| `compute_closeness` 2D defect | Fix in this PR (extends to 3D) | Phase 1 EncounterFlag reads `r_min`; latent 2D bug affects Phase 1 correctness on out-of-plane configurations. Single-line fix; in scope. |
| `requires_deterministic_force` | `false` | Mercurius does not rely on the outer `ctx.force` being deterministic (it computes its own K-weighted forces internally). The inner IAS15's deterministic-force requirement is satisfied by the embedded close-field force model. |
| `controls_own_step_size` | `false` | Outer $\tau$ is consumed exactly per call. The internal IAS15 *does* control its own step within encounter sub-integrations, but that is not visible at the Mercurius trait boundary. |
| `execution_profile` | `Realtime` | Outer per-step wall time bounded by IAS15-on-encountering-subset cost. Adversarial worst case (deep simultaneous encounters at high N) would push into Precision; flagged for empirical re-classification after Tier 4. |

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

Two limits of the changeover function reduce Mercurius to existing first-class behaviours:

- **No-encounters limit** ($\mathrm{dcrit} \to 0$, equivalently $\alpha = 0$): every pair has $L = 1$, the encounter detector flags nothing, the rewind+IAS15 stage is a no-op, and Mercurius reduces to Wisdom-Holman with the same K-weighted-= full kick decomposition. Bound: $|\Delta E / E_0| \le 10^{-12}$ over 200 steps on `quiet_planetary` with planet/Sun mass ratio $10^{-6}$. The bound is set by the difference in jump-step placement between apsis WH (single jump after Kepler) and Mercurius (jump halves on each side of Kepler) — both are 2nd-order symplectic but accumulate rounding differently. A direct trajectory-equality check is *not* part of Tier 1 because the two splits differ at $O(\tau^3 \cdot m_p / m_0)$ per step.
- **All-encounters limit** ($\mathrm{dcrit} \to \infty$): every pair triggers the encounter detector, IAS15 integrates the full system every step, Mercurius reduces to IAS15-with-restart. Bound: per-step trajectory $|\Delta r / r| \le 10^{-10}$ vs an independent `Ias15` instance on `quiet_planetary` over 50 steps (loose because the encounter step uses IAS15-restart per outer step rather than continuous IAS15 state).

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
