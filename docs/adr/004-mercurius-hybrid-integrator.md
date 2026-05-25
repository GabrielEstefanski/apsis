# ADR-004 — Mercurius Close-Encounter Hybrid Integrator

**Status:** Accepted
**Date:** 2026-05-13
**PRs:** [#83](https://github.com/GabrielEstefanski/apsis/pull/83) (implementation),
[#84](https://github.com/GabrielEstefanski/apsis/pull/84) (Kepler-µ fix),
[#85](https://github.com/GabrielEstefanski/apsis/pull/85) (REBOUND parity Tier 1),
[#86](https://github.com/GabrielEstefanski/apsis/pull/86) (perturbation routing fix)
**Lab notebooks:** `docs/experiments/2026-05-13-mercurius-hybrid.md`,
`paper/notebooks/2026-05-13-rebound-parity-mercurius.md`

---

## Context

Apsis shipped two integrators that cover orthogonal regimes:

* **IAS15** — 15th-order Gauss–Radau, adaptive sub-step. Handles close
  encounters by shrinking `dt` globally. Cost: every body in the
  system pays the encounter rate, even when only one pair is
  encountering.
* **Wisdom–Holman** — 2nd-order symplectic. Per-step cost is one
  Kepler step per planet plus pairwise kicks; bounded forever
  regardless of horizon length. Breaks down inside ~3 Hill radii of
  any body because the perturbation expansion's small parameter
  inverts.

Sparse-encounter planetary scenarios (Sun + ≤ 10 planets, Jupiter-
crossing test particle, planetesimal disk) sit between the two: WH
costs are right for the far-field, IAS15 precision is needed for the
close-encounter pair, and neither alone is the answer.

Three regimes drive the case for a hybrid:

1. **Sparse-encounter planetary.** IAS15 alone amplifies encounter
   cost by N. The encountering pair determines `dt`; the rest of the
   system pays it.
2. **Long-horizon hierarchical** (10⁴ – 10⁶ orbits, planet–planetesimal).
   IAS15 saturates at the Brouwer-law floor for energy drift over
   very long horizons; WH preserves secular stability arbitrarily
   long. Apsis had no integrator that handled both regimes within
   one run.
3. **Federated-FPM positioning.** The thesis frames apsis as a
   federation where physics components are first-class. A hybrid
   that *composes two existing first-class integrators into a single
   algorithm* federates the temporal axis the same way perturbation
   crates federate the force axis. Worth implementing on
   positioning grounds alone.

Rein, Hernandez, Tamayo & Brown (2019, *MNRAS* 489, 4632) describe
MERCURIUS: a smooth `K(r/r_crit)` weighting routes each pairwise
interaction between a WH-step "far" half and an IAS15-step "close"
half. REBOUND ships the original implementation (Rein et al. 2019).

---

## Decision

Implement Mercurius as a first-class apsis integrator
(`crates/apsis/src/physics/integrator/mercurius.rs`,
`IntegratorKind::Mercurius`), structurally faithful to Rein et al. 2019
and bit-comparable against REBOUND on a Tier 1 (Kepler) baseline.

### Algorithm structure

* Outer step: WH leapfrog (kick–drift–kick) on the K-attenuated
  far-field force.
* Inner step: IAS15 sub-integration on the (1 − K)-weighted close-
  field force, advanced over the same outer `dt` for any pair whose
  separation crosses `r_crit`.
* Changeover function: smooth polynomial `K(y)` with `K(0) = 0`,
  `K(1) = 1`, `K'(0) = K'(1) = 0`. Default `r_crit = α · max(r_Hill_i,
  r_Hill_j)` with `α = 3` (REBOUND default).

### Validation tiers

| Tier | Scenario | Status | Notebook |
|------|----------|--------|----------|
| 1 — Structural / Kepler | Two-body, no close encounter | PASS | `2026-05-13-mercurius-hybrid.md` §Tier 1 |
| 2 — REBOUND parity | Sun + 4 outer planets + Jupiter-crosser, 10⁴ yr | Tier 1 PASS, Tier 2 chaotic FAIL | `2026-05-13-rebound-parity-mercurius.md` |
| 3 — Cost profile | Sparse-encounter scaling vs IAS15 alone | PASS | `2026-05-13-mercurius-hybrid.md` §Tier 4 |

The Tier 2 chaotic FAIL is documented as expected: at 10⁴ years the
test particle's Lyapunov time is shorter than the integration window,
so trajectory divergence between apsis and REBOUND is dominated by
shadow-Hamiltonian truncation differences both sides accumulate
independently. The integrator is correct; the comparison metric needs
to move from |Δr| to orbital-element drift (per memory entry on
adaptive integrator parity).

### Federation invariants

Mercurius participates in the same `Operator` / `HamiltonianOperator`
contract as the rest of the integrator stack (see ADR-005). The
perturbation-routing fix in #86 closed a gap where the K-weighted
kicks bypassed `ctx.perturbations`; with that fix, registered
operators (1PN, radiation, J2, …) compose with Mercurius the same way
they compose with VV/Y4/WH/IAS15.

The Kepler-µ fix in #84 was a deeper bug surfaced by the parity gate:
the universal Kepler equation kept `r₀ · χ` scaled by √µ when it must
not. The fix landed independently in #84 and lifted both Mercurius
and WH parity at the same time.

---

## Alternatives rejected

| Alternative | Reason rejected |
|---|---|
| IAS15-only with adaptive dt | Encounter cost amplified by N for sparse-encounter scenarios; no long-horizon ceiling. The whole motivation for a hybrid. |
| WH-only with manual close-encounter detection | Re-implements MERCURIUS poorly; loses Rein et al.'s changeover smoothness. |
| External REBOUND binding | Couples apsis to a C dependency for an algorithm we want as a first-class apsis artifact. Federation positioning wants the integrator *in* apsis, not adjacent to it. |
| Higher-order outer step (Yoshida-4 instead of WH) | WH is the published baseline. A higher-order outer would be a deviation from Rein et al. 2019 with no validation literature. Defer to a future Mercurius variant. |

---

## Consequences

**Good:**
- Sparse-encounter planetary scenarios run at WH cost in the far-field
  with IAS15 precision on the encountering pair.
- Long-horizon hierarchical runs become tractable without choosing
  between secular drift (IAS15) and close-encounter blow-up (WH).
- The federation thesis gains a concrete temporal-axis exemplar:
  Mercurius is composed from existing first-class integrators, no
  algorithmic rewrite.

**Neutral:**
- Mercurius is now the recommended integrator for solar-system
  scenarios; IAS15 remains the default for high-precision short-
  horizon work and WH for pure secular runs without close encounters.
- The REBOUND parity Tier 2 metric is documented as a known gap
  pending the orbital-element comparison rewrite.

**Watch out:**
- Close-encounter detection adds a per-step cost proportional to
  pairs near `r_crit`. For dense scenarios (≫ 10 close pairs at any
  step) the overhead approaches IAS15-only. The cost crossover lives
  in the Tier 4 cost notebook.
- Picking `α` larger than ~5 enlarges the IAS15 region until the
  hybrid degrades into IAS15-only with extra overhead. The default
  `α = 3` matches REBOUND and the Rein et al. recommendation.
