# WHFast — paper-baseline symplectic integrator for long-horizon planetary dynamics

**Date:** 2026-05-13

**Subject:** Adds `IntegratorKind::WHFast` to the apsis integrator zoo as a faithful port of Rein & Tamayo 2015 (MNRAS 452, 376–388). Compensated-summation symplectic integrator that pushes the round-off floor from `O(ε · N)` (apsis WH 1991) to `O(ε · √N)`, enabling 10⁹-orbit horizons. Symplectic corrector of order 17 (Wisdom 1996) ON by default for boosted truncation. Does **not** replace `IntegratorKind::WisdomHolman` — adds to zoo.

**Status:** Protocol declared *a priori*, before any code lands. Locks the algorithmic decisions and validation-against-REBOUND scenario before implementation, per [[feedback_research_commit_discipline]].

**Branch:** `feat/whfast-integrator`, branched from `develop` after PR #87 (long-horizon Mercury 1PN) merges.

**Roadmap context:** First locked target in [[project_integrator_zoo_roadmap]] — paper-baseline credential. Subsequent integrators in the zoo (Implicit Midpoint, Gauss-Legendre RK, AR-Chain) follow only after the v0.1 paper.

---

## Abstract

The apsis integrator zoo ships WH (1991) at the foundational "Wisdom-Holman split, leapfrog kick-drift-kick" level, IAS15 at the precision-controlled adaptive level, and Mercurius at the close-encounter-hybrid level. WHFast (Rein & Tamayo 2015) sits between WH and IAS15 in capability profile but is the **market standard** for long-horizon planetary integration — the integrator any reviewer expects to see in a planetary N-body code. Its differentiators over apsis WH (1991) are two:

1. **Compensated summation** (Kahan / Neumaier-style) on every accumulator that grows with `N_steps`. Reduces round-off accumulation from `O(ε · N_steps)` (where the WH 1991 leapfrog saturates around `~10⁶ steps` at the f64 floor for Solar-System orbits) to `O(ε · √N_steps)`. Unlocks `~10⁹` orbit horizons that WH 1991 cannot reach without artefactual energy drift.

2. **Symplectic correctors** (Wisdom 1996, McLachlan 1995) — additional kick/drift compositions that cancel higher-order error terms in the WH split. Order-17 corrector reduces truncation from `O(dt²)` to `O(dt^{17 + 1})` for the corrected slice, at one-time-only cost (correctors compose into a single boundary kick before/after the inner KDK loop). Practical effect at Solar-System cadence: ~3-5× tighter conservation per step, free.

WHFast is paper-baseline. With WHFast in the zoo, apsis covers the four canonical integrator regimes: low-cost real-time (VV / Yoshida4), market-standard long-horizon planetary (WHFast), close-encounter hybrid (Mercurius), and adaptive precision (IAS15). Federation thesis ([[project_thesis_anchor]]): every regime composes with every perturbation through a single contract.

---

## Motivation

Three claims chain into this experiment:

1. **The Mercury 1PN convergence experiment showed f64 round-off becomes visible at ~10⁸ perturbation evaluations** (`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md` §Convergence experiment). That matches the regime where compensated summation pays off — `N_steps ≳ 10⁸` is where WH 1991's `O(ε·N)` floor crosses the science-relevant precision. WHFast's `O(ε·√N)` floor pushes the crossover by ~`√10⁸ ≈ 10⁴`, which is exactly the multi-Gyr Solar System integration regime planetary dynamics papers cite.

2. **The integrator-zoo positioning thesis ([[project_integrator_zoo_roadmap]]) makes WHFast first** — it is the credential reviewers expect before they engage with novel integrators (Mercurius, federation contract). Without WHFast in the zoo, the integrator-of-integrators claim looks under-equipped.

3. **The federated FPM thesis depends on the integrator slot accepting any first-class symplectic-class integrator.** WHFast composes with apsis-1pn and any future perturbation through the same `Mercurius::interaction_step`-style mechanism (PR #86). Validating WHFast + 1PN at the same ~3-7 ppm precision floor as Mercurius + 1PN closes the federation contract demonstration over the entire active-integrator set.

### What this experiment is NOT testing

- **Not 10⁹-orbit endurance.** A multi-Gyr run would saturate at the f64 floor where WHFast's compensated summation matters most, but the wall-time budget per validation run rules it out. Cap at 10⁵ Jupiter orbits (~1.2 Myr) — strong enough to demonstrate the WHFast advantage over WH 1991 without committing to a multi-day run.
- **Not WHFast512 / SIMD-vectorised batch WHFast.** REBOUND ships a separate `whfast512` integrator (Rein, Tamayo & Brown 2024) that 8-way-vectorises the Kepler solver. Out of scope for this PR; would land as a separate `IntegratorKind::WHFast512` if/when SIMD-batched planetary integration becomes a priority.
- **Not Jacobi-coordinate WHFast.** REBOUND WHFast supports both democratic-heliocentric and Jacobi coords; apsis sticks with DH for consistency with WH 1991 and Mercurius. Jacobi support is an opt-in follow-up if a hierarchical-system scenario surfaces a need.
- **Not a replacement for WH 1991.** Both ship in the zoo. WH 1991 stays as the pedagogical / minimal symplectic baseline; WHFast is the production long-horizon path.

---

## Locked design decisions

| Question | Decision | Rationale |
| --- | --- | --- |
| Public name | `IntegratorKind::WHFast`, slug `"whfast"` | Matches REBOUND. Self-documenting. |
| Coordinates | Democratic-heliocentric (DH) | Matches existing `WisdomHolman` and `Mercurius`. Consistent integrator zoo convention; reduces "did the coord transform compose right" as a debugging axis. Jacobi support deferred. |
| Kepler propagator | Existing `apsis::physics::integrator::kepler::kepler_step` (Stumpff universal variable) | Already correct for `μ ≠ 1` post PR #84. Mercurius already uses it; same code path, same numerical behaviour. |
| Compensated summation | Neumaier-style on planet positions, planet velocities, COM position, integrator's own internal accumulators | Rein-Tamayo 2015 §3.2 specifies Kahan/Neumaier on every length-`N_steps` accumulator. Implemented as `(value, compensator)` pairs with `add_cs(value, compensator, increment)` updates. |
| Symplectic corrector | Order 17 (Wisdom 1996), ON by default; opt-out via builder | REBOUND default. ~3-5× truncation reduction at one-time-only boundary cost. Opt-out reserved for studies isolating the compensated-summation contribution from the corrector contribution. |
| Step size | Fixed-step, user-supplied via `with_dt` (matches WH 1991 + Mercurius) | WHFast is fixed-step by construction. Adaptive overlay (a la `DtMode::Adaptive`) breaks symplectic invariants and is silently disabled when `IntegratorKind::WHFast` is selected. |
| Hierarchical-system requirement | Same as WisdomHolman (`HierarchySignal::classify` ≥ Borderline) | WH-class derivation; opens up to non-hierarchical configurations only at the cost of breaking the small-parameter expansion. Same enforcement path as `set_integrator(WisdomHolman)`. |
| Force-model determinism | `requires_deterministic_force = false` | WHFast computes its own K-weighted planet-planet kicks internally (same as Mercurius); does not use `ctx.force` directly. |
| Perturbation handling | `interaction_step` accumulates `ctx.perturbations` after the planet-planet kick, before the velocity update — same Strang-split position as Mercurius post PR #86. The corrector slice does **not** apply perturbations (correctors fold into closed-form Kepler-order error reduction; mixing them with perturbations defeats the construction). | Symmetric `dt` perturbation strength split across the two τ/2 kicks, matching the WH-class convention since 1991. |
| Snapshot codec | Byte 5 in the on-disk codec (continues from `Mercurius = 4`, PR #83) | Standard append. |

### What's NOT a parameter

- **Not WHFast variant choice (Jacobi vs DH).** Locked to DH.
- **Not corrector order beyond 17.** REBOUND offers 11 + 17; we ship 17 only for v0.1. Adding 11 is a configuration-surface follow-up.
- **Not safe-mode / synchronisation flush behaviour.** REBOUND's `safe_mode = 1` flushes the corrector after every step, vs `safe_mode = 0` only at synchronisation boundaries. We use `safe_mode = 1` (flush every step) to keep the public API simple — the corrector is always synchronised when the user reads body state. The lazier `safe_mode = 0` is a perf knob for power users; deferred.

---

## Algorithm

Following Rein & Tamayo 2015 (MNRAS 452, 376) §3.

### Per-step structure (DH coordinates, `safe_mode = 1`)

For each `step(dt)`:

1. **Inertial → DH** (same as Mercurius). Capture COM position / velocity for restoration.
2. **Pre-kick corrector boundary** (only if `with_correctors == true` and integrator is "synchronised", i.e. on its first step or after a flush): apply the order-17 corrector pre-composition. Sequence of small kicks/drifts that absorbs the higher-order error terms the bare KDK leaves on the boundary.
3. **`interaction(τ/2)`** — K-weighted half-kick on planet velocities (same shape as Mercurius). Includes registered perturbations.
4. **`jump(τ/2)`** — DH indirect drift on planet positions.
5. **`com(τ)`** — advance the inertial COM by `τ · v_com`.
6. **`kepler(τ)`** — analytical Kepler drift around the central body for every planet, via `kepler_step`.
7. **`jump(τ/2)`**.
8. **`interaction(τ/2)`**.
9. **Post-kick corrector boundary** (only if `with_correctors == true` and synchronisation requested): inverse of stage 2, restores the un-corrected representation. Required before the user reads body state.
10. **DH → inertial** (same as Mercurius).

Per Wisdom 1996, the corrector is a closed-form composition of `Z_n` operators (each itself a small kick/drift sequence) chosen so the leading-error-term cancellation gives total order 17 in `dt`. The sequence and coefficients are tabulated in Wisdom (1996, AJ 112, 1305) Table 1; reproduced in Rein & Tamayo (2015) §3.2.

### Compensated summation

Every per-step accumulator that grows with the integration horizon is stored as a `(value, compensator)` pair:

- planet positions: `(q_i, c_q_i)` per planet
- planet velocities: `(v_i, c_v_i)` per planet
- COM position: `(com_pos, c_com_pos)`

Updates use Neumaier compensated summation:

```text
add_cs(value, compensator, inc):
    y = inc - compensator
    t = value + y
    compensator = (t - value) - y
    value = t
```

This bounds round-off accumulation to `O(ε · √N_steps)` for length-`N` summations (Higham 2002 §4.5). Compared to the `O(ε · N)` naive accumulator, the floor crossover for Solar-System orbits (`ε = 2.2e-16`, dt ≈ Mercury_period/200) shifts by ~`√(N_critical) = √(1/ε) ≈ 7 × 10⁷` orbits. WH 1991's f64 floor at ~10⁶ steps becomes WHFast's ~7 × 10¹³ steps — well into the multi-Gyr horizon planetary papers cite.

### Symplectic corrector — pedagogical sketch

The bare WH leapfrog has Hamiltonian error `H_err = c_2 · dt² · {H_K, {H_K, H_I}} + O(dt^4)` where `{·, ·}` is the Poisson bracket. The order-17 corrector `Z` is constructed so `Z H_KDK Z^{-1}` has its leading error pushed to `O(dt^{18})`. This is achieved by a closed-form composition of additional Z_n operators with tabulated coefficients (Wisdom 1996 §3, Rein-Tamayo 2015 §3.2 reproduce the coefficients).

The cost is one Z application per integration boundary (begin, end). For a 10⁵-orbit run, the cost is amortised to negligible per-step overhead. The truncation gain is ~3-5× tighter at typical Solar-System cadence, validated empirically against REBOUND WHFast in the §Tier 2 below.

---

## Protocol

### Hypothesis

#### Tier 1 — Federation contract: WHFast + apsis-1pn matches Mercurius + apsis-1pn on the existing Mercury 1PN scenario *(hard gate)*

Re-runs the long-horizon Mercury 1PN scenario (`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`) with WHFast in place of Mercurius, comparing cumulative `Δω` against the GR analytical prediction.

| Metric | Bound | Rationale |
| --- | ---: | --- |
| WHFast + apsis-1pn `\|Δω(end) − Δω_GR(end)\| / \|Δω_GR(end)\|` | ≤ 10⁻⁵ (10 ppm) | Same bound as Mercurius+1PN. WHFast is symplectic same-class; the Mercury IC-precision floor (~3 ppm) dominates regardless of integrator. |
| WHFast vs Mercurius cross-integrator parity | ≤ 5 × 10⁻⁵ (50 ppm) | Same bound as IAS15-vs-Mercurius. Both are symplectic with shared 1PN evaluator; cross-impl drift sits at integrator-truncation level. |

This is a hard gate because it confirms `Mercurius::interaction_step`-style perturbation wiring works for any WH-class integrator, not just Mercurius — the federation thesis claim.

#### Tier 2 — Cross-implementation parity: apsis WHFast vs REBOUND WHFast on Solar System outer 4 planets *(hard gate)*

Sun + Jupiter + Saturn + Uranus + Neptune (no test particle, no encounter), integrated for 10⁵ Jupiter orbits (~1.2 × 10⁶ years) under both apsis WHFast and REBOUND WHFast with corrector order 17 on both sides. Conservation diagnostics compared.

| Metric | Bound | Rationale |
| --- | ---: | --- |
| Cross-impl `\|ΔE/E₀\|` peak | ≤ 5 × 10⁻¹¹ | Mercurius parity (PR #85) hit `3.7 × 10⁻¹¹` at 10⁴ years; WHFast at 10⁵ Jupiter orbits (~10× longer in cumulative cycles) should sit at similar level if compensated summation works. Same compensated-summation algorithm on both sides; the f64 random-walk floors should agree to within the 1PN-formula evaluation noise. |
| Cross-impl `\|ΔLz/Lz₀\|` peak | ≤ 10⁻¹² | Angular momentum conserved exactly by analytical Kepler drift on both sides; cross-impl floor sits at the corrector's f64 noise. |
| Per-side `\|ΔE/E₀\|` peak | ≤ 10⁻¹⁰ each | WHFast's compensated-summation floor on a 10⁵-orbit Solar-System integration. Demonstrates WHFast pushing past WH 1991's `O(ε·N) ≈ 10⁻⁹` saturation for the same horizon. |

Same-scenario REBOUND parity test demonstrates the canonical-reference equivalence (PR #85 pattern). The 50 ppm bound applies to integrators on the same physics; the 5 × 10⁻¹¹ bound here applies to integrators with the same algorithm, which is a tighter claim.

#### Tier 3 — WHFast vs WH 1991 on the same scenario *(reported, no gate)*

Re-runs Tier 2 with `IntegratorKind::WisdomHolman` instead of WHFast on the apsis side. Compares per-side `|ΔE/E₀|` peak. Expected: WH 1991 sits 1-2 orders of magnitude above WHFast on this horizon, demonstrating the compensated-summation advantage. If this difference is **not** observed, either the compensated summation is not actually engaging or the scenario is at too short a horizon to surface the difference; both indicate a follow-up is needed.

### Methodology

Three-side test infrastructure following the existing `validation/rebound-parity/{kepler,figure8,pythagorean,retrograde,mercurius}/` pattern:

1. **apsis WHFast side** (`crates/apsis/examples/rebound_parity_whfast.rs`): instantiates the Solar System outer-4-planets scenario, runs WHFast for 10⁵ Jupiter orbits with corrector ON, writes per-Jupiter-orbit (state, total energy, total Lz) to `validation/rebound-parity/whfast/out/apsis.csv`.
2. **REBOUND WHFast side** (`validation/rebound-parity/whfast/rebound_side.py`): mirrors the apsis side with REBOUND's `whfast` integrator + `safe_mode = 1` + `corrector = 17`.
3. **Comparator** (`validation/rebound-parity/whfast/compare.py`): loads both CSVs, computes Tier 2 metrics, exits 0 iff every gated metric is within tolerance.
4. **Tier 1 (Mercury 1PN)**: extends the existing `validation/mercury-1pn-long-horizon/run.py` orchestrator with a `--include-whfast` flag, which adds a third side (apsis WHFast + apsis-1pn on the same Mercury 1PN scenario).

---

## Results

*Populated post-implementation + run.*

### Tier 1 — WHFast + apsis-1pn vs GR analytical + Mercurius parity

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| WHFast + 1PN Δω rel err vs GR (end) | TBD | ≤ 10⁻⁵ | TBD |
| WHFast vs Mercurius cross-integrator parity | TBD | ≤ 5 × 10⁻⁵ | TBD |

### Tier 2 — apsis WHFast vs REBOUND WHFast (Solar System outer)

| Metric | Observed | Bound | Status |
| --- | ---: | ---: | --- |
| Cross-impl `\|ΔE/E₀\|` peak | TBD | ≤ 5 × 10⁻¹¹ | TBD |
| Cross-impl `\|ΔLz/Lz₀\|` peak | TBD | ≤ 10⁻¹² | TBD |
| Per-side `\|ΔE/E₀\|` peak (apsis) | TBD | ≤ 10⁻¹⁰ | TBD |
| Per-side `\|ΔE/E₀\|` peak (REBOUND) | TBD | ≤ 10⁻¹⁰ | TBD |

### Tier 3 — WHFast vs WH 1991 on the same scenario

| Side | `\|ΔE/E₀\|` peak (Solar outer, 10⁵ Jupiter orbits) |
| --- | ---: |
| apsis WHFast | TBD |
| apsis WH 1991 | TBD |
| Ratio (WH 1991 / WHFast) | TBD (expected ~10²-10³) |

---

## Interpretation

*Populated post-results.*

---

## Decision

*Populated post-interpretation. Possible outcomes:*

- **All gates pass** → WHFast enters the v0.1 paper §Validation table alongside WH 1991, IAS15, and Mercurius. Federation contract validated for the integrator-zoo target set.
- **Tier 1 passes, Tier 2 fails** → WHFast composes with perturbations correctly but apsis WHFast and REBOUND WHFast disagree on the canonical scenario; bisect against the corrector implementation or compensated-summation order.
- **Tier 1 fails** → perturbation wiring through WHFast's interaction_step has a regression vs Mercurius; check `interaction_step` signature compatibility.
- **Tier 3 shows no advantage of WHFast over WH 1991** → either compensated summation is not engaging or the horizon is too short to surface the difference. Diagnostic: enable per-step round-off tracing; look at the f64 floor crossover empirically.

---

## References

- Rein, H., & Tamayo, D. (2015). *WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations.* MNRAS, 452(1), 376–388.
- Wisdom, J. (1996). *Symplectic correctors for canonical Levi-Civita Kustaanheimo-Stiefel regularization.* AJ, 112, 1305.
- Wisdom, J., & Holman, M. (1991). *Symplectic maps for the n-body problem.* AJ, 102, 1528.
- McLachlan, R. I. (1995). *Composition methods in the presence of small parameters.* BIT, 35, 258–268.
- Higham, N. J. (2002). *Accuracy and Stability of Numerical Algorithms* (2nd ed.). SIAM, §4.5 on Kahan / Neumaier compensated summation.
- Existing apsis WH (1991): `crates/apsis/src/physics/integrator/wisdom_holman.rs`.
- Mercurius implementation lab notebook: `docs/experiments/2026-05-13-mercurius-hybrid.md`.
- Long-horizon Mercury 1PN: `docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`.
- Integrator zoo roadmap: [[project_integrator_zoo_roadmap]].
