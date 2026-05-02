# Proposal — close-encounter hybrid integrator

**Status:** Draft. Implementation deferred; this document specifies the work, not its execution.
**Reference:** Rein, H., Hernandez, D. M., Tamayo, D., & Brown, G. (2019). *Hybrid symplectic integrators for planetary dynamics.* MNRAS, 489(4), 4632–4640 — the MERCURIUS paper.
**Out of current scope:** WHFast refactor (tracked separately as TD-008). This proposal assumes one of two starting points: (a) a future WHFast-class implementation in `apsis::physics::integrator`, or (b) Yoshida-4 as the symplectic far-field stand-in until the WH refactor lands.

---

## 1. Motivation

The current close-encounter handling in `apsis` is reactive and implicit: IAS15's adaptive controller shrinks `dt` when the local truncation error grows, which happens during close approaches. This works (Mercury 1PN reproduces GR to ~1 ppm; the Pythagorean parity scenario survives 70 t.u. of repeated close encounters at f64 floor) but has two consequences:

1. **Cost accrual is global, not local.** When two bodies in a 100-body system enter a close approach, the entire integration shrinks `dt` to resolve the encounter, paying the high-order Picard-iteration cost on every body in the system. For a sparse system (one Sun + planets) where encounters are rare, the cost is dominated by the encountering pair, but the cost-per-step is amplified by N.

2. **No symplectic-class option for long-horizon dense scenarios.** A 100-orbit Mercury run under IAS15 is fine; a $10^4$-orbit run remains within the validated regime; a $10^6$-orbit run starts hitting Brouwer-law saturation. Symplectic integrators (WH, WHFast) preserve secular stability over arbitrary horizons but break under hard close encounters. The standard answer is a hybrid: symplectic far-field + non-symplectic high-precision near-field.

The MERCURIUS approach (Rein et al. 2019) is the canonical solution. This proposal adopts it, scaled to apsis's architecture and capability.

## 2. Current state in apsis

A read of the codebase as of `742a64e` (master tip, post v0.1.0-alpha.1):

| Element | Location | Behaviour |
| --- | --- | --- |
| Pairwise minimum distance | `crates/apsis/src/core/system/helpers.rs::compute_closeness` | Computed once per `System::step()`; stored in `System.r_min` for $N < 512$ |
| `System.r_min`, `System.softening_max` | `crates/apsis/src/core/system/mod.rs` | Diagnostic-only fields; surfaced through `Metrics` for UI/logging |
| IAS15 adaptive controller | `crates/apsis/src/physics/integrator/ias15.rs` | Reactive `dt` shrink based on Picard-iteration truncation; floor at `DT_MIN = 1e-12`; `degraded_total` counter on floor-pinned steps |
| `StepResult.degraded` | `crates/apsis/src/physics/integrator/traits.rs` | Boolean flag, passive observation |
| WH → Yoshida-4 fallback | `crates/apsis/src/physics/integrator/wisdom_holman.rs` lines 118–131 | Triggers when central body is not dominant ($M_0 / \sum m_i < 10$); unrelated to close encounters |
| Force-model / integrator pairing | `crates/apsis/src/core/system/config.rs::set_integrator` | Auto-corrects force model on integrator change; no encounter-driven switching |

What does not exist:

- No threshold on `r_min` triggering action.
- No hybrid integrator type that holds two integrator instances.
- No two-zone Hamiltonian decomposition.
- No changeover function $K(d / r_\text{cross})$.
- No Hill-radius computation as an integrator-relevant scale.
- No smooth transition; the only "switching" in the codebase is the WH→Yoshida-4 fallback above.

The terrain is therefore **diagnostic-ready**: `r_min` is plumbed, but no consumer reads it for control.

## 3. Architecture target

Following Rein et al. (2019). Notation: $r_{ij} = |x_i - x_j|$, $r_\text{cross}$ = changeover radius for the pair.

### 3.1 Hamiltonian decomposition

The Hamiltonian is split per pair into a far-field term and a close-field term using a smooth changeover function $K$:

$$
H_{ij}^\text{far} = -G m_i m_j \, K(r_{ij} / r_\text{cross}) / r_{ij}
$$

$$
H_{ij}^\text{close} = -G m_i m_j \, (1 - K(r_{ij} / r_\text{cross})) / r_{ij}
$$

with $H = \sum_{i < j} (H_{ij}^\text{far} + H_{ij}^\text{close})$ identical to Newtonian gravity by construction.

$K$ must be:
- $C^2$ continuous (preserves symplecticity of the far-field WH step).
- $K(0) = 0$ (close-field carries the entire force as $r \to 0$).
- $K(r) = 1$ for $r > r_\text{cross}$ (far-field carries the entire force outside changeover).
- Monotone increasing in $[0, r_\text{cross}]$.

A standard choice (REBOUND default) is the Heaviside-smoothed polynomial:

$$
K(y) = \begin{cases}
0 & y \le 0 \\
y^2 (3 - 2y) & 0 < y < 1 \\
1 & y \ge 1
\end{cases}
$$

with $y = (r - r_\text{inner}) / (r_\text{cross} - r_\text{inner})$. Other choices (smoother polynomials, exponential transitions) appear in the literature.

### 3.2 Changeover scale

$r_\text{cross}$ is set per pair and depends on Hill radius:

$$
r_\text{cross}^{(ij)} = \alpha \cdot r_H^{(ij)} = \alpha \cdot a_{ij} \, \sqrt[3]{(m_i + m_j) / (3 M_\star)}
$$

with $\alpha \in [1, 10]$ a configurable parameter (REBOUND default: $\alpha = 3$). For non-hierarchical systems where there is no obvious "central body" $M_\star$, an alternative scale (closest-pair separation, system-mean Hill radius, etc.) is needed; this is one of the open design questions in §6.

### 3.3 Step integration

Per Rein et al. §2:

1. Drift bodies under $H_\text{Kepler}^\text{far}$ for $\Delta t / 2$ (analytical Kepler step on far-field forces only).
2. Kick bodies under $H^\text{close}$ for $\Delta t$ via IAS15 sub-integration (high-order, adaptive within the Δt window).
3. Drift bodies under $H_\text{Kepler}^\text{far}$ for $\Delta t / 2$.

The IAS15 sub-step in (2) integrates only the close-field forces — typically zero except for pairs near or below $r_\text{cross}$. The cost saving is precisely this: when no pair is in encounter, (2) is a near-no-op; when pairs are encountering, the cost is local to those pairs.

Symplecticity: the far-field steps (1) and (3) are symplectic by construction (analytical Kepler in canonical coordinates). The close-field step (2) is not symplectic, but its accumulated round-off is bounded by IAS15's published behaviour. The smoothness of $K$ ensures the splitting error is $O(\Delta t^3)$ at the boundary, preserving the second-order character of the overall map.

## 4. Phased rollout

Implementation proceeds in four phases, each independently testable and reviewable.

### Phase 1 — Diagnostic surfacing (foundation)

- Add `close_encounter_threshold: Option<f64>` to `PhysicsConfig`.
- Extend `compute_closeness()` to compare `r_min` against the threshold and return an `EncounterFlag` enum (`Far`, `Approaching`, `Close`).
- Emit a `warn_diag!` event on `Close` transitions.
- No behavioural change yet; observability only.
- Tests: encounter detected when threshold crossed; not detected otherwise; threshold default of `None` preserves current behaviour bit-for-bit.

Acceptance: existing test suite passes unchanged; new tests cover the three flag transitions; lab-notebook demo on a Sun + Mercury + scattering body scenario.

### Phase 2 — Hard-switch hybrid integrator

- Create a `HybridIntegrator` type that wraps two integrators (`fast: Box<dyn Integrator>`, `precise: Box<dyn Integrator>`).
- Switch globally between them based on the `EncounterFlag` from Phase 1.
- Hard switch: when `Close` flag transitions on, swap to `precise`; when off, swap to `fast`.
- Document explicitly that hard switching breaks symplecticity at the transition; this is a stepping stone, not the target.
- Tests: switch fires on encounter, energy drift visible at switch boundary (acceptable here because Phase 3 fixes it), no regression on non-encountering scenarios.

Acceptance: hybrid pair `(Yoshida4, IAS15)` demonstrably swaps on a scattering scenario; lab-notebook documents the energy-drift signature at the swap boundary as expected and to-be-removed in Phase 3.

### Phase 3 — Smooth changeover (MERCURIUS proper)

- Replace the hard switch with the two-zone Hamiltonian decomposition.
- Implement the $K$ changeover function.
- Compute Hill radius per pair and derive $r_\text{cross}^{(ij)}$.
- Refactor `HybridIntegrator::step` to:
  1. Compute far-field and close-field decompositions.
  2. Far drift via `fast` integrator on far-field only.
  3. Close kick via `precise` integrator on close-field only.
  4. Far drift on far-field only.
- Tests: swap-boundary energy drift gone (or below f64 floor); long-horizon symplectic behaviour preserved; cross-implementation parity entry against MERCURIUS in REBOUND on a planetary-encounter scenario.

Acceptance: Phase 3 lab-notebook entry against REBOUND `MERCURIUS` integrator on a Solar-System scattering scenario, gating on energy and angular momentum at $\sim 10^{-9}$ over $10^4$ orbits.

### Phase 4 — Configuration surface and integration with existing integrators

- Expose `HybridIntegrator` through `IntegratorKind`.
- Decide on configuration knobs: $\alpha$ (Hill-radius multiplier), changeover-function shape, Hill-radius computation strategy for non-hierarchical systems.
- Update `paper.md` §Design and validation if MERCURIUS-class hybrid is positioned as a v0.1 feature; otherwise ship as v0.2 entry.
- Update `apsis-py` to expose hybrid configuration.
- Update `apsis-app` composition panel (see Issue #30) to surface the hybrid mode.

## 5. Dependencies

- **TD-008 (WH refactor, Issue #28).** Phase 3 ideally targets a WHFast-class symplectic integrator as the `fast` half of the hybrid. With the current WH carrying four documented defects, Phase 3 should either land after the WH refactor or use Yoshida-4 as the symplectic stand-in (lower precision but symplectic). This is a sequencing decision, not a hard block — Phase 1 and Phase 2 can land independently.

- **Long-horizon Mercury 1PN (Issue #29).** Phase 3 results would inform whether the long-horizon Mercury reproducibility can use the hybrid path (lower per-step cost) without losing precision.

- **3D port (already shipped).** Hill-radius and changeover are 3D-aware by construction; the proposal assumes the `Vec3` value type and 3D physics observables already in place.

## 6. Open design questions

Items deliberately not decided in this proposal; flagged for the implementing PR to resolve.

1. **Hill radius for non-hierarchical systems.** What is $M_\star$ when there is no dominant central body? Options:
   - System total mass $\sum m_k$ (overestimates for non-hierarchical).
   - Pair-local: $M_\star = m_i + m_j$ in $r_H^{(ij)}$ formula.
   - Global maximum: largest single mass in the system.
   - Refuse: only support hybrid mode on systems with a clearly dominant body.

2. **Changeover function shape.** Default $C^2$ polynomial works; smoother $C^\infty$ alternatives exist (Rein et al. §2.2 mention) but cost more in evaluation. Choose conservatively or expose as a configuration knob.

3. **Default value of $\alpha$.** REBOUND default is 3 Hill radii. Validate against parity scenarios before pinning.

4. **`fast` integrator selection at Phase 1–2.** If the WH refactor has not landed, use Yoshida-4 as the symplectic stand-in for the fast half. If WH has been refactored, prefer it. Document both pathways.

5. **Per-pair vs global changeover.** Rein et al. compute changeover per pair. Implementing per-pair adds bookkeeping cost; validate that this is necessary vs a single global changeover threshold. Likely necessary for correctness; flagged for explicit verification.

6. **Step-size behaviour.** Does the outer step `dt` need to shrink near close encounters, or does the inner IAS15 sub-integration absorb the cost? Rein et al. argue the inner sub-integration is sufficient; verify on apsis-side scenarios.

## 7. Acceptance criteria for the full proposal

- [ ] Phase 1 tests passing; encounter flag transitions visible on a scattering scenario.
- [ ] Phase 2 hard-switch demonstrable, with lab-notebook documenting the swap-boundary energy drift.
- [ ] Phase 3 smooth changeover landing the swap-boundary drift below the f64 floor; cross-implementation parity entry against REBOUND `MERCURIUS` on a planetary-encounter scenario.
- [ ] Phase 4 configuration surface stable; `paper.md` §Design and validation updated if shipped as v0.1 feature.
- [ ] All open questions in §6 resolved with explicit rationale, recorded in the implementing PR's commits.

## 8. Out of scope for this proposal

- Full WHFast reimplementation (Issue #28, TD-008).
- Adaptive timestep coupled to changeover-function gradient.
- GPU or SIMD acceleration of the close-field IAS15 sub-step.
- Backwards-compatibility with non-hybrid integrator configurations after this lands — those continue to work, just without hybrid behaviour.
- Pulsar-timing-class perturbations or other federation extensions; this proposal modifies the integrator stack, not the perturbation surface.

## 9. References

- Rein, H., Hernandez, D. M., Tamayo, D., & Brown, G. (2019). *Hybrid symplectic integrators for planetary dynamics.* MNRAS, 489(4), 4632–4640.
- Rein, H., & Tamayo, D. (2015). *WHFast: a fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long-term gravitational simulations.* MNRAS, 452(1), 376–388.
- Wisdom, J., & Holman, M. (1991). *Symplectic maps for the n-body problem.* AJ, 102, 1528.
- Chambers, J. E. (1999). *A hybrid symplectic integrator that permits close encounters between massive bodies.* MNRAS, 304(4), 793–799 — original hybrid-integrator-with-close-encounter paper that MERCURIUS refines.
- Issue #28 (TD-008: WH refactor) — sequencing dependency for Phase 3.
- Issue #29 (Long-horizon Mercury 1PN) — informed by hybrid path performance.
- Issue #31 (Barnes-Hut octree) — orthogonal scaling work; hybrid integrator and octree compose at the force-model layer.
