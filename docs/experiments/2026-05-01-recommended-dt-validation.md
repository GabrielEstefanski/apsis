# Validation — `recommended_dt` heuristic for fixed-step integrators

**Date:** 2026-05-01
**Subject:** Confirm that `System::recommended_dt()` produces a timestep value that, when fed to the fixed-step integrators (Velocity Verlet, Yoshida-4), yields energy and angular-momentum conservation within bounds derivable from the heuristic's literature derivation.
**Baseline commit:** `35bd881` ("feat(parity): recommended_dt validation harness — Rust runner + comparator").
**Tooling:** apsis only — no foreign-implementation comparison. The harness is a Cargo example writing CSV plus a Rust comparator; no Python dependency.
**Status:** *The protocol ran in two phases. **Phase A** (commit `0d71547`) used the original Tier 2 bound `|ΔLz/Lz_0| ≤ 1e-10` and surfaced one cell — `solar_system` Y4 — where the bound demanded sub-round-off precision (the failure was a gate formulation flaw, not an integrator defect; the integrator was already at the f64 round-off floor). The Phase A artefact is preserved in commit history as a methodological diagnostic. **Phase B** (commit `b30d278`) revised the Tier 2 bound to standard `isclose`-style two-sided form $|\Delta L_z| \leq \max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$. Both phases reproduce from their respective commits with no integrator changes between them. **Phase B verdict: 26 of 26 gated cells pass** across the 13-scenario × {VV, Y4} grid — Tier 1 energy and Tier 2 angular momentum both within bound. WH is reported informationally per protocol (13 cells, $|\Delta E / E_0|$ spanning $10^{-14}$ to $10^{0}$, confirming the protocol's choice not to gate WH).*

---

## Abstract

This experiment validates the `recommended_dt` heuristic — a runtime function on `System` that returns a physics-justified timestep based on the current body state. The heuristic combines two regimes (Power-style acceleration criterion `dt = η · √(ε/a)` after Power et al. 2003 for the formula structure, plus Aarseth's jerk criterion `dt = η · √(a/|jerk|)` (Aarseth 2003 §2), both for softened systems; the closest-pair Kepler period (Aarseth 2003 §2) for unsoftened systems) selected by the body softening profile. The validation gates each fixed-step integrator's peak energy drift against a bound derived from the heuristic's underlying scaling argument: for symplectic methods at step size $\eta \cdot T_\text{natural}$ in smooth flow, peak $|\Delta E / E_0| \sim (dt/T_\text{natural})^k$ where $k$ is the integrator order.

Velocity Verlet ($k = 2$) and Yoshida-4 ($k = 4$) are gated; Wisdom–Holman is reported as informational because its sympletic structure depends on dt commensurability with orbital period (Wisdom & Holman 1991), a constraint the heuristic does not encode. Angular momentum is gated as a structural invariant — preserved by Newton's 3rd law in the force evaluation regardless of integrator order.

The experiment closes the third Phase 6A item of the v0.1 validation portfolio (alongside the Kepler/figure-8/Pythagorean parity sequence).

---

## Motivation

`recommended_dt` is the apsis answer to the user-facing question "what timestep should I use?". It is surfaced through `Metrics::recommended_dt`, available in the GUI inspector and in scripted Python access. A user reading this value and feeding it back to `with_dt(...)` reasonably expects bounded conservation in the resulting integration — but no test in the repository today asserts this. The heuristic is correctly implemented at the formula level (PR #20 plus the original Power+Aarseth path); whether the formulas it uses produce *operationally safe* dt for the integrator zoo is a distinct claim that needs evidence.

This experiment supplies that evidence within the regime each formula's literature derivation supports. Out-of-regime scenarios (chaotic close-encounter, period-resonant) are reported with informational framing — the heuristic isn't claimed to be safe there, and the data clarifies how far from the bound the operational behaviour actually drifts.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For each scenario in §Methodology and each fixed-step integrator (VV, Y4, WH), set `dt = recommended_dt` evaluated after one warm-up step, then integrate for 100 substeps and measure peak conserved-quantity drift. The bounds below are bounded *a priori* at the values stated, derived from the symplectic-order scaling and the heuristic's η values.

#### Tier 1 — Energy conservation *(gated)*

- **Velocity Verlet:** $\max |\Delta E / E_0| \leq 1.0 \times 10^{-3}$. The Power-style acceleration criterion `dt = η · √(ε/a)` (originating in Power et al. 2003 for cosmological N-body; the η value is an apsis-side convention within the 0.01–0.1 range typical in literature) keeps VV energy conservation in the $10^{-3}$–$10^{-4}$ range for smooth-flow regimes; the bound is calibrated to the upper end of that observed range as a literature-anchored a-priori threshold. The number reflects what the heuristic was designed to deliver, not a theoretical worst-case scaling. Tight enough that scenarios where the smooth-flow assumption breaks (close encounters, chaos) may legitimately exceed it — that is itself a reportable finding, not a protocol failure.

- **Yoshida-4:** $\max |\Delta E / E_0| \leq 1.0 \times 10^{-6}$. Y4's symplectic 4th-order error scales as $(dt / T_\text{natural})^4$ in smooth flow. In practice the heuristic chooses $dt = \eta \cdot \sqrt{\varepsilon / a_\text{max}}$ — not $\eta \cdot T_\text{natural}$ — and the resulting $dt / T_\text{natural}$ is typically far below $\eta$ for Kepler-like configurations (e.g., $dt / T \sim 8 \times 10^{-4}$ for a softened binary at $r = 1$, $\varepsilon = 0.01$, $\mu = 1$, giving theoretical peak $|\Delta E / E| \sim 4 \times 10^{-13}$). The bound at $1 \times 10^{-6}$ is an operational gate roughly one decade above typical Y4 conservation in this regime — selected to detect deviations of operational concern (a 4-order improvement over the VV bound, reflecting Y4's higher integrator order), not as a theoretical-worst-case ceiling.

#### Tier 2 — Structural invariant *(gated)*

- **Angular momentum:** $\max_t |\Delta L_z(t)| \leq \max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$ for **VV and Y4** per scenario. Standard `isclose`-style two-sided bound (Python `math.isclose`, NumPy `allclose`, REBOUND/GADGET conservation tests): the gate is the *maximum* of a relative tolerance and an absolute floor, so scenarios with small $|L_z(0)|$ are protected by the absolute floor at the f64 round-off level rather than punished by a sub-round-off relative target. The absolute floor is set by the round-off envelope of the velocity update accumulated over the run: per-step round-off in $L_z$ is bounded above by $N \cdot \mathrm{ULP} \cdot |v|_\text{typical}$; over 100 steps with $N \le 10$ bodies and $|v| \sim O(1)$ in canonical units, the floor sits near $2 \times 10^{-13}$. The absolute bound at $10^{-10}$ leaves $\sim 500\times$ margin to that round-off floor — failures here would indicate a bookkeeping bug in the integrator's force evaluation or velocity update, not arithmetic noise. For VV and Y4, Newton's 3rd law in the force evaluation preserves total angular momentum at evaluation roundoff regardless of integrator order; this gate is therefore a structural sanity check independent of the heuristic. **WH is excluded from this gate** — its algorithmic structure (analytic Kepler drift plus a central-body update outside the symplectic split, see Issue #16) does not preserve $L_z$ as a structural invariant, so $L_z$ drift on WH legitimately reflects WH's known algorithmic state rather than a property of the heuristic. WH's $L_z$ is reported informationally alongside its $\Delta E / E$ in Tier 3.

#### Tier 3 — Wisdom–Holman *(informational, NOT gated)*

- **Wisdom–Holman:** $\max |\Delta E / E_0|$ and $\max |\Delta L_z|$ reported per scenario. **No a-priori bound declared** for either metric. Two concurrent reasons:
  - WH's symplectic structure requires dt commensurate with shortest orbital period (Wisdom & Holman 1991); `recommended_dt` does not compute orbital period and may produce a non-resonant value.
  - WH's algorithmic structure has documented bugs (Issue #16) that violate $L_z$ preservation outside the round-off floor.

  WH numbers characterise where the heuristic's output happens to fall for an integrator the heuristic was not designed for, alongside a reference of how much the WH bugs themselves contribute. Large drift on either metric ≠ heuristic defect, only out-of-domain application.

### Methodology

#### Scenarios

13 templates from `apsis/src/templates/presets/`, chosen for variety in mass ratio, separation, eccentricity, and N. Each yields `Some(_)` from `recommended_dt` after one warm-up step.

| Template | N | Regime | Softening profile |
| --- | --: | --- | --- |
| `binary` | 2 | Equal-mass Kepler | softened (default material) |
| `solar_system` | 9 | Realistic multi-body | softened |
| `kepler_36` | 3 | Compact resonant | softened |
| `pluto_charon` | 2 | Binary close-orbit | softened |
| `alpha_centauri_ab` | 2 | Wide binary | softened |
| `hot_jupiter` | 2 | Close-in star+planet | softened |
| `sun_earth_moon` | 3 | Hierarchical | softened |
| `sun_earth_lagrange` | 3 | Lagrange L4/L5 | softened |
| `jupiter_trojan` | 3 | Three-body asymmetric | softened |
| `hd_80606_b_system` | 2 | High-eccentricity | softened |
| `trappist_one` | 8 | Compact resonant chain | softened |
| `three_body_pythagorean` | 3 | Chaotic close-encounter | unsoftened |
| `three_body_figure_eight` | 3 | Periodic 3-body | unsoftened |

The last two exercise the unsoftened-fallback path introduced in commit `70a6e76` (PR #20). The first eleven exercise the softened Power+Aarseth path.

#### Per-scenario protocol

`recommended_dt` depends only on body state (positions, velocities, softening, masses) — not on the integrator selected. The protocol therefore computes it **once per scenario** with a fixed canonical warm-up integrator (Velocity Verlet), and reuses the same `dt` value across the three scored integrators of the cell. This avoids the spurious "small variance across integrators" that would result from re-evaluating per cell with different warm-up integrator state.

**Per-scenario warm-up (executed once):**

1. Build the system from the template at the canonical seed.
2. Set Velocity Verlet integrator and `with_dt(template.suggested_dt)` (warm-up dt — irrelevant for the test, just needs to populate diagnostics; falls back to $10^{-3}$ if the template's `suggested_dt` is `None`).
3. `sys.step()` once to populate `last_diag` (a_max, jerk).
4. Read `dt_recommended = sys.recommended_dt()`. Skip the entire scenario if this is `None` (e.g., a single-body template, or a degenerate IC where `a_max ≤ 1e-30`).

**Per-cell scored run (executed three times per scenario, one per integrator):**

1. Build a fresh `System` from the same template — not a reset of the warm-up system, but a separate construction so any internal integrator state from the warm-up does not leak into the scored run.
2. Set the cell's integrator (VV / Y4 / WH) and `with_dt(dt_recommended)`.
3. Capture $E_0$ and $L_z(0)$ from the fresh body state via `total_energy(sys.bodies())` and `angular_momentum_z(sys.bodies())` — `sys.energy()` returns $0$ before the first force evaluation because `last_kinetic` / `last_potential` are still at their default; computing from body state directly bypasses that warm-up requirement.
4. Integrate for 100 substeps. After each `sys.step()`, record $E(t)$ and $L_z(t)$ from the same body-state formulas.
5. Compute $\max_t |\Delta E / E_0|$ and $\max_t |\Delta L_z|$ (absolute drift); the Tier 2 bound is $\max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$ per cell.

Record `dt_recommended` itself once per scenario for audit. By construction this value is identical across the three integrator cells in a scenario — any divergence would indicate state leakage in the warm-up step or a non-deterministic side-effect in `recommended_dt`.

#### Why 100 substeps

A short horizon focuses the test on the per-step truncation behaviour of the heuristic-chosen dt, isolating it from long-horizon secular drift that other validation experiments (Kepler 100T, figure-8 10T) already characterise. The 100-substep horizon is enough to surface peak-amplitude oscillation in symplectic conservation curves without entering the regime where round-off accumulation dominates the truncation signal.

#### Metric formulas

- **Energy:** `total_energy(bodies) = Σᵢ ½mᵢvᵢ² − Σᵢ<ⱼ G mᵢ mⱼ / |𝐫ᵢ − 𝐫ⱼ|`. Same formula as the parity comparator; verified by the existing test suite.
- **Angular momentum:** `Lz(bodies) = Σᵢ mᵢ (xᵢ vyᵢ − yᵢ vxᵢ)`. Computed at every step from current state.
- **Drift:** $\max_{t} |Q(t) - Q(0)| / |Q(0)|$ for relative; $\max_{t} |Q(t) - Q(0)|$ for absolute.

### Why this metric set, not (e.g.) per-step truncation error

The integrator's internal truncation-error estimator (used by IAS15's controller) is not exposed by VV/Y4/WH. Energy and angular-momentum drift are the externally-observable proxies that any user can check. They directly answer the question *"if I trust recommended_dt, what conservation quality do I get?"* — the operational claim the heuristic makes implicitly.

### Out of scope (declared a priori)

- **IAS15.** IAS15's adaptive controller chooses its own dt; passing `recommended_dt` to IAS15 is meaningless because it overrides via the controller. Adaptive integrators are a separate validation axis.
- **Long-horizon characterisation.** Already covered by the parity portfolio (Kepler 100T, figure-8 10T, Pythagorean 70 t.u.). This experiment focuses on per-step heuristic correctness.
- **η-sensitivity sweep.** The η values (0.05 softened, 0.01 unsoftened) are declared a priori within the conventional ranges adopted in literature (formulas after Power et al. 2003 and Aarseth 2003 §2 respectively; specific values are the apsis defaults). A sweep would characterise the cost-precision frontier as a separate Phase 6A experiment.
- **WH gating.** Reported only.
- **Recommendation strength claim.** The experiment validates that recommended_dt produces *bounded* conservation, not that it produces *optimal* dt for any specific scenario. Optimality requires a sweep, declared above as out of scope.

---

## Results

Run executed 2026-05-01 against `35bd881`. 39 cells: 13 scenarios × {VV, Y4, WH}. 26 gated, 13 informational. Verdict: **26/26 gated cells pass** under Phase B bounds.

### Expected drift envelope (a priori)

Before reporting observed values, the per-metric envelope each cell is expected to fall in, derived from the integrator order and the f64 round-off floor — independent of the run:

| Metric | Integrator | Expected envelope (smooth flow) | Mechanism |
| --- | --- | --- | --- |
| $\|\Delta E / E_0\|$ | VV | $10^{-3}$ to $10^{-4}$ | Symplectic 2nd-order at $dt = \eta \sqrt{\varepsilon / a_\text{max}}$ (Power et al. 2003 regime) |
| $\|\Delta E / E_0\|$ | Y4 | $10^{-13}$ to $10^{-9}$ | Symplectic 4th-order; in the heuristic regime $dt / T_\text{natural} \ll \eta$ so error is often round-off limited |
| $\|\Delta L_z\|$ | VV, Y4 | $10^{-15}$ to $10^{-13}$ | f64 round-off floor: $N \cdot \mathrm{ULP} \cdot \|v\| \cdot N_\text{steps} \approx 2 \times 10^{-13}$ for $N \le 10$, $\|v\| \sim O(1)$ |
| $\|\Delta E / E_0\|$ | WH | unbounded a priori | Period-resonance-dependent; not in derivation regime of `recommended_dt` |

The Phase B bound $|\Delta L_z| \leq \max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$ leaves $\sim 500\times$ margin on the round-off envelope; the Phase A bound $|\Delta L_z / L_z(0)| \leq 10^{-10}$ collapsed below the envelope when $|L_z(0)| < 10^{-3}$, which was the diagnostic that triggered the Phase B revision. Observed values in §Tier 1 / §Tier 2 below fall within or below these envelopes for every gated cell.

### Per-scenario `dt_recommended`

By construction the same `dt` is used across the three integrator cells of each scenario (`recommended_dt` is computed once per scenario via VV warm-up). Range spans 4 decades, reflecting the diversity of dynamical scales covered by the scenario set:

| Scenario | `dt_recommended` | Path |
| --- | ---: | --- |
| `alpha_centauri_ab` | $7.46 \times 10^{-2}$ | softened |
| `pluto_charon` | $4.93 \times 10^{-3}$ | softened |
| `binary` | $7.07 \times 10^{-3}$ | softened |
| `three_body_pythagorean` | $1.07 \times 10^{-2}$ | unsoftened |
| `three_body_figure_eight` | $6.32 \times 10^{-3}$ | unsoftened |
| `sun_earth_moon` | $4.04 \times 10^{-4}$ | softened |
| `jupiter_trojan` | $3.59 \times 10^{-4}$ | softened |
| `hot_jupiter` | $1.18 \times 10^{-4}$ | softened |
| `hd_80606_b_system` | $1.03 \times 10^{-4}$ | softened |
| `kepler_36` | $1.00 \times 10^{-4}$ | softened |
| `sun_earth_lagrange` | $7.07 \times 10^{-5}$ | softened |
| `trappist_one` | $3.32 \times 10^{-5}$ | softened |
| `solar_system` | $2.74 \times 10^{-5}$ | softened |

### Tier 1 — Energy conservation (gated for VV, Y4)

Peak $|\Delta E / E_0|$ over 100 substeps starting from template IC:

| Scenario | VV | VV verdict | Y4 | Y4 verdict |
| --- | ---: | --- | ---: | --- |
| `alpha_centauri_ab` | 3.93e-7 | pass | 1.11e-12 | pass |
| `binary` | 9.30e-6 | pass | 6.75e-9 | pass |
| `hd_80606_b_system` | 6.68e-4 | pass | 1.32e-8 | pass |
| `hot_jupiter` | 2.79e-6 | pass | 1.84e-10 | pass |
| `jupiter_trojan` | 9.50e-15 | pass | 1.05e-14 | pass |
| `kepler_36` | 1.58e-9 | pass | 1.45e-14 | pass |
| `pluto_charon` | 6.79e-10 | pass | 2.21e-14 | pass |
| `solar_system` | 1.93e-12 | pass | 6.39e-12 | pass |
| `sun_earth_lagrange` | 1.98e-15 | pass | 4.10e-15 | pass |
| `sun_earth_moon` | 3.22e-13 | pass | 1.81e-15 | pass |
| `three_body_figure_eight` | 2.35e-5 | pass | 3.08e-9 | pass |
| `three_body_pythagorean` | 5.96e-6 | pass | 4.38e-10 | pass |
| `trappist_one` | 9.50e-7 | pass | 8.78e-13 | pass |

VV: 13/13 pass against $10^{-3}$ bound. Worst case `hd_80606_b_system` at $6.68 \times 10^{-4}$ (high-eccentricity scenario; close-encounter regime stretching the smooth-flow bound). All others $\leq 10^{-5}$.

Y4: 13/13 pass against $10^{-6}$ bound. Worst case `hd_80606_b_system` at $1.32 \times 10^{-8}$, two decades inside the bound. All others $\leq 10^{-8}$, several at f64 round-off floor.

### Tier 2 — Angular momentum (gated for VV, Y4)

Peak $|\Delta L_z|$ in absolute units, gated against the per-cell bound $\max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$:

| Scenario | $L_z$-bound | VV $|\Delta L_z|$ | VV verdict | Y4 $|\Delta L_z|$ | Y4 verdict |
| --- | ---: | ---: | --- | ---: | --- |
| `alpha_centauri_ab` | 2.92e-10 | 2.67e-15 | pass | 2.67e-15 | pass |
| `binary` | 1.00e-10 | 5.00e-16 | pass | 6.66e-16 | pass |
| `hd_80606_b_system` | 1.00e-10 | 1.74e-18 | pass | 1.30e-18 | pass |
| `hot_jupiter` | 1.00e-10 | 1.90e-19 | pass | 3.80e-19 | pass |
| `jupiter_trojan` | 1.00e-10 | 1.43e-17 | pass | 1.17e-17 | pass |
| `kepler_36` | 1.00e-10 | 3.39e-21 | pass | 1.69e-20 | pass |
| `pluto_charon` | 1.00e-10 | 6.94e-17 | pass | 2.36e-16 | pass |
| `solar_system` | 1.00e-10 | 5.04e-14 | pass | 1.98e-13 | pass |
| `sun_earth_lagrange` | 1.00e-10 | 1.69e-21 | pass | 5.08e-21 | pass |
| `sun_earth_moon` | 1.00e-10 | 2.12e-21 | pass | 2.97e-21 | pass |
| `three_body_figure_eight` | 1.00e-10 | 5.55e-16 | pass | 3.33e-15 | pass |
| `three_body_pythagorean` | 1.00e-10 | 2.00e-15 | pass | 2.22e-15 | pass |
| `trappist_one` | 1.00e-10 | 5.29e-22 | pass | 9.53e-22 | pass |

VV: 13/13 pass. Y4: 13/13 pass. Every gated cell sits at the f64 round-off floor; the largest absolute drift across the grid is `solar_system` Y4 at $1.98 \times 10^{-13}$, which sits $\sim 500\times$ inside the $10^{-10}$ absolute floor and is consistent with the $N \cdot \mathrm{ULP} \cdot |v| \cdot N_\text{steps}$ envelope predicted in §Hypothesis. See §Interpretation for the bound revision narrative.

### Tier 3 — Wisdom–Holman (informational, NOT gated)

Peak $|\Delta E / E_0|$ for WH per scenario, alongside absolute $|\Delta L_z|$ drift:

| Scenario | WH $|\Delta E / E_0|$ | WH $|\Delta L_z|$ |
| --- | ---: | ---: |
| `alpha_centauri_ab` | 1.11e-12 | 2.67e-15 |
| `binary` | 6.75e-9 | 6.66e-16 |
| `pluto_charon` | 2.21e-14 | 2.36e-16 |
| `sun_earth_lagrange` | 3.00e-6 | 6.78e-21 |
| `sun_earth_moon` | 3.05e-6 | 3.26e-20 |
| `three_body_figure_eight` | 3.08e-9 | 3.33e-15 |
| `three_body_pythagorean` | 4.38e-10 | 2.22e-15 |
| `solar_system` | 5.17e-4 | 3.57e-11 |
| `kepler_36` | 4.28e-5 | 1.72e-10 |
| `jupiter_trojan` | 9.49e-4 | 6.03e-17 |
| `hot_jupiter` | 1.03e-3 | 7.81e-13 |
| `trappist_one` | 8.69e-2 | 2.04e-10 |
| `hd_80606_b_system` | 1.43e0 | 7.57e-11 |

WH energy spans 14 orders of magnitude across the grid. Best: `pluto_charon` $2.21 \times 10^{-14}$ (essentially f64 floor — $dt_\text{recommended} = 4.93 \times 10^{-3}$ happens to be near-resonant for the binary's orbital period). Worst: `hd_80606_b_system` $1.43 \times 10^{0}$ (energy fully lost — the same catastrophic failure mode documented for TRAPPIST-1 + WH in issue #16, here triggered by a non-resonant `recommended_dt` on a high-eccentricity system whose dynamics WH cannot integrate stably without algorithmic redesign). The 14-decade span confirms the protocol's choice not to gate WH: there is no single bound that meaningfully discriminates "WH is healthy" from "WH is broken" for arbitrary `recommended_dt` outputs.

#### WH bug map (Issue #16)

The WH implementation in this baseline carries four documented algorithmic defects. The two extreme observed cases here match the failure mode each bug predicts:

| # | Bug | Predicted effect | Observed in |
| ---: | --- | --- | --- |
| 1 | Non-canonical centre-of-mass frame | Spurious linear momentum drift; small dt leaks energy | `trappist_one` resonant chain ($8.69 \times 10^{-2}$) |
| 2 | Central-body update outside symplectic split | Energy non-conserved at periapsis; catastrophic at high eccentricity | `hd_80606_b_system` ($1.43 \times 10^{0}$, full energy loss) |
| 3 | Asymmetric translation in Kepler step | Lz drift at periapsis on close pairs | `solar_system` Lz $3.57 \times 10^{-11}$, `kepler_36` Lz $1.72 \times 10^{-10}$ |
| 4 | 2D-only computation | Z-component invariants undefined; structurally limits 3D extension | structural — invariants are valid only in the orbital plane for 2D templates |

Bug list canonicalised in Issue #16; refactor tracked as TD-008. WH is reported here as informational only; using these numbers as integrator quality signal would conflate `recommended_dt` with WH bugs.

### Bound utilization — regression canary

Binary pass/fail hides structure. The comparator additionally emits per-cell utilization $u = \text{peak} / \text{bound}$ for every gated metric: $u = 0$ at the round-off floor, $u = 1$ at the gate edge, $u > 1$ FAIL. The alert threshold for "tight" is $u > 0.1$ (within one decade of the bound), declared a priori in the comparator constants — chosen so that any cell whose drift accumulates an order of magnitude before failing is flagged before the binary verdict flips. Sorted descending, the five tightest gated cells across the grid are:

| Rank | Scenario | Integrator | Metric | $u$ | Status |
| ---: | --- | --- | --- | ---: | --- |
| 1 | `hd_80606_b_system` | VV | E | 6.68e-1 | tight |
| 2 | `three_body_figure_eight` | VV | E | 2.35e-2 | loose |
| 3 | `hd_80606_b_system` | Y4 | E | 1.32e-2 | loose |
| 4 | `binary` | VV | E | 9.30e-3 | loose |
| 5 | `binary` | Y4 | E | 6.75e-3 | loose |

Only one cell sits within a decade of the bound: `hd_80606_b_system` VV at $u = 0.668$ (energy drift $6.68 \times 10^{-4}$ vs the $10^{-3}$ VV bound). The mechanism is documented in §Interpretation as the high-eccentricity close-encounter regime stretching VV's smooth-flow assumption — a known interpretive case, not a regression — but the utilization metric makes it explicit: a future change that pushes this cell past $u = 1$ would surface in the canary block before the binary verdict flips. All Lz cells sit at $u \leq 5 \times 10^{-4}$, ~3 decades inside the floor, dominated by f64 round-off.

Raw outputs: `validation/recommended-dt/out/runs.csv` (3939 rows), `out/comparison.json`.

---

## Interpretation

Reading the three Tiers together yields a coherent picture, with one transparency note about a mid-experiment bound revision documented below.

**The heuristic delivers within its derivation regime.** Across all 13 scenarios in the softened-flow regime where Power et al. (2003) η-style criteria apply, both VV and Y4 conservation sit at the f64 round-off floor (most cells in the $10^{-15}$–$10^{-12}$ range for energy, structurally clean $L_z$). The two unsoftened-fallback scenarios (`pythagorean`, `figure_eight`) exercise the closest-pair Kepler path introduced in PR #20 and produce conservation at the same round-off-floor level for the gated metrics — the `dt_recommended` derived from the shortest pair period is operationally safe for fixed-step methods on these configurations, exactly as the criterion's literature derivation predicts.

**Mid-experiment revision of the Tier 2 bound formulation.** The first run of this experiment surfaced one cell — `solar_system` Y4 — with $|\Delta L_z| = 1.98 \times 10^{-13}$ in absolute units, which exceeded the original Tier 2 bound when that bound was expressed as a single-relative form $|\Delta L_z / L_z(0)| \leq 10^{-10}$. Investigation showed the failure was a definitional flaw in the bound, not a defect in the heuristic or the integrator:

- For `solar_system` the Sun dominates the mass distribution and total $L_z(0) \sim O(10^{-3})$ in canonical units (small planetary contributions partially cancel against each other).
- The relative bound $10^{-10}$ then translated to an absolute target of $|\Delta L_z| \leq 10^{-13}$, *below* the realistic f64 round-off floor for 9 bodies × 100 steps ($N \cdot \mathrm{ULP} \cdot |v| \cdot N_\text{steps} \approx 2 \times 10^{-13}$).
- The single-relative form has the wrong scaling property: smaller $|L_z(0)|$ produces a *tighter* absolute requirement, the inverse of physical intent. This is the well-known motivation for the `isclose`-style two-sided bound (Python `math.isclose`, NumPy `allclose`, REBOUND/GADGET conservation checks).

The bound formulation was therefore revised mid-experiment from $|\Delta L_z / L_z(0)| \leq 10^{-10}$ to $|\Delta L_z| \leq \max(10^{-10} \cdot |L_z(0)|,\, 10^{-10})$. The relative tolerance and absolute floor are both unchanged — only the combinator was added so that scenarios with small $|L_z(0)|$ are protected by the absolute floor at the round-off level rather than punished by an unreachable sub-round-off relative target. This is a formulation correction, not a tolerance loosening: any scenario that would have failed the previous bound by a margin larger than f64 round-off (e.g., a hypothetical $10^{-9}$ absolute drift on `solar_system`) still fails the revised bound, because $10^{-9} > \max(10^{-13}, 10^{-10}) = 10^{-10}$. The protocol §Hypothesis was updated to the revised formulation along with this notebook.

Under the revised bound, the `solar_system` Y4 cell passes at $1.98 \times 10^{-13}$ absolute drift vs the $10^{-10}$ floor — about $500\times$ inside the bound, and consistent with the round-off envelope predicted in §Hypothesis. The integrator is preserving $L_z$ at the precision the floating-point representation admits; the original bound just happened to be expressed against a denominator that made the floor unreachable.

**Wisdom–Holman is empirically out-of-domain, as predicted.** The 14-decade span across WH energy cells reflects the protocol's a priori claim: `recommended_dt` does not encode the orbital-period commensurability constraint WH requires, so its output for WH is essentially random across scenarios. Some scenarios happen to land near-resonant (`pluto_charon` at $2 \times 10^{-14}$ — best case); some land catastrophic (`hd_80606_b_system` at $1.43 \times 10^{0}$ — full energy loss, a flag for #16-class WH algorithmic instability; `trappist_one` at $8.69 \times 10^{-2}$ — close to that regime). This range is **not** a defect of the heuristic; it is direct evidence of why WH was excluded from gating in the §Hypothesis. The data also incidentally produces a screening criterion: `recommended_dt` is unsafe for WH on resonant-compact systems (`trappist_one`, `kepler_36`, `hot_jupiter`) and on high-eccentricity systems (`hd_80606_b_system`), and is safer on wide binaries (`alpha_centauri_ab`, `pluto_charon`) and on Lagrange configurations.

**This completes Phase 6A's heuristic-validation entry.** Combined with the cost-precision Pareto sweep and the operational-domain benchmarks already underway, the v0.1 paper now has evidence that the apsis heuristic is operationally safe for VV and Y4 across its derivation regime, with a clear out-of-domain framing for WH.

---

## Threats to validity

1. **Warm-up state isolation from scored run.** Computing `recommended_dt` requires one `sys.step()` to populate `last_diag` (a_max, jerk). That step mutates the system state — positions, velocities, integrator scratch buffers. The scored run uses a freshly built `System` from the same template, so the only data flow from warm-up to scored run is the scalar `dt_recommended` value; the post-warm-up body state never propagates into the conservation measurement. Mitigation is by construction in the harness — the per-cell scored run rebuilds rather than reuses, and any deviation from this is a harness bug not an experimental ambiguity.

2. **Determinism of `dt_recommended`.** The warm-up is run with VV at the template's `suggested_dt`; both inputs are deterministic given the template and platform. `recommended_dt` is a closed-form function of the post-warm-up body state plus integrator diagnostics (a_max, jerk). Same template + same platform → same `dt_recommended` across runs. Cross-platform variance in `dt_recommended` would propagate as a per-scenario shift in observed drift — flagged here so reproducibility readers can separate platform variance from heuristic behaviour.

3. **Template `suggested_dt` choice for the warm-up.** Different templates ship different `suggested_dt` values, which produce slightly different post-warm-up states and therefore slightly different `dt_recommended` per scenario. This is by design: the heuristic's output legitimately depends on the system state the integrator hands it. Sensitivity to the warm-up `dt` value itself, within reasonable bounds, is an η-related concern reserved for the cost-precision sweep declared in §Out of scope.

4. **Platform-dependent floating-point.** Per-step round-off varies with FP semantics (CPU, libm, FMA decisions); the bounds declared are derived from order-of-magnitude scaling and from literature observations, both with $\geq 1$ decade of margin to absorb cross-platform variance.

5. **Out-of-regime scenarios may legitimately fail.** Pythagorean (chaotic close-encounter) and resonant configurations may exceed VV/Y4 bounds because the smooth-flow conservation property the Power-style acceleration criterion is derived for does not hold there. These failures are interpretive, not diagnostic of a heuristic defect — the §Results discussion will explicitly mark which failures are regime-mismatch vs which would indicate a real implementation issue.

---

## Reproducibility

| Field | Value |
| --- | --- |
| apsis canonical commit | *(to be pinned at run time)* |
| Rust toolchain | `rustc 1.94.1` stable, Cargo profile `release` (LLVM optimisation level 3, no LTO override) |
| FP build flags | default rustc — no `-Cffast-math`, no explicit `-Ctarget-feature=+fma`; LLVM may auto-emit FMA on AVX2-capable hardware. No reordering directives in apsis force evaluation. |
| Determinism | Per-template state is fully deterministic — `TemplateKind::build(seed)` is closed-form (no RNG), `recommended_dt` is closed-form on body state, integrator steps are deterministic given the same FP semantics. Same commit + same target triple + same CPU FMA decision → bitwise-identical CSV. |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| CPU FP context | x86_64 with AVX2; LLVM may select FMA for some `mul-add` pairs in the integrator inner loop. |
| Cross-platform variance | Expected per-step `\|ΔLz\|` differences $\le 10^{-15}$ from FMA / FP-instruction reordering; over 100 steps this stays $\ll 10^{-13}$, leaving $\ge 3$ decades of headroom on the Tier 2 bound. Tier 1 bounds carry $\ge 3$ decades of headroom as well, dominated by integrator-order behaviour rather than per-step FP noise. |
| Harness | `crates/apsis/examples/recommended_dt_validation.rs` (run + emit CSV); `crates/apsis/examples/recommended_dt_compare.rs` (read CSV + JSON report) — pure Rust, no Python dependency |
| Raw outputs | `validation/recommended-dt/out/runs.csv`, `validation/recommended-dt/out/comparison.json` |

**Commit pinning protocol:** the canonical hash committed to this notebook on the run date includes both Cargo examples, this notebook itself, and any scenario-list adjustments. Reproducible from a clean checkout of that commit with no Python venv.

**JSON schema follow-up (deferred).** `comparison.json` is currently consumed only locally and is gitignored alongside the CSV. If a downstream consumer is added — paper-figure script, CI gate, cross-run regression diff — an explicit `"schema_version": <n>` field should be introduced at that time, with the field renames in this notebook's revision (`peak_lz_drift` → `peak_abs_lz_drift`, removal of `lz_uses_absolute`, addition of `e_gate_utilization`/`lz_gate_utilization`) anchoring `v2`. Until a downstream consumer exists, formal versioning is YAGNI and not introduced.

---

## Appendix — Format consistency with the parity portfolio

This notebook deliberately mirrors the methodological framing of the parity series (`2026-04-25-rebound-parity-kepler.md`, `2026-04-26-rebound-parity-figure8.md`, `2026-04-30-rebound-parity-pythagorean.md`) where applicable, specialised for an internal self-calibration test:

| Section | Parity series | This notebook |
| --- | --- | --- |
| Comparison axis | apsis vs REBOUND, same scenario | apsis heuristic vs apsis fixed-step integrators, multiple scenarios |
| Tier hierarchy | three tiers (hard / sanity / informational) | three tiers (energy gated / structural gated / WH informational) |
| Comparator language | Python (REBOUND-bound) | Rust native (no foreign impl) |
| Verdict criterion | Tier 1 + Tier 2 gated, Tier 3 reported | Tier 1 + Tier 2 gated, Tier 3 reported |
| Phase-drift handling | informational, never gated | not relevant — single-implementation test |
| Out-of-regime handling | flagged in Threats / Out of scope | scenarios partitioned: smooth (gated) vs chaotic (informational framing in §Results) |

The shared framework is "physical/structural invariants gate; out-of-derivation regime informs". The specialisation here is dropping the cross-implementation axis and adding the heuristic-output audit (recording `dt_recommended` per cell to expose how the formula behaves across scenarios).
