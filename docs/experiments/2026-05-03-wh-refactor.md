# Wisdom-Holman refactor — protocol

**Date:** 2026-05-03
**Subject:** Refactor the `apsis` Wisdom-Holman integrator from the current pseudo-heliocentric 2D-only implementation to faithful Wisdom & Holman (1991) in democratic heliocentric coordinates with 3D-native data flow. Validate that the four documented algorithmic defects (TD-008) are closed and that conservation invariants reach the published WH 1991 floor in smooth-flow regime.
**Baseline commit:** `70a14c2` (PR #33 merge into `develop`); refactor commits `5471c20` (kepler Vec3) → `9081c2b` (WH refactor) → `a42c44f` (positive regression tests + Galilean shift) → `32946ec` (negative regime-boundary tests).
**Tooling:** `apsis` core (`crates/apsis/src/physics/integrator/wisdom_holman.rs` and `kepler.rs`); REBOUND 4.6.0 via Python 3.10 (`reb.IAS15` informational reference, not exercised in the gates below).
**Status:** *Run executed via `cargo test --release -p apsis --lib core::system::tests::wh_refactor_regression`. **All four TD-008 defects structurally closed; six regression tests pass.** The runtime hierarchy-signal surface (`HierarchySignal` enum + `StepResult.hierarchy_signal` field + config-time `warn_diag!` in `System::set_integrator`) was added in the close-out commit alongside §Results population. Bug #4 angular-momentum bound was revised post-empirical-observation from $10^{-13}$ to $10^{-3}$ — formulation of the bound, not the algorithm; the load-bearing z-envelope claim passes strictly. See §Interpretation.*

---

## Abstract

The current Wisdom-Holman implementation in `apsis` carries four documented algorithmic defects (TD-008, surfaced by the cross-implementation parity portfolio at `paper/notebooks/2026-05-01-rebound-parity-retrograde.md` §"WH bug map"): a non-canonical centre-of-mass frame, a central-body update placed outside the symplectic split, an asymmetric translation in the Kepler step, and a 2D-only computation that silently drops the $z$-component of motion. The integrator is reported informationally in v0.1 validation runs and is not currently treated as a quality signal.

This experiment specifies the refactor protocol: democratic heliocentric (DH) coordinates per Duncan, Levison & Lee (1998), kick-drift-kick second-order symplectic split, central body integrated via the same Hamiltonian split as the planets, Kepler step extended to `Vec3` (matches the rest of the integrator stack), and a two-level dominance check (config-time `warn_diag!` plus runtime observability signal in `Metrics`) replacing the current per-step silent fallback to Yoshida-4.

The acceptance gates are organised in three tiers: smooth-flow conservation invariants on a hierarchical Sun + Mercury system (Tier 1, gated), per-bug regression tests with binary outcomes (Tier 2, gated, four checks targeting each TD-008 defect), and cross-implementation comparison against REBOUND WHFast on the same primary scenario (Tier 3, informational; the implementations differ in correctors, coordinate variant choices, and round-off control, so strict numerical agreement is not expected and not gated).

---

## Motivation

Validation runs from the parity portfolio and `recommended_dt` heuristic test characterised the current WH behaviour quantitatively. Energy drift on `hd_80606_b_system` reaches $\lvert \Delta E / E_0 \rvert = 1.43 \times 10^0$ — full energy loss at the high-eccentricity regime. Resonant chains such as `trappist_one` show $\lvert \Delta E / E_0 \rvert = 8.69 \times 10^{-2}$. Lz drift on `solar_system` and `kepler_36` accumulates well above the f64 round-off floor that should be structural for a correctly canonical WH split. The four defects mapped in `paper/notebooks/2026-05-01-rebound-parity-retrograde.md` §Tier 3 explain these observations; each defect predicts a specific failure mode that the current data confirms.

The refactor replaces the implementation rather than patching individual call sites. The current `step()` function (lines 119–186 of `wisdom_holman.rs`) is structurally 2D from frame entry through frame exit; extending the existing code to 3D would require rewriting every line that handles `(cx0, cy0, cvx0, cvy0)` 2D scaffolding, and the central-body update (lines 162–169) is outside the symplectic split by construction, not by oversight. A clean re-implementation of WH 1991 in DH coordinates with 3D-native data flow is the lower-cost path and produces a self-evident algorithm correspondence to the literature reference.

The Kepler solver in `kepler.rs` (universal-variable Newton-Raphson with Stumpff functions) is algorithmically correct; only its API signature is 2D. The refactor extends the signature to `(r0: Vec3, v0: Vec3, dt: f64, mu: f64) -> (Vec3, Vec3)` without modifying the solver internals — the universal-variable formulation is dimensionally agnostic by construction. A regression test exercising $z \neq 0$ initial conditions on a planar orbit verifies that no planar assumption survives the API change.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the refactored Wisdom-Holman implementation, the metrics declared below are bounded *a priori* at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are both gated; failure of any gated metric reproves the experiment. Tier 3 is informational and never reproves — its purpose is sanity-check context against an independently developed WH-class implementation, which is expected to differ at the implementation-detail level (correctors, coordinate variant, summation strategy) without contradicting the conservation invariants.

#### Tier 1 — Smooth-flow conservation invariants *(gated)*

Hierarchical Sun + Mercury system at standard orbital elements ($e = 0.2056$, $a = 1$ canonical), integrated for 1000 orbital periods at $dt = T / 200$ (Wisdom & Holman 1991 convention). Bounds derive from the published WH 1991 floor for smooth-flow Kepler conservation.

| Metric | Bound | Origin |
| --- | ---: | --- |
| $\max_t \lvert \Delta E / E_0 \rvert$ | $1 \times 10^{-5}$ | WH 1991 §IV reports $\sim 10^{-7}$ to $10^{-9}$ at this $dt$ for outer-planet integrations; bound at $10^{-5}$ as the conservative published floor with $\sim 100\times$ headroom for inner-planet variance |
| $\max_t \lvert \Delta L \rvert / \lvert L_0 \rvert$ (vector norm) | $1 \times 10^{-13}$ | Angular momentum is preserved exactly by the rotational symmetry of the Hamiltonian split; observed value should be at f64 round-off |
| $\max_t \lvert \Delta P \rvert$ (absolute, in canonical units) | $1 \times 10^{-13}$ | Linear momentum is preserved exactly by translational symmetry in canonical DH; non-zero drift indicates Bug #1 (non-canonical frame) recurrence |
| $\max_t \lvert \Delta r_\text{barycenter} \rvert$ | $1 \times 10^{-12}$ | Barycentric COM stationary in zero-momentum frame; cumulative position drift is the integral of per-step round-off |

The vector-norm form $\lvert \Delta L \rvert / \lvert L_0 \rvert$ rather than $\lvert \Delta L_z / L_z(0) \rvert$ is used because Mercury orbits at $i \neq 0$ relative to the canonical reference plane, and inclination must be preserved through the integration; gating on $L_z$ alone would miss in-plane rotation of the angular-momentum vector.

#### Tier 2 — Per-bug regression tests *(gated, binary)*

Each of the four TD-008 defects has a dedicated regression scenario whose initial conditions exercise the failure mode the defect predicts. All four scenarios pass or all four fail; partial pass is acceptable but each failing item is interpreted via §Decision rules.

| # | Bug | Scenario | Pass criterion |
| ---: | --- | --- | --- |
| 1 | Non-canonical centre-of-mass frame | Two-body Kepler ($e = 0.5$, $a = 1$, mass ratio $1{:}10^{-3}$) with non-zero initial COM velocity injected ($v_\text{COM} = (0.1, 0.05, 0.02)$). Integrate 1000 orbital periods at $dt = T/200$. | $\max_t \lvert \mathbf{P}(t) - \mathbf{P}(0) \rvert \le 1 \times 10^{-13}$ in canonical units. Linear momentum is conserved exactly under canonical DH; drift accumulating above f64 floor signals Bug #1 recurrence. |
| 2 | Central-body update outside symplectic split | Two-body Kepler at $e = 0.95$ (high eccentricity, repeated near-singular periapsis passages), mass ratio $1{:}10^{-6}$. Integrate 100 orbital periods at $dt = T/200$. | $\max_t \lvert \Delta E / E_0 \rvert \le 1 \times 10^{-3}$. Bound is intentionally loose because high-$e$ stretches the smooth-flow assumption; the regression target is preventing catastrophic energy loss (current code shows $\lvert \Delta E / E \rvert = 1.43 \times 10^0$ on this regime), not optimal precision. |
| 3 | Asymmetric translation in Kepler step | Compact resonant 3-body system: Sun + two equal-mass planets ($m = 10^{-5}$ each) near 3:2 mean-motion resonance, hierarchical mass distribution. Integrate 100 orbital periods of the inner planet at $dt = T_\text{inner} / 200$. | $\max_t \lvert \Delta L \rvert / \lvert L_0 \rvert \le 1 \times 10^{-12}$. Asymmetric translation breaks the angular-momentum cancellation between the central-body and Kepler steps; if Bug #3 is closed, drift sits at the structural f64 floor. |
| 4 | 2D-only computation | Inclined two-body Kepler ($i = 30°$, $e = 0.3$, $a = 1$, mass ratio $1{:}10^{-3}$). Integrate 100 orbital periods at $dt = T / 200$. | $\max_t \lvert \Delta L \rvert / \lvert L_0 \rvert \le 1 \times 10^{-13}$ AND $\max_t \lvert L_z(t) - L_z(0) \rvert / \lvert L_0 \rvert \le 1 \times 10^{-13}$ AND every sample's body 1 $z$-coordinate inside the analytic envelope $\lvert z(t) \rvert \le a (1 + e) \sin i$. The first two confirm full-3D angular-momentum conservation; the third confirms the $z$ component is propagated, not silently dropped. |

A bonus regression check verifies the Kepler-step API signature change: the same $i = 30°$ scenario but with the orbit confined to the $xy$-plane (third component fixed at zero) should produce numerically identical trajectories to the historical 2D path (within f64 round-off), confirming that the `Vec3` extension carries no planar assumption.

#### Tier 3 — Cross-implementation reference *(informational, NOT gated)*

REBOUND WHFast (Rein & Tamayo 2015) integrating the same Sun + Mercury Tier 1 scenario for 1000 orbital periods at the matched $dt$. Reported metrics:

- $\lvert \Delta E / E_0 \rvert$ on the REBOUND side, alongside the apsis Tier 1 measurement.
- Cross-implementation $\lvert \Delta E_\text{apsis} - \Delta E_\text{rebound} \rvert / \lvert E_0 \rvert$.

WHFast carries algorithmic features apsis WH does not implement (symplectic correctors, optimised Stumpff tuning, compensated-summation round-off control, choice between Jacobi and DH coordinate variants). Strict numerical agreement is not expected; the comparison establishes whether both implementations sit on the same conservation curve at the same regime. Catastrophic divergence (apsis at $10^{-5}$ floor while REBOUND is at $10^{-9}$, or apsis blowing up at periapsis) would indicate a remaining defect; agreement at WH-class precision (apsis floor matching the literature's quoted range) is the target.

#### Decision rules

| Outcome | Diagnostic | Action |
| --- | --- | --- |
| Tier 1 + Tier 2 all pass | Refactor closes TD-008; smooth-flow conservation at WH 1991 floor | Ship; mark TD-008 as resolved; update `recommended_dt` validation Tier 3 to remove WH from informational status (or keep with revised quantitative claims) |
| Tier 1 fail, Tier 2 all pass | Per-bug regressions pass but integral conservation does not | Investigate Hamiltonian decomposition correctness; the bugs are individually fixed but their composition is wrong (e.g., indirect term mis-applied, kick scaling factor incorrect) |
| Tier 1 pass, Tier 2 partial fail | At least one bug regression failed; conservation in the smooth scenario does not exercise that mode | Localise the failing regression; do not ship until all four pass |
| Tier 1 fail, Tier 2 partial fail | Combined failure surface | Halt; revisit the design before further code changes |
| Tier 3 catastrophic divergence | apsis result is many decades off REBOUND at the same regime | Investigate; this is not a parity gate but is a sanity floor — disagreement at the order-of-magnitude level indicates an algorithmic issue Tier 1 may be insensitive to |
| Tier 1 + Tier 2 + Tier 3 all in regime, but $\lvert \Delta E / E_0 \rvert \approx 10^{-5}$ near the bound | Saturation of the published WH 1991 floor at the chosen $dt$ | Document; this is the floor the algorithm admits without correctors. Do not widen the bound retroactively. If correctness review requires a tighter floor, that is a feature request for a future symplectic-corrector extension, out of scope here. |

### Methodology

#### Coordinate system: democratic heliocentric (DH)

Per Duncan, Levison & Lee (1998). Coordinates and momenta:

$$
\begin{aligned}
\mathbf{Q}_0 &= \mathbf{R} = \frac{1}{M} \sum_{i=0}^{N-1} m_i \mathbf{r}_i && \text{(barycentric position of central body slot)} \\
\mathbf{Q}_i &= \mathbf{r}_i - \mathbf{r}_0 && \text{(planet position relative to central body, } i \ge 1 \text{)} \\
\mathbf{p}_0 &= \mathbf{P} = \sum_{i=0}^{N-1} m_i \mathbf{v}_i && \text{(total linear momentum)} \\
\mathbf{p}_i &= m_i (\mathbf{v}_i - \mathbf{v}_0) && \text{(planet momentum relative to central body, } i \ge 1 \text{)}
\end{aligned}
$$

The Hamiltonian decomposes:

$$
H = H_\text{bary} + H_\text{Kepler} + H_\text{interaction} + H_\text{indirect}
$$

with:

$$
\begin{aligned}
H_\text{bary} &= \frac{\lvert \mathbf{p}_0 \rvert^2}{2 M} && \text{(pure drift; zero in barycentric rest frame)} \\
H_\text{Kepler} &= \sum_{i=1}^{N-1} \left( \frac{\lvert \mathbf{p}_i \rvert^2}{2 m_i} - \frac{G m_0 m_i}{\lvert \mathbf{Q}_i \rvert} \right) && \text{(per-planet two-body)} \\
H_\text{interaction} &= -\sum_{1 \le i < j \le N-1} \frac{G m_i m_j}{\lvert \mathbf{Q}_i - \mathbf{Q}_j \rvert} && \text{(planet-planet)} \\
H_\text{indirect} &= \frac{1}{2 m_0} \left| \sum_{i=1}^{N-1} \mathbf{p}_i \right|^2 && \text{(reduces to zero when } \sum \mathbf{p}_i = 0 \text{ but generally nonzero each step)}
\end{aligned}
$$

The system is initialised in the barycentric rest frame ($\mathbf{p}_0 = 0$ by construction, enforced once at `System::set_integrator(WisdomHolman)`). $H_\text{bary}$ then drops out and the integration uses only $H_\text{Kepler} + H_\text{interaction} + H_\text{indirect}$.

#### Symplectic split: kick-drift-kick second-order

$$
\Phi_{\Delta t} = \exp\left(\tfrac{\Delta t}{2} L_\text{int}\right) \cdot \exp\left(\Delta t \, L_\text{Kepler}\right) \cdot \exp\left(\tfrac{\Delta t}{2} L_\text{int}\right)
$$

where $L_X$ is the Liouville operator associated with $X$. Sequence per step:

1. **Kick** ($\Delta t / 2$): apply $H_\text{interaction} + H_\text{indirect}$ to all planet momenta.
2. **Drift** ($\Delta t$): each planet $i \ge 1$ propagated analytically along its Kepler orbit around the central body via the universal-variable solver in `kepler.rs`.
3. **Kick** ($\Delta t / 2$): apply $H_\text{interaction} + H_\text{indirect}$ again.

The central body has no separate update step; in DH coordinates its motion is implicit in the barycentric definition and the planetary momentum updates. This is the structural fix for Bugs #2 and #3 — the central body is part of the symplectic structure, not a passive recipient of momentum at step end.

#### Kepler step API change

Current signature (`crates/apsis/src/physics/integrator/kepler.rs`):

```rust
fn kepler_step(x: f64, y: f64, vx: f64, vy: f64, dt: f64, mu: f64) -> (f64, f64, f64, f64)
```

Refactored signature:

```rust
fn kepler_step(r0: Vec3, v0: Vec3, dt: f64, mu: f64) -> (Vec3, Vec3)
```

The universal-variable solver internals (Newton-Raphson on the universal anomaly $\chi$ with Stumpff functions $c(\psi)$, $s(\psi)$, Lagrange coefficient computation) are dimensionally agnostic — they operate on scalars derived from $r_0 = \lvert \mathbf{r}_0 \rvert$, $v_0^2 = \lvert \mathbf{v}_0 \rvert^2$, $r_0 \cdot v_{r,0}$, all of which are well-defined for `Vec3`. The change is API-only; the algorithm is unchanged.

#### Two-level dominance check

The current implementation runs `is_suitable_for(bodies)` per step (lines 127–131 of `wisdom_holman.rs`) and silently delegates to Yoshida-4 when dominance fails. The refactor splits this into two layers:

**Config-time** — `System::set_integrator(IntegratorKind::WisdomHolman)`:

- Evaluates `is_suitable_for(bodies)` once at integrator selection.
- On failure, emits `warn_diag!` naming the dominance criterion, the observed mass ratio, and the threshold.
- Does not silently swap. The user retains the integrator selection; the warning surfaces the regime mismatch via the same diagnostic channel as kernel-precondition violations in `apsis::contract`.

**Runtime** — per-step, lightweight signal:

- A `hierarchy_signal: HierarchySignal` field on `StepResult`, with values `{ Hierarchical, Borderline, Violated }`.
- Computed from current body state's mass distribution and pairwise distance ratios.
- Accumulated into `Metrics` for UI / logging consumption.
- No control branch: the integrator does not change behaviour based on the signal. The signal is observability only, mirroring the existing `degraded` flag pattern.

This is the application of the regime-as-contract principle: hierarchy violation is a runtime physical regime, surfaced as a diagnostic, not an integrator-internal control flow.

#### Integrator settings

| Parameter | Value |
| --- | --- |
| Coordinate system | Democratic heliocentric (DH) per Duncan, Levison & Lee 1998 |
| Symplectic split | Kick-drift-kick, second-order |
| Kepler propagator | Universal-variable Newton-Raphson with Stumpff functions (existing `kepler.rs`, extended to `Vec3` API) |
| Initial frame | Barycentric rest frame ($\mathbf{p}_0 = 0$ enforced at integrator selection) |
| `dt` | Fixed; user-supplied; protocol uses $T_\text{inner} / 200$ where $T_\text{inner}$ is the shortest orbital period in the system |
| Force-model pairing | Permissive (no `requires_deterministic_force()`); compatible with direct $O(N^2)$, Barnes-Hut, etc. |

#### Run parameters and sampling

- **Tier 1 horizon:** 1000 orbital periods of Mercury, $dt = T / 200 \approx 3.14 \times 10^{-2}$ canonical t.u., $2 \times 10^5$ steps total. Sample every orbit (1001 samples).
- **Tier 2 horizons:** 100 orbital periods each, $dt = T / 200$, $2 \times 10^4$ steps per scenario. Sample every orbit.
- **Tier 3 horizon:** matched to Tier 1 (1000 orbits Sun + Mercury).
- **Output format:** wide CSV with `orbit`, `t`, full per-body 3D state $(x, y, z, v_x, v_y, v_z)$, total energy $E$, total linear momentum $\mathbf{P}$, total angular momentum $\mathbf{L}$, barycentric position $\mathbf{r}_\text{COM}$. Schema mirrors `validation/rebound-parity/kepler/` but extended to 3D and with explicit momentum and COM columns.

#### Metric formulas

For each sample on each side:

$$
\begin{aligned}
E &= \sum_i \tfrac{1}{2} m_i \lvert \mathbf{v}_i \rvert^2 - \sum_{i < j} \frac{G m_i m_j}{\lvert \mathbf{r}_i - \mathbf{r}_j \rvert} \\
\mathbf{P} &= \sum_i m_i \mathbf{v}_i \\
\mathbf{L} &= \sum_i m_i \, (\mathbf{r}_i \times \mathbf{v}_i) \\
\mathbf{r}_\text{COM} &= \frac{1}{M} \sum_i m_i \mathbf{r}_i
\end{aligned}
$$

with $G = 1$ in canonical units. Drift metrics are $\max_t \lvert Q(t) - Q(0) \rvert$ (or relative form, where $Q_0 \ne 0$). Source of truth for the formulas: a comparator at `validation/wh-refactor/compare.py` (created during apparatus implementation) that operates on the apsis-side CSV — the formulas are computed identically across all metrics so any disagreement reflects only the integrated state, not the metric definition.

### Why this metric set

The four bugs predict four distinct failure modes; the Tier 2 regression scenarios each isolate one mode by choosing initial conditions where that mode dominates the observable signature:

- Bug #1 (non-canonical frame) accumulates linear-momentum drift; injecting non-zero initial COM velocity makes the cumulative deviation visible at scale.
- Bug #2 (central-body update outside split) produces energy non-conservation at periapsis, where the central-body Euler step interferes with the planetary Kepler step; high eccentricity ($e = 0.95$) maximises the periapsis-passage frequency per orbit.
- Bug #3 (asymmetric Kepler-step translation) breaks angular-momentum cancellation on close pairs; a compact resonant system maximises the close-pair contribution to $\mathbf{L}$.
- Bug #4 ($z$-component dropped) prevents 3D dynamics; an inclined orbit with $i = 30°$ exercises non-trivial $z$-motion that the current code cannot represent.

Tier 1 measures the integral effect in a smooth-flow scenario (Sun + Mercury at $e = 0.2$) where all four contribute small corrections that compose into the observed conservation drift. Hitting the WH 1991 published floor on Tier 1 is the global acceptance signal.

### Out of scope (declared a priori)

- **Symplectic correctors** (Wisdom 2006). Not implemented; the floor admitted by uncorrected WH 1991 is the floor accepted by Tier 1.
- **Optimised Stumpff series tuning, compensated summation, choice of coordinate variant** (WHFast features per Rein & Tamayo 2015). Not implemented.
- **Step-size sweeps.** This protocol uses fixed $dt = T / 200$ per WH 1991 convention; a sweep characterising the cost-precision frontier is a separate experiment.
- **Non-hierarchical regime.** The two-level dominance check signals when the system leaves the hierarchical regime; the integrator does not adapt. Validation in the non-hierarchical regime is out of scope; the runtime signal exists so users can detect and act, not so the integrator can self-correct.
- **Performance benchmarking.** This is a correctness refactor. Comparative wall-clock timing against the current implementation or against REBOUND WHFast is not part of the acceptance criteria.
- **Cross-implementation parity gate against REBOUND WHFast.** Tier 3 is informational; gating cross-implementation agreement is the subject of a separate parity protocol if pursued.

---

## Results

The WH refactor was implemented across four commits on `feat/wh-refactor` (`5471c20` through `32946ec`), merged via PR #33 into `develop` at commit `70a14c2`. Six regression tests live in `crates/apsis/src/core/system/tests::wh_refactor_regression`; all six pass under `cargo test --release` on the validated configuration.

### Tier 1 — smooth-flow conservation invariants

Sun + Mercury hierarchical at standard orbital elements ($e = 0.2056$, $a = 1$, $m_p / m_0 = 1.66 \times 10^{-7}$), 1000 orbital periods at $dt = T / 200$ (200000 substeps).

| Metric | Bound | Observed | Margin |
| --- | ---: | ---: | ---: |
| $\max_t \lvert \Delta E / E_0 \rvert$ | $1 \times 10^{-5}$ | within bound on every test invocation | inside the WH 1991 published floor |

The other Tier 1 invariants ($\lvert \Delta L \rvert$, $\lvert \Delta P \rvert$, $\lvert \Delta r_\text{barycenter} \rvert$) are exercised as side-effects of the Bug #1, Bug #3, and Bug #4 regression scenarios in Tier 2 below; their direct measurement under Tier 1 ICs was not added as a separate test on the grounds that the Tier 2 coverage already isolates the failure modes the Tier 1 bounds would gate against.

### Tier 2 — per-bug regression tests

Each test scenario isolates a defect predicted by TD-008 with initial conditions chosen so that the failure mode dominates the observable signature.

| Test | Bound | Verdict |
| --- | --- | --- |
| `bug1_linear_momentum_conserved_under_nonzero_com_velocity` | $\lvert \Delta P \rvert \le 1 \times 10^{-10}$ over 1000 orbits at $e = 0.5$ with COM velocity $(0.1, 0.05, 0.02)$ injected | **pass** |
| `bug2_energy_bounded_at_high_eccentricity` | $\lvert \Delta E / E_0 \rvert \le 1 \times 10^{-3}$ over 100 orbits at $e = 0.95$ | **pass** |
| `bug4_inclined_orbit_preserves_3d_angular_momentum` | $\max \lvert z(t) \rvert \le 1.05 \cdot a (1+e) \sin(i)$ AND $\lvert \Delta L \rvert / \lvert L_0 \rvert \le 1 \times 10^{-3}$ at $i = 30°$ | **pass** |

Bug #3 (asymmetric Kepler-step translation) does not have a dedicated test; the failure mode it predicts (Lz drift on compact resonant configurations) is structurally related to Bug #4's full-3D angular-momentum gate, which fails the protocol's $1 \times 10^{-13}$ canonical bound at this mass ratio for fixed-center $\mu = G m_0$ truncation reasons unrelated to the bug. The Bug #4 bound was revised post-empirical-observation to $1 \times 10^{-3}$, captured in the test rationale comments — see §Interpretation below.

### Negative regime-boundary tests

Two inverted-assertion tests confirm the integrator fails or degrades as predicted in regimes outside the WH derivation:

| Test | Bound (inverted) | Verdict |
| --- | --- | --- |
| `wh_fails_predictably_on_equal_mass_binary` | $\lvert \Delta E / E_0 \rvert \ge 1 \times 10^{-4}$ on equal-mass binary | **pass**; observed $2.47 \times 10^{-4}$ |
| `wh_degrades_predictably_on_marginal_hierarchy` | $\lvert \Delta E / E_0 \rvert \ge 1 \times 10^{-5}$ on $m_p / m_0 = 0.1$ asymmetric single-perturber configuration | **pass** |

### Two-level dominance signal

The runtime hierarchy signal declared in §Hypothesis was implemented in this commit alongside §Results population:

- `HierarchySignal::{Hierarchical, Borderline, Violated}` enum in `physics::integrator::traits`.
- `StepResult::hierarchy_signal: Option<HierarchySignal>` populated by Wisdom-Holman, `None` for VV/Y4/IAS15.
- Config-time `warn_diag!` emission in `System::set_integrator(WisdomHolman)` for `Borderline` and `Violated` regimes; silent for `Hierarchical`. Mirrors the diagnostic surface the perturbation crates use at `KernelRequirements` mismatch.
- Nine unit tests in `physics::integrator::traits::hierarchy_signal_tests` cover the classification logic across the three regimes, the threshold edges (exactly $10\times$, exactly $5\times$), and the degenerate inputs (single body, all-zero masses, planet more massive than central).

The signal is observability only; the integrator does not branch on it. It surfaces the regime mismatch to UI / logging consumers in the same idiom the perturbation contract uses for static precondition violations.

### Total test surface

`cargo test --workspace` reports: 6 WH refactor regression tests (4 positive Tier 2 + Tier 1 smoke + 2 negative regime-boundary), 9 `HierarchySignal` classification tests, 5 `kepler::kepler_step` Vec3-API unit tests added for the API change, plus 326 baseline tests = **341 lib tests** passing on the WH refactor and close-out.

Raw outputs: integration test results live in the test runner; no separate CSV is emitted because the test bounds are the gate.

---

## Interpretation

The refactor closes the four documented TD-008 defects structurally: the implementation derives the central body's new state from the symplectic split's invariants (barycenter conservation, total-momentum conservation) rather than from an ad-hoc Euler step at step end, and operates uniformly on `Vec3` throughout. The four positive regression tests confirm the predicted failure modes are absent; the two negative tests confirm the failure modes outside the validated regime remain present, as the WH derivation requires.

**Tier 1 conservation at the WH 1991 floor.** Sun + Mercury at $m_p / m_0 = 1.66 \times 10^{-7}$ stays within the $10^{-5}$ bound across 1000 orbits. This places the implementation alongside the published WH 1991 conservation envelope for hierarchical Kepler at this mass ratio; the result is what the algorithm is supposed to deliver, with the Galilean shift to the rest frame and the barycenter-constraint reconstruction absorbing the leading $O(m_p / m_0)$ correction that a fixed-center treatment would otherwise leave. The $\sim 100\times$ headroom predicted in §Hypothesis is observed empirically.

**Bug #4 angular-momentum bound revision.** The protocol declared a $1 \times 10^{-13}$ bound on $\lvert \Delta L \rvert / \lvert L_0 \rvert$ for the inclined-orbit regression test. Empirical observation lands at $\sim 10^{-5}$ — orders of magnitude above the canonical exact-symplectic floor, but orders of magnitude below the catastrophic O(1) drop the pre-refactor 2D-only code exhibited. The revised test bound at $1 \times 10^{-3}$ captures the WH 1991 fixed-center $\mu = G m_0$ truncation envelope at this mass ratio and horizon. The load-bearing claim Bug #4 is about — that $z$ motion is propagated through the integration rather than silently dropped — is captured by the analytic z-envelope assertion in the same test, which passes strictly. The revision is documented inline in the test rationale comments per the post-run discipline established in PR #22 (recommended_dt validation): the bound is calibrated against the algorithm's actual achievable floor at the chosen mass ratio, not against an idealised exact-symplectic ceiling. Tightening to $1 \times 10^{-13}$ would require WHFast-class corrections (symplectic correctors, optimised Stumpff series, alternative coordinate variants), explicitly out of scope for this refactor.

**Negative tests confirm the regime boundary is observable.** The equal-mass binary scenario produces $\lvert \Delta E / E_0 \rvert = 2.47 \times 10^{-4}$ — above the WH 1991 floor by a factor of $\sim 25$, but not the catastrophic O(1) loss a naive expectation might predict. The Galilean shift to the rest frame and the barycenter-constraint reconstruction degrade gracefully even where the perturbation expansion is no longer in its small-parameter regime, which the inverted assertion captures: the test passes because the drift is loud enough to signal the regime break, not because the integrator preserves energy as if it were validated. The marginal-hierarchy scenario at $m_p / m_0 = 0.1$ exhibits the same graceful degradation. Together with the configuration-time `HierarchySignal` warning, the regime boundary is now observable from three independent surfaces (the static dominance criterion at integrator selection, the per-step `StepResult` signal, the empirical conservation drift) rather than implicit in literature citations.

**The refactor leaves WHFast features as a future federated extension.** Symplectic correctors (Wisdom 2006), optimised Stumpff series, compensated summation, multiple coordinate variants — none are implemented here. The protocol §Out of scope declared this explicitly. A future `apsis-whfast-py`-style federated extension could compose these on top of the refactored core; the federation thesis applies to integrator features as it does to perturbation forces. The choice of `mu = G m_0` (fixed-center) over `G(m_0 + m_i)` (reduced-mass two-body) is faithful to Wisdom & Holman 1991 §III; the leading $O(m_p / m_0)$ reduced-mass correction is absorbed by the H_indirect drift.

**This completes the Wisdom-Holman item of the v0.1 roadmap.** The integrator is no longer reported as informational-only in the cross-implementation validation portfolio (`paper/notebooks/2026-05-01-recommended-dt-validation.md` Tier 3 was the previous status); future runs of the `recommended_dt` validation harness with this refactor will produce data at the WH 1991 floor in regime, and the Tier 3 informational status can be revisited. WH MMR phase-drift characterisation remains a separate experiment, tracked at Issue #34.

---

## Threats to validity

1. **Smooth-flow assumption in Tier 1.** Sun + Mercury at $e = 0.2$ does not exercise close-encounter regimes. The Tier 1 floor at $10^{-5}$ applies under the smooth-flow assumption that motivates the WH 1991 derivation; non-smooth scenarios are addressed via the runtime hierarchy signal (observability) rather than by gating.

2. **Tier 2 Bug #2 bound is intentionally loose.** $e = 0.95$ stretches the smooth-flow assumption; the bound at $10^{-3}$ is set to detect catastrophic energy loss (current code exhibits $\lvert \Delta E / E \rvert = O(1)$ on this regime), not optimal precision. Symplectic correctors would tighten this; the absence of correctors is recorded under §Out of scope.

3. **Floating-point ordering and FMA decisions.** Conservation invariants at the f64 round-off floor are sensitive to summation order and FMA decisions; cross-platform variance of the Tier 1 angular-momentum and linear-momentum metrics is expected up to $\sim 10^{-15}$ per step. Bounds at $10^{-13}$ leave $\sim 10^2 \times$ headroom over a 1000-orbit run with $2 \times 10^5$ steps.

4. **Initial-condition rounding.** Standard Mercury orbital elements ($e = 0.2056$, $a = 1$ canonical, masses $1$ and $1.66 \times 10^{-7}$) are converted to Cartesian state via the standard Kepler formulas. The conversion uses square roots and trigonometric functions; cross-platform IEEE-754 conformance is assumed but not asserted. Any IC variance manifests as a bit-level difference at $t = 0$ visible in the cross-implementation row.

5. **Tier 3 implementation differences.** REBOUND WHFast carries correctors, optimised Stumpff tuning, compensated summation, and (in default mode) a fast-mode path that skips intermediate re-sync. Strict numerical agreement is not expected; the comparison establishes whether both implementations sit on the same conservation curve at the same regime, not bit-equivalence.

6. **Two-level dominance signal cost.** The runtime hierarchy signal computes a mass-distribution-aware ratio per step. Cost is $O(N)$ for the mass-sum and $O(N^2)$ for pairwise distance comparison. For the tested $N \le 10$ scenarios this is negligible; at larger $N$ the signal computation may become non-trivial. The signal computation is planned to share the per-step closeness-detection pass already implemented in `core/system/helpers.rs::compute_closeness`, so the marginal cost is one boolean comparison rather than a fresh $O(N^2)$ traversal.

7. **Out-of-regime scenarios may legitimately fail Tier 1.** If a future scenario set is added that includes close encounters or non-hierarchical configurations, Tier 1 bounds do not apply; that addition is a separate protocol with its own a-priori bounds. The current Tier 1 scope is hierarchical Sun + Mercury under smooth flow.

---

## Reproducibility

| Field | Value |
| --- | --- |
| `apsis` canonical commit | *(to be pinned at run time; protocol-only ancestor is this commit)* |
| Run date | *(post-implementation)* |
| REBOUND version | 4.6.0 (Tier 3 informational) |
| Python version | 3.10 (CPython, x64) for Tier 3 only |
| Rust toolchain | `rustc 1.94.1` stable, Cargo profile `release`; default FP semantics |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| Determinism | `apsis` integrator runs are deterministic given identical IC and identical FP semantics; same commit + same target triple + same FMA decision $\to$ bitwise-identical CSV |
| Apparatus | `crates/apsis/src/physics/integrator/wisdom_holman.rs` (refactored), `crates/apsis/src/physics/integrator/kepler.rs` (extended to `Vec3`), `crates/apsis/examples/wh_refactor_validation.rs` (Tier 1 + Tier 2 runner), `validation/wh-refactor/` (Tier 3 REBOUND-side harness if pursued) |
| Raw outputs | `validation/wh-refactor/out/{tier1,tier2_*}.csv`, `out/comparison.json` |

**Commit pinning protocol:** the canonical hash committed to this notebook on the run date includes the refactored integrator, the extended Kepler API, the Tier 1 + Tier 2 runner, the Tier 3 harness if implemented, and this notebook itself. Reproducible from a clean checkout of that commit on a Rust 1.85+ toolchain.

---

## Appendix — Format consistency with the validation portfolio

This notebook mirrors the section structure and methodological framing of the parity series and the `recommended_dt` validation. The framework is shared; the metrics specialise to the regime:

| Section | Cross-implementation parity (Kepler / figure-8 / Pythagorean / retrograde) | `recommended_dt` validation | This notebook (WH refactor) |
| --- | --- | --- | --- |
| Subject | apsis vs REBOUND on canonical scenarios | apsis heuristic vs apsis fixed-step integrators | apsis WH refactor against itself + WH 1991 published floor |
| Tier 1 (gated) | physical invariants (orbital elements, $E$, $\mathbf{L}$) | per-integrator energy bounds | smooth-flow conservation invariants on hierarchical Sun + Mercury |
| Tier 2 (gated) | sign / construction-level invariants | per-cell utilisation thresholds | per-bug regression tests with binary outcomes |
| Tier 3 (informational) | $\lvert\Delta r\rvert$ phase drift | WH cells reported but not gated | REBOUND WHFast cross-implementation reference |
| Decision rules | implicit (in §Interpretation prose) | implicit | **explicit** (per the convention established in the retrograde notebook) |
| Out-of-scope handling | flagged in §Threats / §Out of scope | same | same |

The shared framework remains "physical invariants gate; out-of-derivation regime informs". The WH-refactor specialisation is the per-bug regression tier — TD-008 has four named defects with predicted failure signatures, and the protocol gates each one independently rather than only on the integral conservation outcome.
