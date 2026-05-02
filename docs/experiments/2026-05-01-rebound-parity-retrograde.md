# REBOUND Parity — Retrograde Kepler e=0.5

**Date:** 2026-05-01
**Subject:** Numerical agreement between IAS15 (apsis) and IAS15 (REBOUND) on a canonical retrograde Kepler orbit, completing the (prograde, retrograde) pair for sign-convention coverage of the Kepler limit.
**Baseline commit:** *(to be pinned at run time)*
**Tooling:** apsis IAS15 (`crates/apsis/src/physics/integrator/ias15.rs`), REBOUND 4.6.0 via Python 3.10 (`reb.IAS15`).
**Status:** *Protocol declared a priori. Run pending — §Results, §Interpretation, and the Reproducibility commit hash are deliberately empty until apparatus is implemented and executed.*

---

## Abstract

This experiment extends the Kepler parity result (notebook `2026-04-25-rebound-parity-kepler.md`) to the retrograde orbital orientation. The IC differs from the prograde Kepler test in exactly one sign — the tangential velocity at periapsis is flipped from $+v_\text{peri}$ to $-v_\text{peri}$. Every other component of the IC is held identical: same primary, same secondary, same $(a, e, r_\text{peri}, |v_\text{peri}|)$, same total energy. Only the direction of motion along the orbit is reversed, which inverts the sign of the specific angular momentum $h_z$ while preserving its magnitude.

Where Kepler-prograde validates that apsis IAS15 reproduces the canonical Kepler magnitude invariants at machine precision against REBOUND, the retrograde experiment closes the **sign-convention gap**: any latent bug in cross-product order, in eccentricity-vector orientation, in `atan2` quadrant handling, or in an internal controller assumption that $L_z > 0$ would manifest here as either a Tier 1 magnitude failure (orientation-reversed quantity disagreeing with itself between sides) or a Tier 2 sign-consistency violation (sign flip during the run, or sign mismatch between sides). Such bugs would pass Kepler-prograde silently — the reason this experiment exists.

This is the fourth and final entry of the parity validation portfolio: Kepler-prograde / figure-8 / Pythagorean / Kepler-retrograde, spanning periodic 2-body, periodic 3-body, chaotic 3-body, and sign-flipped 2-body regimes against REBOUND IAS15.

---

## Motivation

The numerical foundation underlying the v0.1 paper claims (Mercury 1PN at 4.4 ppm, perturbation federation) has been validated against REBOUND IAS15 in three regimes already: prograde Kepler at machine precision (`2026-04-25`), the periodic three-body figure-8 at machine precision (`2026-04-26`), and the chaotic three-body Pythagorean at the f64 close-encounter floor (`2026-04-30`). Each regime exercises a distinct facet of integrator behaviour — Keplerian smoothness, three-body coupling, chaotic stiffness — but all three share one structural feature: **non-negative or zero total angular momentum**.

- Kepler-prograde: $L_z > 0$ by IC construction ($v_\text{peri}$ tangential CCW).
- Figure-8: $L_z = 0$ by IC symmetry (Chenciner & Montgomery 2000 ICs).
- Pythagorean: $L_z = 0$ by IC construction (all initial velocities zero).

A reviewer can reasonably ask whether the integrator's correctness depends on this sign convention. Bug classes that would slip through three positive-or-zero-$L_z$ tests but break under retrograde include:

1. Cross-product evaluated in the wrong order in the inner force / angular-momentum loop — sign error invisible while the convention flatters it.
2. Eccentricity-vector formula with a swapped term order — magnitude correct, orientation reversed.
3. `atan2(e_y, e_x)` for $\omega$ used with arguments in the wrong order — quadrant mapping breaks on retrograde half-plane.
4. Internal controller assumption that $L_z > 0$ in a stability check, convergence criterion, or rejection threshold.
5. Sign-dependent intermediate quantities triggering underflow / overflow paths only on retrograde.

Each of these would pass Kepler-prograde silently and fail retrograde explicitly. The experiment is the minimal IC change that exposes them: same orbit, same energy, same magnitude of $L_z$, only the sign of $L_z$ flipped.

The framing remains validation, not competition. Tolerances are identical to Kepler-prograde, because f64 precision and IAS15's order-conditions are sign-agnostic — relaxing them under retrograde would imply a lack of confidence in the symmetry of the underlying physics and method, which is unfounded. The expected outcome is that all gated metrics pass at the same f64 round-off floor as the prograde test; any departure would itself be the scientific finding.

---

## Protocol *(declared a priori, before any code runs)*

### Hypothesis

For the retrograde Kepler two-body system at $e = 0.5$ integrated under IAS15 in both `apsis` and REBOUND with the IC described in §Methodology and matched integrator settings, the metrics declared below are bounded *a priori* at the values stated. Bounds are organised into three tiers reflecting the evidentiary role of each metric.

**Verdict criterion.** Tier 1 and Tier 2 are both gated; failure of any gated metric reproves the experiment. Tier 3 is informational and never reproves — it characterises phase drift, which is not a cross-implementation invariant for adaptive integrators (per `2026-04-25-rebound-parity-kepler.md` §Pilot Interpretation).

#### Tier 1 — Magnitude invariants *(gated)*

These are the constants of motion of pure Kepler dynamics; both per-side conservation and cross-implementation agreement are gated. **Tolerances are identical to Kepler-prograde.** The numerical floor is sign-agnostic — relaxing here would imply distrust of the symmetry of the physics and method, which is unjustified.

The tier deliberately tests only the **magnitudes** of orientation-bearing quantities ($|h|$); the sign of $h$ is gated separately in Tier 2 so the diagnostic distinguishes a magnitude-drift bug from an orientation-flip bug. A bug that preserves $|h|$ but inverts $\mathrm{sign}(h)$ intermittently would pass Tier 1 silently; gating sign separately catches it.

| Metric | Bound | Origin |
| --- | ---: | --- |
| $\max\|\Delta a\|/a$ over both sides | $1 \times 10^{-13}$ | Specific energy $\to a = -\mu/(2\varepsilon)$ |
| $\max\|\Delta e\|$ over both sides | $1 \times 10^{-13}$ | $(E, h) \to e^2 = 1 - h^2/(\mu a)$ |
| $\max\|\Delta\omega\|$ over both sides (rad) | $1 \times 10^{-12}$ | `atan2(e_y, e_x)`; $1/\|e\|$ condition factor justifies one decade of margin |
| $\max\bigl(\bigl\|\,\|h\| - \|h_0\|\,\bigr\| / \|h_0\|\bigr)$ over both sides | $1 \times 10^{-13}$ | Specific angular momentum **magnitude only** (sign gated in Tier 2) |
| $\max\|\Delta E/E_0\|$ apsis | $1 \times 10^{-13}$ | IAS15 machine-precision conservation (Rein & Spiegel 2015 §4) |
| $\max\|\Delta E/E_0\|$ rebound | $1 \times 10^{-13}$ | Same |
| Cross-impl $\max\|E_\text{apsis} - E_\text{rebound}\|/\|E_0\|$ | $1 \times 10^{-13}$ | Two correct IAS15 sides agree on $E$ to f64 floor |

Total Tier-1 gated metrics: **7**, mirroring Kepler-prograde row-for-row.

#### Tier 2 — Orientation invariants *(gated, binary)*

Tier 2 gates the **sign of $h$** as a discrete invariant. These checks are binary (pass/fail) and capture the bug classes Tier 1 magnitude bounds cannot — `sign(h)` flipping mid-run, $|h|$ collapsing toward zero, or the two implementations disagreeing on orientation.

Let $h_0 := h(t = 0)$ be the initial specific angular momentum. For the retrograde IC, $h_0 < 0$ by construction.

| # | Check | Definition | Pass criterion |
| ---: | --- | --- | --- |
| 1 | apsis sign consistency | $\mathrm{sign}(h_\text{apsis}(t_k)) = \mathrm{sign}(h_0)$ at every sample $t_k$, AND $\|h_\text{apsis}(t_k)\| > \varepsilon_\text{floor}$ at every $t_k$ | All samples satisfy both |
| 2 | rebound sign consistency | Same on the REBOUND side | All samples satisfy both |
| 3 | cross-impl sign agreement | $\mathrm{sign}(h_\text{apsis}(t_k)) = \mathrm{sign}(h_\text{rebound}(t_k))$ at every sample $t_k$ | All samples agree |

The near-zero floor is set at $\varepsilon_\text{floor} = 1 \times 10^{-10}$. For the IC declared in §Methodology, $|h_0| = \sqrt{\mu \, a \, (1 - e^2)} = \sqrt{1 \cdot 1 \cdot 0.75} \approx 0.866$ in canonical units, so $\varepsilon_\text{floor}$ sits about 10 decades below $|h_0|$. This places it well above any plausible f64 round-off accumulation over 100 orbits ($\sim 10^{-12}$ absolute ceiling for IAS15-class methods on Kepler smooth flow) while still firing on any pathological collapse toward zero. It is a defensive guard, not a routine threshold.

Tier 2 has no continuous numerical bound — these are exact sign checks. If the integrator is correct, every sample passes by construction; if any bug class enumerated in §Motivation is present, it manifests as a binary failure at the first affected sample.

#### Tier 3 — Geometric coherence *(informational, NOT gated)*

- **Per-body position drift** — $\max_t |\mathbf{r}_{1,\text{apsis}}(t) - \mathbf{r}_{1,\text{rebound}}(t)|$ over all sample times. Reported as context. The expected magnitude is the same as Kepler-prograde ($\sim 10^{-9}$ peak around orbit 81; see prograde §Pilot Interpretation). No tolerance is declared because phase drift is not a cross-implementation invariant under adaptive controllers — gating on it would conflate physical disagreement with controller-level ULP divergence, an error the prograde notebook diagnosed and corrected.

### Methodology

#### Initial conditions

Two-body system in canonical units ($G = 1$):

| Body | Mass | $x$ | $y$ | $v_x$ | $v_y$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 (primary) | $1$ | $0$ | $0$ | $0$ | $0$ |
| 2 (secondary) | $10^{-6}$ | $r_\text{peri} = 0.5$ | $0$ | $0$ | **$-v_\text{peri}$** |

with $a = 1$, $e = 0.5$, $r_\text{peri} = a(1 - e) = 0.5$, and $v_\text{peri} = \sqrt{(1 + e) / (a(1 - e))} \approx 1.732051$.

The single difference vs Kepler-prograde is the sign of $v_y$ on the secondary: $+v_\text{peri}$ → $-v_\text{peri}$. Every other component, mass, length, and physical scale is held identical to the prograde test. This isolates sign-convention coverage as the only experimental variable.

By construction:

- $\mathbf{P}(0) = (0, -10^{-6} \cdot v_\text{peri}, 0) \approx (0, -1.73 \times 10^{-6}, 0)$ — non-zero, but reflected from the primary by the centre-of-mass-shift step (see below).
- $h_0 := (x_2 - x_1) \cdot (v_{y,2} - v_{y,1}) - (y_2 - y_1) \cdot (v_{x,2} - v_{x,1}) = 0.5 \cdot (-v_\text{peri}) - 0 = -0.866$ — explicitly negative, as the experimental design requires.
- $\varepsilon = \tfrac{1}{2} v_\text{peri}^2 - \mu / r_\text{peri} = 1.5 - 2 = -0.5$, identical to prograde.
- $|h_0| \approx 0.866$, identical to prograde.

The centre of mass is shifted to the origin and zeroed in momentum on both sides before integration starts (same convention as Kepler-prograde and the figure-8 / Pythagorean tests), so the sign on $\mathbf{P}_0$ above is corrected to zero before any integration step. The sign of $h_0$ is preserved because the COM-shift is a Galilean transform that does not affect angular momentum about the COM.

#### Integrator settings

| Parameter | apsis IAS15 | REBOUND IAS15 |
| --- | --- | --- |
| Initial timestep | $T/1000 \approx 6.28 \times 10^{-3}$ | $T/1000 \approx 6.28 \times 10^{-3}$ |
| Adaptive control | active, default tolerance | active, default `epsilon` |
| Force model | direct $O(N^2)$ (forced via pairing rule, ADR-003) | direct (REBOUND default) |
| Exact finish time | not enforced | `sim.exact_finish_time = 1` |

Identical to Kepler-prograde. The orbital period $T = 2\pi \sqrt{a^3 / \mu} = 2\pi$ in canonical units does not depend on the sign of $v$, so the initial timestep and horizon scale are the same.

#### Run parameters and sampling

- **Total integration:** 100 orbital periods ($T = 2\pi$ in canonical units). Same horizon as Kepler-prograde.
- **Output cadence:** 1 sample per orbital period plus initial state — 101 samples per body per side. Schema mirrors `validation/rebound-parity/kepler/` byte-for-byte to preserve cross-experiment comparability of the parity portfolio.
- **Output format:** wide CSV with `sample`, `t`, full per-body state $(x, y, v_x, v_y)$ for both bodies, and total energy $E$.

#### Metric formulas

Tier 1 quantities are computed identically on both sides — same formula, same constants, same per-sample state vectors. Only the integrated state itself differs.

For each sample on each side, with $\vec{r} = \vec{r}_2 - \vec{r}_1$, $\vec{v} = \vec{v}_2 - \vec{v}_1$, and $\mu = G(m_1 + m_2)$:

$$
\begin{aligned}
\varepsilon &= \tfrac{1}{2}\,v^2 - \mu/r &&\text{(specific energy)} \\
a &= -\mu / (2\varepsilon) &&\text{(semi-major axis)} \\
h &= x\,v_y - y\,v_x &&\text{(specific angular momentum, $z$-component, signed)} \\
e^2 &= 1 - h^2 / (\mu \, a) &&\text{(eccentricity squared)} \\
\vec{e}_\text{vec} &= \big((v^2 - \mu/r)\,\vec{r} - (\vec{r}\cdot\vec{v})\,\vec{v}\big) / \mu &&\text{(eccentricity vector)} \\
\omega &= \mathrm{atan2}(e_y, e_x) &&\text{(argument of periapsis, rad)}
\end{aligned}
$$

For Tier 1, $|h|$ is reported as $\sqrt{h^2}$ (i.e., the magnitude); the **signed** $h$ is preserved separately for Tier 2 sign checks. This separation is the load-bearing methodological choice that distinguishes a magnitude-drift bug from an orientation-flip bug.

### Why this metric set, and why the sign / magnitude separation

The justification for gating on orbital invariants rather than $|\Delta r|$ is the same as Kepler-prograde (`2026-04-25` §Pilot Interpretation): adaptive integrators with non-deterministic `dt` sequences accumulate orbital phase drift that conflates with geometric drift in the $|\Delta r|$ metric, even when the two trajectories live on the same Kepler ellipse. The invariant set $(a, e, \omega, |h|, E)$ measures what the physics actually conserves.

The retrograde-specific addition is the **sign / magnitude separation of $h$** in Tier 1 vs Tier 2. A naive `|Δh|/|h_0|` test conflates two distinct bug modes:

1. *Magnitude-drift bug:* $|h|$ drifts secularly, sign preserved. Caught by Tier 1.
2. *Orientation-flip bug:* $|h|$ preserved, $\mathrm{sign}(h)$ inverted intermittently or globally. Tier 1 (which uses $|h|$) **does not see this** — the absolute-value strips orientation. Caught only by Tier 2.

Reporting them as separate gates makes the diagnostic unambiguous: a Tier-1-only failure points at energy or radial bookkeeping; a Tier-2-only failure points at cross-product order, eccentricity-vector composition, or `atan2` argument order; a joint failure points at a deeper inner-loop defect.

### Out of scope (declared a priori)

- **Trajectory-level prediction.** Same as Kepler-prograde — phase-level $|\Delta r|$ is informational only, since adaptive controllers produce ULP-different `dt` sequences that accumulate as orbital phase drift.
- **3D / inclined retrograde.** This experiment is planar 2D. Out-of-plane retrograde (orbit inclination > 90°) requires the 3D port (Phase 0 of the doc roadmap) and is reserved for v0.2.
- **Eccentricity sweep.** Only $e = 0.5$ is tested, matching Kepler-prograde for direct comparability. Sensitivity to eccentricity is a separate study and not relevant to the sign-convention claim of this experiment.
- **Substep-economy comparison.** Same justification as in Pythagorean (`2026-04-30`) — making cross-implementation substep-count comparison scientifically meaningful requires standardised parity telemetry that does not exist today; deferred to a parity-portfolio-wide enhancement issue.
- **Time-reversal test.** A symmetric retrograde-vs-time-reversed-prograde comparison would be a separate axis; this experiment fixes attention on the IC sign-flip, not on integrator time-reversal symmetry.

---

## Results

*Pending run. §Results, §Interpretation, and the Reproducibility canonical-commit hash will be populated post-run, in a separate commit, against the apparatus commit hash. The protocol declared above is frozen at this commit and will not be retroactively altered to match observed values; any post-run protocol change will be recorded as a separate commit with explicit two-phase framing and rationale, per the discipline established in PR #22 (recommended_dt validation).*

---

## Interpretation

*Pending run.*

---

## Threats to validity

1. **Floating-point ordering.** The two IAS15 implementations sum forces in different orders, producing different ULP-level rounding. This is the dominant source of any residual differences observed; the orbital-invariant metrics measured at 1–10 ULP confirm the floor is at f64 precision and not above it. Sign of $h$ does not enter the floating-point sum order, so this threat does not affect Tier 2.

2. **FMA usage.** apsis is built with default Rust FP semantics; REBOUND is C with potential FMA via the compiler. Different FMA decisions produce small but systematic deviations within the same ULP envelope. No evidence (in the prograde experiment) of FMA-induced bias above the round-off floor. FMA decisions are not sign-dependent under IEEE-754, so this threat is symmetric vs prograde.

3. **Adaptive controller details.** Both implementations follow Rein & Spiegel 2015 for the Picard predictor-corrector loop and the $(\varepsilon/\mathrm{err})^{1/7}$ controller, but micro-decisions in the controller (when to grow `dt`, marginal-convergence handling) propagate ULP-level differences in `err` into ULP-level differences in `dt`. Through the orbit those differences accumulate as orbital phase drift, not invariant drift. Tier 1 metrics are insensitive to this divergence by construction; Tier 3 reports the resulting $|\Delta r|$ as informational context.

4. **Initial-condition rounding.** $r_\text{peri} = a(1-e) = 0.5$ is exact in f64. $v_\text{peri} = \sqrt{(1+e)/(a(1-e))} = \sqrt{3}$ involves a square root and may round differently between Rust's and Python's `sqrt`; both should produce the same IEEE-754 result on x86-64 with default rounding because $\sqrt{}$ is a correctly-rounded operation in IEEE-754 (754-2008 §5.4.1). Sign-flipping a correctly-rounded f64 is an exact bit-level operation (just toggles the sign bit), so the $v_y$ on each side is bit-identical to the negation of its prograde counterpart.

5. **Sign of $h_0$ representation.** A negative-zero $h_0$ at $t = 0$ is theoretically possible if $v_\text{peri}$ on the apsis side is positive-zero by f64 representation and on the REBOUND side is negative-zero (or vice versa) — IEEE-754 distinguishes $+0$ from $-0$. By construction $v_\text{peri} = \sqrt{3} > 0$, so $-v_\text{peri}$ is a strictly negative finite number on both sides, ruling out the negative-zero pathology.

6. **Centre-of-mass shift convention.** Both sides apply COM-shift before integration starts. The shift is a Galilean transform: it preserves $h$ exactly (in particular, preserves its sign). Any residual difference after the shift would manifest as an IC drift at $t = 0$, which the cross-implementation $|h_\text{apsis} - h_\text{rebound}|$ at $t = 0$ would surface immediately. Expected $\|\Delta h\|(t=0) = 0$ to bit precision.

7. **Out-of-regime scenarios may legitimately fail.** Same caveat as Kepler-prograde: the smooth-flow assumption underlying the bound derivation holds for $e = 0.5$ Kepler. At higher eccentricity or in close-encounter regimes the bound's derivation does not apply (cf. Pythagorean `2026-04-30` §Energy drift). $e = 0.5$ is comfortably within the regime; this caveat is for context, not for hedging this experiment's claim.

---

## Reproducibility

| Field | Value |
| --- | --- |
| apsis canonical commit | *(to be pinned at run time)* |
| REBOUND version | 4.6.0 |
| Python version | 3.10 (CPython, x64) |
| Rust toolchain | `rustc 1.94.1` stable, Cargo profile `release`; default FP semantics (no `-Cffast-math`-equivalent) |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| Determinism | Same convention as Kepler-prograde and Pythagorean — single-threaded IAS15 on both sides; same commit + same target triple + same FMA decision → bitwise-identical CSV. |
| Harness | `validation/rebound-parity/retrograde/run.py` (cargo example → REBOUND side → comparator) |
| Apsis side | `crates/apsis/examples/rebound_parity_retrograde.rs` |
| Raw outputs | `validation/rebound-parity/retrograde/out/{apsis,rebound}.csv`, `out/comparison.json` |

**Commit pinning protocol:** the canonical hash committed to this notebook on the run date will include the apsis-side cargo example, the Python harness under `validation/rebound-parity/retrograde/`, and this notebook itself. The harness will be reproducible on a clean checkout of that commit with the dependencies pinned in `validation/rebound-parity/retrograde/requirements.txt` (identical to the Kepler / figure-8 / Pythagorean set: `numpy`, `rebound==4.6.0`).

---

## Appendix — Format consistency with the parity portfolio

This notebook deliberately mirrors the section structure and methodological framing of the three preceding parity notebooks. The framework is shared; the metrics specialise to the regime:

| Section | Kepler-prograde | Figure-8 | Pythagorean | Kepler-retrograde |
| --- | --- | --- | --- | --- |
| Regime | periodic 2-body, $L_z > 0$ | periodic 3-body, $L_z = 0$ | chaotic 3-body, $L_z = 0$ | periodic 2-body, $L_z < 0$ |
| Tier 1 | orbital elements + energy | $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ | $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ | orbital elements + energy (mag only on $h$) |
| Tier 2 | (none — magnitude tier covers it) | (none — same) | (none — same) | **sign($h$) consistency, binary** |
| Tier 3 | informational $|\Delta r|$ | informational $|\Delta r|$ | informational $|\Delta r|$, expected $O(1)$ | informational $|\Delta r|$ |
| Sign-convention coverage | $L_z > 0$ only | $L_z = 0$ | $L_z = 0$ | $L_z < 0$ — closes the gap |
| Horizon | $100\,T$ | $10\,T$ + $50\,T$ | $70$ canonical t.u. | $100\,T$ |

The shared framework remains "physical invariants gate; geometric coherence informs". The retrograde specialisation is the explicit Tier 2 sign-consistency gate, which other parity notebooks do not need (their $L_z$ either is positive by construction or is exactly zero by symmetry, neither of which admits the orientation-flip bug class). The Tier 1 magnitude-only treatment of $h$ in this notebook is the dual of that addition: it isolates magnitude-drift from orientation-flip into separate diagnostic channels.
