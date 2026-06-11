# REBOUND Parity — Plummer Cluster (N = 10³, Softened Kernel)

**Date:** 2026-06-09

**Subject:** Moderate-N validation of apsis IAS15 in the softened-kernel regime: a Plummer sphere of N = 10³ equal-mass bodies in virial equilibrium, integrated with Plummer softening on both apsis and REBOUND from bit-identical initial conditions. Per-side conservation gates, cross-implementation invariant parity, and statistical-equilibrium observables.

**Status:** *Protocol declared a priori, before any apparatus code. Phase 0 (pilot, N = 256) executed 2026-06-10: conventions and floor models settled — see §Phase 0 results and gate freeze. Phase 1 (N = 10³, gated) executed 2026-06-10 against the frozen gates: **12/12 metrics pass** — see §Phase 1 results.*

---

## Abstract

The existing parity portfolio (Kepler, figure-8, Pythagorean, retrograde, MERCURIUS) validates apsis at N ≤ 6. This experiment extends the portfolio to N = 10³ — the regime the softened kernel exists for, and the regime in which the library's Exactness precondition mechanism has its physically motivated use case. At cluster N there is no trajectory oracle: N-body systems are exponentially unstable on a timescale of order the crossing time (Miller 1964; Goodman, Heggie & Hut 1993), and energy conservation alone does not certify a converged trajectory (Boekholt & Portegies Zwart 2015). The protocol therefore gates on the conserved integrals of the softened Hamiltonian — E, **L**, **P**, **r**_COM — per side and cross-implementation, and reports statistical-equilibrium observables (virial ratio, Lagrangian radii) as informational context. A registration-only contract assertion confirms that attaching the exactness-requiring `apsis-1pn` operator to the softened system emits the Exactness diagnostic.

---

## Motivation

Three gaps close at once. First, the validated-regime claim: the paper scopes apsis to N ≤ 10³, but no validation exercises N above 6; this experiment supplies the missing point at the boundary. Second, the softening narrative: the Exactness counter-test demonstrates softening as a *violation* (cluster-scale ε on a solar-system problem); this experiment demonstrates softening as *legitimate physics* (the regime it belongs to), completing both directions of the precondition mechanism. Third, the contract tie-in: the same softened system that is correct physics for the cluster is a precondition violation for `apsis-1pn`, so the registration warning is asserted in its native habitat rather than constructed ad hoc.

The framing remains validation, not competition: no performance comparison against REBOUND, no claim of cluster-dynamics capability against collisional codes (NBODY6-class regularisation is out of scope by design). The claim is that apsis IAS15 and REBOUND IAS15 integrate the same softened Hamiltonian at N = 10³ to the precision f64 admits.

---

## Protocol *(declared a priori)*

### Phase 0 — pilot (N = 256, informational)

A full dry-run of the harness at N = 256: IC generation, both sides, comparator, all self-tests. Purpose: (a) measure wall-time and accepted-step count per side to size the N = 10³ run; (b) validate the floor models for the L/P/COM gates against a real run before freezing them; (c) shake out the 3D comparator. No gated claims at N = 256. **Gates for phase 1 are frozen when phase 0 closes; phase 1 results never feed back into the gates.**

### Initial conditions

Plummer model in N-body (Hénon) units: G = 1, M = 1, E = −1/4, crossing time t_cr = GM^{5/2}/(2|E|)^{3/2} = 2√2 (Heggie & Mathieu 1986). The Plummer scale radius follows from the virial theorem in these units: W = 2E = −1/2 and the Plummer potential energy is W = −3πGM²/(32a), so a = 3π/16 ≈ 0.589; the half-mass radius is r_h = a/√(2^{2/3} − 1) ≈ 1.305 a ≈ 0.77. N = 10³ equal masses m = 1/N. Positions and isotropic velocities sampled from the Plummer distribution function by the standard rejection method (Aarseth, Hénon & Wielen 1974), fixed integer seed.

The generator — not either side — performs exactly one centring pass: subtract the mass-weighted mean position and velocity, so **r**_COM(0) = **P**(0) = 0 to generator round-off (residual recorded in the IC file header). **L**(0) is *not* zeroed; it is small (O(N^{−1/2}) statistical) and the L gates measure drift from L(0), not from zero. The ICs are written once to a committed CSV at shortest round-trip f64 precision and read by both sides; on IEEE-754 platforms the parsed doubles are bit-identical, eliminating IC construction as a confound (the same property the Pythagorean integer ICs provided by accident, here provided by protocol).

### Softening

Plummer softening, identical ε on both sides: apsis via the softened kernel (K(r) = 1/√(r² + ε²)), REBOUND via its gravitational softening parameter (Rein & Liu 2012). The value is the optimal softening for an N-particle Plummer sphere, ε_opt = 0.98 N^(−0.26) in Plummer scale lengths (Athanassoula, Fady, Lambert & Bosma 2000, §3.2; fit valid for 30 ≤ N ≤ 3 × 10⁵, which spans both phases here): ε = 0.1626 a ≈ 0.0958 length units for N = 10³, and ε = 0.2318 a ≈ 0.1365 for the N = 256 pilot, each from the same formula. Athanassoula et al.'s force expression places ε identically to the apsis kernel and to REBOUND, so all three state the same convention.

Two guards remain in force. The convention is verified empirically by a single-pair smoke test — the acceleration of one synthetic pair at fixed separation, both sides, protocol ε, agreement asserted at the f64 floor (executed once during protocol drafting: agreement at 2.6 × 10⁻¹³ relative; re-asserted by the harness before every gated run). And one bookkeeping caveat is recorded: REBOUND's reported total energy does not include the softening term in the potential (verified empirically against the closed form, v4.6.0), which is why no per-side energy bookkeeping enters the metrics — the comparator computes every invariant, including the softened PE, from the snapshots alone.

### Integrator settings

| Parameter | apsis IAS15 | REBOUND IAS15 |
| --- | --- | --- |
| Initial timestep | 10⁻³ | 10⁻³ |
| Adaptive control | active, default tolerance | active, default `epsilon` |
| Softening | kernel ε (protocol value) | `sim.softening` (same value) |
| Force model | direct O(N²) | direct |
| Exact finish time | not enforced | `sim.exact_finish_time = 1` |

### Horizon

t_final = 10 time units ≈ 3.5 t_cr. Rationale: the half-mass relaxation time t_rh = 0.138 N √(r_h³/GM) / ln(0.4N) ≈ 16 t.u. for N = 10³ with the unsoftened Coulomb logarithm; at the protocol ε the logarithm shrinks to ln(r_h/ε) ≈ 2.1, giving t_rh ≈ 45 t.u. The horizon therefore sits at ≈ 0.2 t_rh (≤ 0.6 t_rh even on the unsoftened estimate), so virial stationarity is the statistical expectation, while ~3.5 crossing times exercise the adaptive controller through the full range of cluster-internal timescales. Relaxation-driven structural drift over this window is slow but not zero; the Lagrangian-radius expectation below is stated accordingly.

### Sampling and output

10 snapshots per time unit → 101 snapshots. Long-format CSV per side (`sample, t, body, x, y, z, vx, vy, vz`), ~10⁵ rows; the comparator computes every invariant itself from the snapshots — neither side's internal energy bookkeeping enters the metrics.

### Hypothesis — gated metrics

**Verdict criterion.** Tier 1 and Tier 2 gate; failure of any reproves the experiment. Tier 3 is informational and never reproves. The comparator computes all invariants of the **softened** Hamiltonian: PE uses K(r) = 1/√(r² + ε²) with the protocol ε. Computing the unsoftened PE here would manufacture a spurious drift; this is a definitional requirement, not a tolerance choice.

**Tier 1 — conserved integrals of the softened flow (gated).**

- **Energy**, per side and cross-implementation. The softened cluster has no close-encounter stiffness (that is what ε is for), so the smooth-flow round-off model applies: per-side gate = h · 13ε_f64 · √N_steps (round-off random walk over accepted steps, the same model as the Kepler/retrograde gates), with N_steps taken from each side's integrator telemetry, plus the measurement floor ε_f64 · √n_pair for the comparator's 5 × 10⁵-term PE sum. Cross-implementation gate: √2 × the larger side. Headroom h = 10.
- **Angular momentum** |Δ**L**| from L(0), per side and cross-implementation, full 3D vector norm. Floor: Wilkinson cancellation scale 10ε_f64 · max_t Σᵢ mᵢ(|xᵢ||v_yᵢ| + …) per component; cross ×√2. One term of this floor is decided at phase 0 and frozen before the gated run: whether a √N_steps accumulation factor multiplies the scale. The planar comparators' scale-only form passed with wide margins at N = 3, where encounter spikes inflate the cancellation scale; at N = 10³ the scale stays O(1) and per-step round-off may dominate instead. The choice is made against the phase-0 measurement with a written justification, before any N = 10³ data exists.

**Tier 2 — construction-level sanity (gated, weak).**

- **Linear momentum** |Δ**P**| from 0, per side and cross. Pair-force antisymmetry cancels each pair's contribution identically; floor = 10ε_f64 · max_t Σᵢ mᵢ|**v**ᵢ| (component-wise), cross ×√2.
- **Centre of mass** |Δ**r**_COM| from 0, per side and cross. Drift model 1.5 ε_f64 · P_scale · t_final / M, as in the Pythagorean comparator; cross ×2.

**Tier 3 — statistical equilibrium and geometric coherence (informational, never gated).**

- **Virial ratio** Q(t) = −T/W (W the softened potential energy). Expectation: an initial transient of a few crossing times (the distribution function is sampled from the unsoftened Plummer model and then evolved under softened forces, so the system starts out of equilibrium by O(ε²/a²) ≈ 3 % at the protocol ε), settling to a stationary Q ≈ 0.5 with a percent-level ε-dependent offset (the softened potential is not homogeneous of degree −1, so the virial theorem holds only approximately). Reported per side and overlaid cross-implementation.
- **Lagrangian radii** (10 %, 50 %, 90 % mass fractions about the COM). Expectation: approximately stationary, with slow relaxation-driven drift bounded by the t_final/t_rh ≈ 0.2 window; cross-implementation overlay expected to agree closely even after trajectory-level divergence, because the radii are statistics of the same realisation.
- **Per-body position drift** max |**r**_apsis − **r**_rebound|. Expected to reach O(1) within a few crossing times (per-star exponential instability; Miller 1964, Goodman, Heggie & Hut 1993). Reported for shape; the contrast "trajectories diverge, statistics and invariants agree" is the cluster-scale form of the Pythagorean result.

### Comparator self-tests

Before loading any run data, the comparator asserts itself against synthetic states with hand-computed references: a two-body softened state (PE checked against the closed form at the protocol ε), and a 3-body 3D state with nonzero L_x, L_y (full vector L, P, COM checked). The comparator imports nothing from apsis; the no-import rule is enforced by test.

### Contract assertion (registration-only)

The apsis-side harness registers `PostNewtonian1PN` on the softened system before integration and asserts exactly one Exactness diagnostic (required `exact`, provided `softened`); the operator is then dropped and the integration proceeds without it. No claim involves integrating 1PN dynamics in cluster units; the assertion is that the precondition mechanism fires in the regime softening is legitimately used.

### Out of scope *(declared a priori)*

- **Trajectory-level claims** — no reference trajectory exists at this N (see §Abstract citations).
- **Dynamical-outcome claims** (binary formation rates, core collapse, mass segregation) — these require ensemble statistics over many realisations; a single seeded run supports no such claim.
- **Performance comparison** between implementations.
- **Softening-parameter studies** (ε sweeps, optimal-softening analysis) — ε is a fixed protocol constant here.
- **N = 10⁴** — a possible informational follow-up if phase-0 wall-time extrapolation admits it; not part of the gated claim.

---

## Phase 0 results and gate freeze (2026-06-10)

The pilot ran the full harness at N = 256 (ε = 0.1365, t_final = 10, 101
snapshots): convention smoke test, registration gate, both sides, comparator.
Per-side wall time: apsis 14.6 s (1693 accepted steps), REBOUND 4.2 s (854
steps); the N = 10³ run extrapolates to minutes. The single-pair convention
check agreed with the closed form at 1.5 × 10⁻¹² (apsis) and 3.2 × 10⁻¹³
(REBOUND); REBOUND's controller never drove dt below 2.8 × 10⁻³ — the
softening removes the close-encounter stiffness, as intended.

| Metric (max over samples) | apsis | REBOUND | cross |
| --- | ---: | ---: | ---: |
| \|ΔE/E₀\| | 2.8 × 10⁻¹⁵ | 3.1 × 10⁻¹⁵ | 1.7 × 10⁻¹⁵ |
| \|Δ**L**\| | 4.4 × 10⁻¹⁷ | 5.0 × 10⁻¹⁷ | 3.7 × 10⁻¹⁷ |
| \|Δ**P**\| | 5.6 × 10⁻¹⁷ | 5.5 × 10⁻¹⁷ | 6.1 × 10⁻¹⁷ |
| \|Δ**r**_COM\| | 3.4 × 10⁻¹⁶ | 3.6 × 10⁻¹⁶ | 9.7 × 10⁻¹⁷ |

**Gate freeze.** The accumulation question left open in §Hypothesis is
settled by the pilot: the L/P floors are frozen at the **scale-only**
Wilkinson form, with no √N_steps factor. The observed drifts sit 23–78×
below the static cancellation floors (\|Δ**L**\| ≈ 5 × 10⁻¹⁷ against
2.0 × 10⁻¹⁵; \|Δ**P**\| ≈ 6 × 10⁻¹⁷ against 1.3 × 10⁻¹⁵) — pair-force
antisymmetry and the central-force torque cancellation hold at every step,
so per-step round-off does not random-walk into these invariants. The
√N_steps variant would sit three decades above the observable, trading
detection power for slack. Energy keeps the round-off-walk model of
§Hypothesis (observed 500× under it on the pilot); COM keeps the drift
model. All twelve gated metrics pass at the frozen floors on the pilot data.

**Statistical observables.** Q(0) = 0.521 — the +4 % offset of measuring the
exactly-virialised unsoftened realisation with the softened potential, at the
predicted O(ε²/a²) magnitude — settling to Q = 0.469 ± 0.011 over the second
half. Lagrangian radii drift slowly (r₅₀: 0.821 → 0.865) within the
sub-relaxation expectation.

**Tier-3 expectation revised.** The declared expectation of O(1) per-body
cross-implementation drift does not apply at the protocol softening: the
pilot measured max \|Δ**r**\| = 3.2 × 10⁻¹⁴ at 3.5 crossing times. The
per-star exponential instability that drives the unsoftened e-folding
(Miller 1964; Goodman, Heggie & Hut 1993) is fed by close encounters, which
the optimal-softening choice deliberately suppresses; the realised growth
rate from the round-off floor is correspondingly slow. The Tier-3 quantity
remains informational; for phase 1 the expectation is drift far below O(1),
with the contrast "trajectories and statistics both agree" replacing the
unsoftened "trajectories diverge, statistics agree".

**Convention check, recorded form.** The smoke test gates each side at
1 × 10⁻⁹ relative to the closed form — three decades above the measured
agreement, seven below the O(ε²) ≈ 10⁻² signature of a misplaced ε — rather
than at the f64 floor, which a finite-difference acceleration estimate does
not reach. The L floor uses the component-summed Wilkinson scale, an upper
bound of the per-component form declared in §Hypothesis.

---

## Phase 1 results (2026-06-10) — gated run, N = 10³

All twelve gated metrics pass at the floors frozen in §Phase 0. Per-side
wall time: apsis 214 s (2412 accepted steps), REBOUND 75 s (1165 steps);
REBOUND's minimum accepted dt was 2.7 × 10⁻³ — no close-encounter stiffness,
matching the pilot. The step-count ratio (≈ 2.1×) is consistent across both
phases.

| Metric (max over samples) | apsis | REBOUND | cross | gate |
| --- | ---: | ---: | ---: | ---: |
| \|ΔE/E₀\| | 5.9 × 10⁻¹⁵ | 6.4 × 10⁻¹⁵ | 2.2 × 10⁻¹⁵ | 2.6–4.2 × 10⁻¹² |
| \|Δ**L**\| | 5.4 × 10⁻¹⁷ | 5.2 × 10⁻¹⁷ | 3.0 × 10⁻¹⁷ | 1.9–2.8 × 10⁻¹⁵ |
| \|Δ**P**\| | 3.8 × 10⁻¹⁷ | 4.2 × 10⁻¹⁷ | 4.1 × 10⁻¹⁷ | 1.3–1.8 × 10⁻¹⁵ |
| \|Δ**r**_COM\| | 6.5 × 10⁻¹⁶ | 6.5 × 10⁻¹⁶ | 9.7 × 10⁻¹⁷ | 1.9–3.8 × 10⁻¹⁵ |

Energy sits ~450× under its gate on both sides; the conserved integrals of
the softened Hamiltonian agree between implementations at the f64 round-off
floor through 3.5 crossing times at N = 10³.

**Statistical observables.** Q(0) = 0.512 — the softened-measurement offset
again at its predicted magnitude, O(ε²/a²) ≈ 2.6 % at this ε, scaling down
from the pilot's 4 % exactly as the ε ratio implies — settling to
Q = 0.481 ± 0.004. The t = 0 half-mass radius is 0.778, on the Plummer
r_h ≈ 0.77; Lagrangian radii drift slowly (r₅₀: 0.778 → 0.814) within the
sub-relaxation window, and the two implementations' radii agree to ~10⁻¹⁵
throughout. Per-body cross-implementation drift peaks at 2.0 × 10⁻¹⁴ —
trajectories and statistics both agree, per the revised Tier-3 expectation.

---

## Threats to validity

1. **Softening-convention mismatch.** A different placement of ε between the two implementations would produce slow, systematic invariant divergence easily mistaken for integrator disagreement. Mitigated twice: the published-convention check and the single-pair smoke test, both before any gated run.
2. **Initial virial transient misread as drift.** The unsoftened-DF/softened-force mismatch produces a real, expected transient in Q(t). It is declared above with its mechanism; it does not enter any gate.
3. **Relaxation within the horizon.** At ≈ 0.2 t_rh, Lagrangian radii drift slowly for physical reasons. Stationarity is therefore an *approximate* expectation, stated with its bound, and informational only.
4. **Floor models extrapolated from N ≤ 6 scenarios.** The Wilkinson and round-off-walk forms were derived and exercised at small N; the step-accumulation behaviour at N = 10³ is exactly what phase 0 exists to check before gates freeze.
5. **Exponential per-star instability.** As in the Pythagorean scenario, chaos amplifies ULP-level controller differences into O(1) trajectory separation; all gated quantities are invariant under this divergence by construction.

---

## Reproducibility

| Field | Value |
| --- | --- |
| apsis canonical commit | `cbf63a5` (gate freeze; protocol-only ancestor `24a0a91`) |
| REBOUND version | 4.6.0 (recorded in the run CSV header; requirements bound `>=4.0,<5.0` as in the existing portfolio) |
| Python | 3.10 (CPython, repo venv) |
| Operating system | Microsoft Windows 11 Pro for Workstations, x64 |
| IC files | `ics_n256.csv`, `ics_n1000.csv` (committed; seed 20260609 in headers) |
| Harness | `validation/rebound-parity/plummer-cluster/run.py` (smoke → registration gate → apsis side → REBOUND side → comparator) |
| Raw outputs | `out/{apsis,rebound}.csv`, `out/{apsis,rebound}_stats.json`, `out/comparison.json` |
