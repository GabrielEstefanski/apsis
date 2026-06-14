# Mercury 1PN Error Budget — Floor Decomposition of the Precession Gate

**Date:** 2026-06-10

**Subject:** Decomposition of the residual between the integrated 1PN
Sun–Mercury precession (500 orbits, IAS15) and the closed-form
first-order prediction, into derived and measured floors: the
second-order secular term of the 1PN equation of motion itself, the
two-body O(m/M) correction, the c-convention offset, and the round-off
realisation noise of the adaptive integration. Replaces the
single-number framing of the cross-platform section with a budget in
which every component is either derived in this notebook or measured by
a declared protocol.

**Status:** *Protocol declared a priori 2026-06-10. **Phase A closed
2026-06-11**: A1 — second-order coefficient in closed form, verified by
independent extended-precision integration, parameterization-dependent
at $O(\varepsilon^2)$; A2 — two-body coefficient $C(e)$ measured,
$C(0.20563) = 8.0617$, Mercury floor $1.34\times10^{-6}$; A3 — exact.
**Phase B closed 2026-06-13**: the residual decomposes as an
integrator endpoint-sampling term $Q(\nu_\text{end})$ plus the A1+A2
floors (H1 confirmed, H2 closes to $+4.8\times10^{-12}$, H3 round-off
class); the endpoint artefact is fixed by exact-finish-time integration
(ADR-015) and the gate tightens to $9.2\times10^{-6}$.*

---

## Framing

The gate compares the integration of the test-particle 1PN equation of
motion (Anderson et al. 1975)

$$\mathbf{a}_\text{1PN} = \frac{G M}{c^2 r^2} \left[ \left( \frac{4 G M}{r} - v^2 \right) \hat{\mathbf{n}} + 4 \, (\hat{\mathbf{n}} \cdot \mathbf{v}) \, \mathbf{v} \right]$$

against the first-order secular result of the *same* equation,
$\Delta\omega_1 = 6\pi G M / (c^2 a (1 - e^2))$ per orbit. Two
consequences fix the framing. First, physics absent from the equation
of motion (the true 2PN term) contributes nothing to this residual — it
separates the model from nature, not the integration from its oracle.
Second, the comparison has an irreducible derivation floor: the closed
form is first-order in $\varepsilon \equiv G M / (c^2 a (1 - e^2))$,
while the integration carries the full secular content of the equation,
so a second-order secular term $k \varepsilon^2$ ($k$ = O(1–10), derived
in Phase A) bounds the agreement from below at

$$\varepsilon_\text{Mercury} = \frac{1}{10065.130^2 \times 0.387098 \times (1 - 0.20563^2)} \approx 2.66 \times 10^{-8},$$

with $c = 10065.130$ the gate's `C_SOLAR_UNITS` (IAU julian-year value;
the `for_units` Gaussian convention is $10065.32$, 18.9 ppm higher).

Between that floor and the observed residual sit, in expected order of
size: the c-convention offset (exact: $\Delta\omega \propto c^{-2}$, so
a relative gap $\delta$ in $c$ shifts the prediction by $-2\delta$; the
19 ppm Gaussian-vs-IAU gap gives 38 ppm), the round-off realisation
noise of the adaptive integration (measured by ensemble in Phase B),
and the two-body $O(m/M) \approx 1.7 \times 10^{-7}$ correction to the
test-particle formula.

## Hypotheses *(declared a priori)*

- **H1 (non-additivity).** The libm-vs-UCRT transcendental difference
  does not add an independent floor: a 1-ULP difference in the
  controller's `pow` re-seeds the round-off random walk rather than
  biasing it. Test: the cross-implementation Δω shift must fall within
  the distribution of Δω over an ensemble of 1-ULP initial-condition
  twin runs. H1 false (shift outside the ensemble spread) is a
  reportable systematic, not a failed experiment.
- **H2 (budget closure).** The measured residual decomposes as
  *convention offset* (deterministic) $\pm\,\sigma_\omega$ (ensemble)
  $+$ *derivation floors* ($k\varepsilon^2$, $O(m/M)$), with no
  unexplained remainder larger than $\sigma_\omega$. H2 false is
  likewise a finding to report.
- **H3 (growth regime).** $\sigma_\omega(N) \propto N^\alpha$ with
  $\alpha \approx 1/2$ — the round-off random-walk regime (Brouwer 1937
  phase-error class). The exponent is fitted, never assumed: the same
  measurement discriminates the three known regimes —
  $\alpha \approx 1$ (coherent, truncation-dominated),
  $\alpha \approx 1/2$ (random walk), $\alpha \approx 0$
  (bounded/quasi-periodic) — so an off-hypothesis exponent identifies
  the regime rather than merely failing the test. The
  controller-tolerance sweep (B4) cross-checks the truncation-vs-round-off
  discrimination from an independent axis; the two must agree.

## Phase A — derivations *(in this notebook, before apparatus)*

Working rule for every item: the derivation is carried out
mathematically here — no coefficient is taken from memory — and each
result is verified numerically by an independent high-precision
integration before it anchors a gate.

- **A1 — second-order secular advance of the 1PN equation of motion.**
  Orbit-averaged perturbation theory on the Anderson acceleration
  carried to second order in $\varepsilon$; output is the coefficient
  $k$ in $\Delta\omega = \Delta\omega_1 (1 + k\varepsilon +
  O(\varepsilon^2))$ ($k$ dimensionless, expected O(1–10)). Numerical
  verification: integrate the same equation of motion in extended
  precision (mpmath, two-body, test-particle mass) with $c$
  artificially reduced so $\varepsilon$ spans $[10^{-4}, 10^{-2}]$; fit
  the residual against the first-order formula; the fitted coefficient
  must match the derived $k$ within the fit's own confidence interval.
  Literature anchor: the second-order term of the exact Schwarzschild
  geodesic advance (Will 1993) sets the expected order of magnitude
  only — the gate's equation of motion is the truncated 1PN force,
  whose second-order secular coefficient need not coincide with the
  geodesic's; the geodesic value is never used as the oracle.
- **A2 — two-body correction.** Derivation of the O(m/M) shift of the
  apsidal advance when the primary moves (the gate integrates two
  bodies; the formula assumes a fixed source). Verification: runs at
  m, m/10, m/100; the residual component must scale linearly with m/M.
- **A3 — convention offset.** $\Delta\omega \propto c^{-2}$ exactly, so
  two conventions differing by relative $\delta$ in $c$ predict advances
  differing by $-2\delta$ to first order in $\delta$; with
  $\delta = 18.9$ ppm (ADR-014), the offset is $37.8$ ppm. Stated, not
  gated — it is exact at this order.

## Phase A results (2026-06-11)

### A1 — second-order secular advance: closed form, two parameterizations

The azimuthal equation admits an exact first integral,
$h_0 = h\,e^{4\mu u/c^2}$ with $h = r^2\dot\varphi$ and $u = 1/r$ (the
exponent is determined, not assumed, by the symbolic gate G0). Reducing
to the orbit equation $u(\varphi)$ with $h(u)$ exact and carrying
Lindstedt–Poincaré to second order in $\delta = 1/c^2$ gives, with
$s = \Omega\varphi$ and $u_0 = (\mu/h_0^2)(1 + e\cos s)$:

$$\omega_1 = -\frac{3\mu^2}{h_0^2}, \qquad
  \omega_2 = -\frac{\mu^4\,(19 + 2e^2)}{2 h_0^4}, \qquad
  u_1 = \frac{\mu^3 (e^2+5)}{h_0^4},$$

with $u_2$ a finite harmonic polynomial (no secular terms; the
order-by-order residuals vanish identically). The advance per radial
period, $\Delta\omega = 2\pi(\Omega^{-1} - 1)$, reproduces
$6\pi\varepsilon$ exactly at first order (gate G1) and yields

$$k_\text{inv}(e) = \frac{37 + 2e^2}{6}
  \qquad \text{with} \quad \varepsilon = \frac{\mu}{c^2 p},
  \quad p = h_0^2/\mu .$$

The independent verification integrates the same force in extended
precision (40 significant digits, Taylor IVP; periapsis-to-periapsis
angle by root-finding on $\dot r$) over an amplified ladder
$\varepsilon \in [10^{-6}, 10^{-4}]$ from Newtonian periapsis initial
conditions, and measures

$$k(e{=}0.2) = -3.4200004 \pm 3\times10^{-7}, \qquad
  k(e{=}0.4) = -4.9800003 \pm 3\times10^{-7}$$

— opposite in sign to $k_\text{inv}$. The two are reconciled exactly by
the parameterization: **the coefficient is convention-dependent at
$O(\varepsilon^2)$.** The integration's $\varepsilon$ is built from the
osculating elements of the periapsis initial condition, whose
semi-latus rectum differs from the invariant one:
$p_\text{inv} = p_\text{osc}\,e^{8\mu u_0 \delta}$ at the initial point,
and with $\mu u_0 \delta = \varepsilon_\text{osc}(1+e)$ at periapsis the
first-order term $6\pi\varepsilon_\text{inv}$ contributes $-8(1+e)$ to
the osculating-convention coefficient:

$$k_\text{osc}(e) = k_\text{inv}(e) - 8(1+e)
  = -\frac{11}{6} - 8e + \frac{e^2}{3},$$

predicting $-3.42000$ and $-4.98000$ at the two eccentricities —
agreement with the measured values at their jackknife uncertainty. The
linear-in-$e$ term is legitimate in this convention: starting at
periapsis breaks the $e \to -e$ parity that the invariant
parameterization respects ($k_\text{inv}$ is even in $e$).

The gate's closed-form prediction uses the osculating elements of the
initial condition, so the budget takes the osculating convention: for
Mercury ($e = 0.20563$), $k_\text{osc} \approx -3.46$ and the derivation
floor is

$$\left|k_\text{osc}\right|\varepsilon \approx 9.22\times10^{-8}
  \ \text{(relative)},$$

sharpening the order-of-magnitude estimate of §Framing (the floor sits
three decades below the observed residual, as anticipated). Derivation
and verification artefacts:
`paper/notebooks/scripts/error_budget_k_symbolic.py` (exact-integral
gate, Lindstedt orders, harmonic residuals — all asserted) and
`paper/notebooks/scripts/error_budget_k_numerical.py` (ladder, gates:
first-order recovery, 40-vs-55-digit agreement ≥ 25 significant digits,
Newtonian null at $c \to 10^{30}$).

### A3 — convention offset

Closed as stated in the protocol: $\Delta\omega \propto c^{-2}$ gives a
$-2\delta = -37.8$ ppm offset between the two $c$ conventions, exact at
first order in $\delta = 18.9$ ppm.

### A2 — two-body $O(m/M)$ correction

Two-body integration of the implemented pairwise force in the gate's
own conventions (primary at rest at the origin, secondary at Newtonian
periapsis with $\mu = 1$ elements; relative-orbit geometric advance;
same extended-precision machinery as A1) over a mass-ratio ladder
$q \in \{0, 10^{-5}, 10^{-4}, 10^{-3}\}$ at $\varepsilon = 10^{-5}$,
$e = 0.2$:

$$\frac{\Delta\omega(q) - \Delta\omega(0)}{\Delta\omega(0)} = C\,q + O(q^2),
\qquad C = 8.291,$$

with the quadratic term contributing $0.5\,\%$ at $q = 10^{-3}$, $C$
independent of $\varepsilon$ to $0.36\,\%$ across a decade, the
Newtonian null at $3\times10^{-37}$ rad (a Keplerian relative orbit
does not precess at any $q$ — the effect is purely 1PN-relative), and
the $q \to 0$ limit matching the single-body A1 integration to 40
significant digits.

$C$ was first measured at $e = 0.2$; its eccentricity dependence was
then probed by re-measuring on a grid
(`paper/notebooks/scripts/error_budget_a2_eccentricity.py`). $C$ falls
steeply with $e$ — $C(0.15) = 11.08$, $C(0.2) = 8.291$ (reproducing the
single-$e$ value, cross-check), $C(0.25) = 6.61$, $C(0.30) = 5.49$, a
local slope $dC/de \approx -33.6$ — so the coefficient must be taken at
Mercury's own eccentricity, not at $0.2$:

$$C(e_\text{Mercury} = 0.20563) = 8.0617 \quad
  (\text{direct measurement; pre-registered at} \approx 8.074
  \ \text{before the run}).$$

For Mercury, $q = 1.66\times10^{-7}$ then gives a two-body floor of

$$C(0.20563)\,q \approx 1.34\times10^{-6}\ \text{(relative)},$$

the dominant derivation floor — roughly $8\times$ the naive $m/M$ ceiling
that treats the coefficient as unity.

---

## Phase B — ensemble measurement *(gates frozen at end of Phase A)*

All runs: the gate scenario unchanged (IAS15, initial dt 10⁻⁴,
500 orbits unless stated), both constructors (`for_units` and
`from_raw_c`), current mainline.

- **B1 — realisation noise.** $K = 25$ twin runs per convention, each
  perturbing Mercury's initial $x$-position by 1 ULP (the induced
  change in the predicted advance is $O(\text{ULP}/a)$, negligible).
  Output: $\sigma_\omega$ and the full endpoint distribution.
- **B2 — H1 test.** Position of the libm-vs-UCRT shift (0.002
  arcsec/century, controller-pow notebook) within the B1 distribution.
- **B3 — H3 test.** $\sigma_\omega$ at $N \in \{100, 250, 500, 1000,
  2000\}$ orbits; the growth exponent $\alpha$ fitted from the series
  (regime identification per H3).
- **B4 — controller-tolerance sweep.** Residual at
  $\epsilon_b \in \{10^{-7} \ldots 10^{-11}\}$: movement with
  $\epsilon_b$ indicates a truncation component; a plateau indicates
  round-off domination.
- **B5 — re-measurement.** The two conventions' central residuals with
  the post-ADR-014 constants — the numbers that replace the
  cross-platform section's current values, each then placed in the
  budget table of §Verdict.

## Phase B results (2026-06-13)

### B1 — raw residual

$K = 25$ single-ULP twins per convention at $N = 500$ give a signed
residual that is *not* centred on zero:

| Constructor | mean | $\sigma_\omega$ | central (ulp = 0) |
| --- | --- | --- | --- |
| `from_raw_c` | $-5.53\times10^{-5}$ | $4.10\times10^{-5}$ | $-7.03\times10^{-5}$ |
| `for_units`  | $-5.89\times10^{-5}$ | $3.70\times10^{-5}$ | $-1.16\times10^{-4}$ |

In angle units both the ensemble mean ($\approx -1.5\times10^{-8}$ rad)
and $\sigma_\omega$ ($\approx 1.0\times10^{-8}$ rad) are constant in $N$
over $N \in [100, 2000]$ (B3; raw growth exponent
$\alpha_\text{angle} = -0.01 \pm 0.05$). A constant-in-$N$ bias *and*
spread are incompatible with accumulation noise — both live at the
endpoint, not in the integration.

### B′ — endpoint sampling *(amendment, declared 2026-06-12)*

`System::integrate_until` exited at the first accepted step with
$t \ge t_\text{end}$, sampling the endpoint up to one adaptive sub-step
past $t_\text{end} = N P_0$, at a small osculating true anomaly
$\nu_\text{end} \ne 0$. The osculating $\omega$ (Newtonian e-vector, the
gate's measurement) carries O($\varepsilon$) short-period structure;
differentiating the osculating definitions along the 1PN flow — inputs
the perturbed Binet equation and $h' = -4\mu\delta h u'$, both from the
force — gives

$$\frac{d\omega}{d\nu} = \varepsilon \left[ -\frac{3}{e}\cos\nu + 3
  - 5\cos 2\nu + e\cos\nu \right],$$

whose secular part reproduces $6\pi\varepsilon$ per orbit (an independent
recovery of the A1 first-order result) and whose periodic part is odd in
$\nu$ — invisible to periapsis-to-periapsis measurements. A fixed-time
endpoint at anomaly $\nu$ past the N-th periapsis instead sees

$$Q(\nu) = \varepsilon \left[ 3\nu - \left(\tfrac{3}{e} - e\right)\sin\nu
  - \tfrac{5}{2}\sin 2\nu \right], \qquad
  Q'(0) = -\varepsilon\,\frac{(3-e)(1+e)}{e}.$$

The $3/e$ amplification is the conditioning of the e-vector direction at
small $e$; for Mercury $Q'(0) = -16.38\,\varepsilon \approx
-4.36\times10^{-7}$ rad per radian, so a half-sub-step overshoot
($\nu \approx 0.03$) produces the observed $-1.4\times10^{-8}$ rad. All
closed forms are asserted symbolically
(`paper/notebooks/scripts/error_budget_endpoint_symbolic.py`; gates
GB-h, GB0–GB5) and verified in extended precision
(`error_budget_endpoint_numerical.py`): at the *exact* endpoint the
fitted O($\varepsilon$) coefficient is $-2.4\times10^{-6}$ — against the
f64 ensembles' $-0.56$, confirming the offset is an integrator-endpoint
artefact, not a derivation floor — and $Q(\nu)$ reproduces displaced
endpoints parameter-free to second order in $\varepsilon$.

### Closure

With the endpoint state $(t_\text{overshoot}, \nu_\text{end})$ emitted
per run, the measured angle residual regresses on $Q(\nu_\text{end})$
with slope $1.000\,000\,6 \pm 0.000\,000\,2$, and subtracting
$Q(\nu_\text{end})$ collapses the ULP-twin spread from $10^{-8}$ to
$10^{-14}$ rad. The corrected mean is the floor sum: with the e-correct
floors (§A2, $C(0.20563) = 8.0617$, giving $A2 = 1.34\times10^{-6}$ and
$A1 = -9.22\times10^{-8}$) the prediction $+1.246\times10^{-6}$ relative
matches the measured corrected floor to $+4.8\times10^{-12}$. The entire
Phase-B "realisation noise" was endpoint sampling; the gate residual
decomposes per run as $Q(\nu_\text{end})$ (deterministic) plus the
A1 + A2 derivation floors.

### Verdicts

- **H1** — confirmed: the libm-vs-UCRT shift sits at $1.1\,\sigma_\omega$
  of the B1 distribution. The cross-implementation difference is an
  endpoint re-draw, not accumulated arithmetic.
- **H2** — closes: the predicted floors match the corrected residual to
  $+4.8\times10^{-12}$ relative, with no remainder above the
  $10^{-14}$ rad ensemble noise. Phase C is not required.
- **H3** — the raw $\alpha_\text{angle} = -0.01 \pm 0.05$ identified the
  bounded regime, but the bounded quantity was the endpoint-sampling
  error; on corrected residuals the integration noise grows as
  $\alpha_\text{angle} = +0.36 \pm 0.01$ at $\sigma \sim 10^{-14}$ rad
  (round-off class). The B4 tolerance sweep agrees from an independent
  axis — the corrected residual is $\epsilon_b$-independent across
  $[10^{-7}, 10^{-11}]$ while the raw offset tracks the sub-step size.
- **B5** — the single-run central residuals are endpoint-sampling draws
  under overshoot semantics, not characteristic numbers of the
  constructors; the budget's reproducible content is the floor sum and
  the endpoint mechanism.

### Resolution

Exact-finish-time integration (ADR-015) clips the final step onto
$t_\text{end}$: the overshoot is then zero, the gate residual is the
floor sum $+4.581\times10^{-6}$ relative (pre-registered from the pre-fix
endpoint states, reproduced to seven significant digits), and the
Mercury gate tightens from $10^{-4}$ to $9.2\times10^{-6}$.

### Protocol note

The endpoint floor was identified from Phase-B data, then derived blind:
closed-form in $e$, evaluated at Mercury's $e$ only afterwards, verified
in extended precision (no f64 integrator in the loop) and against the
ensemble with the overshoot pre-registered ($\langle\delta t\rangle =
5.3\times10^{-3}$ predicted before the instrumented re-run measured
$5.35\times10^{-3}$). The §Verdict rule "no bound adjusted after gated
data exists" is met by declaration and blind derivation, not by tuning.

## Phase C — arithmetic/cadence decomposition *(conditional)*

Run only if H2 fails: replay the accepted-dt sequence of a reference
run with the alternative `pow` implementation, separating
per-step-arithmetic from substep-cadence contributions. Requires a
test-only dt-replay hook in the IAS15 controller — an integrator change
scoped and reviewed separately before any Phase C work.

## Verdict criteria

Phase A closes when each derived coefficient passes its numerical
verification. Phase B's gates: H1 and H3 verdicts as declared; the
budget table must account for the measured residuals within
$\pm\sigma_\omega$ (H2). Any unexplained remainder is reported as such.
No bound is adjusted after gated data exists.

## Out of scope *(declared a priori)*

- Implementing 2PN or EIH physics (separate work).
- Claims against observation (the budget characterises the gate, not
  nature).
- Cross-platform re-verification (the existing portfolio covers it).
- Performance.

## Reproducibility *(to be completed at run time)*

| Field | Value |
| --- | --- |
| apsis canonical commit | *(apparatus commit; this protocol is its ancestor)* |
| Scenario | `crates/apsis-1pn/tests/mercury_precession_gate.rs` constants |
| Pow-oracle inputs | controller-pow notebook capture (42,662 inputs), reused |
| Harness | `validation/mercury-1pn-error-budget/` |
