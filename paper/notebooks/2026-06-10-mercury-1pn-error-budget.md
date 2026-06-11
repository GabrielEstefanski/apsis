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

**Status:** *Protocol draft — phases and hypotheses declared a priori.
Phase A derivations to be carried out in this notebook before any
apparatus code; phase-B gates freeze when phase A closes.*

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

$$\varepsilon_\text{Mercury} = \frac{1}{10065.32^2 \times 0.387098 \times (1 - 0.20563^2)} \approx 2.7 \times 10^{-8}.$$

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
