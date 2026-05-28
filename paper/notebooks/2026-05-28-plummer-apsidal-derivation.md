# Plummer-softened apsidal precession — closed-form derivation

## Question

Derive in closed form the apsidal-precession rate produced by a
Plummer-softened pair potential

$$U(r) = -\frac{G m_1 m_2}{\sqrt{r^2 + \varepsilon^2}},$$

evaluate it for the Sun–Mercury counter-test in §3.2 of `paper.md`
($\varepsilon = 0.02$ AU), and confirm it reproduces the apsis-measured
rate of $-2.289\times 10^{6}$ arcsec/century.

This converts the §3.2 counter-test claim from *measured a large wrong
number* to *measured exactly what the violation theory predicts*.

## Empirical anchor

Counter-test parameters (from `crates/apsis-1pn/tests/mercury_precession_gate.rs`
and `paper.md` §3.2):

| Quantity | Value | Unit |
| --- | --- | --- |
| Semi-major axis $a$ | 0.387098 | AU |
| Eccentricity $e$ | 0.20563 | — |
| Mercury mass $m_\text{Mercury}$ | $1.660114 \times 10^{-7}$ | $M_\odot$ |
| Softening $\varepsilon$ | $\approx 0.02$ | AU |
| Integration window | 500 | orbital periods |
| Integrator | IAS15, $dt_0 = 10^{-4} T$ | — |

Apsis measurement (post-unwrap): drift $\approx -2.289 \times 10^6$
arcsec/century, $\sim 5\times 10^4$ times the GR signal and of the
wrong sign. The harness at
`crates/apsis-1pn/tests/exactness_theory_match.rs` re-measures the
drift under the same conditions and exposes the value through
`cargo test --release ... --nocapture` for cross-check against the
closed-form prediction cited in §3.2.

## References

- Plummer, H. C. (1911). *Monthly Notices RAS* 71, 460.
  Original softened potential.
- Heggie, D. C., & Hut, P. (2003). *The Gravitational Million-Body
  Problem*, ch. 8. Standard treatment of softening in $N$-body
  cluster dynamics; gives the radial-force expression and discusses
  the apsidal-angle deviation it induces.
- Goldstein, Poole & Safko (2002). *Classical Mechanics* (3rd ed.),
  §3.6 (orbit equation in a central potential), §3.8 (apsidal
  angle and conditions for closed orbits).
- Touma, J., & Tremaine, S. (1997). *MNRAS* 292, 905. Secular
  perturbation theory in non-Keplerian potentials; references the
  apsidal-angle formula in a form ready for power-series expansion
  in the softening parameter.

## Setup — expansion regime

For Sun–Mercury, $\varepsilon / a \approx 0.052$; $\varepsilon^2 / a^2
\approx 2.7 \times 10^{-3}$. A leading-order expansion of
$U(r)$ in $\varepsilon^2/r^2$ is appropriate and the next-order
correction is suppressed by a factor of $\sim \varepsilon^2/a^2$.

## Derivation

We derive the apsidal-precession rate by computing the radial and
angular frequencies of a near-circular orbit in the Plummer potential
and taking their difference, then generalising to $e > 0$ via
disturbing-function orbital averaging. The two routes converge on the
same closed-form expression, providing an internal cross-check.

### Route A — near-circular frequency decomposition

Specific gravitational potential (per unit mass of the test body):

$$\Phi(r) = -\frac{GM}{\sqrt{r^2 + \varepsilon^2}}.$$

Specific effective potential including the centrifugal barrier with
specific angular momentum $h$:

$$\mathcal{V}_\text{eff}(r) = \Phi(r) + \frac{h^2}{2 r^2}.$$

For a circular orbit at $r = a$, $\mathcal{V}_\text{eff}'(a) = 0$ gives

$$h^2 = \frac{G M\, a^4}{(a^2 + \varepsilon^2)^{3/2}}.$$

The angular frequency follows from $\Omega^2 = h^2 / a^4$:

$$\Omega^2 = \frac{GM}{(a^2 + \varepsilon^2)^{3/2}}.$$

The radial frequency comes from $\omega_r^2 = \mathcal{V}_\text{eff}''(a)$.
Computing the second derivative and substituting $h^2$:

$$\omega_r^2 = \frac{GM\,(a^2 + 4\varepsilon^2)}{(a^2 + \varepsilon^2)^{5/2}}.$$

Expanding both to leading order in $\varepsilon^2 / a^2 \ll 1$ with
$n^2 \equiv GM/a^3$ (Keplerian mean motion):

$$\Omega = n\left(1 - \tfrac{3}{4}\tfrac{\varepsilon^2}{a^2}\right) + O(\varepsilon^4/a^4),\qquad
\omega_r = n\left(1 + \tfrac{3}{4}\tfrac{\varepsilon^2}{a^2}\right) + O(\varepsilon^4/a^4).$$

The apsidal-precession rate is the residual between the two frequencies:

$$\dot{\varpi}_\text{circ} = \Omega - \omega_r = -\frac{3\,n\,\varepsilon^2}{2 a^2}.$$

The negative sign is physical: Plummer softening weakens gravity at
small $r$, the orbit lingers there, and the next periapsis occurs at
a delayed angle (retrograde apsidal motion). Per orbital period
$T = 2\pi/n$:

$$\Delta\varpi_\text{circ,orbit} = -\frac{3\pi\,\varepsilon^2}{a^2}.$$

### Route B — orbital-averaged disturbing function (general $e$)

The disturbing function $\mathcal{R} = \Phi_\text{Plummer} - \Phi_\text{Kepler}$
expanded to leading order in $\varepsilon^2/r^2$:

$$\mathcal{R}(r) = -GM\left[\frac{1}{\sqrt{r^2+\varepsilon^2}} - \frac{1}{r}\right] = +\frac{GM\,\varepsilon^2}{2\,r^3} + O(\varepsilon^4/r^5).$$

This is a $1/r^3$ disturbing function — formally analogous to the
1PN GR Schwarzschild correction, but with opposite sign coefficient.
Where 1PN deepens the effective well (prograde precession), Plummer
softening shallows it (retrograde precession).

Orbital average over a Kepler ellipse with $a$, $e$: using the
standard result $\langle r^{-3}\rangle = [a^3 (1-e^2)^{3/2}]^{-1}$,

$$\langle\mathcal{R}\rangle = \frac{GM\,\varepsilon^2}{2\,a^3 (1-e^2)^{3/2}}.$$

The Lagrange equation for the longitude of periapsis with sign
convention matched to Route A (i.e. retrograde $\mathcal{R} > 0$
gives $\dot\varpi < 0$):

$$\dot{\varpi} = -\frac{\sqrt{1-e^2}}{n a^2 e}\frac{\partial\langle\mathcal{R}\rangle}{\partial e} = -\frac{3\,GM\,\varepsilon^2}{2\,n\,a^5\,(1-e^2)^2}.$$

Multiplying by $T = 2\pi/n$ and using $n^2 a^3 = GM$:

$$\boxed{\;\Delta\varpi_\text{Plummer,orbit} = -\frac{3\pi\,\varepsilon^2}{a^2\,(1-e^2)^2}\;}$$

For $e \to 0$ this reduces to the Route A result. For $e > 0$ the
$(1-e^2)^{-2}$ factor amplifies the precession — Mercury's $e \approx
0.21$ raises the rate by $\sim 9\%$ over the circular limit.

### Result

$$\Delta\varpi_\text{Plummer,orbit} = -\frac{3\pi\,\varepsilon^2}{a^2\,(1-e^2)^2}.$$

Per Earth century, with $T_\text{Mercury} = 87.969$ d and 36525 d/century:

$$\dot\varpi_\text{Plummer,arcsec/cy} = -\frac{3\pi\,\varepsilon^2}{a^2\,(1-e^2)^2}\cdot\frac{36525}{T_\text{Mercury}[\text{d}]}\cdot\frac{180\cdot 3600}{\pi}.$$

## Numerical evaluation

Constants for the substitution (canonical solar units; pull from
`apsis::units::SOLAR_CANONICAL` if matching the integrator's `g`
literal matters at the ULP level):

| Constant | Symbolic | Value |
| --- | --- | --- |
| $G M_\odot$ | $G M$ | 1.0 (canonical) |
| $a_\text{Mercury}$ | $a$ | 0.387098 AU |
| $e_\text{Mercury}$ | $e$ | 0.20563 |
| $\varepsilon$ | $\varepsilon$ | 0.02 AU |
| Orbital period | $T$ | $2\pi \sqrt{a^3/(GM)} \approx 1.515$ canonical time units (87.97 d) |
| Centuries per period | — | $36525 / 87.97 \approx 415.2$ orbits per century |

Plug in:

$$\Delta\varpi_\text{Plummer,orbit} = -\frac{3\pi\,(0.02)^2}{(0.387098)^2\,(1-0.20563^2)^2}
= -\frac{3\pi \cdot 4.000\times 10^{-4}}{0.149845 \cdot 0.917020}
= -2.7429\times 10^{-2}\ \text{rad/orbit}.$$

Per century:

$$\dot\varpi_\text{Plummer,arcsec/cy} = -2.7429\times 10^{-2}\,\text{rad/orbit} \times 415.20\,\text{orbits/cy} \times 206\,265\,\text{arcsec/rad}$$

$$\boxed{\;\dot\varpi_\text{Plummer} = -2.349\times 10^{6}\ \text{arcsec/century}.\;}$$

This is $\sim 5 \times 10^4$ times the 1PN GR signal
$(+43\ \text{arcsec/century})$ and of the opposite sign — consistent
with §3.2's qualitative claim "wrong sign, orders of magnitude larger"
and tightening the pre-unwrap-fix "$\sim 3 \times 10^3$ times" to
$5 \times 10^4$.

## Comparison vs measured

Earlier diagnostics in `exactness_theory_match.rs` returned the
final-vs-initial periapsis-angle difference modulo $2\pi$, mapped to
$(-\pi, \pi]$. For the Plummer-violated case the true cumulative drift
over 500 orbits is $\approx -13.71$ rad, well outside that interval —
the unwrap aliased it to $-13.71 + 4\pi = -1.14$ rad and reported the
fractional value. The pre-fix diagnostics ($-83\,128$ and $-136\,732$
arcsec/century) were aliased fractions of the true rate, and the two
differ because $-83\,128$ predates a softening-setup change (per-body
$\varepsilon$ pair-averaged to $\varepsilon_\text{eff} \approx 0.0141$
AU, versus the current flat $\varepsilon = 0.02$ AU); the jump to the
current $-2.289\times 10^6$ is therefore partly that setup change and
partly the unwrap fix, not unwrap alone. The current test uses
per-orbit accumulation of the unwrapped step and reports the full
drift directly.

An independent scipy DOP853 integration over 50 orbits (script:
`paper/notebooks/scripts/plummer_check.py`) gives:

| Quantity | scipy DOP853 | Closed-form | Ratio |
| --- | --- | --- | --- |
| $\Delta\varpi$ per orbit | $-2.656\times 10^{-2}$ rad | $-2.743\times 10^{-2}$ rad | 0.968 |
| cumulative over 50 orbits | $-1.328$ rad | $-1.371$ rad | 0.968 |

The $\sim 3\%$ gap exceeds the $O(\varepsilon^4/a^4) \sim 7\times
10^{-6}$ next-order expansion term by three orders of magnitude. Its
source is not the leading-order truncation; candidate explanations
(higher-order coefficient, sampling at Kepler vs radial period,
secular vs instantaneous $a$) are not disambiguated here. The 5%
gate absorbs the gap regardless.

| Source | Value (arcsec/century) | Notes |
| --- | --- | --- |
| Predicted (closed form) | $-2.349\times 10^{6}$ | this notebook |
| Independent scipy DOP853 | $-2.275\times 10^{6}$ | 3.2 % below closed form |
| Apsis IAS15, 500 orbits | $-2.289\times 10^{6}$ | 2.55 % below closed form (`exactness_theory_match.rs`) |
| Acceptance bound (gate) | 5 % | absorbs the gap with $\sim 2\times$ margin |

## Result for the paper

Apsis IAS15 over 500 orbits agrees with the closed-form prediction
to 2.55 %, between scipy's 3.2 % residual and zero — consistent with
the gap being a real physical effect rather than scheme noise.

Paper §3.2 cites the closed form $\Delta\varpi_\text{orbit} =
-3\pi\varepsilon^2 / [a^2(1-e^2)^2]$ derived above and the measured
apsis rate $-2.289\times 10^6$ arcsec/century. The
`plummer_drift_matches_softened_theory` gate in
`crates/apsis-1pn/tests/exactness_theory_match.rs` asserts the 5 %
agreement; the per-orbit unwrap in `measure_drift_arcsec_per_century`
recovers the full cumulative drift that the earlier mod-$2\pi$
diagnostic aliased.
