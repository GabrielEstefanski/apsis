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
wrong sign. (The paper quotes $-2.288\times 10^6$, the
full-potential quadrature value; the $\sim 0.04\,\%$ offset between
the osculating post-unwrap measurement here and the quadrature value
is the osculating-vs-geometric definition difference, not a
discrepancy.) The harness at
`crates/apsis-1pn/tests/softened_plummer_precession_validation.rs` re-measures the
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
  §3.5 (orbit equation and integrable power-law potentials), §3.6
  (apsidal angle and conditions for closed orbits).
- Landau, L. D., & Lifshitz, E. M. (1976). *Mechanics* (3rd ed.),
  §15. Perihelion-precession formula for a small perturbation to the
  Kepler potential; used for the first-order $\varepsilon^4$ next-order
  term below.
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

## Full-potential apsidal-angle quadrature (rigorous oracle)

The closed form above is the leading order in $\varepsilon^2$. The
apsidal precession for the *full* Plummer potential needs no expansion:
for a central potential the angle swept periapse-to-periapse is the
quadrature (Goldstein §3.5)

$$\Delta\varpi = 2\int_{r_\text{min}}^{r_\text{max}}\frac{L/r^2}{\sqrt{2\,(E-\Phi(r)) - L^2/r^2}}\,\mathrm{d}r - 2\pi,
\qquad \Phi(r) = -\frac{GM}{\sqrt{r^2+\varepsilon^2}},$$

with $E$, $L$ fixed by the §3.2 initial condition (periapsis
$r_0 = a(1-e)$, tangential vis-viva speed under the Kepler $\mu$). The
companion script `plummer_apsidal_quadrature.py` evaluates it (numpy
Gauss–Legendre + bisection; no scipy). At $\varepsilon = 0$ it
returns $\Delta\varpi = 1.7\times10^{-11}$ rad/orbit — the Kepler orbit
closes, so the quadrature is trustworthy to $\sim10^{-11}$ rad, eight
orders below the signal.

Being a spatial quadrature of the full potential, it is independent of
the LRL-vector $\omega$ definition, of Kepler-vs-radial sampling, and of
any time-integrator. It is the oracle the closed form and the
time-integrators are compared against below.

## Comparison vs measured

Earlier diagnostics in `softened_plummer_precession_validation.rs` returned the
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

The cleanest comparison is the precession *rate* (rad per unit time):
dividing each source by its own period removes the per-radial-vs-per-
Kepler "which orbit" ambiguity (the radial period exceeds the Kepler
period by $0.80\%$ here, an $O(\varepsilon^2)$ effect that would bias a
naive per-orbit comparison).

| Source | Rate (rad/time) | vs quadrature |
| --- | --- | --- |
| Quadrature (full Plummer) | $-1.76559\times10^{-2}$ | — (oracle) |
| Closed form (leading order) | $-1.81261\times10^{-2}$ | $+2.66\%$ |
| Apsis IAS15, 500 orbits | $-1.76624\times10^{-2}$ | $+0.04\%$ |
| scipy DOP853, 50 orbits | $-1.75543\times10^{-2}$ | $-0.58\%$ |

The integrator entries are their reported drift of the osculating
$\omega$ per Kepler period, converted to a rate (sources
`softened_plummer_precession_validation.rs` and `plummer_check.py`). The apsis run
carries Plummer + 1PN by design — the counter-test is precisely the
exactness-requiring 1PN operator on a softened kernel; the softening
artifact dominates the relativistic signal by $\sim 5\times10^4$, so the
1PN term is $\approx 2\times10^{-5}$ of the measured drift (below the
gate) and the comparison to the pure-Plummer quadrature holds. Apsis
reproduces the quadrature apsidal-precession rate to $0.04\%$ for this
orbit. The earlier "$\sim3\%$ gap vs the
closed form" is the closed form's own leading-order truncation
($+2.66\%$) plus scipy's residual ($-0.58\%$, over a shorter 50-orbit
window) — not an apsis error.

**Next-order structure.** An earlier draft mis-stated the next-order
term as $O(\varepsilon^4/a^4)\approx7\times10^{-6}$ — that is the
*square* of the expansion parameter. The correct relative next-order is
$O(\varepsilon^2/a^2)$: the rate sweep gives a clean $\varepsilon^2$ law
(closed/quadrature $-\,1 = +66.8\,\varepsilon^2 + \ldots$, flat across
$\varepsilon\in[0.002,\,0.02]$). The $O(\varepsilon^4)$ per-orbit
correction splits into two pieces:

- the first-order contribution of the $\varepsilon^4$ potential term
  $V_4 = -3GM\varepsilon^4/8r^5$, via the Landau–Lifshitz §15 precession
  formula $\Delta\varpi = \partial_L\!\oint \delta U\,\mathrm{d}t$ (which
  reproduces the leading term exactly):
  $$\Delta\varpi_a = \frac{15\pi\,(4+3e^2)\,\varepsilon^4}{8\,a^4(1-e^2)^4},$$
  confirmed against the quadrature — subtracting it cuts the
  $\varepsilon=0.02$ per-radial-period gap from $-1.81\%$ to $-1.06\%$;
- the second-order contribution of $V_2 = GM\varepsilon^2/2r^3$, isolated
  numerically as the residual (a clean $\varepsilon^4$ term), whose
  closed form requires second-order secular theory and is not derived
  here.

The quadrature captures the full precession (all orders in
$\varepsilon$) to its $\sim10^{-11}$ floor, so it — not the truncated
series — is the value the gate asserts against.

## Result for the paper

§3.2 reports the full-potential apsidal-angle quadrature as the oracle —
$-1.766\times10^{-2}$ rad/time for the full Plummer potential on this
orbit — against which apsis agrees to $0.04\%$. The leading-order closed
form $\Delta\varpi_\text{orbit} = -3\pi\varepsilon^2/[a^2(1-e^2)^2]$ sits
$+2.66\%$ above it, a quantified $O(\varepsilon^2)$ next-order effect
(the piece $\Delta\varpi_a$ derived above, plus a clean $\varepsilon^4$
residual) — a converging approximation, not a discrepancy.

The gate `plummer_drift_matches_quadrature` in
`crates/apsis-1pn/tests/softened_plummer_precession_validation.rs` asserts apsis against
the quadrature oracle (reference value pinned from the
`plummer_apsidal_quadrature.py` script), so the bound reflects the
measurement's own precision rather than absorbing the closed form's
truncation. The per-orbit unwrap in `measure_drift_arcsec_per_century`
recovers the full cumulative drift the earlier mod-$2\pi$ diagnostic
aliased.
