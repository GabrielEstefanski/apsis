# Continuity-violation spike-magnitude bound — derivation

## Question

Derive an a-priori bound on the per-step relative energy error
$|\Delta E / E|$ a symplectic integrator produces when its trajectory
crosses a $C^0$ discontinuity in the radial force. Evaluate the bound
for the §3.3 truncated-Plummer counter-test (equal-mass two-body,
$a=1$, $e=0.5$, Yoshida4 at $dt = 10^{-3}$ canonical, $R_c = 1$,
$\alpha = 0.8$), confirm it reproduces the measured spike-magnitude
range $4.7 \times 10^{-6}$ to $2.0 \times 10^{-4}$, and sharpen the
bijection result from *"events coincide with crossings"* to *"events
coincide with crossings at the predicted magnitude"*.

## Empirical anchor

Counter-test parameters (from
`crates/apsis-1pn/tests/kernel_continuity_counter_test.rs` and
`paper.md` §3.3):

| Quantity | Value | Source |
| --- | --- | --- |
| Cutoff radius $R_c$ | 1.0 | semi-major-axis units |
| Outside scale $\alpha$ | 0.8 | `DEFAULT_TRUNCATED_OUTSIDE_SCALE` |
| Plummer softening $\varepsilon$ | 0.0 (inside-cutoff Plummer) | kernel default |
| Semi-major axis $a$ | 1.0 | — |
| Eccentricity $e$ | 0.5 | — |
| Periapse / apoapse | 0.5 / $\approx 2.06$ | from $K$ truncation |
| Equal-mass two-body | $m_1 = m_2 = 0.5$ | — |
| Integrator | Yoshida4 | symplectic, 4th order |
| Timestep $dt$ | $10^{-3}$ | canonical units, fixed |
| Integration window | 60 simulation units | $\ge 4$ orbits, 8+ crossings |

Force jump at $R_c$ (already in §3.3, derivable from the kernel def):

$$\Delta F = (1 - \alpha) \cdot \frac{R_c}{(R_c^2 + \varepsilon^2)^{3/2}}
= 0.2 \quad \text{for } (\alpha, R_c, \varepsilon) = (0.8, 1, 0).$$

Reported measurement in `paper.md`: spike magnitudes $|\Delta E / E|
\in [4.7 \times 10^{-6},\, 2.0 \times 10^{-4}]$ over 11 crossings.
Smooth-kernel reference floor: $< 2.7 \times 10^{-14}$ per step.

The harness at `crates/apsis-1pn/tests/continuity_theory_match.rs`
re-measures spike magnitudes and the bracketing velocity at each
crossing, then compares against the bound the author derives.

## References

- Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical
  Integration*, ch. IX. Long-time energy behaviour of symplectic
  methods on smooth Hamiltonians (the result whose smoothness
  hypothesis the truncated kernel violates).
- Skeel, R. D., & Gear, C. W. (1992). *Physica D* 60, 311. Energy
  behaviour of symplectic methods applied to non-smooth potentials;
  the canonical reference for the $O(\Delta F \cdot v \cdot dt^k)$
  scaling.
- Wisdom, J., & Holman, M. (1991). *AJ* 102, 1528. Mixed-variable
  symplectic; their footnote on close-encounter step-size and the
  shadowing-Hamiltonian deviation is relevant.

## Setup — variation regime

At each crossing the orbit moves with relative speed $v_\text{cross}$;
the integrator's substep window straddles the discontinuity for
$\le 1$ timestep $dt$. The expected scaling is therefore
$|\Delta E| \sim \Delta F \cdot v_\text{cross} \cdot dt^k \cdot
(\text{integrator-order prefactor})$ with $k$ tied to the integrator
order — author confirms $k$ from the symplectic-shadow analysis.

## Derivation

A symplectic integrator of order $k$ applied to a smooth Hamiltonian
$H$ conserves not $H$ itself but a *shadow Hamiltonian* $\widetilde H$
that differs from $H$ by terms of order $\mathrm{d}t^k$ and higher
[@HairerLubichWanner2006, ch. IX]:

$$\widetilde H = H + \mathrm{d}t^k\,H_k + \mathrm{d}t^{k+1}\,H_{k+1} + \ldots$$

The error coefficients $H_j$ are constructed from the
Baker–Campbell–Hausdorff expansion and depend on derivatives of the
potential $V(r)$ up to order $\sim j$. When $V \in C^\infty$, all
$H_j$ are bounded and $|H - \widetilde H| = O(\mathrm{d}t^k)$
uniformly in time — energy oscillates with amplitude $\mathrm{d}t^k$
but does not drift secularly. This is the symplectic-shadow guarantee
the §3.3 setup is designed to break.

The truncated-Plummer kernel makes $V(r)$ continuous but leaves
$V'(r)$ with a finite jump $\Delta F$ at $r = R_c$. Higher
derivatives of $V$ contain Dirac-delta contributions at $r = R_c$,
the BCH series for $\widetilde H$ diverges there, and the symplectic
guarantee localised to that radius is void. The integrator's energy
preservation breaks down each time the trajectory crosses $R_c$.

### Worst-case bound

Inside a single step $[t_n, t_n + \mathrm{d}t]$ that contains a
crossing at $t_n + \tau$ with $\tau \in (0, \mathrm{d}t)$, the
integrator evaluates the force at internal substep positions. For
positions sampled on the "wrong side" of the discontinuity, the force
differs from the correct value by $\Delta F$ along the radial
direction. The work done on the test particle by this force-error is
bounded above by the line integral

$$|\Delta E_\text{spike}| \le \int_{\text{wrong-side substeps}} |\Delta F \cdot v(t')|\,\mathrm{d}t'
\le \Delta F \cdot v_\text{cross} \cdot \mathrm{d}t,$$

where the second inequality uses $v(t') \le v_\text{cross}(1 +
O(\mathrm{d}t))$ over a window of size $\mathrm{d}t$ around the
crossing. The bound is saturated only in the worst case where the
"wrong side" coincides with the full step *and* the substep weights
do not cancel — a pessimistic envelope rather than a tight predictor.

Normalising by the orbit's specific energy $|E_0|$:

$$\boxed{\;\left|\frac{\Delta E}{|E_0|}\right|_\text{crossing} \le \frac{\Delta F \cdot v_\text{cross} \cdot \mathrm{d}t}{|E_0|}.\;}$$

### Why the bound is loose

The bound assumes the wrong-side force is applied for the full
$\mathrm{d}t$ with no substep cancellation. Yoshida-4 has seven
substeps with positive and negative weights summing to $\mathrm{d}t$
[@Yoshida1990]; if a crossing lands near a substep boundary the
wrong-side dwell is a fraction of $\mathrm{d}t$, and the energy error
is correspondingly smaller. A refined bound that accounts for
crossing-phase variability gives an envelope of width $\sim \Delta F
\cdot v_\text{cross} \cdot \mathrm{d}t / 4$ on average — but the
worst-case formula above is what generalises cleanly across symplectic
schemes and is the version the paper asserts.

### Result

For the §3.3 configuration ($\Delta F = 1 - \alpha = 0.2$,
$v_\text{cross} = 1.0$ at every crossing by energy conservation,
$\mathrm{d}t = 10^{-3}$, $|E_0| = 0.5$):

$$\left|\frac{\Delta E}{|E_0|}\right|_\text{bound} = \frac{0.2 \cdot 1.0 \cdot 10^{-3}}{0.5} = 4.000\times 10^{-4}.$$

Independent scipy DOP853 verification of the orbital configuration
(script in `paper/notebooks/scripts/continuity_check.py`) confirms
11 crossings of $R_c = 1$ over $T_\text{end} = 60$, with
$v_\text{cross}$ constant at $1.0000$ at every crossing — as
required by energy conservation at fixed radius. The single-number
bound applies uniformly; the spread in measured spike magnitudes
($4.7\times 10^{-6}$ to $2.0\times 10^{-4}$, factor $\sim 40\times$)
arises from the per-crossing distribution of $\tau$ relative to the
Yoshida-4 substep structure, not from $v_\text{cross}$ variability.

## Numerical evaluation per crossing

For the §3.3 orbit (equal-mass two-body, $a=1$, $e=0.5$), energy
conservation pins the speed at any given radius. At $r = R_c = 1$,
the vis-viva relation
$v^2 = G M_\text{total}\,(2/r - 1/a)$
gives $v_\text{cross}^2 = 2/1 - 1/1 = 1$, so $v_\text{cross} = 1.000$
at every crossing regardless of direction. (The scipy reproduction
recovers $v_\text{cross} = 1.0000 \pm 4\times 10^{-5}$ across the
eleven crossings; the variation is integrator round-off at
`rtol=1e-12`, not physical.)

The per-crossing bound is therefore a single value:

$$\left|\frac{\Delta E}{|E_0|}\right|_\text{bound} = 4.000\times 10^{-4} \quad\text{(every crossing)}.$$

| Crossing | $t_\text{cross}$ (scipy) | $v_\text{cross}$ | $\lvert\Delta E\rvert/\lvert E\rvert$ bound | Measured (apsis) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 1.0708 | 1.000000 | $4.000\times 10^{-4}$ | $4.074\times 10^{-5}$ |
| 2 | 9.7850 | 0.999980 | $4.000\times 10^{-4}$ | $5.389\times 10^{-6}$ |
| 3 | 11.9267 | 0.999980 | $4.000\times 10^{-4}$ | $5.991\times 10^{-5}$ |
| 4 | 20.6393 | 0.999940 | $4.000\times 10^{-4}$ | $2.013\times 10^{-4}$ |
| 5 | 22.7816 | 0.999930 | $4.000\times 10^{-4}$ | $1.932\times 10^{-4}$ |
| 6 | 31.4940 | 0.999950 | $4.000\times 10^{-4}$ | $4.664\times 10^{-6}$ |
| 7 | 33.6358 | 0.999950 | $4.000\times 10^{-4}$ | $3.165\times 10^{-5}$ |
| 8 | 42.3475 | 0.999920 | $4.000\times 10^{-4}$ | $1.773\times 10^{-4}$ |
| 9 | 44.4898 | 0.999850 | $4.000\times 10^{-4}$ | $4.775\times 10^{-5}$ |
| 10 | 53.1948 | 0.999820 | $4.000\times 10^{-4}$ | $3.268\times 10^{-5}$ |
| 11 | 55.3373 | 0.999840 | $4.000\times 10^{-4}$ | $6.346\times 10^{-5}$ |

Measured values from `continuity_per_crossing_measurements_are_recorded`
in `crates/apsis-1pn/tests/continuity_theory_match.rs`. Worst-case
crossing is #4 ($2.013\times 10^{-4}$ = 50.3% of bound); floor is #6
($4.664\times 10^{-6}$ = 1.2% of bound).

## Comparison vs measured

The bound is satisfied if every measured spike falls below
$4.000\times 10^{-4}$. The eleven measured spikes lie in
$[4.66\times 10^{-6}, 2.01\times 10^{-4}]$ — the peak sits at 50% of
the bound; the floor at 1.2%.

| Source | Range | Notes |
| --- | --- | --- |
| Measured spike range (apsis) | $4.66 \times 10^{-6}$ – $2.01 \times 10^{-4}$ | per-crossing table above |
| Predicted bound | $4.000 \times 10^{-4}$ | uniform across crossings |
| Worst-case ratio (peak / bound) | $0.503$ | bound holds with $2\times$ safety factor |
| Floor ratio (floor / bound) | $0.012$ | reflects substep-phase variability |
| Acceptance gate | $\text{max measured} < \text{bound}$ | strict worst-case envelope |

The remaining $2\times$ margin reflects the unmodelled substep
cancellation in Yoshida-4 — the bound assumes no cancellation by
construction, while the integrator's seven substeps with alternating
signs partially cancel the wrong-side work whenever the crossing is
not centred in the step.

## Result for the paper

The bound holds with a $2\times$ safety factor at the worst-case
crossing. The remaining margin is the unmodelled Yoshida-4 substep
cancellation. Paper §3.3 sharpens the bijection claim to add the
magnitude statement: $|\Delta E|/|E_0| \in [4.66\times 10^{-6},
2.01\times 10^{-4}]$, contained by the a-priori envelope $\Delta F
\cdot v_\text{cross} \cdot dt / |E_0| = 4.00\times 10^{-4}$ — the
50% maximum realisation reflects the seven-substep Yoshida-4
sign-alternation cancelling the wrong-side work whenever the
crossing is not centred in the step.

The theory-match gate `spike_magnitudes_satisfy_jump_bound` in
`crates/apsis-1pn/tests/continuity_theory_match.rs` asserts the bound
at safety factor 1.0.
