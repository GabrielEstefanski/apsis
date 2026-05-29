# REBOUNDx parity — Sun–Mercury 1PN (apsis-1pn vs REBOUNDx `gr`)

## Question

Does apsis-1pn's first post-Newtonian operator reproduce the perihelion
precession of an independent, peer-reviewed code? This is the first parity
check against a **REBOUNDx effect** (the existing portfolio compares only the
base REBOUND integrator). The harness lives at
`validation/reboundx-parity/gr-mercury/`.

## The two formulations (why this is not a bit-parity check)

- **apsis-1pn:** Schwarzschild test-particle 1PN, applied pairwise in the
  **inertial** frame, explicit:
  $\mathbf{a}_{1\text{PN}}(i\leftarrow j) = \frac{G m_j}{c^2 r^2}\left[(4 G m_j/r - v_i^2)\,\hat{\mathbf{n}} + 4(\hat{\mathbf{n}}\cdot\mathbf{v}_i)\,\mathbf{v}_i\right]$
  (Anderson et al. 1975).
- **REBOUNDx `gr`:** single-dominant-mass 1PN in **Jacobi** coordinates with an
  iterative velocity solve, referenced to $\mu = G m_0$.

They differ in gauge, coordinate system and solve method. For Sun–Mercury
(mass ratio $\sim 1.7\times10^{-7}$) both reduce to the test-particle limit, so
the comparison **measures the formulation/gauge difference**, which we cannot
bound analytically (it would need the PN gauge transformation between the two
forms — not derived here). The gate is therefore set from the *measured* floor,
and the analytic 43″/century is the physics anchor.

## Setup

| Quantity | Value |
| --- | --- |
| Bodies | Sun ($m=1$), Mercury ($m=1.660114\times10^{-7}$), COM frame |
| $a$, $e$ | 0.387098 AU, 0.20563 |
| Units | solar-canonical (AU, yr/2π, M⊙, G=1) |
| $c$ | $1.006513002441656681\times10^4$ (C_SOLAR_UNITS; **bit-matched** both sides) |
| Integrator | IAS15 (both), $dt_0 = 10^{-3}\,T$ |
| Window | 500 orbital periods |

The Sun-side and reboundx-side use the same f64 IC expressions; `c` is parsed
from the apsis CSV header so the two match bit-for-bit. The reboundx side lands
at apsis's *actual* sample times (`exact_finish_time=1`), removing "sampled at
slightly different physical times" as a spurious source of difference.

## Metrics and a-priori tolerances

Per the rebound-parity lesson
(`validation/rebound-parity/kepler/compare.py`): adaptive IAS15 is not
bit-deterministic across implementations, so $|\Delta\mathbf{r}|$ at fixed times
conflates orbital **phase** drift with geometric drift. The gated metrics are
therefore orbital invariants and the **secular precession rate** (the
gauge-invariant 1PN observable); $|\Delta\mathbf{r}|$ is informational.

| Gated metric | Tolerance | Rationale |
| --- | --- | --- |
| $\lvert\Delta a\rvert/a$, $\lvert\Delta e\rvert$, $\lvert\Delta h\rvert/h$ cross | $10^{-13}$ | ULP floor, as the rebound-parity scenarios |
| $\lvert\Delta E\rvert/E_0$ cross | $10^{-13}$ | idem |
| precession apsis-vs-reboundx (gr) | $2\times10^{-5}$ | empirical formulation floor $\sim7\times10^{-7}$, $\sim$30× headroom |

Reported (not gated here): precession vs analytic Schwarzschild (accuracy is
owned by `crates/apsis-1pn/tests/mercury_precession_gate.rs`, 28 ppm), per-side
Newtonian-energy drift (1PN does not conserve it; both drift identically), and
$\lvert\Delta\mathbf{r}\rvert$.

The precession tolerance is **empirical** ($7\times10^{-7}\to2\times10^{-5}$);
tighten once cross-platform data exists.

## Measured results

**1PN-off control** (validates the harness): all cross-implementation
invariants and energy at the ULP floor ($\sim5\times10^{-15}$); precession
$\approx 0$ (Kepler does not precess). The harness, ICs and sampling are
correct independent of the 1PN physics.

**gr (1PN on):**

| | apsis | REBOUNDx | analytic |
| --- | --- | --- | --- |
| apsidal precession (″/century) | $+42.9783$ | $+42.9784$ | $+42.9824$ |

- apsis vs REBOUNDx precession: $7.2\times10^{-7}$ (gate $2\times10^{-5}$).
- apsis vs analytic: $-9.3\times10^{-5}$; REBOUNDx vs analytic: $-9.3\times10^{-5}$
  (both share the same $dt$-limited IAS15 floor at this $dt_0$).
- Cross-implementation invariants ($a,e,h$) and energy at the ULP floor
  ($\sim7\times10^{-15}$).
- $\lvert\Delta\mathbf{r}\rvert \sim 10^{-10}$ (phase-contaminated, informational).

## Conclusion

apsis-1pn and REBOUNDx `gr` — independent codes with different formulation,
gauge, coordinates and solve — agree on Mercury's perihelion precession to
$7\times10^{-7}$, and both agree with the analytic Schwarzschild advance to
$\sim10^{-4}$. The formulation/gauge difference is empirically below
$10^{-6}$ for this orbit, dominated by the mass ratio and the integrator
floor; it vanishes in the test-particle limit, which is apsis-1pn's valid
regime.

## Scope

This confirms correctness in the **test-particle regime**. It does not stress
comparable-mass post-Newtonian dynamics — the apsis-1pn pairwise form is not
valid there and warns at registration. A comparable-mass check would compare
against REBOUNDx `gr_full`, out of scope here.

## Reproducibility

Linux only (REBOUNDx does not build on Windows/MSVC). From the scenario
directory, with `cargo` and a reboundx venv:

```bash
pip install --no-cache-dir -r requirements.txt
python run.py
```

Frozen reference CSVs and `comparison_{gr,control}.json` are committed under
`out/`. rebound 4.6.0, reboundx 4.6.2.
