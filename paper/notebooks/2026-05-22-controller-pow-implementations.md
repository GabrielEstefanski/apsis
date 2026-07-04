# 2026-05-22 — Controller `pow` implementations: oracle comparison on the Mercury 1PN workload

**Subject:** Routing the IAS15 adaptive step-size ratio through `libm::pow` shifted the Mercury 1PN perihelion advance measured by `apsis-1pn` over 500 orbits from $4.4 \times 10^{-6}$ relative error (Windows UCRT) to $2.8 \times 10^{-5}$ (cross-platform `libm`/glibc). This notebook compares the two `pow` implementations against an IEEE-754 correctly-rounded oracle on the controller's actual input distribution and traces the result into the integrated trajectory.

---

## Motivation

The IAS15 adaptive step-size controller (Rein & Spiegel 2015 §2.3) chooses the next step size from
$$
\Delta t_\text{next} = \Delta t_\text{current} \cdot 0.9 \cdot \left(\frac{\varepsilon}{\text{err}}\right)^{1/7}
$$
implemented in `crates/apsis/src/physics/integrator/ias15.rs::optimal_dt` as a single transcendental call. Routing that call through the pure-Rust `libm` crate obtains bit-identical trajectory output across heterogeneous x86_64 hosts; the same change moved the Mercury 1PN error reported by `crates/apsis-1pn/examples/mercury_perihelion.rs` by 0.002 arcsec/century.

Either UCRT `pow` is non-conformant with the IEEE-754 specification for transcendentals, in which case the `libm` adoption corrected an implementation defect, or both implementations satisfy the specification (correctly-rounded transcendentals are recommended but not required by the standard), in which case the trajectory shift reflects 1-ULP rounding distributions propagating through the adaptive controller. The cases imply different user-facing claims. The experiment below discriminates by direct measurement on the controller's actual input distribution.

---

## Method

### Workload capture

The `optimal_dt` call site was temporarily instrumented to emit each `(eps/err)` argument value to stderr during a release-mode run of `cargo run --release --example mercury_perihelion -p apsis-1pn` on Windows AMD Zen 4 (Ryzen 5 7600X), Rust 1.94.1. The integration is Sun + Mercury at canonical orbital elements ($a = 0.387098$ AU, $e = 0.20563$) under IAS15 + `PostNewtonian1PN::from_raw_c(C_SOLAR_UNITS)` for 500 Mercury orbits with the initial step seed $\Delta t = 10^{-3}$. The capture yielded 42,662 unique controller-input values, each accepted-substep `(eps/err)` produced by the controller's error estimate over the integration.

### Comparison primitives

For each captured input $x_i$, three computations:

1. `f64::powf(x_i, 1.0/7.0)` — Rust standard, routes to UCRT `pow` on Windows.
2. `libm::pow(x_i, 1.0/7.0)` — `libm` crate v0.2 (pure-Rust port of MUSL libm).
3. Oracle: `mpmath.power(mpf(x_i), mpf(1)/mpf(7))` at 60-digit decimal precision, then `float()` to f64 — produces the IEEE-754 round-to-nearest-even f64 result, the unique correctly-rounded value the standard recommends.

All three values compared by bit pattern (`f64::to_bits`); ULP distance measured as signed difference of bit patterns.

### Statistical aggregation

Counts of bit-exact matches and ULP distance distributions aggregated over the full 42,662-input set produced by the Mercury 1PN integration.

---

## Results

### Bit-exact match against IEEE-754 oracle

| implementation | exact match | $-1$ ULP | $+1$ ULP | total off-by-one |
| --- | --- | --- | --- | --- |
| Windows UCRT (`f64::powf`) | 41,367 (96.97 %) | 1,295 (3.03 %) | 0 (0.00 %) | 3.03 % |
| `libm` crate (`libm::pow`) | 40,652 (95.29 %) | 1,667 (3.91 %) | 343 (0.80 %) | 4.71 % |

### Mutual agreement

UCRT and `libm` produce bit-identical output on 40,747 of 42,662 inputs (95.51 %). On the 1,915 disagreements:

| disagreement category | count |
| --- | --- |
| UCRT matches oracle, `libm` does not | 1,315 |
| `libm` matches oracle, UCRT does not | 600 |
| Both off by the same ULP direction | 695 |
| Both off by different ULP directions | 0 |

### Mercury trajectory under each implementation

Same example, two `optimal_dt` variants:

| controller `pow` | rate (arcsec/century) | relative error vs analytical $6\pi GM_\odot N / [c^2 a (1-e^2)]$ | absolute error |
| --- | --- | --- | --- |
| Windows UCRT (`f64::powf`) | 42.993 | $+4.439 \times 10^{-6}$ | $+1.5 \times 10^{-3}$ arcsec |
| `libm::pow` | 42.991 | $-2.802 \times 10^{-5}$ | $-1.5 \times 10^{-3}$ arcsec |

---

## Interpretation

Both `pow` implementations are within the IEEE-754 specification for transcendentals (the standard permits but does not require correctly-rounded results). UCRT errs only toward zero (1,295 cases at $-1$ ULP, 0 at $+1$ ULP); `libm` errs in both directions with residual $-1$ ULP bias (1,667 at $-1$ ULP, 343 at $+1$ ULP). UCRT matches the oracle on 1.7 points more inputs than `libm` on this distribution.

Neither implementation is defective. The Windows UCRT result of $4.4 \times 10^{-6}$ reflected UCRT's specific rounding distribution interacting with the IAS15 controller's substep selection over 42,662 accepted substeps; the accumulated trajectory residual against the closed-form analytical formula evaluated in f64 was smaller for UCRT than for `libm` on this scenario.

The Mercury divergence between the two `pow` implementations is propagation of 1-ULP rounding differences through 42,662 IAS15 controller decisions. The accumulated shift is 0.002 arcsec/century: per-call rounding within IEEE-754 tolerance for transcendentals propagates to $\mathcal{O}(10^{-5})$ cumulative trajectory error over 500 orbits. This bounds the controller's sensitivity to last-ULP rounding noise at this scenario's natural step cadence.

The cross-platform-consistent value, reproduced bit-identically by `libm` (this experiment) and glibc (Phase 1 of `paper/notebooks/2026-05-20-cross-platform-determinism.md`), is $2.8 \times 10^{-5}$. Both $4.4 \times 10^{-6}$ and $2.8 \times 10^{-5}$ are several orders of magnitude below the current observational precision of Mercury's perihelion advance (Verma & Fienga 2014) and below the experimental constraint on the PPN $\gamma$ parameter from Cassini Doppler tracking (Bertotti et al. 2003).

---

## Reproducibility

### Workload capture

Apply the following diff to `crates/apsis/src/physics/integrator/ias15.rs::optimal_dt`:

```diff
-        let ratio = libm::pow(self.epsilon / err, 1.0 / 7.0);
+        let x = self.epsilon / err;
+        let ratio = libm::pow(x, 1.0 / 7.0);
+        eprintln!("RATIO_CAPTURE x={:.18e} libm={:.18e}", x, ratio);
```

Then:

```bash
cargo run --release --example mercury_perihelion -p apsis-1pn 2> capture.txt
grep RATIO_CAPTURE capture.txt | awk '{print $2}' | sed 's/x=//' > x_inputs.txt
```

### `pow` comparison

A standalone Rust binary reads each captured `x` from stdin and prints both `f64::powf(x, 1.0/7.0)` and `libm::pow(x, 1.0/7.0)` with their f64 bit patterns. Oracle and aggregation in Python:

```python
import csv, struct
from mpmath import mpf, mp, power
mp.dps = 60
ONE_SEVENTH = mpf(1) / mpf(7)

def f64_bits(x: float) -> int:
    return struct.unpack('<Q', struct.pack('<d', x))[0]

with open('pow_results.csv') as f:
    for row in csv.DictReader(f):
        x = float(row['x_dec'])
        oracle_bits = f64_bits(float(power(mpf(x), ONE_SEVENTH)))
        ucrt_bits = int(row['ucrt_hex'], 16)
        libm_bits = int(row['libm_hex'], 16)
        # ucrt_bits - oracle_bits, libm_bits - oracle_bits give signed ULP distances
```

### Mercury runs under each `pow`

The `libm` result reproduces with the adaptive step-size ratio routed through `libm::pow` (Mercury rel = $-2.802 \times 10^{-5}$). The UCRT result reproduces with the ratio routed through `f64::powf` on Windows (Mercury rel = $+4.439 \times 10^{-6}$ on Windows AMD Zen 4). Both values are deterministic across runs on identical inputs.

---

## Implications

The Mercury 1PN agreement cited in v0.1 user-facing documentation updates from 4.4 ppm to 28 ppm against the closed-form GR prediction. Cross-platform reproducibility (`paper/notebooks/2026-05-20-cross-platform-determinism.md`) is preserved. The previously-cited value is attributable to the Windows UCRT operating point.

> **Superseded (2026-06).** ADR-014 moved the `for_units` $c$, and the exact-finish-time fix (ADR-015) removed the endpoint-sampling term that dominated both values above. The documented residual is now $4.6\times10^{-6}$, CI-gated below $9.2\times10^{-6}$; derivation in `2026-06-10-mercury-1pn-error-budget.md`. The ULP methodology and the UCRT-vs-libm comparison in this notebook are unaffected.

The method (workload capture, comparison against an arbitrary-precision oracle rounded to f64, ULP distribution aggregation) applies to any adaptive integrator or operator using transcendentals on an integration-critical path.

---

## References

- Rein, H. & Spiegel, D. S. (2015). *IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits.* MNRAS **446**, 1424. §2.3 for the step-size controller formula.
- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics.* Cambridge University Press.
- Bertotti, B., Iess, L., & Tortora, P. (2003). *A test of general relativity using radio links with the Cassini spacecraft.* Nature **425**, 374. PPN $\gamma$ constraint.
- Verma, A. K., & Fienga, A. (2014). *Use of MESSENGER radioscience data to improve planetary ephemeris and to test general relativity.* A&A **561**, A115. Modern Mercury perihelion observational precision.
- IEEE 754-2019. *IEEE Standard for Floating-Point Arithmetic.* §9.2 on recommended (not required) correctly-rounded transcendentals.
- Previous lab notebook: `paper/notebooks/2026-05-20-cross-platform-determinism.md` (Phase 1 Linux glibc result; cross-platform bit-equality demonstration).
- `crates/apsis/src/physics/integrator/ias15.rs::optimal_dt`.
- `crates/apsis-1pn/examples/mercury_perihelion.rs`.
- `crates/apsis-1pn/tests/mercury_precession_gate.rs` — CI gate at 100 ppm tolerance, passes for both `pow` variants.
