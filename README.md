# APSIS

*A Federated Model for Composable N-Body Force Artifacts*

[![License: Apache 2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

APSIS treats N-body force perturbations as first-class scientific
artifacts. The model is *federated* in Rust: each force is an
independent Cargo crate — developed, versioned, and cited
separately, without central integration into a monolithic codebase.
Python users see a single distribution: `pip install apsis` brings
the core simulator and official operator submodules (`apsis.gr`,
`apsis.radiation`, `apsis.central`, …) under one import. A simulation's
physical model is in its dependency graph: `Cargo.toml` declares
the forces, `Cargo.lock` pins them bit-precisely. The simulator is
infrastructure for composing those artifacts.

The numerical core ships seven integrators; IAS15 (Rein & Spiegel,
2015) is the default. Cross-implementation parity is tested at **1
ULP** of f64 (for the validated reference configuration) across four
canonical regimes (Kepler $e = 0.5$, Chenciner–Montgomery figure-8,
Pythagorean, $10^4$-orbit retrograde Kepler); under
[`apsis.gr`](crates/apsis-1pn/) (post-Newtonian 1PN), Mercury's
perihelion precession matches the analytic 1PN GR prediction to
within **28 ppm**, reproduced bit-identically across Windows and
Linux on x86_64.

> **Status.** Pre-release (`v0.1.0` alpha). 3D-aware physics core
> (Vec3, inclined orbits, 3D observables). Multiple published operator
> crates exercise the federation contract end-to-end. Public API
> stabilised but not yet tagged; citation DOI pending first Zenodo
> release.

---

## Statement of need

A force perturbation in a published N-body simulation lives in the
methods section of a paper and, sometimes, in a fork of an established
framework. The fork is not citable; the prose drifts; the next group
reimplements the same effect from scratch. The framework — REBOUND
(Rein & Liu, 2012), REBOUNDx (Tamayo et al., 2020), MERCURIUS
(Rein et al., 2019), NBODY6/7 (Aarseth, 2003) — is mature, citable,
and validated, but it absorbs every extension into a single binary
with one citation covering everything.

APSIS enables a publication path where perturbations become
independently versioned scientific artifacts. A force is a Cargo crate
declaring its physical preconditions on the gravitational kernel via
the `KernelRequirements` type — the 1PN crate declares
`exact_and_smooth()`; future crates declare a different combination
of exactness and continuity invariants depending on the physics. The
core matches the declared requirements against the active kernel at
`System::add_hamiltonian_perturbation` (or the non-conservative
counterpart) and emits a structured diagnostic for each violated
invariant. Forgetting a precondition surfaces during model
registration rather than during post hoc numerical validation.

Operationally: `Cargo.toml` declares the forces a paper uses,
`Cargo.lock` pins them bit-precisely, and a follow-up paper extending
the model adds one line. This is reproducibility at the
force-composition level, distinct from script-level reproducibility —
the latter captures the configuration but not the physics
implementation.

The perturbation itself becomes an independently validated, citable,
and reusable scientific unit.

> APSIS does not attempt to replace mature integrators or optimize
> numerical performance. Its contribution is orthogonal: defining how
> physical models are structured, published, and composed.

## Kernel invariants

The APSIS core guarantees, independently of any registered perturbation:

- **deterministic integration** given identical initial conditions,
  integrator, and seed;
- **bitwise-consistent Newtonian force evaluation** —
  `compute(bodies)` returns the same accelerations to f64 ULP across
  calls with identical state;
- **additive-only perturbation composition** — a registered
  perturbation accumulates into a scratch buffer; it cannot read or
  mutate the base force evaluation.

These are the invariants `KernelRequirements` declarations are
matched against (§ Statement of need); they hold across every
integrator the workspace ships for the entire lifetime of a
`System`.

The guarantees are formalised as executable specification in
[`apsis::contract`](crates/apsis/src/contract.rs) — twelve CI tests
covering kernel invariants, composition rules, and the failure
model, co-located with the prose statement of each guarantee. See
the §Design and validation section of [`paper.md`](paper.md) for
the formal treatment.

## Quickstart

### Python

APSIS is a runtime for composing physics distributed as crates.
Internal forces are submodules of the apsis distribution
(`apsis.gr`, `apsis.radiation`, `apsis.central`, …); external
forces ship as `apsis-plugin-X` packages with the same registration
contract. The simulation script is a *composition* of physics, not
a *configuration* of a monolith.

`pip install apsis` will work after the first PyPI release. Today,
build from source via [`maturin`](https://github.com/PyO3/maturin):

```bash
git clone https://github.com/GabrielEstefanski/apsis && cd apsis
pip install maturin
maturin develop --release
```

Then:

```python
import apsis
from apsis.gr import PostNewtonian1PN

sun = apsis.Body.star(mass=1.0)
mercury = apsis.Body.rocky(
    mass=1.66e-7, position=(0.387, 0.0), velocity=(0.0, 1.61),
)

sys = apsis.System(
    bodies=[sun, mercury], units=apsis.units.SOLAR,
    integrator="ias15", dt=1e-3,
)
sys.add_hamiltonian_perturbation(
    PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
)
sys.integrate_for(100.0)

print(sys)
```

Adding an internal force is one Cargo crate + one feature flag in
`apsis-python`. Adding an external force is `pip install apsis-plugin-yourforce`
and one extra import. Reproducing a paper's physical model is reading
its `requirements.txt` (Python) or `Cargo.toml` (Rust). The runtime
never changes; the composition does.

### Rust

The same Mercury 1PN scenario as a runnable example, with the GR
comparison printed inline (Rust 1.85+):

```bash
git clone https://github.com/gabrielbragaestefanski/apsis
cd apsis
cargo run --release --example mercury_perihelion -p apsis-1pn
```

Expected output (abridged):

```text
Mercury + Sun + 1PN @ IAS15
  T_mercury      = 1.513251 sim units
  integrating    = 500 orbits  →  t = 756.63
  ...
── GR comparison over 500 orbits ──
  predicted Δω      = +2.509976e-04 rad  (+51.7720 arcsec)
  measured  Δω      = +2.509906e-04 rad  (+51.7705 arcsec)
  relative error    = -2.802e-05
  rate              = 42.991 arcsec/century  (GR expects 43)
```

The same number is asserted in CI, gate-style:

```bash
cargo test --release -p apsis-1pn --tests -- --ignored
```

## A researcher-first API

A script to integrate a preset system with explicit integrator choice reads
in the terms a scientist uses to think about the simulation:

```rust
use apsis::core::system::System;
use apsis::templates::TemplateKind;

fn main() {
    let mut sys = System::from_template(TemplateKind::SolarSystem)
        .with_dt(1e-3);

    sys.integrate_for(100.0);

    println!("dE/E = {:.3e}", sys.energy_delta());
    println!("dLz  = {:.3e}", sys.lz_delta());
}
```

Bodies are built the same way:

```rust
use apsis::domain::body::Body;

let sun     = Body::star(1.0);
let mercury = Body::rocky(3e-6)
    .at(0.307, 0.0)
    .with_velocity(0.0, 1.98);
```

See [`crates/apsis/examples/`](crates/apsis/examples/)
and [`crates/apsis-1pn/examples/`](crates/apsis-1pn/examples/)
for ten runnable examples covering the Kepler 2-body problem, the solar
system integrated long, the three-body figure-eight, the Pythagorean
close-encounter problem, the Mercury perihelion test, preset
enumeration, scaling benchmarks, and the apsis side of each
cross-implementation parity scenario.

## Architecture: federation, library-first

The Rust workspace is a UI-free physics core plus a federation of
independently citable force crates. Python users see a single
`apsis` distribution that bundles every internal operator behind
one import. The core does not depend on bindings or downstream
consumers; CI enforces the separation.

| crate | role | dependencies |
|---|---|---|
| [`apsis`](crates/apsis/) | The library. Physics, integrators, public extension API. | Pure library; no UI or rendering dependencies. |
| [`apsis-1pn`](crates/apsis-1pn/) | First downstream force crate: 1PN Schwarzschild correction. Reference implementation of the federation contract. | **Only** `apsis`. |
| [`apsis-radiation`](crates/apsis-radiation/) | Radiation pressure + Poynting–Robertson drag (Burns 1979). | **Only** `apsis`. |
| [`apsis-central`](crates/apsis-central/) | Central-potential perturbations (observable-inversion exemplar, Tamayo 2019). | **Only** `apsis`. |
| [`apsis-py-core`](crates/apsis-py-core/) | Capsule transport + extractors (rlib) — used by the apsis Python distribution and any external `apsis-plugin-X` cdylib. | `apsis`, `pyo3`. |
| [`apsis-python`](crates/apsis-python/) | PyO3 cdylib backing the `apsis` Python distribution. Bundles every internal operator behind feature flags. | `apsis`, `apsis-py-core`, operator crates. |

Python source lives at the repository root in [`apsis/`](apsis/);
maturin builds the cdylib in `crates/apsis-python` via the root
`pyproject.toml`. Adding an internal force is adding a Cargo crate
and a feature flag; adding an external force is publishing
`apsis-plugin-X` against the public extension API.

The interactive visualisation shell lives in a separate repository
at [`GabrielEstefanski/apsis-app`](https://github.com/GabrielEstefanski/apsis-app);
it is a downstream consumer of the public `apsis` API and is not
part of the library's validated surface or release cadence.

## Fine-physics guardrail

A Plummer-softened kernel with `ε ≈ 0.02 AU` (a typical cluster-scale
ε for a solar-mass body) introduces a numerical apsidal precession
**≈ 5 × 10⁴ larger** than the 43 arcsec/century GR effect for Mercury,
with the *opposite* sign — the leading-order closed form
$\dot\varpi_\text{Plummer} = -3 n \varepsilon^2 / [2 a^2 (1 - e^2)^2]$
sets the scale; the exact value from the full-potential apsidal-angle
quadrature is $\approx -2.29 \times 10^6$ arcsec/century, which apsis
reproduces to $0.04\,\%$ (derivation in [`paper.md`](paper.md) §3.2).
It is invisible at the integrator level — energy still conserves to
machine precision.

This class of error is otherwise difficult to detect: the standard
conservation invariants (energy, angular momentum) remain satisfied bit-for-bit
while the precession measurement diverges from the physics being modelled.
The only upstream signal is a quantitative comparison against an analytic
prediction — which is exactly the step a researcher is likely to skip when
the simulator *looks* correct under every usual check.

The library surfaces the trap at the type level. The default kernel
is `NewtonKernel::exact()` (ε = 0); fine-physics scenarios stay safe
without opting in. Perturbations whose derivation depends on the
bit-exact `1/r` kernel (GR, tidal dissipation) or on a smooth
Hamiltonian (any symplectic-splitting derivation) override

```rust
fn kernel_requirements(&self) -> KernelRequirements {
    KernelRequirements::exact_and_smooth()
}
```

on the `Operator` trait. Registering such an operator on top of a
softened kernel — opted into via `System::with_kernel(Arc::new(NewtonKernel::new(ε > 0)))`,
typically for cluster work — emits a `warn_diag!` diagnostic at
registration time naming the violated invariant.

## What this library is NOT

For research where the simulator is the primary tool — solar-system
integration, hybrid close-encounter regimes, collisionless large-N,
stellar evolution, hydrodynamics — use REBOUND, MERCURIUS, NBODY6/7,
or GADGET. APSIS is the tool when the perturbation itself is the
contribution.

Out of current scope: stellar evolution, hydrodynamics, collisionless
large-N. Published validation runs at $N \le 10^3$ bodies (horizons
up to $10^4$ orbits for the long-horizon parity gate); larger-$N$
behaviour is the subject of the v0.2 scaling notebook.

## Validation

What is verified in CI:

- **Unit + integration test suite** across the workspace covering
  energy conservation on canonical scenarios (Kepler circular,
  Pythagorean three-body, figure-eight), IAS15 determinism on seeded
  close encounters, and conservation-contract assertions on the
  public API.
- **Per-operator validation gates**: 1PN (Mercury precession),
  radiation (Burns 1979 β-table dust decay), central-force
  (Tamayo 2019 round-trip).
- **Release-mode Mercury gate**: `cargo test --release -p apsis-1pn
  -- --ignored` asserts Mercury's precession within 100 ppm of the
  analytic 1PN GR prediction over 500 orbits, with the achieved
  figure 28 ppm bit-identical across Windows and Linux on x86_64.
- **Cross-implementation parity portfolio**: against REBOUND's IAS15
  on four canonical scenarios spanning periodic 2-body, periodic
  3-body, chaotic 3-body, and sign-flipped 2-body regimes. All gated
  invariant metrics (energy, angular momentum, orbital elements where
  defined, linear momentum and centre-of-mass for the three-body
  scenarios) agree at 1 ULP in regime. Each scenario carries an
  *a priori* protocol notebook (initial conditions, integrator
  settings, and tolerances declared before the run) and a
  self-contained Python harness:
  - **Kepler-prograde ($e = 0.5$, 100 orbits):** seven gated metrics
    at 1–3 ULP; informational $\lvert\Delta r\rvert$ at
    $1.57 \times 10^{-9}$ (peak orbit 81).
    Notebook: [`paper/notebooks/2026-04-25-rebound-parity-kepler.md`](paper/notebooks/2026-04-25-rebound-parity-kepler.md).
  - **Figure-8 (Chenciner–Montgomery, 10 periods gated + 50 periods
    informational):** twelve gated metrics organised in three
    evidentiary tiers (hard physical invariants, construction-level
    sanity, geometric coherence) at 1 ULP.
    Notebook: [`paper/notebooks/2026-04-26-rebound-parity-figure8.md`](paper/notebooks/2026-04-26-rebound-parity-figure8.md).
  - **Pythagorean (Burrau 1913, 70 canonical t.u.):** structural
    invariants ($\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$)
    at f64 round-off floor on both sides; energy bound exceeded
    symmetrically by both implementations in the chaotic
    close-encounter regime — a documented regime mismatch with the
    smooth-flow bound's derivation, not a parity defect. 98 % event
    alignment (44/45 close-encounter peaks matched within
    $3 \times 10^{-2}$ t.u.).
    Notebook: [`paper/notebooks/2026-04-30-rebound-parity-pythagorean.md`](paper/notebooks/2026-04-30-rebound-parity-pythagorean.md).
  - **Kepler-retrograde ($L_z < 0$; long-horizon $10^4$ orbits +
    100-orbit checkpoint):** ten gated metrics × two horizons; all
    twenty pass at 1–10 ULP. Closes the sign-convention coverage
    gap and the long-horizon stability gate identified during the
    GR-readiness review. Brouwer-law growth confirmed
    ($\sim 8\times$ across $100\times$ horizon, slightly below
    $\sqrt{N}$, consistent with IAS15's near-symplectic structure).
    Notebook: [`paper/notebooks/2026-05-01-rebound-parity-retrograde.md`](paper/notebooks/2026-05-01-rebound-parity-retrograde.md).
- **`recommended_dt` heuristic** (utility for fixed-step integrators):
  13 scenarios × 3 integrators × 100 substeps; 18 / 18 gated cells
  pass, 21 informational cells (WH + four out-of-regime scenarios:
  high-eccentricity Kepler, periodic and chaotic 3-body). Per-cell
  utilization ($u = \text{peak}/\text{bound}$) emitted as a
  regression canary.
  Note: [`docs/experiments/2026-05-01-recommended-dt-heuristic.md`](docs/experiments/2026-05-01-recommended-dt-heuristic.md).
## Further reading

The repository carries the full methodological record a software paper
normally cites only in passing — reviewers and users who want the
decisions and the failed experiments behind a number can follow the
trail directly:

- [`docs/overview.md`](docs/overview.md) — the project's scope, physical
  model, architecture, and known limitations in prose.
- [`docs/integrator.md`](docs/integrator.md) — per-integrator contract:
  execution profile, determinism requirement, selection rubric.
- [`docs/softening.md`](docs/softening.md) — Plummer softening derivation,
  per-body scaling rule, and the regime in which it is trustworthy.
- [`docs/adr/`](docs/adr/) — architectural decision records.
- [`docs/experiments/`](docs/experiments/) — lab notebooks pairing
  *a priori* protocols with executed runs and post-mortem analysis.
- [`paper/notebooks/`](paper/notebooks/) — reviewer-facing parity
  and validation protocols cited from the paper.
- [`validation/`](validation/) — runnable cross-implementation
  harnesses (one directory per reference tool).

## License

Licensed under the [Apache License, Version 2.0](LICENSE).

## Citing this work

*(A Zenodo DOI and JOSS reference will appear here after the first tagged
release.)*

## References

- Rein, H., & Spiegel, D. S. (2015). *IAS15: a fast, adaptive,
  high-order integrator for gravitational dynamics, accurate to machine
  precision over a billion orbits.* **MNRAS, 446(2), 1424–1437.**
- Everhart, E. (1985). *An efficient integrator that uses Gauss–Radau
  spacings.* In A. Carusi & G. B. Valsecchi (Eds.), *Dynamics of Comets:
  Their Origin and Evolution*, **Astrophysics and Space Science Library
  115**, 185–202. Springer.
- Rein, H., & Liu, S.-F. (2012). *REBOUND: an open-source multi-purpose
  N-body code for collisional dynamics.* **A&A, 537, A128.**
- Tamayo, D., Rein, H., Shi, P., & Hernandez, D. M. (2020). *REBOUNDx: a
  library for adding conservative and dissipative forces to otherwise
  symplectic N-body integrations.* **MNRAS, 491(2), 2885–2901.**
- Chenciner, A., & Montgomery, R. (2000). *A remarkable periodic
  solution of the three-body problem in the case of equal masses.*
  **Annals of Mathematics, 152(3), 881–901.**
- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics.*
  Cambridge University Press.
- Einstein, A. (1915). *Explanation of the perihelion motion of Mercury
  from the general theory of relativity.* **Preussische Akademie der
  Wissenschaften, Sitzungsberichte, 831–839.**
