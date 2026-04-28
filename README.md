# APSIS

*A Federated Model for Composable N-Body Force Artifacts*

APSIS treats N-body force perturbations as first-class scientific
artifacts. The model is *federated*: each force is an independent
Rust crate — developed, versioned, and cited separately, without
central integration into a monolithic codebase. Every force ships
with a Python binding through `apsis-py-core`'s cross-extension
transport, preserving the same contract across the FFI without
reimplementing physics. A simulation's physical model is therefore
not embedded in code, but in its dependency graph: `Cargo.toml`
declares the forces, `Cargo.lock` pins them bit-precisely. The
simulator is infrastructure for composing those artifacts.

The core integrator is IAS15 (Rein & Spiegel, 2015), audited against
the algorithmic specification §2–3 and validated against REBOUND's
IAS15 on Kepler ($e = 0.5$, 100 orbits) and the Chenciner–Montgomery
figure-8 (10 periods); all gated invariant metrics agree at **1 ULP**
of f64 machine epsilon. The first downstream artifact,
[`apsis-1pn`](crates/apsis-1pn/), reproduces Mercury's perihelion
precession to **4.4 ppm** of the GR prediction over 500 orbits, gated
in CI.

> **Status.** Pre-release (`v0.1.0` alpha). The integrator and contract
> machinery are 2D; the 3D port is the next breaking-change milestone
> (v0.2). Public API stabilised but not yet tagged; citation DOI pending
> first Zenodo release.

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

APSIS replaces that publication path. A force is a Cargo crate
declaring its physical preconditions on the gravitational kernel via
the `KernelRequirements` type — the 1PN crate declares
`exact_and_smooth()`; future crates declare a different combination
of exactness and continuity invariants depending on the physics. The
core matches the declared requirements against the active kernel at
`System::add_perturbation` and emits a structured diagnostic for each
violated invariant. Forgetting a precondition surfaces as a
registration warning, not as a wrong number in a paper.

Operationally: `Cargo.toml` declares the forces a paper uses,
`Cargo.lock` pins them bit-precisely, and a follow-up paper extending
the model adds one line. This is reproducibility at the
force-composition level, distinct from script-level reproducibility —
the latter captures the configuration but not the physics
implementation.

> APSIS does not attempt to replace mature integrators or optimize
> numerical performance. Its contribution is orthogonal: defining how
> physical models are structured, published, and composed.

The IAS15 integrator and the Mercury 4.4 ppm result are evidence that
the contract machinery operates against numerics meeting the field's
precision floor — not the headline claim. Use REBOUND/REBOUNDx when
the simulator is the primary tool; use APSIS when the perturbation is.

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
integrator (Velocity Verlet, Yoshida-4, Wisdom-Holman, IAS15) for
the entire lifetime of a `System`. Phase 1 of the v0.2 milestone
turns these guarantees into typed contracts with CI tests asserting
each one directly.

## Quickstart

### Python

APSIS is not a simulation library — it is a runtime for composing
physics distributed as crates. Each force lives in its own
pip-installable package and is registered with a system at runtime;
the simulation script is a *composition* of physics, not a
*configuration* of a monolith. The example below composes two crates:
the `apsis` runtime and `apsis-1pn`, a force crate implementing the
1PN relativistic correction — the effect responsible for Mercury's
perihelion precession.

`pip install apsis apsis-1pn` will work from v0.2.0. Today, build
from source via [`maturin`](https://github.com/PyO3/maturin):

```bash
git clone https://github.com/gabrielbragaestefanski/apsis && cd apsis
pip install maturin
maturin develop --release --manifest-path crates/apsis-py/Cargo.toml
maturin develop --release --manifest-path crates/apsis-1pn-py/Cargo.toml
```

Then:

```python
import apsis
import apsis_1pn  # the 1PN force, distributed as an independent package

sun = apsis.Body.star(mass=1.0).unsoftened()
mercury = apsis.Body.rocky(
    mass=1.66e-7, position=(0.387, 0.0), velocity=(0.0, 1.61),
).unsoftened()

sys = apsis.System(
    bodies=[sun, mercury], units=apsis.units.SOLAR,
    integrator="ias15", dt=1e-3, exact_gravity=True,
)
sys.add_perturbation(apsis_1pn.PostNewtonian1PN.solar_units())
sys.integrate_for(100.0)

print(sys)
```

Adding a force is `pip install apsis-yourforce` and one extra import.
Reproducing a paper's physical model is reading its
`requirements.txt`. The runtime never changes; the composition does.

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
  predicted Δω      = +2.509427e-04 rad  (+51.7606 arcsec)
  measured  Δω      = +2.509438e-04 rad  (+51.7609 arcsec)
  relative error    = +4.449e-06
  rate              = 42.983 arcsec/century  (GR expects 43)
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
use apsis::physics::integrator::IntegratorKind;
use apsis::templates::TemplateKind;

fn main() {
    let mut sys = System::from_template(TemplateKind::SolarSystem)
        .with_integrator(IntegratorKind::Ias15)
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
    .with_velocity(0.0, 1.98)
    .unsoftened();                    // see "Fine-physics" below
```

See [`crates/apsis/examples/`](crates/apsis/examples/)
and [`crates/apsis-1pn/examples/`](crates/apsis-1pn/examples/)
for ten runnable examples covering the Kepler 2-body problem, the solar
system integrated long, the three-body figure-eight, the Pythagorean
close-encounter problem, the Mercury perihelion test, preset
enumeration, scaling benchmarks, and the apsis side of each
cross-implementation parity scenario.

## Architecture: federation, library-first

The workspace is six crates split by role: a UI-free physics core, a
Python façade, and a federation of independently citable extension
points. The core does not know the app or the bindings exist; CI
enforces the separation.

| crate | role | dependencies |
|---|---|---|
| [`apsis`](crates/apsis/) | The library. Physics, integrators, public extension API. | Zero UI: `cargo tree -p apsis` resolves no `egui`/`wgpu`/`eframe`. |
| [`apsis-py`](crates/apsis-py/) | Python binding (PyO3, abi3-py39). Façade-only. | `apsis`, `pyo3`, `numpy`. |
| [`apsis-py-core`](crates/apsis-py-core/) | Cross-extension transport (rlib): `Box<dyn PerturbationForce>` ↔ `PyCapsule`. | `apsis`, `pyo3`. |
| [`apsis-1pn`](crates/apsis-1pn/) | First downstream force crate: 1PN Schwarzschild correction. Reference implementation of the federation contract. | **Only** `apsis`. |
| [`apsis-1pn-py`](crates/apsis-1pn-py/) | Python binding for `apsis-1pn`. Reference implementation of the contract at the Rust/Python boundary. | `apsis-1pn`, `apsis-py-core`. |
| [`apsis-app`](crates/apsis-app/) | Optional interactive egui/wgpu shell. Not part of the library's validated surface. | `apsis`, `egui`, `wgpu`, `eframe`. |

Direction: every binding and every force crate depends on `apsis`
through the public extension API only — never `pub(crate)`, never
core internals. Adding a force is adding a crate.

## Fine-physics guardrail

The default material-scaled Plummer softening (`ε ≈ 0.02 AU` on the Sun)
introduces a numerical apsidal precession that is **≈ 2 × 10³ larger** than
the 43 arcsec/century GR effect for Mercury. It is invisible at the
integrator level — energy still conserves to machine precision.

This class of error is otherwise difficult to detect: the standard
conservation invariants (energy, angular momentum) remain satisfied bit-for-bit
while the precession measurement diverges from the physics being modelled.
The only upstream signal is a quantitative comparison against an analytic
prediction — which is exactly the step a researcher is likely to skip when
the simulator *looks* correct under every usual check.

The library surfaces the trap at the type level. Perturbations whose
derivation depends on the bit-exact `1/r` kernel (GR, tidal dissipation)
or on a smooth Hamiltonian (any symplectic-splitting derivation)
override

```rust
fn kernel_requirements(&self) -> KernelRequirements {
    KernelRequirements::exact_and_smooth()
}
```

on the `PerturbationForce` trait. Registering such a perturbation into
a system with softened bodies emits a `warn_diag!` diagnostic at
registration time, with per-body softening statistics naming the
violated invariant. Dismiss by

```rust
let sun = Body::star(1.0).unsoftened();             // per body
let sys = System::from_template(..).with_exact_gravity();   // whole system
```

Both helpers exist because a reviewer reading fluent-builder code should read
the *intent*, not the field assignment.

## What this library is NOT

APSIS occupies a different *category* of system from the established
N-body codes — not a feature-thin alternative to them. The
[orthogonality declaration](#statement-of-need) makes this concrete:
APSIS does not replace mature integrators or chase numerical performance.

For research where the simulator is the primary tool — solar-system
integration with extra forces, hybrid close-encounter regimes,
collisionless large-N, stellar evolution, hydrodynamics — use REBOUND,
MERCURIUS, NBODY6/7, or GADGET. APSIS is the tool when **the perturbation
is the scientific contribution** and the question is how to publish,
version, and compose it. APSIS trades ecosystem maturity for
composability and publication clarity; the choice between the two is
a property of the research question, not of the codebase.

Out of scope at v0.1: 3D integration (planned v0.2), symplectic
compositions beyond Yoshida-4, MERCURIUS-style close-encounter switching,
stellar evolution, hydrodynamics, collisionless large-N.

## Validation

What is verified in CI:

- **241 unit tests** in the core covering energy conservation on canonical
  scenarios (Kepler circular, Pythagorean three-body, figure-eight),
  IAS15 determinism on seeded close encounters, conservation-contract
  assertions on the public API, and direct unit tests pinning the IAS15
  warmstart against the analytical Pascal-triangle transformation
  derived in Everhart (1985).
- **13 tests in the 1PN plugin**: 7 unit (sign convention, magnitude,
  additivity, speed-of-light limit), 4 in the Mercury-precession gate,
  and 2 debug-mode contract (softened-system-warns, unsoftened-system-silent).
- **Release-mode Phase-3 gate**: `cargo test --release -p apsis-1pn
  -- --ignored` asserts Mercury's precession within 1 % of GR over 300
  orbits. 4.4 ppm is the achieved figure.
- **Cross-implementation parity portfolio**: against REBOUND's IAS15
  on two canonical scenarios, with all gated invariant metrics
  (energy, angular momentum, orbital elements where defined, linear
  momentum and centre-of-mass for the figure-8) agreeing at 1 ULP.
  Each scenario carries an *a priori* protocol notebook
  (initial conditions, integrator settings, and tolerances declared
  before the run) and a self-contained Python harness:
  - **Kepler ($e = 0.5$, 100 orbits):** seven gated metrics at 1–3 ULP;
    informational $\lvert\Delta r\rvert$ at $2.18 \times 10^{-12}$.
    Notebook: [`docs/experiments/2026-04-25-rebound-parity-kepler.md`](docs/experiments/2026-04-25-rebound-parity-kepler.md).
  - **Figure-8 (Chenciner–Montgomery, 10 periods):** twelve gated metrics
    organised in three evidentiary tiers (hard physical invariants,
    construction-level sanity, geometric coherence) at 1 ULP;
    informational $\lvert\Delta r\rvert$ at $9.44 \times 10^{-13}$.
    Notebook: [`docs/experiments/2026-04-26-rebound-parity-figure8.md`](docs/experiments/2026-04-26-rebound-parity-figure8.md).
- **Workspace isolation**: `cargo build -p apsis` resolves no
  UI dependency.

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
  `001-wall-time-budget.md` on the interactive timestep model,
  `002-sim-rate-target.md` on frame-pacing, and
  `003-integrator-execution-profile.md` on why the default is
  Yoshida-4 rather than IAS15 for render-loop contexts.
- [`docs/experiments/`](docs/experiments/) — lab-notebook entries for
  reproducible experiments run during development. Each entry pairs an
  *a priori* protocol with the executed run and a post-mortem analysis;
  the directory currently records:
  - the IAS15 phase-profile breakdown,
  - a null result on the Picard noise floor,
  - the operational-domain benchmark suite that motivated the
    versioned baseline harness,
  - the Kepler and figure-8 cross-implementation parity protocols
    (the notebooks the validation portfolio is anchored to), and
  - the IAS15 controller architecture audit
    ([`2026-04-26-ias15-warmstart-bug.md`](docs/experiments/2026-04-26-ias15-warmstart-bug.md)),
    which documents the three controller divergences from
    Rein & Spiegel (2015) that the figure-8 parity scenario
    surfaced and the line-by-line resolution.
- [`validation/`](validation/) — runnable cross-implementation harnesses
  one directory per reference tool (currently REBOUND), each with its
  own Python `run.py` orchestrator and a comparator that emits a
  structured JSON report alongside the CSV outputs.

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
