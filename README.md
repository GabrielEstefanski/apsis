# APSIS

*Verified Extension Contracts for N-Body Simulation in Rust*

A Rust N-body gravitational simulation library with an adaptive IAS15
integrator (Rein & Spiegel, 2015), audited at the controller level
against the algorithmic specification, and a compiler-enforced public
extension API. Validated against two independent reference signals: an
out-of-tree companion crate reproducing Mercury's perihelion precession
to **4.4 parts per million** of the General-Relativistic prediction, and
a cross-implementation parity portfolio — Kepler ($e = 0.5$, 100 orbits)
and the Chenciner–Montgomery figure-8 (10 periods) — against REBOUND's
IAS15 implementation, with all gated invariant metrics agreeing at
**1 ULP** of f64 machine epsilon.

> **Scope.** The solver is currently 2D. 3D is a planned, deliberately
> breaking change — the current surface is frozen at 2D so the API-contract
> machinery can be exercised end-to-end against a real physical result
> before the coordinate dimension changes.
>
> *Status: pre-release (`v0.1.0` alpha). Public API is stabilised but not yet
> tagged; citation DOI pending first Zenodo release.*

---

## Statement of need

N-body gravitational simulation in solar-system-scale physics is dominated by
a small number of mature C/Fortran codes — REBOUND (Rein & Spiegel, 2012),
MERCURY (Chambers, 1999), NBODY6/7 (Aarseth, 2003) — each with decades of
community validation. This library does not seek to replace them.

It fills a narrower niche: **a Rust-native N-body library providing an
adaptive IAS15-style integrator behind a public API whose invariants are
promoted to type-level, CI-enforced contracts.** To the authors' knowledge,
the specific combination — Rust, a validated IAS15 implementation, and
extension contracts enforced by compilation rather than convention — is
not currently available elsewhere. Concretely, the claim means:

- Physical preconditions (exact `1/r` gravity, determinism seed, softening
  contracts) are declared in code at the type of each extension point, not
  left to prose in a methods section. Forgetting them surfaces as a build-time
  warning, not a silently-wrong result at publication time.
- Third-party physics extensions compose against the core through a
  `PerturbationForce` trait in an out-of-tree crate, with nothing in the
  core reaching for `pub(crate)` or internals. The contract is
  **compilation**, not convention.
- Validation uses the canonical test physicists have reached for a century:
  the perihelion precession of Mercury. The out-of-tree
  [`apsis-1pn`](crates/apsis-1pn/) crate reproduces the textbook
  43 arcsec/century result at 4.4 ppm relative error, on an isolated build
  that never touches the core's sources.

This is a **software-methods** contribution, not a new-physics contribution.

## Quickstart

Prerequisites: Rust 1.85+ (`rustup install stable`).

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

## Architecture: library-first, app-as-side

The workspace is three crates deliberately split by role:

| crate | role | dependencies |
|---|---|---|
| [`apsis`](crates/apsis/) | The library. Physics, integrators, public API. | Zero UI: `cargo tree -p apsis` resolves no `egui`/`wgpu`/`eframe`. |
| [`apsis-1pn`](crates/apsis-1pn/) | Out-of-tree companion crate: 1PN correction via `PerturbationForce`. | **Only** `apsis`. Reviewed as the paper's Phase-3 gate. |
| [`apsis-app`](crates/apsis-app/) | Optional interactive egui/wgpu shell. **Not** part of the library's validated surface. | `egui`, `wgpu`, `eframe`. |

The direction is `apsis-app` → `apsis`; the core does not
know the app exists. CI enforces the separation.

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

The library surfaces the trap at the type level. Perturbations whose signal
measures a deviation from `1/r` (GR, J2 oblateness, tidal dissipation)
override

```rust
fn requires_exact_gravity(&self) -> bool { true }
```

on the `PerturbationForce` trait. Registering such a perturbation into a
system with softened bodies emits a `warn_diag!` diagnostic at registration
time, with per-body softening statistics. Dismiss by

```rust
let sun = Body::star(1.0).unsoftened();             // per body
let sys = System::from_template(..).with_exact_gravity();   // whole system
```

Both helpers exist because a reviewer reading fluent-builder code should read
the *intent*, not the field assignment.

## What this library is NOT

Honest scope for reviewers. This library **does not** provide:

- Symplectic composition integrators beyond 4th order Yoshida (no SABA, no
  higher-order splittings).
- A hybrid close-encounter regime switcher (no MERCURIUS equivalent).
- Stellar evolution, hydrodynamics, or collisionless large-N (no GADGET
  equivalent).
- 3D integration (see the Scope note at the top — a planned breaking
  change, not a regression).
- Python bindings (possible via `pyo3` as future work; out of scope for v1).

For any of the above, use REBOUND, MERCURY, or NBODY6. This library's
positioning is narrow and deliberate.

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
