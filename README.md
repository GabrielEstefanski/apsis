# gravity-sim

A Rust N-body gravitational simulation library with an adaptive IAS15-style
integrator (in the sense of Rein & Spiegel, 2015) and a compiler-enforced
public extension API. Validated by an out-of-tree plugin crate reproducing
Mercury's perihelion precession to **4.4 parts per million** of the
General-Relativistic prediction.

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
the specific combination — Rust, a validated IAS15 implementation, and a
plugin contract enforced by compilation rather than convention — is not
currently available elsewhere. Concretely, the claim means:

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
  [`gravity-sim-1pn`](crates/gravity-sim-1pn/) crate reproduces the textbook
  43 arcsec/century result at 4.4 ppm relative error, on an isolated build
  that never touches the core's sources.

This is a **software-methods** contribution, not a new-physics contribution.

## Quickstart

Prerequisites: Rust 1.85+ (`rustup install stable`).

```bash
git clone https://github.com/gabrielbragaestefanski/gravity-sim
cd gravity-sim
cargo run --release --example mercury_perihelion -p gravity-sim-1pn
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
cargo test --release -p gravity-sim-1pn --tests -- --ignored
```

## A researcher-first API

A script to integrate a preset system with explicit integrator choice reads
in the terms a scientist uses to think about the simulation:

```rust
use gravity_sim_core::core::system::System;
use gravity_sim_core::physics::integrator::IntegratorKind;
use gravity_sim_core::templates::TemplateKind;

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
use gravity_sim_core::domain::body::Body;

let sun     = Body::star(1.0);
let mercury = Body::rocky(3e-6)
    .at(0.307, 0.0)
    .with_velocity(0.0, 1.98)
    .unsoftened();                    // see "Fine-physics" below
```

See [`crates/gravity-sim-core/examples/`](crates/gravity-sim-core/examples/)
and [`crates/gravity-sim-1pn/examples/`](crates/gravity-sim-1pn/examples/)
for seven runnable examples covering Kepler 2-body, the solar system
integrated long, the three-body figure-eight, the Pythagorean problem, the
Mercury perihelion test, and preset enumeration.

## Architecture: library-first, app-as-side

The workspace is three crates deliberately split by role:

| crate | role | dependencies |
|---|---|---|
| [`gravity-sim-core`](crates/gravity-sim-core/) | The library. Physics, integrators, public API. | Zero UI: `cargo tree -p gravity-sim-core` resolves no `egui`/`wgpu`/`eframe`. |
| [`gravity-sim-1pn`](crates/gravity-sim-1pn/) | The out-of-tree plugin demonstration. 1PN correction via `PerturbationForce`. | **Only** `gravity-sim-core`. Reviewed as the paper's Phase-3 gate. |
| [`gravity-sim-app`](crates/gravity-sim-app/) | Optional interactive egui/wgpu shell. **Not** part of the library's validated surface. | `egui`, `wgpu`, `eframe`. |

The direction is `gravity-sim-app` → `gravity-sim-core`; the core does not
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
system with softened bodies emits a `warn_diag!` diagnostic at plugin
registration, with per-body softening statistics. Dismiss by

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

- **200 unit tests** in the core covering energy conservation on canonical
  scenarios (Kepler circular, Pythagorean three-body, figure-eight), IAS15
  determinism on seeded close encounters, and conservation-contract
  assertions on the public API.
- **11 tests in the 1PN plugin**: 7 unit (sign convention, magnitude,
  additivity, speed-of-light limit), 2 release-mode integration
  (`mercury_precession_matches_gr_within_one_percent`,
  `baseline_newtonian_kepler_is_closed`), and 2 debug-mode contract
  (softened-system-warns, unsoftened-system-silent).
- **Release-mode Phase-3 gate**: `cargo test --release -p gravity-sim-1pn
  -- --ignored` asserts Mercury's precession within 1 % of GR over 300
  orbits. 4.4 ppm is the achieved figure.
- **Workspace isolation**: `cargo build -p gravity-sim-core` resolves no
  UI dependency.

## Further reading

The repository carries the full methodological record a software paper
normally cites only in passing — reviewers and users who want the
decisions and the failed experiments behind a number can follow the
trail directly:

- [`docs/overview.md`](docs/overview.md), [`docs/integrator.md`](docs/integrator.md),
  [`docs/forces.md`](docs/forces.md), [`docs/softening.md`](docs/softening.md),
  [`docs/stability.md`](docs/stability.md) — domain-level documentation of
  the physics stack and its regime of applicability.
- [`docs/adr/`](docs/adr/) — architectural decision records.
  `001-wall-time-budget.md` on the interactive timestep model,
  `002-sim-rate-target.md` on frame-pacing, and
  `003-integrator-execution-profile.md` on why the default is
  Yoshida-4 rather than IAS15 for render-loop contexts.
- [`docs/experiments/`](docs/experiments/) — lab-notebook entries for
  reproducible experiments run during development, including the
  IAS15 phase-profile breakdown, a null-result on the Picard noise
  floor, and the solar-system stutter diagnosis that motivated the
  versioned baseline harness.
- [`docs/references.md`](docs/references.md) — consolidated bibliography.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).

## Citing this work

*(A Zenodo DOI and JOSS reference will appear here after the first tagged
release.)*

## References

- Rein, H., & Spiegel, D. S. (2015). *IAS15: a fast, adaptive,
  high-order integrator for gravitational dynamics, accurate to machine
  precision over a billion orbits.* **MNRAS, 446(2), 1424–1437.**
- Rein, H., & Liu, S.-F. (2012). *REBOUND: an open-source multi-purpose
  N-body code for collisional dynamics.* **A&A, 537, A128.**
- Tamayo, D., Rein, H., Shi, P., & Hernandez, D. M. (2020). *REBOUNDx: a
  library for adding conservative and dissipative forces to otherwise
  symplectic N-body integrations.* **MNRAS, 491(2), 2885–2901.**
- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics.*
  Cambridge University Press.
- Einstein, A. (1915). *Explanation of the perihelion motion of Mercury
  from the general theory of relativity.* **Preussische Akademie der
  Wissenschaften, Sitzungsberichte, 831–839.**
