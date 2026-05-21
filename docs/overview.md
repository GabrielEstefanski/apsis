# System Overview

## 1. Purpose and scope

`apsis` is an open-source N-body gravitational simulator implemented
in Rust. The project's published contribution is the `apsis`
library; the interactive visualisation shell lives in a separate
repository ([`GabrielEstefanski/apsis-app`](https://github.com/GabrielEstefanski/apsis-app))
as a downstream consumer of the public API, not part of the
library's validated surface (see the workspace
[`README.md`](../README.md) for the paper-level positioning).

The solver is 3D: `Vec3` value type, per-body `z` / `vz`, 3D direct
N-body kernel, IAS15 dense interpolation in 3D, orbital elements
with inclination $i$ and node $\Omega$, and 3D energy /
angular-momentum / centre-of-mass observables.

---

## 2. Physical model and assumptions

The core physics follows classical Newtonian mechanics, with optional
perturbation operators (post-Newtonian corrections, radiation pressure,
central forces) supplied by independent operator crates published
against the public extension API.

**Modelled:**

- Universal gravitation (Newton's law of gravitation). The
  gravitational constant $G$ is derived from the active
  [`UnitSystem`](../crates/apsis/src/units.rs); the canonical
  (Hénon) unit system uses $G = 1$ by convention.
- Isolated system boundary: no external forces beyond those registered
  explicitly via the operator extension point — see
  [`crates/apsis/src/physics/integrator/operator.rs`](../crates/apsis/src/physics/integrator/operator.rs)
  for the three-trait split (`Operator`, `HamiltonianOperator`,
  `NonConservativeOperator`).
- Plummer softening as a property of the gravitational kernel
  (`Kernel::Newton { softening }`), not of individual bodies — see
  [`softening.md`](softening.md) for derivation and trustworthiness
  regime. The default kernel is `Kernel::Newton::exact()` ($\epsilon = 0$);
  softening is opt-in for cluster-scale work.
- Operator-published perturbations: post-Newtonian corrections
  ([`apsis-1pn`](../crates/apsis-1pn/), Schwarzschild test-particle 1PN),
  radiation pressure and Poynting–Robertson drag
  ([`apsis-radiation`](../crates/apsis-radiation/), Burns 1979),
  central-force perturbations
  ([`apsis-central`](../crates/apsis-central/), Tamayo 2019).

**Not modelled:**

- Collisions, mergers, or fragmentation (bodies are point masses even
  when the plugin layer measures their physical radius).
- Hydrodynamics, stellar structure, or collisionless large-N.
- Large-$N$ scaling: published validation runs at $N \le 10^3$
  bodies (horizons up to $10^4$ orbits for the long-horizon parity
  gate); larger-$N$ behaviour is the subject of the v0.2 scaling
  notebook.

---

## 3. Workspace architecture

The project is a Cargo workspace split by role:

| crate | role |
|---|---|
| [`apsis`](../crates/apsis/) | The library: domain types, physics, integrators, `apsis::contract` formalisation, public API. Pure library; no UI or rendering dependencies. |
| [`apsis-1pn`](../crates/apsis-1pn/) | First downstream force crate: 1PN Schwarzschild correction (Anderson 1975). Reference implementation of the federation contract. Depends only on `apsis` through the public API. |
| [`apsis-radiation`](../crates/apsis-radiation/) | Radiation pressure + Poynting–Robertson drag (Burns 1979). |
| [`apsis-central`](../crates/apsis-central/) | Central-potential perturbations (Tamayo 2019, observable-inversion exemplar). |
| [`apsis-py-core`](../crates/apsis-py-core/) | Capsule transport + extractors (rlib). Used by the `apsis` Python distribution and any external `apsis-plugin-X` cdylib. |
| [`apsis-python`](../crates/apsis-python/) | PyO3 cdylib backing the `apsis` Python distribution. Bundles every internal operator behind feature flags. |

The interactive visualisation shell lives in a separate
repository at [`GabrielEstefanski/apsis-app`](https://github.com/GabrielEstefanski/apsis-app);
it is a downstream consumer of the public `apsis` API.

The Python package source lives at the repository root in
[`apsis/`](../apsis/); maturin builds the cdylib via the root
`pyproject.toml`. Users `pip install apsis` and write
`from apsis.gr import PostNewtonian1PN` — every internal operator
ships under one import.

The dependency direction is monotone: every binding and every force
crate depends on `apsis` through the public extension API only —
never `pub(crate)`, never core internals. Adding a force is adding a
crate.

### 3.1 Core subsystems

Inside `apsis`, the code is organised by responsibility:

- [`domain/`](../crates/apsis/src/domain/) — `Body`,
  `NamedBody`, and the body-array layout used by the force kernels.
- [`physics/`](../crates/apsis/src/physics/) — force
  models (direct $O(N^2)$ + Barnes-Hut with a Newton kernel that
  optionally carries Plummer softening), the integrator stack
  (Velocity Verlet, Yoshida 4, Wisdom-Holman, WHFast, Mercurius,
  Implicit Midpoint, IAS15), orbital elements, energy and
  angular-momentum diagnostics, and the operator extension traits
  (`HamiltonianOperator`, `NonConservativeOperator`).
- [`core/`](../crates/apsis/src/core/) — the `System`
  orchestrator, adaptive dt/θ controllers, hook registry,
  structured diagnostic event bus, and metrics assembly.
- [`io/`](../crates/apsis/src/io/) — CSV recorder for headless runs.
- [`records/`](../crates/apsis/src/records/) — the apsis-record
  binary trajectory format: TOML provenance header, Snapshot /
  Diagnostic / Event / ResumeState / Trailer frames, BLAKE3-covered
  frame stream, byte-equal reproducibility gate.
- [`templates/`](../crates/apsis/src/templates/) — the
  `TemplateKind` enum and builders for built-in presets (Solar
  System, TRAPPIST-1, figure-eight, etc.), consumed uniformly by
  the app, the headless runner, and test scripts.

---

## 4. Numerical integration

Each integrator declares its per-step cost class via
[`IntegratorKind::is_adaptive`](../crates/apsis/src/physics/integrator/traits.rs):

| integrator | order | per-step cost | primary use |
|---|---|---|---|
| IAS15 | 15th (adaptive Gauss–Radau) | Adaptive | **Default**; paper-grade trajectory |
| Mercurius | hybrid (WH + IAS15) | Bounded outer | Planetary systems with close encounters |
| WHFast | 2nd (Keplerian, compensated) | Bounded | Long-horizon planetary integration |
| Wisdom–Holman | mixed-order | Bounded | Hierarchical analytical Kepler drift |
| Implicit Midpoint | 2nd (A-stable Gauss–Legendre) | Bounded | BH binaries / equal-mass / particle clouds |
| Yoshida 4 | 4th (symplectic) | Bounded | Higher-order symplectic, fixed cost |
| Velocity Verlet | 2nd (symplectic) | Bounded | Quick exploration, educational runs |

See [`integrator.md`](integrator.md) for the per-integrator contract
and the force-model determinism rule (IAS15 + direct summation is
the coherent default; Barnes-Hut is opt-in via
`set_exact_threshold(N)`). See [ADR-013](adr/013-default-integrator-ias15.md)
for the default-choice rationale.

---

## 5. Validation

Validation rests on two independent reference signals: an
analytic-physics test (Mercury's perihelion precession against General
Relativity) and a cross-implementation parity portfolio (REBOUND IAS15
on canonical scenarios). Together they bound the library's correctness
both physically and numerically.

### 5.1 Analytic physics — Mercury perihelion (Phase-3 gate)

The library's headline physical validation is the Mercury perihelion
precession test shipped by the `apsis-1pn` crate:

- 500 Mercury orbits under IAS15 + 1PN + unsoftened gravity.
- Measured perihelion drift: $+42.983$ arcsec/century.
- General-Relativistic prediction: $+43$ arcsec/century.
- Relative error: $\approx 1$ part per million on developer hardware
  (at the f64 noise floor of the test-particle 1PN approximation;
  the prior `9caaef2` IAS15 controller refactor exposed a latent
  velocity-prediction flaw that, once fixed, moved the residual
  from a 4.4 ppm systematic bias to ~1 ppm stochastic round-off —
  see [`experiments/2026-04-28-ias15-velocity-prediction-bug.md`](experiments/2026-04-28-ias15-velocity-prediction-bug.md)).

This figure is asserted in CI via
`cargo test --release -p apsis-1pn -- --ignored` at a 100-ppm
threshold that absorbs cross-platform LLVM / libm variance, so any
regression that would invalidate the paper's headline figure fails
the build.

### 5.2 Cross-implementation parity — REBOUND IAS15

A reproducible parity portfolio compares `apsis` IAS15 against
REBOUND's IAS15 across four canonical scenarios spanning periodic
2-body, periodic 3-body, chaotic 3-body, and sign-flipped 2-body
regimes. Each scenario carries an *a priori* protocol notebook
(initial conditions, integrator settings, metrics, and tolerances
declared before the run) and a self-contained Python harness under
[`../validation/rebound-parity/`](../../validation/rebound-parity/).
All four pass at f64 machine epsilon in regime:

| Scenario | Regime | Horizon | Gated metrics | Worst observed | Protocol |
| --- | --- | --- | --- | --- | --- |
| Kepler-prograde ($e = 0.5$) | periodic 2-body, $L_z > 0$ | 100 orbits | 7 (orbital elements + energy) | 1–3 ULP | [2026-04-25](experiments/2026-04-25-rebound-parity-kepler.md) |
| Figure-8 (Chenciner–Montgomery) | periodic 3-body, $L_z = 0$ | 10 periods + 50 informational | 12 (3 evidentiary tiers) | 1 ULP | [2026-04-26](experiments/2026-04-26-rebound-parity-figure8.md) |
| Pythagorean (Burrau 1913) | chaotic 3-body | 70 canonical t.u. | structural invariants ($\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$) at f64 floor; energy bound exceeded symmetrically by both implementations in the close-encounter regime (regime mismatch with the smooth-flow bound, not a parity defect) | 98 % event alignment | [2026-04-30](experiments/2026-04-30-rebound-parity-pythagorean.md) |
| Kepler-retrograde ($L_z < 0$) | sign-flipped 2-body, long horizon | $10^4$ orbits + 100-orbit checkpoint | 10 (7 magnitude + 3 sign-consistency) at both horizons | 1–10 ULP | [2026-05-01](experiments/2026-05-01-rebound-parity-retrograde.md) |

The figure-8 parity scenario in particular drove an architecture audit
of the IAS15 controller against the algorithmic specification in Rein
& Spiegel (2015); the audit is documented in
[`experiments/2026-04-26-ias15-warmstart-bug.md`](experiments/2026-04-26-ias15-warmstart-bug.md)
and is the canonical reference for the controller-level invariants
the implementation now upholds. The retrograde long-horizon scenario
closes the long-horizon stability gate identified during the
GR-readiness review as a precondition for 1PN-class extensions in
the federation thesis.

### 5.3 Unit-test surface

Beyond the reference signals, the library ships a unit-test surface
covering: energy conservation on canonical scenarios (Kepler,
Pythagorean three-body, figure-eight); IAS15 determinism on seeded
close encounters; conservation-contract assertions on the public API;
the IAS15 warmstart against the Pascal-triangle transformation
(Everhart 1985); the twelve `apsis::contract` tests covering kernel
invariants, composition rules, and the failure model. The 1PN plugin
ships its own suite, including the Mercury-precession gate and
softening-violation contract diagnostics.

---

## 6. Known limitations

- **Wisdom-Holman carries documented algorithmic defects (TD-008)**
  and is not treated as a quality signal in validation runs.
- **Large-N is not the design centre.** The solver handles up to
  $\sim 10^3$ bodies comfortably in the currently validated regime;
  for collisionless $10^6+$-body problems, use GADGET or PKDGrav.
- **Per-step cost of IAS15 is unbounded in stiff regimes.** Intrinsic
  to an adaptive high-order integrator. The soft warn at
  `N > ADAPTIVE_BODY_SOFT_WARN` surfaces this at integrator selection
  time; opt into a fixed-step integrator if a bounded per-step cost
  is required. See [ADR-013](adr/013-default-integrator-ias15.md).

---

## 7. Further reading

| Topic | Document |
| --- | --- |
| Integrator contracts & selection rubric | [`integrator.md`](integrator.md) |
| Plummer softening derivation | [`softening.md`](softening.md) |
| Architectural decisions | [`adr/`](adr/) |
| Reproducible experiments (lab notebooks) | [`experiments/`](experiments/) |
