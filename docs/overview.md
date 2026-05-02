# System Overview

## 1. Purpose and scope

`apsis` is an open-source N-body gravitational simulator implemented
in Rust. The project's published contribution is the `apsis`
library; the interactive visualisation shell `apsis-app` is an
optional side consumer, not part of the library's validated public
surface (see the workspace [`README.md`](../README.md) for the paper-level
positioning).

The solver is **3D-aware** as of the v0.1 alpha: `Vec3` value type,
per-body `z` / `vz`, 3D-aware direct N-body kernel, IAS15 dense
interpolation in 3D, orbital elements with inclination $i$ and node
$\Omega$, and 3D-aware energy / angular-momentum / centre-of-mass
observables. The previous 2D-only surface was frozen until the
API-contract machinery had been exercised end-to-end against the
Mercury perihelion test; the 3D port followed without breaking the
contract surface.

---

## 2. Physical model and assumptions

The core physics follows classical Newtonian mechanics with an optional
first-post-Newtonian correction supplied by an out-of-tree plugin.

**Modelled:**

- Universal gravitation (Newton's law of gravitation, `G = 1` in
  simulation units).
- Isolated system boundary: no external forces beyond those registered
  explicitly via the [`PerturbationForce`](../crates/apsis/src/physics/integrator/perturbation.rs)
  extension point.
- Plummer-softened pairwise forces with per-body material-scaled
  softening length — see [`softening.md`](softening.md) for the
  derivation and the published convention this follows.
- Relativistic corrections as an opt-in plugin (see
  [`apsis-1pn`](../crates/apsis-1pn/)) — Schwarzschild
  test-particle 1PN, applied pairwise.

**Not modelled:**

- Collisions, mergers, or fragmentation (bodies are point masses even
  when the plugin layer measures their physical radius).
- Radiation pressure beyond the experimental support in
  `physics/radiation` (not part of the validated surface).
- Hydrodynamics, stellar structure, or collisionless large-N.
- Large-$N$ scaling beyond the currently validated regime
  ($N \le 10^3$).

---

## 3. Workspace architecture

The project is a Cargo workspace of six crates, split by role:

| crate | role |
|---|---|
| [`apsis`](../crates/apsis/) | The library: domain types, physics, integrators, `apsis::contract` formalisation, public API. Resolves no UI dependency — enforced in CI. |
| [`apsis-1pn`](../crates/apsis-1pn/) | First downstream force crate: Einstein-Infeld-Hoffmann (Schwarzschild limit) 1PN correction. Reference implementation of the federation contract. Depends only on `apsis` through the public API. |
| [`apsis-py`](../crates/apsis-py/) | Python binding (PyO3, abi3-py39). Researcher-first kwargs API exposing `Body`, `IntegratorKind`, `System`, `Trajectory`, and adaptive-controller diagnostics. |
| [`apsis-py-core`](../crates/apsis-py-core/) | Cross-extension transport (rlib): `Box<dyn PerturbationForce>` ↔ `PyCapsule`. Allows out-of-tree force crates to expose Python bindings without re-implementing physics. |
| [`apsis-1pn-py`](../crates/apsis-1pn-py/) | Python binding for `apsis-1pn`. Reference implementation of the federation contract at the Rust/Python boundary. |
| [`apsis-app`](../crates/apsis-app/) | Optional interactive shell: egui/wgpu event loop, camera, panels, and GPU-side rendering. Not part of the validated surface. |

The dependency direction is monotone: every binding and every force
crate depends on `apsis` through the public extension API only —
never `pub(crate)`, never core internals. Adding a force is adding a
crate. `cargo tree -p apsis` resolves zero UI dependencies. A CI job
asserts this on every push.

### 3.1 Core subsystems

Inside `apsis`, the code is organised by responsibility:

- [`domain/`](../crates/apsis/src/domain/) — `Body`,
  `NamedBody`, `Material`, softening primitives, and the
  `BodyField` plugin trait used by the UI to surface scalar fields
  per body.
- [`physics/`](../crates/apsis/src/physics/) — force
  models (Barnes-Hut + Plummer kernel), integrators (Velocity
  Verlet, Yoshida 4, Wisdom-Holman, IAS15), orbital elements,
  energy and angular-momentum diagnostics, and the
  `PerturbationForce` extension trait.
- [`core/`](../crates/apsis/src/core/) — the `System`
  orchestrator, adaptive dt/θ controllers, hook registry,
  structured diagnostic event bus, trail ring buffer, and
  metrics assembly.
- [`io/`](../crates/apsis/src/io/) — snapshot
  serialisation (`.grav` binary format), CSV recorder, headless
  run configuration (`.toml`).
- [`templates/`](../crates/apsis/src/templates/) — the
  `TemplateKind` enum and builders for built-in presets (Solar
  System, TRAPPIST-1, figure-eight, etc.), consumed uniformly by
  the app, the headless runner, and test scripts.

---

## 4. Numerical integration

Four integrators are available, each declared with an
[`ExecutionProfile`](../crates/apsis/src/physics/integrator/traits.rs)
that downstream code reads to decide how to drive them:

| integrator | order | execution profile | primary use |
|---|---|---|---|
| Velocity Verlet | 2nd (symplectic) | Real-time | Quick exploration, educational runs |
| Yoshida 4 | 4th (symplectic) | Real-time | **Default**; interactive playback at any N |
| Wisdom-Holman | mixed-order | Real-time | Informational only — carries four documented algorithmic defects (TD-008) |
| IAS15 | 15th (adaptive Gauss-Radau) | Precision | Machine-precision off-line integration |

See [`integrator.md`](integrator.md) for the detailed contract of each
integrator, including the force-model determinism requirement that
pairs IAS15 with direct-O(N²) gravity rather than Barnes-Hut; see
[ADR-003](adr/003-integrator-execution-profile.md) for the rationale.

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

Beyond the reference signals, the library ships an extensive unit-test
surface covering: energy conservation on canonical scenarios (Kepler,
Pythagorean three-body, figure-eight); IAS15 determinism on seeded
close encounters; conservation-contract assertions on the public API;
direct unit tests pinning the IAS15 warmstart against the analytical
Pascal-triangle transformation derived in Everhart (1985); the
twelve `apsis::contract` tests covering kernel invariants, composition
rules, and the failure model as executable specification; and 3D
validation tests for orbital algebraic correctness and integrator
behaviour. The 1PN plugin ships a further 13 tests including the
Mercury-precession gate and softening-violation contract diagnostics.

---

## 6. Known limitations

- **Wisdom-Holman carries documented algorithmic defects (TD-008)**
  and is not treated as a quality signal in validation runs.
- **Single-precision trails.** The interactive trail buffer stores
  positions as `f32` for GPU efficiency. This affects only the
  visual history, not the physics — the integrator state is `f64`
  throughout.
- **Large-N is not the design centre.** The solver handles up to
  $\sim 10^3$ bodies comfortably in the currently validated regime;
  for collisionless $10^6+$-body problems, use GADGET or PKDGrav.
- **Integrator zoo is deliberately narrow.** Only the four listed
  above are provided. Higher-order composition methods (SABA),
  fourth-order Hermite, and hybrid close-encounter switchers
  (MERCURIUS) are out of current scope.
- **Per-step cost of IAS15 is unbounded in stiff regimes.** This
  is intrinsic to an adaptive high-order integrator; driving
  IAS15 from a render loop at high N is not recommended. See
  [ADR-003](adr/003-integrator-execution-profile.md) for the
  rationale that led to the execution-profile contract.

---

## 7. Further reading

| Topic | Document |
| --- | --- |
| Integrator contracts & execution profiles | [`integrator.md`](integrator.md) |
| Plummer softening derivation | [`softening.md`](softening.md) |
| Wall-time budget vs steps-per-frame | [`adr/001-wall-time-budget.md`](adr/001-wall-time-budget.md) |
| Sim-rate as primary speed control | [`adr/002-sim-rate-target.md`](adr/002-sim-rate-target.md) |
| Integrator execution profile & force-model compatibility | [`adr/003-integrator-execution-profile.md`](adr/003-integrator-execution-profile.md) |
| IAS15 per-phase wall-time breakdown (experiment) | [`experiments/2026-04-22-ias15-phase-profile.md`](experiments/2026-04-22-ias15-phase-profile.md) |
| Picard noise-floor null result (experiment) | [`experiments/2026-04-22-picard-noise-floor.md`](experiments/2026-04-22-picard-noise-floor.md) |
| Operational-domain benchmarks (experiment) | [`experiments/2026-04-24-operational-domain-benchmarks.md`](experiments/2026-04-24-operational-domain-benchmarks.md) |
| Kepler cross-implementation parity protocol | [`experiments/2026-04-25-rebound-parity-kepler.md`](experiments/2026-04-25-rebound-parity-kepler.md) |
| Figure-8 cross-implementation parity protocol | [`experiments/2026-04-26-rebound-parity-figure8.md`](experiments/2026-04-26-rebound-parity-figure8.md) |
| IAS15 controller architecture audit | [`experiments/2026-04-26-ias15-warmstart-bug.md`](experiments/2026-04-26-ias15-warmstart-bug.md) |
