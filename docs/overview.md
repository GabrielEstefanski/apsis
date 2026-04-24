# System Overview

## 1. Purpose and scope

`gravity-sim` is an open-source N-body gravitational simulator implemented
in Rust. The project's published contribution is the `gravity-sim-core`
library; the interactive visualisation shell `gravity-sim-app` is an
optional side consumer, not part of the library's validated public
surface (see the workspace [`README.md`](../README.md) for the paper-level
positioning).

The solver is currently **2D**. A 3D extension is planned and will be a
deliberately breaking change; the current 2D surface was frozen so the
API-contract machinery could be exercised end-to-end against the Mercury
perihelion test before the coordinate dimension changes.

---

## 2. Physical model and assumptions

The core physics follows classical Newtonian mechanics with an optional
first-post-Newtonian correction supplied by an out-of-tree plugin.

**Modelled:**

- Universal gravitation (Newton's law of gravitation, `G = 1` in
  simulation units).
- Isolated system boundary: no external forces beyond those registered
  explicitly via the [`PerturbationForce`](../crates/gravity-sim-core/src/physics/integrator/perturbation.rs)
  extension point.
- Plummer-softened pairwise forces with per-body material-scaled
  softening length — see [`softening.md`](softening.md) for the
  derivation and the published convention this follows.
- Relativistic corrections as an opt-in plugin (see
  [`gravity-sim-1pn`](../crates/gravity-sim-1pn/)) — Schwarzschild
  test-particle 1PN, applied pairwise.

**Not modelled:**

- Collisions, mergers, or fragmentation (bodies are point masses even
  when the plugin layer measures their physical radius).
- Radiation pressure beyond the experimental support in
  `physics/radiation` (not part of the validated surface).
- Hydrodynamics, stellar structure, or collisionless large-N.
- 3D coordinates; see the scope note above.

---

## 3. Workspace architecture

The project is a Cargo workspace of three crates, split by role:

| crate | role |
|---|---|
| [`gravity-sim-core`](../crates/gravity-sim-core/) | The library: domain types, physics, integrators, public API. Resolves no UI dependency — enforced in CI. |
| [`gravity-sim-1pn`](../crates/gravity-sim-1pn/) | Out-of-tree plugin demonstration: Einstein-Infeld-Hoffmann (Schwarzschild limit) 1PN correction. Depends only on `gravity-sim-core` through the public API. |
| [`gravity-sim-app`](../crates/gravity-sim-app/) | Optional interactive shell: egui/wgpu event loop, camera, panels, and GPU-side rendering. Not part of the validated surface. |

The read direction is `gravity-sim-app → gravity-sim-core`. The core
crate has no awareness of the app, and
`cargo tree -p gravity-sim-core` resolves zero UI dependencies. A
CI job asserts this on every push.

### 3.1 Core subsystems

Inside `gravity-sim-core`, the code is organised by responsibility:

- [`domain/`](../crates/gravity-sim-core/src/domain/) — `Body`,
  `NamedBody`, `Material`, softening primitives, and the
  `BodyField` plugin trait used by the UI to surface scalar fields
  per body.
- [`physics/`](../crates/gravity-sim-core/src/physics/) — force
  models (Barnes-Hut + Plummer kernel), integrators (Velocity
  Verlet, Yoshida 4, Wisdom-Holman, IAS15), orbital elements,
  energy and angular-momentum diagnostics, and the
  `PerturbationForce` extension trait.
- [`core/`](../crates/gravity-sim-core/src/core/) — the `System`
  orchestrator, adaptive dt/θ controllers, hook registry,
  structured diagnostic event bus, trail ring buffer, and
  metrics assembly.
- [`io/`](../crates/gravity-sim-core/src/io/) — snapshot
  serialisation (`.grav` binary format), CSV recorder, headless
  run configuration (`.toml`).
- [`templates/`](../crates/gravity-sim-core/src/templates/) — the
  `TemplateKind` enum and builders for built-in presets (Solar
  System, TRAPPIST-1, figure-eight, etc.), consumed uniformly by
  the app, the headless runner, and test scripts.

---

## 4. Numerical integration

Four integrators are available, each declared with an
[`ExecutionProfile`](../crates/gravity-sim-core/src/physics/integrator/traits.rs)
that downstream code reads to decide how to drive them:

| integrator | order | execution profile | primary use |
|---|---|---|---|
| Velocity Verlet | 2nd (symplectic) | Real-time | Quick exploration, educational runs |
| Yoshida 4 | 4th (symplectic) | Real-time | **Default**; interactive playback at any N |
| Wisdom-Holman | mixed-order | Real-time | Dominant-primary Keplerian hierarchy |
| IAS15 | 15th (adaptive Gauss-Radau) | Precision | Machine-precision off-line integration |

See [`integrator.md`](integrator.md) for the detailed contract of each
integrator, including the force-model determinism requirement that
pairs IAS15 with direct-O(N²) gravity rather than Barnes-Hut; see
[ADR-003](adr/003-integrator-execution-profile.md) for the rationale.

---

## 5. Validation

The library's top-line validation is the Mercury perihelion precession
test shipped by the `gravity-sim-1pn` crate:

- 500 Mercury orbits under IAS15 + 1PN + unsoftened gravity.
- Measured perihelion drift: `+42.983 arcsec / century`.
- General-Relativistic prediction: `+43 arcsec / century`.
- Relative error: `4.4 × 10⁻⁶` (≈ 4 parts per million).

This number is asserted in CI via
`cargo test --release -p gravity-sim-1pn -- --ignored`, so any
regression that would invalidate the paper's headline figure fails
the build.

Beyond Mercury, the library ships ~200 unit tests covering energy
conservation on canonical scenarios (Kepler, Pythagorean three-body,
figure-eight), IAS15 determinism on seeded close encounters, and
conservation-contract assertions on the public API.

---

## 6. Known limitations

- **2D only.** See the scope note above.
- **Single-precision trails.** The interactive trail buffer stores
  positions as `f32` for GPU efficiency. This affects only the
  visual history, not the physics — the integrator state is `f64`
  throughout.
- **Large-N is not the design centre.** The solver handles up to
  ~10³ bodies comfortably; for collisionless 10⁶+ body problems,
  use GADGET or PKDGrav.
- **Integrator zoo is deliberately narrow.** Only the four listed
  above are provided. Higher-order composition methods (SABA),
  fourth-order Hermite, and hybrid close-encounter switchers
  (MERCURIUS) are out of scope for v1.
- **Per-step cost of IAS15 is unbounded in stiff regimes.** This
  is intrinsic to an adaptive high-order integrator; driving
  IAS15 from a render loop at high N is not recommended. See
  [ADR-003](adr/003-integrator-execution-profile.md) for the
  stutter-diagnosis that led to the execution-profile contract.

---

## 7. Further reading

| Topic | Document |
|---|---|
| Integrator contracts & execution profiles | [`integrator.md`](integrator.md) |
| Plummer softening derivation | [`softening.md`](softening.md) |
| Wall-time budget vs steps-per-frame | [`adr/001-wall-time-budget.md`](adr/001-wall-time-budget.md) |
| Sim-rate as primary speed control | [`adr/002-sim-rate-target.md`](adr/002-sim-rate-target.md) |
| Integrator execution profile & force-model compatibility | [`adr/003-integrator-execution-profile.md`](adr/003-integrator-execution-profile.md) |
| IAS15 per-phase wall-time breakdown (experiment) | [`experiments/2026-04-22-ias15-phase-profile.md`](experiments/2026-04-22-ias15-phase-profile.md) |
| Picard noise-floor null result (experiment) | [`experiments/2026-04-22-picard-noise-floor.md`](experiments/2026-04-22-picard-noise-floor.md) |
| Solar-system stutter diagnosis (experiment) | [`experiments/2026-04-22-solar-system-stutter-diagnosis.md`](experiments/2026-04-22-solar-system-stutter-diagnosis.md) |
