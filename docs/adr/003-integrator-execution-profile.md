# ADR-003 — Integrator Execution Profile and Force-Model Compatibility

**Status:** Accepted
**Date:** 2026-04-23
**Branch:** bench/ias15-structured-rings-quality

---

## Context

The simulator shipped with IAS15 as the default integrator, paired
with Barnes-Hut for force evaluation above `EXACT_THRESHOLD = 64`.
That default was inherited from the position "IAS15 is what REBOUND
uses; match REBOUND." It was wrong for us.

Two facts, taken together, make the pairing untenable:

1. **IAS15 requires bit-reproducible forces.** The 15th-order method
   solves an implicit system by Picard predictor-corrector within each
   adaptive sub-step. Convergence of that iteration requires
   `f(x, v, t)` to be a deterministic function of state —
   bit-identical across calls with the same `(x, v, t)` to within
   f64 ULP.

2. **Barnes-Hut is not bit-reproducible across Picard iterations.**
   Between iterations, body positions drift by small amounts (that is
   what Picard is for). The quadtree is rebuilt on each `compute`
   call. A body near a cell boundary crosses to a different leaf in
   response to sub-ULP position drift; the multipole approximation
   for that body's far-field changes discretely; the force on that
   body changes discretely across iterations. Picard interprets this
   as non-convergence, the controller rejects the attempt, shrinks
   `dt`, and cascades toward `DT_MIN`.

This was not a theoretical concern. Evidence:

* **`solar_n641` stress test**: at N=641 with BH at θ=0.5, IAS15
  produced a 194% rejection rate (pre RMS-norm fix) and degenerated
  into tens of thousands of rejected evaluate calls per duration
  unit. The RMS-norm fix dropped wall time 23.5% but did not
  eliminate the cascade.

* **`structured_rings_n641` diagnostic**: a scenario designed to be
  maximally friendly to IAS15 (N=641, regular geometry, pair
  separations bounded well above the softening length) still
  cascaded when paired with BH. 19 minutes of single-thread CPU for
  `duration = 6.0` and it did not complete. Forcing the engine to
  direct O(N²) (`set_exact_threshold(usize::MAX)`) did not eliminate
  the cascade either — in the time budget of that diagnostic, we
  could not exclude that the scenario itself has stiffness we did
  not characterise. But the pairing incompatibility is independent
  of scenario: it is a property of BH vs. Picard's fixed-point
  requirement.

* **Reference implementation**: REBOUND pairs IAS15 exclusively with
  direct O(N²) summation. It is not an implementation shortcut but a
  consequence of the same determinism analysis.

Further: interactive testing of the `solar_system` preset (641 bodies)
confirmed the pathology under realistic loads. The `fine` rendering
preset failed to maintain frame cadence, and `medium` exhibited
sustained periodic hitches. REBOUND deliberately positions IAS15 as
an offline integrator driven by `reb_integrate(sim, tmax)`, not as
a render-loop integrator. The original `apsis` configuration imposed
a constraint the algorithm does not satisfy at any scale: bounded
per-step wall time.

The original framing — "tune the scenario so IAS15 is happy" — was
optimisation-of-symptom. The architectural frame — "IAS15 is a
precision method whose role is offline; the default must be a
real-time-compatible integrator" — is the honest one.

---

## Decision

The simulator adopts a **single structural contract** for integrator
and force-model compatibility, implemented as two trait-level
properties and enforced in exactly one place.

### Contract

* **Each integrator declares two properties**:
  * `execution_profile() -> ExecutionProfile` — `Realtime` (bounded
    per-step wall time, safe for the render loop) or `Precision`
    (unbounded per-step wall time; must run off-thread to
    completion). Defaults `Realtime`.
  * `requires_deterministic_force() -> bool` — whether the
    integrator's mathematical construction requires the force
    function to be bit-reproducible across calls with identical
    state. Defaults `false`.

* **Each force model declares one property**:
  * `is_deterministic() -> bool` — whether the current configuration
    of the force model produces bit-reproducible forces (a property
    of configuration, not only of type). Defaults `true`.

* **`System::set_integrator` is the single enforcement point**. When
  the new integrator requires deterministic forces and the current
  force model is not deterministic, the force model is
  auto-reconfigured (exact threshold raised so the BH branch is
  unreachable) and a structured `warn_diag!` event is emitted on
  stderr. No silent
  behaviour change: the user sees exactly what was corrected and
  why.

### Integrator overrides

| Integrator      | Execution profile | Requires deterministic force |
| --------------- | ----------------- | ---------------------------- |
| Velocity Verlet | Realtime          | false                        |
| Yoshida 4       | Realtime          | false                        |
| Wisdom–Holman   | Realtime          | false                        |
| IAS15           | **Precision**     | **true**                     |

### Default integrator

`System::new` now constructs with `IntegratorKind::Yoshida4` (4th-
order symplectic, bounded per-step cost, publication-quality
orbital conservation at real-time cadence). IAS15 remains available
via `set_integrator` but is opt-in, with the user understanding it
is a precision tool.

`PhysicsConfig::default` matches. The UI's integrator combo box
retains IAS15 as a selectable option with its existing label and
description.

### DT_MIN as scenario-stiffness signal

`DT_MIN = 1e-12` stays in IAS15 as the adaptive-shrink floor. When
the controller saturates the floor (`dt_try <= DT_MIN`) and accepts a
degraded step on that branch, IAS15 now emits a per-occurrence
`eprintln!` WARN with the current `dt`, the floor value, and the
running degraded count. This reframes the floor from a silent
rubber-stamp into a diagnostic signal that surfaces scenario-level
stiffness (close-encounter geometry below softening, N too high for
the controller to resolve). Deadline-driven degraded accepts — the
cooperative budget exhaustion path used in interactive runs — stay
silent by design; they are not a scenario indictment.

### What this is NOT

* **Not a mode flag.** There is no `ExecutionMode { Interactive,
  Offline }` top-level state. The profile lives on the integrator
  and is consumed by callers (physics thread, UI) as a derived
  value. One source of truth.

* **Not a force-model hierarchy.** The boolean `is_deterministic` is
  sufficient today. Both trait doc-comments flag the evolution path
  to a `DeterminismLevel { Strict, Approximate { bound }, Nondeterministic }`
  enum for when a second non-trivial force model (FMM, GPU with
  reduction noise) makes the `Strict` / `Approximate` distinction
  load-bearing. Until then the boolean does not encode spurious
  precision.

* **Not a user-visible lockdown.** The UI will expose force-model
  overrides in an Advanced panel (follow-up work), gated by a
  warning. Users who want to reintroduce the cascade for
  debug / experimentation can do so consciously.

---

## Consequences

### Good

* Interactive playback regains its contract: the default integrator
  has bounded per-step cost. `fine` preset at N=1400+ becomes
  tractable again because it no longer runs an unbounded-cost
  algorithm.
* The pairing rule is centralised. Adding a new integrator
  (e.g. WHFast with its own determinism requirements) plugs into the
  existing enforcement without touching physics-thread, UI, or
  bench code.
* IAS15 remains a first-class citizen for what it is good at: paper-
  quality trajectories, long-term Kepler / Pythagorean reference
  runs, determinism-audit work. The bench harness continues to
  cover it with bit-exact validation against `ias15.toml`.
* The `DT_MIN` floor becomes informative rather than mute. A user
  running IAS15 on a pathological scenario now sees it as a log
  event the moment it happens.

### Neutral

* The UI combo box still lists IAS15. Users who select it
  consciously accept the stutter; callers of `set_integrator` see
  the auto-correction warning.
* Existing benchmark scenarios retain their IAS15 coverage. The
  `solar_n641` stress test and `structured_rings_n641` diagnostic
  are unaffected — the pairing is now automatic via the enforcement
  in `set_integrator`.

### Watch out

* Precision-run UI chrome (banner, progress bar, disable
  interactive perturbations) is follow-up work and not part of this
  ADR's scope. Today, selecting IAS15 in the UI drops the user into
  the current real-time loop; the warn log tells them the pairing
  adjusted, but the loop itself will stutter at large N. The
  precision-run chrome makes that experience coherent.
* `PhysicsConfig::default().integrator` now disagrees with any
  serialised `run.toml` that pinned IAS15 as the integrator. Load
  paths already drive integrator selection explicitly (see
  `save_modal.rs`), so this is cosmetic rather than breaking.
* The boolean determinism property will need to become an enum when
  the next force model lands. Both doc-comments call this out; the
  call-site in `System::set_integrator` will be the single point of
  update.

---

## Implementation

Delivered in commit 8d655d0:

* `ExecutionProfile` enum in `physics::integrator::traits`.
* `Integrator::{execution_profile, requires_deterministic_force}`
  trait methods with defaults; IAS15 overrides.
* `ForceModel::is_deterministic` trait method with default;
  `GravityForceModel` returns `true` only when
  `exact_threshold() >= 10_000` (engine's clamp ceiling — the BH
  branch is unreachable for any practical N at that threshold).
* `System::set_integrator` enforcement + `eprintln!` WARN.
* `System::new` and `PhysicsConfig::default` switched to Yoshida 4.
* IAS15 per-occurrence WARN on `dt_try <= DT_MIN` degraded accepts.
* `ScenarioSpec::force_exact` diagnostic flag removed (pairing is
  now automatic).
* `structured_rings_n641` kept in `scenarios.rs` as a named diagnostic
  but excluded from `all()` — its doc-comment explains the status.

Tests: 179 library tests pass; no behaviour regressions.

Follow-up (separate branches):

* Precision-run UI chrome: banner + progress bar + perturbation lock
  during IAS15 runs. Consumer: `ExecutionProfile::Precision`.
* Advanced UI panel exposing force-model overrides with warning.
* `DeterminismLevel` enum upgrade when the next force model arrives.
* ADR consumer documentation: update `docs/integrator.md` to reflect
  the new default and the profile contract.
