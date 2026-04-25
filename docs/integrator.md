# Integrators

The simulator supports four integration schemes, chosen at runtime via
`System::set_integrator` or the UI combo box. All share the same
[`ForceModel`](../crates/apsis/src/physics/integrator/force_model.rs) interface
and the same `IntegratorContext` plumbing — changing integrator
never touches the force, perturbation, or diagnostic code.

This document captures the contract (what each integrator is for,
what it requires) rather than the numerical derivation (which lives
in the source-level doc-comments).

## Selection rubric

| If you want...                               | Pick                         |
| -------------------------------------------- | ---------------------------- |
| Interactive playback at any N                | **Yoshida 4** (default)      |
| Cheaper playback (lower accuracy OK)         | Velocity Verlet              |
| Paper-grade Keplerian hierarchy              | Wisdom–Holman                |
| Paper-grade trajectory with close encounters | **IAS15** (precision mode)   |

`Yoshida 4` is the default in both `System::new` and
`PhysicsConfig::default`. It is 4th-order symplectic, has bounded
per-step wall time at any realistic N, and conserves orbital energy
at publication quality for bound orbits.

## Execution profile — real-time vs precision

Every integrator declares an
[`ExecutionProfile`](../crates/apsis/src/physics/integrator/traits.rs) that
downstream code (physics thread, UI) reads to decide how to drive
it.

| Integrator      | Profile     | Why                                                                 |
| --------------- | ----------- | ------------------------------------------------------------------- |
| Velocity Verlet | `Realtime`  | Fixed per-step cost, O(N²) or O(N log N) per step.                  |
| Yoshida 4       | `Realtime`  | Same — 4 evals per step but still bounded.                          |
| Wisdom–Holman   | `Realtime`  | Analytic Kepler + perturbation; no adaptation.                      |
| IAS15           | `Precision` | Adaptive Gauss-Radau; `dt` can shrink unboundedly in stiff regimes. |

`Precision` means the caller must expect unbounded per-step wall
time. In practice this means running the integrator off-thread to
completion with a progress indicator, not inside a 60 Hz render
loop. REBOUND pairs IAS15 with a scripted `reb_integrate(sim, tmax)`
entry point for the same reason: IAS15 is an offline precision tool,
not a real-time one.

## Force-model determinism contract

Some integrators require the force function `f(x, v, t)` to be a
deterministic function of state — bit-reproducible across calls with
the same inputs to within f64 ULP.

| Integrator      | `requires_deterministic_force` | Why                                                                 |
| --------------- | ------------------------------ | ------------------------------------------------------------------- |
| Velocity Verlet | `false`                        | 2nd-order; per-step error absorbs tree-noise at O(dt²).             |
| Yoshida 4       | `false`                        | Same argument at O(dt⁴).                                            |
| Wisdom–Holman   | `false`                        | Kepler is analytic; perturbation is stepped at low order.           |
| IAS15           | **`true`**                     | Picard predictor–corrector diverges under non-deterministic forces. |

Force models declare the dual property via
[`ForceModel::is_deterministic`](../crates/apsis/src/physics/integrator/force_model.rs).
The default gravity engine returns `true` iff its
`exact_threshold >= DIRECT_MODE_THRESHOLD` — i.e. the Barnes-Hut
branch is unreachable for any practical N.

**Barnes-Hut is NOT deterministic in the Picard sense.** The
position-dependent quadtree rebuild means a body near a cell
boundary can cross leaves in response to sub-ULP position drift
between Picard iterations, which changes the multipole approximation
for that body's far-field discretely. IAS15's fixed-point iteration
reads this as non-convergence and cascades `dt` toward `DT_MIN`. See
[ADR-003](adr/003-integrator-execution-profile.md) for the full
derivation and the production-scale evidence that drove this
constraint.

## Enforcement

`System::set_integrator` is the **single enforcement point** for
the pairing rule. When the new integrator requires deterministic
forces and the current force model is not deterministic, the force
model is auto-reconfigured (exact threshold raised so BH is
bypassed) and a `warn_diag!` event is emitted with structured
fields (`integrator`, `exact_threshold_before`, `exact_threshold_after`).

Downstream code (physics thread, UI, benchmark runner) does not
re-check the pairing. The invariant holds by construction after
each `set_integrator` call.

A second advisory fires at the same call site if the new integrator
is `Precision` and the current body count is above
`PRECISION_BODY_SOFT_WARN` (200). This is a hint, not a block: the
user may proceed with IAS15 at large N, but the soft warn surfaces
the stutter expectation early rather than when the first frame drop
arrives.

## Scenario stiffness signal (IAS15 only)

IAS15's adaptive controller has a hard floor at `DT_MIN = 1e-12`.
When the controller wants to shrink below the floor (typically due
to close-encounter geometry below the softening length or N too
high for the controller to resolve at the configured ε), it accepts
a degraded step and increments `AdaptiveStats.degraded`.

Each such floor-hit emits a `warn_diag!` event with the current
`dt`, the running count, and a hint. The log rate is thinned:
first three occurrences verbatim, then every power of two
(4, 8, 16, ...). This keeps the signal visible at low frequency
without drowning stderr on pathological scenes.

Degraded accepts triggered by the cooperative deadline (physics
thread budget exhausted) do not emit a log — they are expected in
interactive precision runs and not a scenario indictment. Both
causes still accumulate into the unified `degraded` counter.

## Further reading

* Source: [`crates/apsis/src/physics/integrator/`](../crates/apsis/src/physics/integrator/)
* ADR-003: [`adr/003-integrator-execution-profile.md`](adr/003-integrator-execution-profile.md)
  — full rationale for the execution-profile and determinism
  contract.
