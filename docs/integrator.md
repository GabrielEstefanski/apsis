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

## Per-step cost: bounded vs adaptive

`IntegratorKind::is_adaptive()` distinguishes the cost classes the
zoo exposes:

| Integrator       | `is_adaptive` | Per-step cost                                              |
| ---------------- | ------------- | ---------------------------------------------------------- |
| Velocity Verlet  | `false`       | Bounded by force evaluation.                               |
| Yoshida 4        | `false`       | Same — 4 evals per step, bounded.                          |
| Wisdom-Holman    | `false`       | Analytic Kepler + perturbation; no adaptation.             |
| WHFast           | `false`       | Same, with compensated summation.                          |
| Mercurius        | `false`       | Outer step is bounded; IAS15 only fires on close encounter.|
| Implicit Midpoint| `false`       | Fixed-step with bounded inner-iteration cap.               |
| IAS15            | **`true`**    | Adaptive controller can shrink `dt` arbitrarily.           |

Adaptive integrators (IAS15) advertise unbounded per-step wall time:
a stiff regime can spend seconds on one logical step while the
controller cascades toward `DT_MIN`. REBOUND pairs IAS15 with a
scripted `reb_integrate(sim, tmax)` entry point for the same reason.

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
bypassed). The correction is silent — inspect
`sys.force_model().exact_threshold()` to audit the post-state.

Downstream code does not re-check the pairing. The invariant holds
by construction after each `set_integrator` call.

A second advisory fires at the same call site if the new integrator
is adaptive and the current body count is above
`ADAPTIVE_BODY_SOFT_WARN` (200). Hint, not a block: the user may
proceed with IAS15 at large N, but the soft warn surfaces the
per-step-cost expectation early rather than when the first stall
arrives.

## Scenario stiffness signal (IAS15 only)

IAS15's adaptive controller has a hard floor at `DT_MIN` $= 10^{-12}$.
When the controller wants to shrink below the floor (typically due to
close-encounter geometry below the softening length or $N$ too high
for the controller to resolve at the configured $\epsilon$), it
accepts a degraded step and increments `AdaptiveStats.degraded`.

Each such floor-hit emits a `warn_diag!` event with the current `dt`,
the running count, and a hint. The log rate is thinned: first three
occurrences verbatim, then every power of two ($4, 8, 16, \ldots$).
This keeps the signal visible at low frequency without drowning stderr
on pathological scenes.

### Diagnostic counters in `AdaptiveStats`

`AdaptiveStats` (returned by `Integrator::adaptive_stats()`) carries
four cumulative counters in addition to the substep tally; together
they are the cheapest signal of controller health and are surfaced
unconditionally (no feature flag, single `saturating_add` per accept):

| Field | Healthy regime | What an elevated value means |
| --- | --- | --- |
| `rejections` | $\ll$ `substeps` | Controller is rejecting at a rate the spec's halving cannot keep up with |
| `picard_iters` / `attempts` | $\sim 2$–$3$ | Predictor–corrector is starting too far from the converged $b$ |
| `picard_stagnations` | $\ll$ `substeps` | Picard residual saturating above `PICARD_TOL` — typically a sign of warmstart bias against the true $b$ |
| `shrink_grow_cycles` | $\ll$ `substeps` | Controller chatter; the dt proposal alternates between shrinking and growing on consecutive accepts |
| `degraded` | $0$ | `DT_MIN` floor escape clause fired |

A run that disagrees on any of these by orders of magnitude from the
healthy regime is a controller-health regression, even when the gated
energy and angular-momentum metrics still pass at machine precision.
This is the methodological observation that motivated the figure-8
parity scenario in the validation portfolio: invariant-passing alone
is insufficient evidence of an honest cross-implementation match for
adaptive integrators (see
[`experiments/2026-04-26-ias15-warmstart-bug.md`](experiments/2026-04-26-ias15-warmstart-bug.md)
for the full argument).

For investigations that need finer resolution than these counters
provide, the optional `ias15-diag` Cargo feature compiles in a
per-attempt trace — one tab-separated line per controller decision,
gated at runtime by the env var `APSIS_IAS15_TRACE=1` and throttled
by `APSIS_IAS15_TRACE_CAP` (default 2000 events) so cascade scenarios
do not bury the diagnostic in $10^{8}$ duplicate lines.

## First-Same-As-Last (FSAL)

The IAS15 sub-step ends with a force evaluation at the accepted
end-of-step body positions; that result is also the next sub-step's
start-of-step $a_0$. Re-evaluating forces a second time at the same
positions is unnecessary work, so the integrator caches the post-accept
acceleration buffer and reuses it on the next call (the canonical FSAL
property of any explicit/implicit method whose stage-0 evaluation
coincides with the previous step's stage-end). The cache is invalidated
whenever any precondition that ties the cached $a_0$ to the next call's
parameters can change: capacity resize, `recenter_bodies` translation,
buffer-length mismatch, or a difference in `g_factor` or in the
perturbation count. The optimisation is bit-identical to the
non-cached path on the canonical parity scenarios (Kepler, figure-8) —
it elides the duplicate evaluation, not the result.

## Further reading

* Source: [`crates/apsis/src/physics/integrator/`](../crates/apsis/src/physics/integrator/)
* ADR-003: [`adr/003-integrator-execution-profile.md`](adr/003-integrator-execution-profile.md)
  — full rationale for the execution-profile and determinism
  contract.
