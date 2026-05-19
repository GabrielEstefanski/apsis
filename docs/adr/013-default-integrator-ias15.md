# ADR-013 — Coherent Defaults: IAS15 + Direct Summation

**Status:** Accepted
**Date:** 2026-05-19

**Supersedes (in part):** [ADR-003](003-integrator-execution-profile.md) §Decision (`Default integrator`)

**Depends on:** [ADR-010](010-extract-apsis-app.md) (apsis-app extraction)

---

## Context

[ADR-003](003-integrator-execution-profile.md) (2026-04-23) switched
the `System::new` default from IAS15 to Yoshida 4 to keep per-step
wall time bounded for the in-tree interactive shell. The decision was
correct for the constraint of the time — a render-loop driver cannot
absorb IAS15's unbounded per-step cost in stiff regimes.

[ADR-010](010-extract-apsis-app.md) (2026-05-15) moved the shell to
a separate repository. With the shell gone, the apsis core no longer
ships a driver that imposes a bounded-per-step-cost contract. The
remaining drivers — headless runs via `apsis-record`, Python `Sim`,
benchmark harness, reviewer scripts cited from the validation
portfolio — all run to completion off any frame budget.

The default therefore now serves a different audience: a reader who
opens the library via the paper or `apsis-record`'s reproducibility
claim and constructs `System::new(...)` to reproduce a published run.
That audience expects the integrator that the validation portfolio
uses for paper-grade trajectories. The portfolio uses IAS15
throughout (Mercury 1PN long-horizon, REBOUND parity, federation
gates). A Yoshida-4 default would silently hand the reader a
4th-order symplectic integrator whose energy oscillates at $O(dt^4)$
— adequate for most cases but not the integrator the paper cites.

Switching the integrator default to IAS15 surfaces a second
inconsistency that ADR-003 inherited: `GravityForceModel::default()`
ships with `exact_threshold = 64`, so Barnes-Hut activates above
that body count. IAS15 requires a deterministic force across Picard
iterations; Barnes-Hut's position-dependent tree rebuild violates
that contract (see [ADR-003] §Context). With both defaults active,
`System::new(many_bodies, ...)` followed by `step()` would silently
trigger a Picard cascade toward `DT_MIN`. Patching the constructor
to auto-raise the threshold would resolve the symptom, not the
incoherence; it would also mutate user-visible state silently in the
single API entry point where silent magic is most surprising.

The principle: **defaults across components should be internally
consistent for correctness; optimisations that alter mathematical
properties belong opt-in.** Barnes-Hut is a legitimate optimisation
for throughput, visualisation, and large-N exploration, but its
spatial discretisation and tree-topology dependence make it an
algorithmic approximation — not a transparent micro-optimisation.

## Decision

Two coupled changes that together yield coherent defaults:

1. `System::new(...)` constructs with `IntegratorKind::Ias15`.
2. `GravityForceModel::default()` initialises with
   `exact_threshold = DIRECT_MODE_THRESHOLD` (always-direct). The
   `EXACT_THRESHOLD = 64` constant is removed; Barnes-Hut is no
   longer the default code path at any body count.

Together they produce two coherent user profiles:

**Scientific / precision (default).** `System::new(bodies, units)`
yields IAS15 + direct summation: paper-grade trajectory, no
approximations, reproducible across runs and platforms.

**Performance / scaling (opt-in).** The user explicitly enables
Barnes-Hut via `set_exact_threshold(N)` (or the equivalent builder
method) and accepts the trade: approximation, possible tree-topology
discontinuities, accuracy↔throughput. If they also want IAS15 with
this configuration, the existing silent-pairing enforcement in
`set_integrator` raises the threshold back — the only path that
triggers silent reconfiguration is one the user invoked explicitly.

Additional notes:

- The Python `Sim(...)` constructor is unaffected — it has no
  default; the integrator is a required keyword argument.
- The adaptive-integrator scale advisory at `System::set_integrator`
  (soft warn when `is_adaptive() && N > ADAPTIVE_BODY_SOFT_WARN`)
  continues to fire and now reaches the default path.

## Consequences

### Good

- `System::new` produces a state that is safe by construction
  rather than safe by procedural patching. A reader auditing the
  defaults can verify coherence statically — no runtime mutation in
  the constructor.
- The default matches REBOUND's default (`IAS15`) and the integrator
  the validation portfolio uses. A reader who reproduces a portfolio
  run via `System::new` gets the same integrator the run was
  validated against, without having to read the source.
- The accuracy↔performance decision becomes explicit. Users who
  want Barnes-Hut throughput now confront its mathematical
  approximations as an active choice rather than inheriting them
  silently at moderate N.
- Removes a coupling that was load-bearing only while the in-tree
  shell existed — the design has caught up with the post-ADR-010
  topology.
- The silent-pairing enforcement at `set_integrator` rarely fires
  on real workloads now; defaults don't conflict, so it activates
  only when the user explicitly enabled Barnes-Hut and then asked
  for IAS15.

### Neutral

- Tests that explicitly set an integrator (every benchmark, every
  parity scenario, every controller test) are unaffected. The
  single test that depended on default Barnes-Hut at N=80
  (`ias15_selection_forces_deterministic_force_model`) is updated
  to call `set_exact_threshold(64)` explicitly in its baseline.
- The Python surface stays identical; user code does not change.

### Watch out

- A newcomer who does `System::new(...)` on a many-body scene and
  calls `step()` in a tight loop will hit IAS15's per-step variance
  (no Picard cascade — defaults are coherent — but adaptive cost
  itself is now the default behaviour). The soft warn at
  `N > ADAPTIVE_BODY_SOFT_WARN = 200` exists to surface this; the
  newcomer can opt into a fixed-step integrator via
  `set_integrator(...)` if a bounded per-step cost is needed.
- Users who relied on default Barnes-Hut for throughput (no
  external evidence anyone did) now need to opt in via
  `set_exact_threshold(N)`. The transition is a one-line change
  and the resulting API is more honest about what they're asking
  for.

## Implementation

- `System::new` default: `IntegratorKind::Ias15`.
- `BarnesHutEngine` default: `exact_threshold = DIRECT_MODE_THRESHOLD`.
  The `EXACT_THRESHOLD = 64` constant is removed; comments and tests
  that referenced it as a meaningful boundary are updated.
- `tests::ias15_selection_forces_deterministic_force_model` baseline
  configures `set_exact_threshold(64)` explicitly before asserting
  Barnes-Hut is active.
- `docs/integrator.md`, `docs/overview.md`: selection rubric, zoo
  table, and Enforcement section updated.
- `docs/softening.md`: "switch to Yoshida 4" recommendation kept as
  one of several stiff-regime mitigations, no longer framed as the
  default integrator.
- ADR-003 supersession note updated to reference this ADR.
