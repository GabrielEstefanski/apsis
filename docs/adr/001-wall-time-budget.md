# ADR-001 — Replace `steps_per_frame` with a wall-time budget model

**Date:** 2026-04-21
**Status:** Accepted — implementation on `feat/wall-budget`

---

## Context

`steps_per_frame: u32` controls how many physics steps run per frame. It
works for fixed-cost integrators but breaks with IAS15 (adaptive):

| Integrator      | Force evals / step | × 100k steps/frame |
|-----------------|--------------------|--------------------|
| Velocity Verlet | 1                  | 100 k evals        |
| Yoshida 4       | 4                  | 400 k evals        |
| IAS15           | ~14–23             | ~2.3 M evals       |

Symptoms: the UI freezes when pausing or deleting a body under IAS15 +
high `steps_per_frame`; the slider requires the user to know the
integrator's internal cost model.

**Quick-win already in `develop`:** `MAX_BATCH_WALL_MS = 33` in
`physics_thread.rs` — the batch breaks after 33 ms of wall-clock even
if `steps_per_frame` has not been reached. This fixes the immediate
freeze but does not fix the abstraction.

---

## Decision

Replace `steps_per_frame` with a **per-batch wall-time budget**
(`batch_budget_ms: u32`).

The physics thread runs steps until the wall-clock budget is consumed.
The integrator decides how many steps fit — not the user.

```rust
// Before
while i < steps_per_frame {
    system.step();
    i += 1;
}

// After
let deadline = Instant::now() + Duration::from_millis(batch_budget_ms);
while Instant::now() < deadline {
    system.step();
}
```

The existing `yr/s` display becomes the primary speed feedback.

---

## Alternatives rejected

| Alternative                      | Reason                                          |
|----------------------------------|-------------------------------------------------|
| Keep the quick-win as permanent  | Does not fix the wrong abstraction              |
| Target sim rate (`yr/s`)         | More complex; requires per-integrator cost estimates |

---

## Implementation plan (`feat/wall-budget`)

**Step 1 — `physics_thread.rs`**

- `PhysicsCmd::SetStepsPerFrame(u32)` → `SetBatchBudgetMs(u32)`
- `steps_per_frame: u32` → `batch_budget_ms: u32` in the inner loop
- Inner loop: `while i < steps_per_frame` → `while Instant::now() < deadline`
- Keep `MAX_BATCH_WALL_MS` as a hard safety cap above the maximum user budget

**Step 2 — `PhysicsHandle`**

- `set_steps_per_frame` → `set_batch_budget_ms`
- Update every call-site in `ui.rs`

**Step 3 — `SimulationApp`**

- `steps_per_frame: u32` → `batch_budget_ms: u32` (suggested default: 8 ms)

**Step 4 — `playbar.rs`**

- Slider `× N steps` → `X ms` (range: 1–100 ms)
- Remove the "↑ dt for speed" hint under IAS15 (becomes unnecessary)

**Step 5 — Snapshot / config**

- Verify whether `steps_per_frame` is persisted in `.grav` — if so, the field
  is silently ignored on load so existing saves keep working.

**Step 6 — Tests + smoke run**

- `cargo test` — no test references `steps_per_frame` directly
- Velocity Verlet: equal or better throughput
- IAS15: responsive UI at any budget

**Step 7 — Merge and cleanup**

- Remove `MAX_BATCH_WALL_MS` from `develop` (redundant after merge) or keep
  it as defence in depth.

---

## Consequences

**Positive.** The UI never freezes; the slider is integrator-agnostic; the
`yr/s` display becomes the natural speed feedback.

**Watch out.** `steps_per_frame` disappears from the API — a call-site
checklist is required before starting. Throughput behaviour changes under
Velocity Verlet on lightweight scenes (more steps fit in the same budget
than before).

---

## References

- REBOUND: `reb_integrate(sim, tmax)` — integrates until a simulation time,
  not until N steps.
- [`crates/gravity-sim-core/src/core/physics_thread.rs`](../../crates/gravity-sim-core/src/core/physics_thread.rs) — current batch loop + `MAX_BATCH_WALL_MS`.
- [`crates/gravity-sim-app/src/app/panel/playbar.rs`](../../crates/gravity-sim-app/src/app/panel/playbar.rs) — current speed slider.
