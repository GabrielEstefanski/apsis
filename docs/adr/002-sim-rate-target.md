# ADR-002 — Sim-Rate Target as Primary Speed Control

**Status:** Accepted  
**Date:** 2026-04-21  
**Branch:** feat/sim-rate

---

## Context

ADR-001 replaced `steps_per_frame` with `batch_budget_ms` — a wall-clock budget the
physics thread may consume per frame.  This is an improvement over the raw step count,
but it still exposes the CPU mindset to the user: "how many milliseconds of CPU may the
sim eat?" is the wrong question.

The right question — the one every serious N-body code (REBOUND, Gadget, Universe Sandbox)
answers — is: **"how fast should simulated time advance?"**

REBOUND's canonical API is `reb_integrate(sim, tmax)`: the caller sets a sim-time target
and the integrator fills it however it can.  Universe Sandbox exposes a "Simulation Speed"
slider in yr/s that the user can override; the underlying CPU throttle is invisible.

`batch_budget_ms` leaks two implementation details:

1. It assumes the user knows (or cares) how much CPU a physics step costs.
2. It couples simulation speed to hardware — the same slider produces different yr/s on
   different machines.

A `sim_rate_target` (sim units/s, displayed as yr/s) is hardware-agnostic: the user says
"simulate 1 yr/s" and the physics thread does whatever it takes to achieve that, up to the
invisible CPU safety cap.

---

## Decision

Replace the user-facing `batch_budget_ms` slider with a `sim_rate_target` control
(sim units/s, displayed as yr/s).

### Core loop change

Each batch, compute the sim-time the frame should advance:

```
t_target = system.t() + sim_rate_target × wall_delta
```

where `wall_delta` is the wall time since the previous batch started (clamped to avoid
spiral-of-death on slow frames).  The batch loop runs `while system.t() < t_target` with
the existing `MAX_BATCH_WALL_MS` (500 ms) as the hard CPU cap.

### What changes

| Before (ADR-001) | After (ADR-002) |
|------------------|-----------------|
| `batch_budget_ms: u32` — CPU budget | `sim_rate_target: f64` — sim units/s target |
| Slider: 1–200 ms (log) | Slider: 0.1 yr/s – 10 000 yr/s (log) |
| Label: `100ms` | Label: `1.0 yr/s` |
| Physics: run until deadline | Physics: run until `t_target` or CPU cap |
| Actual yr/s: derived, read-only | Actual yr/s: primary feedback |

### Dense-output bonus

With sim-rate as the primary control, the renderer's `t_render` advance is trivial:

```
t_render += sim_rate_target × wall_delta_render
```

No exponential moving average needed (that was a workaround for missing sim-rate target).

### What stays invisible

`MAX_BATCH_WALL_MS = 500` remains as a hard CPU cap.  If the physics thread cannot keep up
(e.g. large N, small dt, high yr/s target), the actual yr/s drops below target — the
shortfall is displayed in the playbar so the user knows why the simulation feels slow.

---

## Consequences

**Good:**
- Hardware-agnostic speed control — same slider value means the same simulated yr/s on
  every machine.
- Natural "fast-forward" semantics: drag to 1 000 yr/s → the sim advances 1 000 years per
  real second (or as fast as it can).
- Dense-output `t_render` arithmetic becomes a single multiply — no EMA, no teleport risk.
- Aligns with REBOUND / Universe Sandbox mental model.

**Neutral:**
- Physics thread must measure elapsed wall time between batches (`batch_start` already
  exists; add `prev_batch_start`).
- Paused state: `sim_rate_target` is irrelevant; thread sleeps as before.

**Watch out:**
- Very high targets with IAS15's adaptive dt: the integrator may take large steps that
  overshoot `t_target` slightly.  Acceptable — overshoot by at most one step.
- dt-driven integrators (VV, Y4, WH) with a fixed small dt and a high yr/s target will
  hit `MAX_BATCH_WALL_MS` before reaching `t_target`.  The displayed shortfall communicates
  this clearly.

---

## Implementation Plan

1. Add `sim_rate_target: f64` to `SimulationApp`; default 2π (≈ 1 yr/s).
2. Replace playbar budget slider with yr/s slider (log, 0.01 – 100 000 sim-units/s).
3. Add `PhysicsCmd::SetSimRateTarget(f64)` + `PhysicsHandle::set_sim_rate_target`.
4. Physics loop: replace deadline-based loop with t_target-based loop.
5. Display: show `(actual / target) %` tint when actual < 80% of target.
6. Smoke test: `cargo test` + manual check at 1 yr/s, 1 000 yr/s, max yr/s.
