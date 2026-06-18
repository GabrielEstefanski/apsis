# ADR-015 — Exact Finish Time for `integrate_until`

**Status:** Accepted
**Date:** 2026-06-12

**PR:** [#207](https://github.com/GabrielEstefanski/apsis/pull/207) — shipped the
exact-finish-time fix and tightened the Mercury 1PN gate to `9.2×10⁻⁶`;
reconciled into the paper in #211.

**Amends:** the run-loop semantics documented in ADR-004 (sub-step
granularity); `System::integrate_until` / `integrate_for`.

---

## Context

`System::integrate_until(t_end)` looped `while t < t_end { step() }`,
exiting at the first accepted step with `t ≥ t_end`. The endpoint state
therefore sat up to one step past the requested time — for IAS15 one
adaptive sub-step (the configured `dt` is only the controller's seed,
ADR-004), about `5×10⁻³` time units on the Mercury 1PN scenario.

The Mercury error-budget analysis
(`paper/notebooks/2026-06-10-mercury-1pn-error-budget.md`) measured the
consequence. The osculating ω of a 1PN orbit carries an O(ε)
short-period term whose slope at periapsis is `−ε(3−e)(1+e)/e` per
radian of true anomaly — the `3/e` factor is the conditioning of the
eccentricity-vector direction. Sampling the endpoint half a sub-step
late moved the gate measurement by `−1.5×10⁻⁸` rad (`−6×10⁻⁵`
relative), and the sub-step-sized endpoint jitter across single-ULP
initial-condition twins produced the entire observed ensemble spread
(`σ ≈ 10⁻⁸` rad). Subtracting the derived endpoint term collapsed the
spread to `10⁻¹⁴` rad. The gate's `10⁻⁴` tolerance was nearly consumed
by the worst endpoint draw (`−1.16×10⁻⁴` observed on a twin); a
different platform's step sequence could fail it spuriously.

REBOUND integrates to the exact requested time by default
(`exact_finish_time=1`, Rein & Liu 2012; the final timestep is
reduced), so cross-code comparisons at "the same t" also inherited the
offset.

## Decision

1. `integrate_until` lands exactly on `t_end` **by default**: when the
   next step would cross the target, the step is clipped to the
   remainder. After the loop, `t` is set to `t_end` exactly,
   collapsing the final `t += dt` round-off.
2. Fixed-step integrators are clipped by the orchestrator through
   their `dt_hint` (they consume it verbatim). Self-adaptive
   integrators expose `Integrator::cap_next_step(max_dt)` — a one-shot
   cap honoured for a single `step()` call.
3. A clipped IAS15 step must not disturb the controller: the step
   samples the trajectory, not the error landscape. The pre-clip
   `dt_next` is restored after a capped accept, and the warmstart
   record is dropped (its `b` history is scaled to the clipped dt; a
   large-ratio extrapolation gives Picard a meaningless predictor).
   Cost per sampling boundary: the clipped step plus one cold-Picard step.
4. Opt-out: `System::set_exact_finish_time(false)` restores the
   run-whole-steps semantics, for callers that prefer an undisturbed
   fixed-step symplectic rhythm over endpoint accuracy (the same
   trade-off REBOUND documents for WHFast).

## Consequences

- The Mercury precession gate residual becomes the derived budget: the
  derivation floors plus the exact-endpoint phase deficit, `+4.581×10⁻⁶`
  relative, pre-registered from pre-fix data and reproduced to seven
  significant digits. The gate tightens from `10⁻⁴` to `9.2×10⁻⁶`
  (2× the deterministic central).
- ULP-twin endpoint spread drops from `10⁻⁸` rad to `10⁻¹⁴` rad over
  500 Mercury orbits; cross-platform step-sequence differences no
  longer move fixed-time measurements.
- A fixed-step integrator's final clipped step breaks its uniform-dt
  rhythm once per `integrate_until` call — a one-time, non-secular
  perturbation, opt-out available.
- Sampling loops should use absolute targets (`integrate_until(t_k)`),
  not accumulated relative durations, to avoid summing fp round-off in
  the targets themselves.
- The parity notebooks' tabulated endpoint deltas were measured under
  overshoot semantics and are queued for refresh.

## Verification

`crates/apsis/tests/exact_finish_time.rs` (exact landing for all seven
integrators, opt-out, energy floor, controller-rhythm bound);
`crates/apsis-1pn/examples/error_budget_run.rs` post-fix:
`t_overshoot = 0`, `rel_err` within the pre-registered band
(`validation/audit/ledger.md`, 2026-06-12 entry).
