# Mercury 1PN error budget — Phase B ensemble

Implementation of the Phase-B protocol declared in
[`paper/notebooks/2026-06-10-mercury-1pn-error-budget.md`](../../paper/notebooks/2026-06-10-mercury-1pn-error-budget.md),
plus the Phase-B' endpoint-sampling correction (derived in
`paper/notebooks/scripts/error_budget_endpoint_symbolic.py`, verified in
`paper/notebooks/scripts/error_budget_endpoint_numerical.py`).

## Files

| File           | Purpose                                                                                          |
| -------------- | ------------------------------------------------------------------------------------------------ |
| `ensemble.py`  | Orchestrator: runs the Rust example across ULP/N/eps_b/constructor grids; phases `smoke`, `b1`, `b3`, `b4` |
| `out/`         | Generated artefacts (`b1.csv`, `b3.csv`, `b4.csv`, logs). Git-ignored.                           |

The Rust measurement binary lives at
[`crates/apsis-1pn/examples/error_budget_run.rs`](../../crates/apsis-1pn/examples/error_budget_run.rs).
Constants (`A`, `E`, `M_MERCURY`, IC construction, osculating-ω measurement) mirror
`crates/apsis-1pn/tests/mercury_precession_gate.rs` exactly.

## CSV schema

```text
orbits,ulp,constructor,eps_b,measured_rad,predicted_rad,rel_err,t_overshoot,nu_end
```

`rel_err` is signed: `(measured − predicted) / predicted`. `t_overshoot`
is the time by which `integrate_until` exceeded `N · el0.period` (up to
one adaptive IAS15 sub-step). `nu_end` is the osculating true anomaly at
the endpoint; the Phase-B' offset function
`Q(ν) = ε(3ν − (3/e − e)sinν − (5/2)sin2ν)` converts it into the
deterministic endpoint-sampling part of the residual.

## Quick start

```powershell
# From the workspace root:
python validation/mercury-1pn-error-budget/ensemble.py --phase smoke
```

Requires Python 3.9+ (stdlib only — no pip install needed).

## Phases

| Phase   | Runs          | Output         | Reports                                                   |
| ------- | ------------- | -------------- | --------------------------------------------------------- |
| `smoke` | 3 × 2 = 6    | stdout table   | raw\_c ulp=0 \|rel\_err\| < 1e-4                          |
| `b1`    | 25 × 2 = 50  | `out/b1.csv`   | raw + Q-corrected mean/σ, ulp=0 centrals, floor check     |
| `b3`    | 25 × 5 = 125 | `out/b3.csv`   | σ-vs-N growth exponent, raw vs Q-corrected                |
| `b4`    | 12 × 5 = 60  | `out/b4.csv`   | raw offset tracks endpoint step; corrected plateaus       |
