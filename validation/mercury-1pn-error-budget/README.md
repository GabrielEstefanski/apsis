# Mercury 1PN error budget — Phase B ensemble

Implementation of the Phase-B protocol declared in
[`paper/notebooks/2026-06-10-mercury-1pn-error-budget.md`](../../paper/notebooks/2026-06-10-mercury-1pn-error-budget.md).

## Files

| File           | Purpose                                                                                               |
| -------------- | ----------------------------------------------------------------------------------------------------- |
| `ensemble.py`  | Orchestrator: runs the Rust example across ULP/N/constructor grids; phases `smoke`, `b1`, `b3`       |
| `out/`         | Generated artefacts (`b1.csv`, `b3.csv`). Git-ignored.                                               |

The Rust measurement binary lives at
[`crates/apsis-1pn/examples/error_budget_run.rs`](../../crates/apsis-1pn/examples/error_budget_run.rs).
Constants (`A`, `E`, `M_MERCURY`, IC construction, osculating-ω measurement) mirror
`crates/apsis-1pn/tests/mercury_precession_gate.rs` exactly.

## Quick start

```powershell
# From the workspace root:
python validation/mercury-1pn-error-budget/ensemble.py --phase smoke
```

Requires Python 3.9+ (stdlib only — no pip install needed).

## Phases

| Phase   | Runs          | Output         | Gates                                          |
| ------- | ------------- | -------------- | ---------------------------------------------- |
| `smoke` | 3 × 2 = 6    | stdout table   | raw\_c ulp=0 rel\_err < 1e-4                  |
| `b1`    | 25 × 2 = 50  | `out/b1.csv`   | mean, sigma\_omega, B5 central per constructor |
| `b3`    | 25 × 5 = 125 | `out/b3.csv`   | fitted alpha with std error (H3)               |
