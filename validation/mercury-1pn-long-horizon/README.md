# Long-horizon Mercury 1PN

Implementation of the protocol declared in [`docs/experiments/2026-05-13-mercury-1pn-long-horizon.md`](../../docs/experiments/2026-05-13-mercury-1pn-long-horizon.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                              |
| ------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `compare.py`        | Loads `out/ias15.csv` (and optionally `out/mercurius.csv`), unwraps the periapsis-orientation trajectory, gates against the closed-form GR prediction. Pure stdlib, no external dependency. |
| `run.py`            | Orchestrator: cargo-runs the IAS15 (always) and Mercurius (`--include-mercurius`) sides, then `compare.py`. |
| `out/`              | Generated artefacts (CSVs, JSON report). Git-ignored.                                                                |

The `apsis` sides live as Cargo examples under [`crates/apsis-1pn/examples/`](../../crates/apsis-1pn/examples/) (apsis-1pn is the natural home — it depends on apsis and ships the perturbation that the example registers).

## Quick start

```text
python run.py                       # Tier 1 only (IAS15 + apsis-1pn)
python run.py --include-mercurius   # Tier 1 + 2 + 3 (post PR #86)
```

The orchestrator exits non-zero if any of the gated metrics defined in the protocol notebook §Hypothesis exceeds its *a priori* tolerance.

### Tier 1 — IAS15 + apsis-1pn

- `|Δω_measured(end) − Δω_GR(end)| / |Δω_GR(end)|` ≤ 1 × 10⁻⁵ (10 ppm)
- per-orbit linearity `R²` ≥ 0.99999

### Tier 2 — Mercurius + apsis-1pn (post PR #86)

Same bounds as Tier 1, applied independently to the Mercurius side.

### Tier 3 — Cross-integrator parity (post PR #86)

- `|Δω_IAS15(end) − Δω_Mercurius(end)| / |Δω_GR(end)|` ≤ 5 × 10⁻⁵ (50 ppm)

## Output

`out/comparison.json` carries the structured per-metric verdict (observed value, tolerance, pass/fail flag, optional context including the per-orbit GR rate and absolute measured + predicted Δω at the end of the run). The console output prints the same data in tabular form. A non-zero exit code from `run.py` indicates at least one gated metric failed; the JSON report enumerates which.

## Why canonical Hénon units (G = 1)

`apsis-1pn`'s `PostNewtonian1PN::solar_units()` constructor carries `c = 10065.13` AU per (year/2π) — `c_SI · (year_s / 2π) / AU_SI`. The matching unit system has time = year/(2π), which is exactly the canonical Hénon `UnitSystem::canonical()` time unit. Re-using these conventions matches the existing 500-orbit gate (`crates/apsis-1pn/tests/mercury_precession_gate.rs`) and removes "did the units rescale correctly" as a source of error. Mercury orbital elements (`a = 0.387098 AU`, `e = 0.20563`) numerically transcribe one-to-one into the canonical system.
