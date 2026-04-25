# REBOUND parity — Kepler e=0.5

Implementation of the parity protocol declared in [`docs/experiments/2026-04-25-rebound-parity-kepler.md`](../../../docs/experiments/2026-04-25-rebound-parity-kepler.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                                  |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `rebound_side.py`   | Runs REBOUND IAS15 with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                        |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes the three metrics from §Hypothesis, exits 0 iff within tolerance    |
| `run.py`            | Orchestrator: runs the apsis side via `cargo run`, then `rebound_side.py`, then `compare.py`                              |
| `requirements.txt`  | Pinned Python dependencies                                                                                                |
| `out/`              | Generated artefacts (CSVs, plots). Git-ignored.                                                                           |

The `apsis` side lives at [`crates/apsis/examples/rebound_parity_kepler.rs`](../../../crates/apsis/examples/rebound_parity_kepler.rs).

## Quick start

```text
python -m venv .venv
# Windows:
.venv\Scripts\activate
# Linux/macOS:
source .venv/bin/activate
pip install -r requirements.txt
python run.py
```

The orchestrator exits non-zero if any of the three metrics defined in the protocol notebook (`max|Δr|`, `|ΔE/E_0|` per side, cross-implementation energy drift) exceeds the *a priori* tolerance.

## Output

| Path                 | Content                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------- |
| `out/apsis.csv`      | Wide-format CSV with 101 rows: orbit, time, body 0/1 state, total energy. Produced by the apsis Cargo example.   |
| `out/rebound.csv`    | Same schema, produced by `rebound_side.py`.                                                                   |
| `out/comparison.json` | Per-orbit metrics + summary (max values, pass/fail per metric).                                               |
| `out/report.md`      | Human-readable summary that pastes into the protocol notebook's §Results section at experimental-run time.    |

`out/` is git-ignored. Final paper-anchoring numbers land in the protocol notebook (committed) or via Zenodo (separate archive).
