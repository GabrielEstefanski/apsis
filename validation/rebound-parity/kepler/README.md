# REBOUND parity — Kepler e=0.5

Implementation of the parity protocol declared in [`docs/experiments/2026-04-25-rebound-parity-kepler.md`](../../../docs/experiments/2026-04-25-rebound-parity-kepler.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                                  |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `rebound_side.py`   | Runs REBOUND IAS15 with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                        |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes the seven gated metrics from §Revised Methodology, exits 0 iff within tolerance |
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

The orchestrator exits non-zero if any of the seven gated metrics defined in the protocol notebook's §Revised Methodology exceeds its *a priori* tolerance:

- $\lvert \Delta a \rvert / a$ — fractional drift in semi-major axis
- $\lvert \Delta e \rvert$ — eccentricity drift
- $\lvert \Delta \omega \rvert$ — periapsis-orientation drift (radians)
- $\lvert \Delta h \rvert / h$ — fractional angular-momentum drift
- $\lvert \Delta E / E_0 \rvert$ per side (apsis, REBOUND) — energy conservation per implementation
- cross-implementation $\lvert \Delta E \rvert / \lvert E_0 \rvert$ — energy agreement between implementations

Point-wise $\max \lvert \Delta r \rvert$ is preserved in the output as **informational context only**, not as a gate. See §Pilot Interpretation in the notebook for why phase drift makes $\lvert \Delta r \rvert$ an inadequate parity criterion for adaptive integrators.

## Output

| Path                 | Content                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------- |
| `out/apsis.csv`      | Wide-format CSV with 101 rows: orbit, time, body 0/1 state, total energy. Produced by the apsis Cargo example.   |
| `out/rebound.csv`    | Same schema, produced by `rebound_side.py`.                                                                   |
| `out/comparison.json` | Structured report: gated metrics (observed, tolerance, pass/fail) plus the informational $\lvert \Delta r \rvert$ value. |

`out/` is git-ignored. Final paper-anchoring numbers land in the protocol notebook (committed) or via Zenodo (separate archive).
