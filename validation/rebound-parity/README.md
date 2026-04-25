# REBOUND parity

Cross-implementation parity between `apsis` and [REBOUND](https://github.com/hannorein/rebound) on canonical scenarios.

The framing is **validation** — establishing that `apsis`'s numerical foundation produces trajectories consistent with the reference implementation that the literature treats as standard. The framing is *not* competitive comparison: agreement within a documented tolerance is the goal, not numerical superiority. See [`memory/project_paper_positioning.md`](https://github.com/gabrielbragaestefanski/apsis) for the principle.

## Scenarios

| Scenario             | Status                | Protocol notebook                                                                                                       |
| -------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `kepler/`            | implemented           | [`docs/experiments/2026-04-25-rebound-parity-kepler.md`](../../docs/experiments/2026-04-25-rebound-parity-kepler.md)    |
| `pythagorean/`       | planned (Phase 6A)    | TBD                                                                                                                     |
| `retrograde/`        | planned (Phase 6A)    | TBD                                                                                                                     |

## Running a scenario

Each scenario has a Python orchestrator (`run.py`) that runs the `apsis` side (via `cargo run --release --example`), the REBOUND side (via Python), and the comparator. Run from the scenario directory:

```text
cd validation/rebound-parity/<scenario>
python -m venv .venv
# Windows:
.venv\Scripts\activate
# Linux/macOS:
source .venv/bin/activate
pip install -r requirements.txt
python run.py
```

The orchestrator exits non-zero if any metric falls outside the tolerance declared *a priori* in the protocol notebook.

## Adding a new scenario

1. Open the protocol notebook in `docs/experiments/YYYY-MM-DD-rebound-parity-<scenario>.md` with ICs, integrator settings, metrics, and tolerances declared *a priori*.
2. Implement the `apsis` side as a Cargo example in `crates/apsis/examples/rebound_parity_<scenario>.rs`.
3. Mirror the directory structure of `kepler/` for the Python side.
4. Update this README's scenarios table.
