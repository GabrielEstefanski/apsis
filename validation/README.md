# Validation

External validation of `apsis` results against reference implementations and analytic solutions. This directory holds reproducible cross-implementation comparisons; each subdirectory targets one reference tool, and each scenario within is a self-contained, runnable harness.

## Layout

| Path                | Purpose                                                  |
| ------------------- | -------------------------------------------------------- |
| `rebound-parity/`   | Cross-implementation comparison against REBOUND          |

Future siblings (e.g. `analytic-kepler/`, `dashboard/`) live at this level.

## Conventions

- **Each scenario directory** contains: a Python orchestrator (`run.py`), individual scripts per side (e.g. `rebound_side.py`), a comparator (`compare.py`), pinned dependencies (`requirements.txt`), and an `out/` directory (git-ignored) for generated artefacts.
- **Scenario protocols** — initial conditions, integrator settings, metrics, and tolerances declared *a priori* — live as lab notebooks under [`../docs/experiments/`](../docs/experiments/) and are referenced from each scenario's README.
- **Constants in scripts** must mirror the protocol notebook in lockstep. A change here is a protocol change; the notebook updates with it.
- **`apsis` side** of any comparison is implemented as a Cargo example under [`../crates/apsis/examples/`](../crates/apsis/examples/) so it builds with the rest of the workspace and is exercised by CI.

## Reproducibility

Each scenario is intended to be runnable on a clean checkout with a Rust toolchain (≥ 1.85) plus a Python environment satisfying that scenario's `requirements.txt`. Outputs land in the local `out/` and are not committed; the protocol notebook records the run-time hashes (apsis commit, REBOUND version, OS/CPU, FMA flags) at experimental run time.
