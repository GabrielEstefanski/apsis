# Validation

External validation of `apsis` results against reference implementations, analytic solutions, and cross-host reproducibility. Each subdirectory is a self-contained, runnable harness.

## Layout

| Path                          | Purpose                                                                          |
| ----------------------------- | -------------------------------------------------------------------------------- |
| `rebound-parity/`             | Cross-implementation comparison against REBOUND on canonical scenarios           |
| `mercury-1pn-long-horizon/`   | Mercury 1PN integrated over 1000 years (~4150 orbits) vs the closed-form GR rate |
| `recommended-dt/`             | Step-size heuristic validation for fixed-step integrators (VV, Y4, WH)           |
| `cross-platform/`             | Bit-identical trajectory reproduction across Windows and Linux on x86_64         |

## Conventions

- **Scenario directories** (under `rebound-parity/` and `mercury-1pn-long-horizon/`) follow a common shape: a Python orchestrator (`run.py`), individual scripts per side (e.g. `rebound_side.py`), a comparator (`compare.py`), pinned dependencies (`requirements.txt`), and an `out/` directory (git-ignored) for generated artefacts. `cross-platform/` uses a different shape — see its own README.
- **Scenario protocols** — initial conditions, integrator settings, metrics, and tolerances declared *a priori* — live as lab notebooks under [`../docs/experiments/`](../docs/experiments/) or [`../paper/notebooks/`](../paper/notebooks/) and are referenced from each scenario's README.
- **Constants in scripts** mirror the protocol notebook in lockstep. A change here is a protocol change; the notebook updates with it.
- **`apsis` side** of any comparison is implemented as a Cargo example under [`../crates/apsis/examples/`](../crates/apsis/examples/) or [`../crates/apsis-1pn/examples/`](../crates/apsis-1pn/examples/) so it builds with the rest of the workspace and is exercised by CI.

## Reproducibility

Each scenario runs on a clean checkout with the Rust toolchain pinned by `rust-toolchain.toml` plus a Python environment satisfying that scenario's `requirements.txt`. Outputs land in the local `out/` and are not committed. The cross-platform scenario commits its captured outputs (`cross-platform/windows/` and `cross-platform/linux/`) because the bit-equal claim is verifiable from the captured files.
