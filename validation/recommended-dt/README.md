# Recommended-dt heuristic

Validation of the `recommended_dt` step-size heuristic for fixed-step integrators (Velocity Verlet, Yoshida 4, Wisdom–Holman) across 13 template scenarios. Heuristic shape, gates, and current verdict are recorded in [`docs/experiments/2026-05-01-recommended-dt-heuristic.md`](../../docs/experiments/2026-05-01-recommended-dt-heuristic.md); constants in the Cargo examples below mirror that note in lockstep.

## Layout

| Path                                                                                                            | Purpose                                                                                                                  |
| --------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| [`crates/apsis/examples/recommended_dt_validation.rs`](../../crates/apsis/examples/recommended_dt_validation.rs) | Runs the 13 templates at `recommended_dt` under VV/Y4/WH; writes per-step energy + angular-momentum to `out/runs.csv`    |
| [`crates/apsis/examples/recommended_dt_compare.rs`](../../crates/apsis/examples/recommended_dt_compare.rs)       | Loads `out/runs.csv`, computes the gated metrics from the protocol §Methodology, emits `out/comparison.json`             |
| `out/`                                                                                                          | Generated artefacts (CSV + JSON report). Git-ignored.                                                                    |

The cargo examples live in `apsis`'s example tree rather than in this directory so they build with the rest of the workspace and are exercised by CI. Default output path is `validation/recommended-dt/out/` relative to the workspace root.

## Quick start

```bash
cargo run --release --example recommended_dt_validation -p apsis
cargo run --release --example recommended_dt_compare    -p apsis
```

The comparator exits non-zero if any of the gated metrics defined in the protocol notebook §Methodology exceeds its *a priori* tolerance.

## Output

`out/runs.csv` holds the per-step energy and angular-momentum trace across all 39 (template × integrator) combinations. `out/comparison.json` carries the structured per-metric verdict (observed value, tolerance, pass/fail flag). Both files are git-ignored; the protocol notebook records the paper-anchoring numbers.
