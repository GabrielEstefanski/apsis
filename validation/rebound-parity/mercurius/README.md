# REBOUND parity — Mercurius

Implementation of the parity protocol declared in [`docs/experiments/2026-05-13-rebound-parity-mercurius.md`](../../../docs/experiments/2026-05-13-rebound-parity-mercurius.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                                                  |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `rebound_side.py`   | Runs REBOUND MERCURIUS with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                                   |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes the seven gated metrics from §Hypothesis, exits 0 iff within tolerance             |
| `run.py`            | Orchestrator: runs the apsis side via `cargo run`, then `rebound_side.py`, then `compare.py`                                             |
| `requirements.txt`  | Pinned Python dependencies                                                                                                                |
| `out/`              | Generated artefacts (CSVs, JSON report). Git-ignored.                                                                                     |

The `apsis` side lives at [`crates/apsis/examples/rebound_parity_mercurius.rs`](../../../crates/apsis/examples/rebound_parity_mercurius.rs).

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

The orchestrator exits non-zero if any of the seven gated metrics defined in the protocol notebook's §Hypothesis exceeds its *a priori* tolerance.

### Tier 1 — conservation parity

- ΔE/E₀ peak per side (apsis, REBOUND) — 2nd-order method floor (≤ 10⁻⁸ each)
- cross-impl ΔE/E₀ peak — independent-implementations agreement (≤ 5×10⁻⁹)
- cross-impl ΔLz/Lz₀ peak — angular momentum agreement (≤ 10⁻¹⁰)

### Tier 2 — test-particle orbital element parity (end of run)

- Δa/a — semi-major axis (≤ 10⁻⁵)
- Δe/e — eccentricity (≤ 10⁻⁵)
- Δi/i — inclination (≤ 10⁻⁵)

Point-wise max |Δr| on the test particle is preserved in the output as **informational context only**, not as a gate. See §Methodology in the notebook for why phase drift through the IAS15-driven encounter step makes |Δr| an inadequate parity criterion at the end of a 10⁴-year run.

## Output

`out/comparison.json` carries the structured per-metric verdict (observed value, tolerance, pass/fail flag, optional context). The console output prints the same data in tabular form. A non-zero exit code from `run.py` indicates at least one gated metric failed; the JSON report enumerates which.

## Why solar AU-year units (G ≈ 4π²)

Mercurius is a planetary-scenario integrator; running it in canonical units (G = 1) would obscure the connection to physical timescales (Jupiter orbital period, encounter durations) that the protocol's α = 3 default and 10⁴-yr horizon are calibrated against. The Python side computes G from the same SI constants the apsis side uses (`G_SI · M_sun_kg · yr_s² / AU_m³`) so the two sides land on bit-equivalent G values.
