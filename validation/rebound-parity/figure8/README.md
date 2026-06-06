# REBOUND parity — Figure-8 choreography

Implementation of the parity protocol declared in [`paper/notebooks/2026-04-26-rebound-parity-figure8.md`](../../../paper/notebooks/2026-04-26-rebound-parity-figure8.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                                  |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `rebound_side.py`   | Runs REBOUND IAS15 with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                        |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes the ten gated metrics from §Hypothesis (Tier 1 + Tier 2), exits 0 iff within tolerance |
| `run.py`            | Orchestrator: runs the apsis side via `cargo run`, then `rebound_side.py`, then `compare.py`                              |
| `requirements.txt`  | Pinned Python dependencies                                                                                                |
| `out/`              | Generated artefacts (CSVs, plots). Git-ignored.                                                                           |

The `apsis` side lives at [`crates/apsis/examples/rebound_parity_figure8.rs`](../../../crates/apsis/examples/rebound_parity_figure8.rs).

## Quick start

<details>
<summary>POSIX (Linux / macOS)</summary>

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python run.py
```

</details>

<details>
<summary>Windows (PowerShell)</summary>

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python run.py
```

</details>

For the Appendix-A informational $50\,T$ sanity run (not part of the gated parity claim):

```text
python run.py --periods 50 --output-dir ./out-50T
```

## Gated metrics (Tier 1 + Tier 2)

The orchestrator exits non-zero if any of the ten gated metrics defined in the protocol notebook's §Hypothesis exceeds its *a priori* tolerance.

**Tier 1 — hard physical invariants** (energy and angular momentum; the load-bearing parity statement):

- $\lvert \Delta E / E_0 \rvert$ per side (apsis, REBOUND) — energy conservation per implementation
- cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$ — energy agreement between implementations
- $\lvert \Delta \mathbf{L} \rvert$ per side (apsis, REBOUND) — angular-momentum conservation, full vector norm, absolute units ($\mathbf{L}_0 \approx \mathbf{0}$)
- cross-impl $\lvert \Delta \mathbf{L} \rvert$ — angular-momentum agreement, absolute

**Tier 2 — construction-level sanity** (linear momentum and centre of mass; weak gates that catch book-keeping bugs):

- $\lvert \Delta \mathbf{P} \rvert$ per side, cross-impl $\lvert \Delta \mathbf{P} \rvert$ — linear-momentum conservation and agreement, absolute ($\mathbf{P}_0 \approx \mathbf{0}$)
- $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ per side, cross-impl $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ — centre-of-mass conservation and agreement, absolute ($\mathbf{r}_\text{COM}(0) \approx \mathbf{0}$)

**Tier 3 — geometric coherence** is reported as informational only:

- $\max \lvert \mathbf{r}_\text{apsis}(t) - \mathbf{r}_\text{rebound}(t) \rvert$ per body — phase-drift contaminated, not a parity criterion. See the protocol notebook §"Why this metric set, not $\lvert \Delta r \rvert$" and the Kepler notebook §Pilot Interpretation for why point-wise position metrics are inadequate gates for adaptive high-order integrators.

## Output

| Path                  | Content                                                                                                       |
| --------------------- | ------------------------------------------------------------------------------------------------------------- |
| `out/apsis.csv`       | Wide-format CSV with 2001 rows ($10\,T$ baseline at 200 samples/period): sample index, time, three-body state, total energy. Produced by the apsis Cargo example. |
| `out/rebound.csv`     | Same schema, produced by `rebound_side.py`.                                                                   |
| `out/comparison.json` | Structured report: Tier-1, Tier-2, Tier-3 metrics with observed values, tolerances, and pass/fail per gated metric. |

`out/` is git-ignored. Final paper-anchoring numbers land in the protocol notebook (committed) or via Zenodo (separate archive).

## Sampling cadence — analysis vs report

Two cadences are decoupled by design:

- **Analysis cadence (dense)**: 200 samples/period $\times$ 10 periods $=$ 2001 samples. Both sides emit at this resolution; the comparator's $\max(\cdot)$ aggregates run over the dense set. The gate is evaluated against the worst-case sample over the full horizon.
- **Report cadence (sparse)**: 4 samples/period (every $T/4$) $\to$ 41 representative rows published in the protocol notebook's §Results evolution table. Sparse cadence is for human readability and does not enter gating.

The CSV always contains the dense set; the report cadence is a slice taken at notebook-writing time.

## Cross-reference — Kepler scenario

The Kepler parity scenario at [`../kepler/`](../kepler/) follows the same conceptual framework — a-priori protocol notebook $\to$ cargo example $\to$ REBOUND side $\to$ comparator with structured JSON output — specialised to a Keplerian regime. Differences:

| Aspect | Kepler | Figure-8 |
| --- | --- | --- |
| Gated invariants | orbital elements $(a, e, \omega, h)$ + energy | global integrals $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ |
| Per-side check | conservation against analytic Kepler invariants | conservation against IC values |
| Metric tier hierarchy | flat (seven gated) | three tiers (hard / sanity / informational) |
| Horizon | $100\,T$ | $10\,T$ (gated) + $50\,T$ (informational, Appendix A) |
| Sampling | 1/period | 200/period (analysis), 4/period (report) |
