# REBOUND parity — Pythagorean three-body (Burrau 1913)

Implementation of the parity protocol declared in [`docs/experiments/2026-04-30-rebound-parity-pythagorean.md`](../../../docs/experiments/2026-04-30-rebound-parity-pythagorean.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File                | Purpose                                                                                                                   |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `rebound_side.py`   | Runs REBOUND IAS15 with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                         |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes the twelve gated metrics from §Hypothesis (Tier 1 + Tier 2), exits 0 iff within tolerance |
| `run.py`            | Orchestrator: runs the apsis side via `cargo run`, then `rebound_side.py`, then `compare.py`                               |
| `requirements.txt`  | Pinned Python dependencies                                                                                                 |
| `out/`              | Generated artefacts (CSVs, JSON report). Git-ignored.                                                                      |

The `apsis` side lives at [`crates/apsis/examples/rebound_parity_pythagorean.rs`](../../../crates/apsis/examples/rebound_parity_pythagorean.rs).

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

For an informational longer-horizon stress run (not part of the gated parity claim — see protocol notebook §Out of scope):

```text
python run.py --horizon 200 --output-dir ./out-200tu
```

## Gated metrics (Tier 1 + Tier 2)

The orchestrator exits non-zero if any of the twelve gated metrics defined in the protocol notebook's §Hypothesis exceeds its *a priori* tolerance.

**Tier 1 — hard physical invariants** (energy and angular momentum; the load-bearing parity statement):

- $\lvert \Delta E / E_0 \rvert$ per side (apsis, REBOUND) — energy conservation per implementation
- cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$ — energy agreement between implementations
- $\lvert \Delta \mathbf{L} \rvert$ per side (apsis, REBOUND) — angular-momentum conservation, full vector norm, absolute units ($\mathbf{L}_0 = \mathbf{0}$ by IC)
- cross-impl $\lvert \Delta \mathbf{L} \rvert$ — angular-momentum agreement, absolute

**Tier 2 — construction-level sanity** (linear momentum and centre of mass; weak gates that catch book-keeping bugs):

- $\lvert \Delta \mathbf{P} \rvert$ per side, cross-impl $\lvert \Delta \mathbf{P} \rvert$ — linear-momentum conservation and agreement, absolute ($\mathbf{P}_0 = \mathbf{0}$ by IC)
- $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ per side, cross-impl $\lvert \Delta \mathbf{r}_\text{COM} \rvert$ — centre-of-mass conservation and agreement, absolute ($\mathbf{r}_\text{COM}(0) = \mathbf{0}$ by IC, verified algebraically in the protocol notebook §Methodology)

**Tier 3 — geometric coherence** is reported per-sample, never aggregated into a pass/fail criterion:

- $\max \lvert \mathbf{r}_\text{apsis}(t) - \mathbf{r}_\text{rebound}(t) \rvert$ per body — phase-drift contaminated and Lyapunov-amplified for the Pythagorean dynamics; expected to reach $O(1)$ before the horizon. See the protocol notebook §"Why this metric set, not $\lvert \Delta r \rvert$" — and §Verdict criterion for the explicit statement that no Tier-3 quantity ever fails the experiment.

## Output

| Path                  | Content                                                                                                                |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `out/apsis.csv`       | Wide-format CSV with 2101 rows (70 t.u. baseline at 30 samples/t.u.): sample index, time, three-body state, total energy. Produced by the apsis Cargo example. |
| `out/rebound.csv`     | Same schema, produced by `rebound_side.py`.                                                                            |
| `out/comparison.json` | Structured report: Tier-1, Tier-2, Tier-3 metrics with observed values, tolerances, and pass/fail per gated metric.    |

`out/` is git-ignored. Final paper-anchoring numbers land in the protocol notebook (committed) or via Zenodo (separate archive).

## Sampling cadence — analysis vs report

Two cadences are decoupled by design:

- **Analysis cadence (dense)**: 30 samples/t.u. $\times$ 70 t.u. $=$ 2101 samples. Both sides emit at this resolution; the comparator's $\max(\cdot)$ aggregates run over the dense set. The gate is evaluated against the worst-case sample over the full horizon.
- **Report cadence (sparse)**: 4 samples/t.u. (every $0.25$ t.u.) $\to$ 281 representative rows published in the protocol notebook's §Results evolution table. Sparse cadence is for human readability and does not enter gating.

The CSV always contains the dense set; the report cadence is a slice taken at notebook-writing time. Tier-1 and Tier-2 metrics are integration-level invariants and are insensitive to analysis cadence above a density threshold; Tier-3 $\lvert \Delta \mathbf{r} \rvert$ is cadence-dependent in principle but is informational only — see protocol notebook §Sample-density sensitivity for the full caveat.

## Cross-reference — figure-8 and Kepler scenarios

This scenario shares the conceptual framework — a-priori protocol notebook $\to$ cargo example $\to$ REBOUND side $\to$ comparator with structured JSON output — with the Kepler ([`../kepler/`](../kepler/)) and figure-8 ([`../figure8/`](../figure8/)) scenarios, specialised to the chaotic three-body regime. Differences:

| Aspect | Kepler | Figure-8 | Pythagorean |
| --- | --- | --- | --- |
| Regime | periodic 2-body | periodic 3-body | chaotic 3-body |
| Gated invariants | orbital elements + energy | global integrals $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ | global integrals $E$, $\mathbf{L}$, $\mathbf{P}$, $\mathbf{r}_\text{COM}$ |
| Per-side check | conservation against analytic Kepler invariants | conservation against IC values | conservation against IC values |
| Metric tier hierarchy | flat (seven gated) | three tiers (hard / sanity / informational) | three tiers (hard / sanity / informational) |
| Horizon | $100\,T$ | $10\,T$ (gated) + $50\,T$ (informational) | $70$ canonical t.u. (gated; no period concept) |
| Tier-3 expected magnitude | $\sim 10^{-12}$ (figure-8 post-fix) | $\sim 10^{-12}$ | $O(1)$ — Lyapunov-amplified |

The shared framework is "physical invariants gate; geometric coherence informs". The Pythagorean specialisation makes the geometric-coherence framing load-bearing: in a chaotic regime, the only honest cross-implementation comparison is on the conserved quantities.
