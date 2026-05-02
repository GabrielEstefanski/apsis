# REBOUND parity — Retrograde Kepler e=0.5

Implementation of the parity protocol declared in [`docs/experiments/2026-05-01-rebound-parity-retrograde.md`](../../../docs/experiments/2026-05-01-rebound-parity-retrograde.md). Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

The single IC difference vs the Kepler-prograde scenario is a sign-flip on the secondary's tangential periapsis velocity ($+v_\text{peri} \to -v_\text{peri}$). Every other component, mass, length, and physical scale is held identical. This isolates sign-convention coverage as the only experimental variable; the experiment closes the $L_z < 0$ gap left by the prograde / figure-8 / Pythagorean tests, which all sit at $L_z \ge 0$.

## Files

| File                | Purpose                                                                                                                       |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `rebound_side.py`   | Runs REBOUND IAS15 with the protocol ICs, exports `out/rebound.csv` matching the apsis-side schema                            |
| `compare.py`        | Loads `out/apsis.csv` and `out/rebound.csv`, computes Tier 1 + Tier 2 + Tier 3 metrics at both horizons, emits decision rule  |
| `run.py`            | Orchestrator: runs the apsis side via `cargo run`, then `rebound_side.py`, then `compare.py`                                  |
| `requirements.txt`  | Pinned Python dependencies                                                                                                    |
| `out/`              | Generated artefacts (CSVs, plots). Git-ignored.                                                                               |

The `apsis` side lives at [`crates/apsis/examples/rebound_parity_retrograde.rs`](../../../crates/apsis/examples/rebound_parity_retrograde.rs).

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

Wall-time estimate: $\sim 10$–$30$ seconds total for the apsis side (10000 orbits at IAS15) plus a similar duration for the REBOUND side, dominated by IAS15's adaptive substep selection. CSV output is $\sim 1$ MB per side.

## Two-horizon design

The orchestrator runs a single $10^4$-orbit integration on each side and the comparator evaluates **two horizons** from the same data:

- **Long-horizon gate** ($10^4$ orbits) — primary verdict; closes the long-horizon stability gate identified during the GR-readiness review for the federation thesis.
- **Short-horizon checkpoint** (100 orbits) — direct comparability with Kepler-prograde at the matched horizon.

Both horizons must pass for the experiment to verdict `pass`. The comparator emits a per-horizon decision-rule outcome (PASS / TIER1-FAIL / TIER2-FAIL / DEEP-FAIL / BROUWER-SATURATION) so post-run interpretation does not depend on prose.

## Tier structure (from the notebook §Hypothesis)

**Tier 1 — magnitude invariants (gated, 7 metrics)** — identical to Kepler-prograde:

- $\lvert \Delta a \rvert / a$
- $\lvert \Delta e \rvert$
- $\lvert \Delta \omega \rvert$ (rad)
- $\bigl\| \lvert h \rvert - \lvert h_0 \rvert \bigr\| / \lvert h_0 \rvert$ (cross-impl, magnitude only)
- $\lvert \Delta E / E_0 \rvert$ apsis
- $\lvert \Delta E / E_0 \rvert$ rebound
- cross-impl $\lvert \Delta E \rvert / \lvert E_0 \rvert$

**Tier 2 — sign(h) consistency (gated, 3 binary checks)** — retrograde-specific:

- apsis: $\mathrm{sign}(h(t)) = \mathrm{sign}(h_0)$ at every sample $\wedge$ $\lvert h(t) \rvert > \varepsilon_\text{floor}$ at every sample
- rebound: same on the REBOUND side
- cross-impl: $\mathrm{sign}(h_\text{apsis}(t)) = \mathrm{sign}(h_\text{rebound}(t))$ at every sample

$\varepsilon_\text{floor} = 10^{-10}$ — defensive guard quantified in the notebook §Hypothesis.

**Tier 3 — geometric coherence (informational, NOT gated)**:

- $\max \lvert \Delta r \rvert$ — phase-drift contaminated, expected to saturate at $O(1)$ at the long horizon.

## Output

| Path                 | Content                                                                                                                                  |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `out/apsis.csv`      | Wide-format CSV with 10001 rows: orbit, time, body 0/1 state, total energy. Produced by the apsis Cargo example.                          |
| `out/rebound.csv`    | Same schema, produced by `rebound_side.py`.                                                                                              |
| `out/comparison.json` | Per-horizon structured report: Tier 1 + Tier 2 + Tier 3 + decision-rule verdict.                                                          |

`out/` is git-ignored. Final paper-anchoring numbers land in the protocol notebook (committed) or via Zenodo (separate archive).
