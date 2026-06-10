# REBOUND parity — Plummer cluster (N = 10³, softened kernel)

Implementation of the parity protocol declared in
[`paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md`](../../../paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md).
Constants in this directory's scripts mirror the notebook in lockstep — changes here are protocol changes.

## Files

| File               | Purpose                                                                                          |
| ------------------ | ------------------------------------------------------------------------------------------------ |
| `generate_ics.py`  | Seeded Plummer sampler (Aarseth, Hénon & Wielen 1974); centring pass; per-realisation rescale to E = −1/4; writes the committed IC CSVs |
| `ics_n256.csv`     | Committed phase-0 (pilot) initial conditions — the canonical artefact; the generator is provenance |
| `ics_n1000.csv`    | Committed phase-1 (gated) initial conditions                                                     |
| `smoke_pair.py`    | Single-pair softening-convention check: both sides vs the closed form, gate 1e-9                 |
| `rebound_side.py`  | REBOUND IAS15 with `sim.softening = ε`, sampled at apsis's actual times                          |
| `compare.py`       | Loads both long-format CSVs, computes all invariants (softened PE) itself; embedded self-tests   |
| `run.py`           | Orchestrator: smoke test → registration gate → apsis side → REBOUND side → comparator; reports wall-time per step |
| `out/`             | Generated artefacts. Git-ignored.                                                                |

The apsis side lives at
[`crates/apsis/examples/rebound_parity_plummer_cluster.rs`](../../../crates/apsis/examples/rebound_parity_plummer_cluster.rs);
the registration-only Exactness assertion at
[`crates/apsis-1pn/tests/plummer_cluster_registration_gate.rs`](../../../crates/apsis-1pn/tests/plummer_cluster_registration_gate.rs).

## Quick start (phase 0, N = 256)

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python run.py --n 256
```

Phase 0 is informational: the comparator reports every metric and both candidate
L/P floor models but never fails the run. The gated phase-1 run (`--n 1000`)
asserts the floors frozen at phase-0 close (protocol notebook, §Phase 0 results
and gate freeze). To keep per-phase artefacts side by side, pass
`--output-dir out-n256` / `--output-dir out-n1000` (the default `out/` is
overwritten by each run).

## Single source of truth

ε is computed once by the generator (`0.98 · N^(−0.26)` Plummer scale lengths;
Athanassoula et al. 2000) and embedded in the IC header (`# eps=…`). The apsis
example, the REBOUND side, and the comparator all parse it from there — no
second copy exists.
