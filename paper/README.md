# paper/

Manuscript-supplementary material for the apsis paper submission.

This directory currently holds:

- `notebooks/` — lab notebooks cited in or supporting the paper's
  validation portfolio (REBOUND parity scenarios, Mercury 1PN long-
  horizon convergence, recommended-dt validation).

Pending under Issue #137 (arxiv-grade reproducibility pipeline):

- `Makefile` — `make paper` / `make figures` / `make data` / `make tables`
- `scripts/` — Python scripts that produce figures and tables
- `data/` — cached simulation outputs (CSV / parquet) feeding the scripts
- `figures/` — manuscript-ready PDF figures
- `tables/` — auto-generated benchmark / parity tables

For the manuscript draft itself, see `paper.md` at the repository
root.

For lab notebooks that document internal investigations (bug
forensics, performance experiments, design spikes, architectural
transitions), see `docs/experiments/`. The split is by audience:
`paper/notebooks/` is reviewer-facing, `docs/experiments/` is
maintainer-facing.
