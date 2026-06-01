# Top-level reviewer-facing entry points for the paper pipeline.
#
# Figures are generated from frozen data snapshots under
# `paper/figures/data/`. `make figures` re-renders the committed PDFs
# from those snapshots and is deterministic. To refresh a data
# snapshot, run the corresponding validation harness and copy its
# output into `paper/figures/data/`.

.PHONY: help figures paper validation clean data lint-paper
.PHONY: figures-mercury-1pn figures-rebound-parity-trajectories figures-rebound-parity-brouwer
.PHONY: figures-plummer-apsidal-convergence data-plummer-convergence
.PHONY: validation-mercury-1pn validation-recommended-dt

FIGURES_DIR  := paper/figures
SCRIPTS_DIR  := $(FIGURES_DIR)/scripts
DATA_DIR     := $(FIGURES_DIR)/data
NOTEBOOK_DIR := paper/notebooks/scripts

PAPER_FIGURES := \
	$(FIGURES_DIR)/mercury_1pn_long_horizon.pdf \
	$(FIGURES_DIR)/rebound_parity_trajectories.pdf \
	$(FIGURES_DIR)/rebound_parity_brouwer.pdf \
	$(FIGURES_DIR)/plummer_apsidal_convergence.pdf

help:
	@echo "Reviewer targets:"
	@echo "  make figures               Re-render all paper figures from frozen data"
	@echo "  make paper                 Compile paper.md to paper.pdf (requires pandoc)"
	@echo "  make validation            Run every validation harness"
	@echo "  make clean                 Remove generated figures, paper.pdf, harness outputs"
	@echo
	@echo "  make data                  Regenerate figure data snapshots (slow; runs harnesses)"
	@echo "  make lint-paper            Check paper.md for hyphen-wrap render artifacts"
	@echo
	@echo "Per-figure / per-harness:"
	@echo "  make figures-mercury-1pn"
	@echo "  make figures-rebound-parity-trajectories"
	@echo "  make figures-rebound-parity-brouwer"
	@echo "  make figures-plummer-apsidal-convergence"
	@echo "  make data-plummer-convergence"
	@echo "  make validation-mercury-1pn"
	@echo "  make validation-recommended-dt"

# ── Figures ──────────────────────────────────────────────────────────── #

figures: $(PAPER_FIGURES)

$(FIGURES_DIR)/mercury_1pn_long_horizon.pdf: \
		$(SCRIPTS_DIR)/mercury_1pn_long_horizon.py \
		$(DATA_DIR)/mercury_1pn_long_horizon_ias15.csv
	python $(SCRIPTS_DIR)/mercury_1pn_long_horizon.py

figures-mercury-1pn: $(FIGURES_DIR)/mercury_1pn_long_horizon.pdf

REBOUND_PARITY_DATA := \
	$(DATA_DIR)/rebound_parity_kepler_apsis.csv \
	$(DATA_DIR)/rebound_parity_kepler_rebound.csv \
	$(DATA_DIR)/rebound_parity_figure8_apsis.csv \
	$(DATA_DIR)/rebound_parity_figure8_rebound.csv \
	$(DATA_DIR)/rebound_parity_pythagorean_apsis.csv \
	$(DATA_DIR)/rebound_parity_pythagorean_rebound.csv \
	$(DATA_DIR)/rebound_parity_retrograde_apsis.csv \
	$(DATA_DIR)/rebound_parity_retrograde_rebound.csv

$(FIGURES_DIR)/rebound_parity_trajectories.pdf: \
		$(SCRIPTS_DIR)/rebound_parity_trajectories.py \
		$(REBOUND_PARITY_DATA)
	python $(SCRIPTS_DIR)/rebound_parity_trajectories.py

$(FIGURES_DIR)/rebound_parity_brouwer.pdf: \
		$(SCRIPTS_DIR)/rebound_parity_brouwer.py \
		$(DATA_DIR)/rebound_parity_retrograde_apsis.csv \
		$(DATA_DIR)/rebound_parity_retrograde_rebound.csv
	python $(SCRIPTS_DIR)/rebound_parity_brouwer.py

figures-rebound-parity-trajectories: $(FIGURES_DIR)/rebound_parity_trajectories.pdf
figures-rebound-parity-brouwer: $(FIGURES_DIR)/rebound_parity_brouwer.pdf

$(FIGURES_DIR)/plummer_apsidal_convergence.pdf: \
		$(SCRIPTS_DIR)/plummer_apsidal_convergence.py \
		$(DATA_DIR)/plummer_apsidal_convergence.csv
	python $(SCRIPTS_DIR)/plummer_apsidal_convergence.py

figures-plummer-apsidal-convergence: $(FIGURES_DIR)/plummer_apsidal_convergence.pdf

# ── Data snapshots ───────────────────────────────────────────────────── #
# Slow regeneration of figure data (the figures themselves render from the
# frozen CSVs above and stay deterministic). Refresh a snapshot only when the
# underlying physics or protocol changes, then commit the updated CSV.

data: data-plummer-convergence

# apsis geometric sweep (Rust) -> apsis CSV; quadrature oracle + closed form
# (numpy) merge -> the convergence CSV the figure reads.
data-plummer-convergence:
	cargo run --release --example softened_plummer_sweep -p apsis -- \
		--output $(DATA_DIR)/apsis_softened_sweep.csv
	python $(NOTEBOOK_DIR)/plummer_apsidal_convergence_data.py

# ── Paper ────────────────────────────────────────────────────────────── #

paper: paper.pdf

paper.pdf: paper.md paper.bib $(PAPER_FIGURES)
	pandoc paper.md \
		--bibliography=paper.bib \
		--citeproc \
		--number-sections \
		--highlight-style=tango \
		--pdf-engine=xelatex \
		-o paper.pdf

# Guard against the hard-wrap-after-hyphen artifact: a markdown soft line
# break renders as a space, so a line ending in a word hyphen (`bit-` with
# `identically` on the next line) prints as "bit- identically". Fails if any
# paper.md line ends in a hyphen other than the YAML `---` fences.
lint-paper:
	@if grep -nP -- '-$$' paper.md | grep -vP ':---+$$'; then echo "ERROR: paper.md has a line ending in a hyphen (renders as 'word- word'); join the compound onto one line." && exit 1; else echo "lint-paper OK: no trailing-hyphen wraps in paper.md"; fi

# ── Validation ───────────────────────────────────────────────────────── #

validation: validation-mercury-1pn validation-recommended-dt
	@echo "REBOUND parity scenarios run individually: validation/rebound-parity/<scenario>/run.py"
	@echo "Cross-platform workflow: validation/cross-platform/README.md"

validation-mercury-1pn:
	cd validation/mercury-1pn-long-horizon && python run.py

validation-recommended-dt:
	cargo run --release --example recommended_dt_validation -p apsis --quiet
	cargo run --release --example recommended_dt_compare -p apsis --quiet

# ── Clean ────────────────────────────────────────────────────────────── #

clean:
	rm -f $(PAPER_FIGURES) paper.pdf
	rm -rf validation/mercury-1pn-long-horizon/out
	rm -rf validation/recommended-dt/out
