# Contributing to apsis

Thanks for your interest in apsis. This document describes the workflow,
conventions, and quality gates used in this repository.

## Quick start

```bash
git clone https://github.com/GabrielEstefanski/apsis
cd apsis
cargo build --workspace
cargo nextest run --workspace
```

Python distribution (requires `maturin`):

```bash
python -m venv .venv
. .venv/bin/activate    # or .venv\Scripts\activate on Windows
pip install maturin
maturin develop --release
pytest tests/
```

## Branch flow

- `develop` is the integration branch. Open feature PRs against `develop`.
- `master` receives release / hotfix only.
- Branch naming: `feat/<short-name>`, `fix/<short-name>`, `chore/<short-name>`,
  `refactor/<short-name>`, `docs/<short-name>`.

## Pre-commit checks

Before every commit on Rust changes, run the trio:

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo nextest run --workspace
```

CI enforces these. Skipping locally just costs an extra round-trip.

## Commit messages

Conventional Commits:

```text
type(scope): one-line summary

Optional body explaining motivation, trade-offs, and what changed.
Reference issues and PRs by number.
```

Types in use: `feat`, `fix`, `refactor`, `chore`, `docs`, `ci`, `perf`.
Scopes follow workspace structure (`python`, `kernel`, `examples`,
`apsis-1pn`, `apsis-radiation`, `apsis-central`, `ci`, …).

## Pull requests

PR descriptions should include:

- **Summary**: what changed and why
- **Test plan**: checkboxes for verification steps (`cargo test`, CI gates,
  manual reproduction)
- **References**: linked issues, ADRs that apply, relevant lab notebooks

Fix-and-merge cycles are normal; review feedback should be addressed in
follow-up commits on the same branch, squashed before merge if the
intermediate history is noise.

## Architectural decisions

Significant architectural changes ship with an ADR in
[`docs/adr/`](docs/adr/). The format is short: context, decision,
consequences, alternatives considered. Existing ADRs document the
federated operator model, kernel-as-system-parameter, citation
provenance, and the consolidated Python distribution.

When proposing a refactor that touches the public extension API or the
operator contract, draft the ADR alongside the implementing PR.

## Lab notebooks

Numerical experiments (validation portfolios, performance benchmarks,
new integrator characterisations) ship with a notebook entry under
[`docs/experiments/`](docs/experiments/). The convention is

1. Protocol declared *a priori* (initial conditions, integrator,
   metrics, tolerances, expected result)
2. Run executed against the protocol
3. Post-mortem analysis added below

Notebooks freeze after the run. Don't edit them retroactively;
file a follow-up notebook if the analysis changes.

## Testing

- Unit tests live next to the code (`mod tests` in `.rs` files).
- Integration tests live in `tests/` directories per crate.
- Cross-implementation parity tests against REBOUND live in
  `crates/apsis/examples/rebound_parity_*.rs` and run as `cargo
  run --release --example rebound_parity_<name>`.
- Python smoke tests live in `tests/` at the repo root and run via
  `pytest tests/` after `maturin develop`.

Tests assert behaviour and physics, not implementation details.
Refactor-survivability is the criterion: a test that breaks when a
private buffer changes shape but the integrator still produces correct
trajectories is over-coupled to implementation.

## Filing issues

Use the GitHub issue tracker. Templates aren't required; useful
ingredients in any issue:

- What you expected and what you observed
- Minimal reproduction (commit hash, OS, integrator, scenario)
- Impact (gate / regression / paper-relevant)

For numerical claims (energy drift, orbital element bounds), include
the actual measurement with units and the comparison reference.

## Licensing

apsis is licensed under the Apache License, Version 2.0 (see
[LICENSE](LICENSE)). Contributions are accepted under the same license
without separate Contributor License Agreement.

## Code of conduct

This project follows the [Contributor Covenant Code of
Conduct](CODE_OF_CONDUCT.md). By participating you agree to abide by
its terms.

## Questions

Open a GitHub Discussion or an issue. For security-sensitive reports,
email `gabrielbragaestefanski@gmail.com` directly rather than filing
publicly.
