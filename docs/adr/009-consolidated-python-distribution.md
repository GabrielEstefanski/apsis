# ADR-009 — Consolidated Python Distribution

**Status:** Accepted
**Date:** 2026-05-15
**Supersedes (in part):** ADR-005 §"Per-operator Python wrappers"

---

## Context

Pre-consolidation, the Python surface was distributed across three
crates that produced two PyPI packages:

- `apsis-py` (cdylib) → `pip install apsis`
- `apsis-py-core` (rlib) — capsule transport, linked by every cdylib
- `apsis-1pn-py` (cdylib) → `pip install apsis-1pn`

The intent (ADR-005) was that every operator crate would ship its own
`-py` cdylib alongside it, giving each force a separately-citable
Python package. After the first concrete consumer (`apsis-1pn`) and
two more operator crates (`apsis-radiation`, `apsis-central`) the
cost-benefit ratio of that approach inverted:

- Each `-py` crate carried ~150 LOC of duplicated boilerplate
  (Cargo.toml, lib.rs, __init__.py, _native.pyi, README, CI config,
  maturin config). Two crates was tolerable; five would not have been.
- PyO3 cross-cdylib `#[pyclass]` identity is fragile, so each `-py`
  cdylib needed its own `apsis.Perturbation` workaround via the
  pure-Python wrapper class. The workaround was load-bearing for
  external plugins but offered no benefit to internally-shipped
  operators.
- Researcher UX was worse than REBOUND/REBOUNDx: REBOUND ships one
  Python package containing every effect, with submodule structure
  (`reboundx.<effect>`) preserving conceptual separation. The
  apsis approach forced `pip install apsis apsis-1pn` and two
  `import` lines.

## Decision

Consolidate the Python surface into a single `apsis` distribution
backed by one cdylib (`apsis-python`), with each internal operator
exposed as a submodule of the `apsis` package:

```python
import apsis
from apsis.gr import PostNewtonian1PN
```

Repository layout:

```text
apsis/                       # Python package source (root level)
  __init__.py
  gr/
    __init__.py
pyproject.toml               # distribution metadata + maturin config
crates/
  apsis/                     # core (unchanged name)
  apsis-1pn/                 # operator crate (Rust only)
  apsis-radiation/           # operator crate (Rust only)
  apsis-central/             # operator crate (Rust only)
  apsis-py-core/             # capsule transport + extractors (rlib)
  apsis-python/              # PyO3 cdylib (manifest-path target)
```

The Python package source lives at the repository root in `apsis/`,
not inside any Cargo crate. Maturin reaches into `crates/apsis-python`
via the root `pyproject.toml`'s `manifest-path` setting, decoupling
the Python distribution layout from the Rust workspace layout.

Internal operators are gated behind Cargo features in
`apsis-python/Cargo.toml` (`gr` enables `apsis-1pn`, future operators
add their own features). The published wheel uses `--all-features`;
local dev / CI matrix can build subsets.

## Architectural rule

`apsis-py-core` **never imports the apsis Python package by name**
and **never defines a `#[pyclass]`**. Both rules exist for concrete
reasons:

- Importing the package would couple the binding kit to the canonical
  distribution layout, breaking reuse from any context that does not
  have the apsis Python package installed (a Rust-only downstream, an
  alternative binding, an extracted GUI shell).
- PyO3 cross-cdylib `#[pyclass]` identity is fragile: a class defined
  in an `rlib` and linked into multiple `cdylib`s registers distinct
  Python classes per cdylib, breaking `isinstance` at the boundary.
  The `apsis-python` cdylib owns every `#[pyclass]`; external plugin
  cdylibs do too, for their own classes.

The bright-line review heuristic: anything that touches the apsis
Python package by name (`py.import("apsis")` or `apsis.<attr>`)
belongs in `apsis-python`, never in `apsis-py-core`. Anything that
defines `#[pyclass]` or `#[pymodule]` belongs in a cdylib, never in
an rlib.

## Plugin contract preserved

External `apsis-plugin-X` crates retain the same surface:

- Build a `Box<dyn HamiltonianOperator>`
- Wrap via `apsis_py_core::box_into_capsule`
- Return as `apsis.Perturbation` (constructed in plugin's own Python
  glue using the wrap template documented in `apsis-py-core`)
- User attaches via `sys.add_hamiltonian_perturbation(plugin.MyForce.for_units(...))`

The wrap-dance template lives in the `apsis-py-core` README;
plugin authors copy it into their own cdylib's source, where the
`"apsis"` import name actually belongs.

## Consequences

**Architectural wins:**

- Plugin author cognitive load drops: one binding kit (`apsis-py-core`),
  one wrap pattern (template), one place to look (the README)
- Internal operators add via Cargo feature flag, not new crate +
  `pip install` + `import` line
- JOSS reviewer reads one Python distribution, one architecture
- `apsis-py-core` becomes testable in isolation (no apsis Python
  package required for its tests)

**Migration cost:**

- `apsis-py` → `apsis-python` (Cargo crate name) — no published wheel
  consumers exist yet (alpha-only); internal workspace deps update in
  one commit
- `apsis-1pn-py` → consolidated into `apsis.gr` submodule
- v0.1.0-alpha.1 release notes carry a forward-pointer; install
  instructions on root README updated to `pip install apsis` only

**Backward compatibility for external plugins:**

- `apsis-py-core::box_into_capsule` / `take_box_from_capsule` API
  unchanged (capsule version still `_v3`)
- `apsis.Perturbation` Python class unchanged
- Plugin authors who wrote against the v0.1.0-alpha.1 surface need
  zero code changes

**Future operators:**

When a new federated operator gets a Python adapter, the steps are:

1. Add a Cargo feature in `apsis-python/Cargo.toml`
2. Add `crates/apsis-python/src/operators/<name>.rs` with the
   `#[pyclass]` and `register(parent)` function
3. Add `apsis/<name>/__init__.py` re-exporting from `apsis._native.<name>`
4. Add `tests/test_<name>_smoke.py`

No new Cargo crate, no new PyPI package, no CI matrix entry.

## Anti-pattern this rules out

A future PR proposing "create `apsis-newforce-py` cdylib for the next
operator" should be rejected. Internal operators consolidate; only
out-of-tree third-party plugins ship as separate cdylibs.

## References

- ADR-005 (federated perturbation operators) — the original per-operator-
  wrapper design, superseded in part by this ADR
- ADR-006 (operator preconditions) — the kernel-requirement contract
  preserved across the consolidation
- ADR-007 (citation provenance) — the per-operator Citation infrastructure
  preserved across the consolidation
- REBOUND/REBOUNDx single-package precedent
