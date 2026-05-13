# ADR-007 — Citation Provenance Contract

**Status:** Accepted
**Date:** 2026-05-13
**PRs:** [#93](https://github.com/GabrielEstefanski/apsis/pull/93) (contract + Python bindings),
[#94](https://github.com/GabrielEstefanski/apsis/pull/94) (second publisher: `apsis-radiation`)

---

## Context

The federation thesis ([[project_thesis_anchor]]) treats every
perturbation crate as a citable scientific artifact: it has a paper,
a DOI, a versioned Cargo dependency, and (when the build is captured
from a git checkout) a commit hash. A paper that cites apsis with
1PN, radiation pressure, J2, and a custom researcher-supplied
operator should be able to enumerate the full reference list directly
from the operator stack — *the dependency graph is the references
list*.

Pre-#93 there was no machinery for this. An author building a
multi-operator simulation had to manually maintain a parallel
bibliography listing each crate they registered, with no guard against
the bibliography drifting away from the actual operator versions in
use. Two runs that reported "1PN + radiation" might have been built
against different crate versions or different commits and produce
materially different trajectories — the standard reproducibility
envelope was incomplete.

Beyond reproducibility, the federation thesis needs the citation
surface to be *opt-out for the federation, opt-in for individual
operators*. A new perturbation crate that ships without a citation
silently breaks the federation's claim to be a citable system. The
trait surface needs to make the omission visible.

---

## Decision

Add a `citation()` method to the `Operator` trait that returns an
`Option<Citation>`, aggregate across the registered stack at the
`System` level, and capture the implementing crate's source state at
the operator's crate compile site (not apsis core's).

### `Citation` struct

```rust
pub struct Citation {
    pub bibtex: &'static str,           // Full entry / entries
    pub doi: Option<&'static str>,       // Bare suffix, e.g. "10.1086/153180"
    pub crate_name: &'static str,        // env!("CARGO_PKG_NAME")
    pub crate_version: &'static str,     // env!("CARGO_PKG_VERSION")
    pub commit_hash: Option<&'static str>, // option_env!("APSIS_<CRATE>_GIT_COMMIT")
}
```

All fields `&'static str` so the citation is zero-cost at runtime
and trivially embeddable in `Box<dyn Operator>`. Rich BibTeX text
goes in `bibtex` rather than splitting into structured fields —
paper.md and supplementary-material consumers want the raw entry,
and operators with multiple references concatenate them in one
string.

### Build-time SHA capture

Each publishing crate ships a `build.rs` that runs `git rev-parse
HEAD` and emits `cargo:rustc-env=APSIS_<CRATE>_GIT_COMMIT=<sha>`.
The runtime treats an empty string as "no commit known" and renders
the citation without a commit line. Re-runs only when `HEAD`
changes (`cargo:rerun-if-changed=../../.git/HEAD`), so the
incremental build cost is one `git rev-parse` per branch switch /
commit.

The capture happens in the **operator's crate**, not apsis core, so
`crate_name` and `crate_version` reflect the implementing source —
1PN's citation pins `apsis-1pn` 0.1.0 even when apsis core is at a
different version, and a downstream researcher's custom operator
crate pins itself.

### Aggregation

```rust
impl System {
    pub fn citations(&self) -> Vec<Citation>;   // dispatch order
    pub fn provenance(&self) -> String;          // human-readable block
}
```

`citations()` walks Hamiltonian → non-conservative → observer stacks
in registration order, skipping operators that return default
`None`. The order is stable so consumers can diff two `provenance()`
outputs across runs to confirm the operator stack stayed bit-equal.

`provenance()` renders the standard layout:

```text
Provenance (1 operator):

  apsis-1pn 0.1.0 (commit 66a8683)
    DOI: 10.1086/153180
    @article{anderson1975, ...}
    @book{will1993, ...}
```

Suitable for paper supplementary material or for embedding into a
snapshot file.

### Python bindings

`System.citations()` returns a `list[dict]` with keys
`crate_name` / `crate_version` / `commit_hash` / `doi` / `bibtex`.
`System.provenance()` returns the rendered string. The Python layout
is identical to Rust — the same diff between two runs works across
language bindings.

### Federation invariant

`Operator::citation` defaults to `None`. Operators with no canonical
reference (test fakes, internal tooling) inherit the default and
are silently skipped. Published perturbation crates (`apsis-1pn`,
`apsis-radiation`, future `apsis-j2`, …) **must** override.
Verified by per-crate tests
(`apsis_1pn::tests::citation_pins_anderson_1975_and_will_1993`,
`apsis_radiation::tests::both_operators_cite_burns_1979`).

The omission is visible: a federation member that ships without a
citation breaks no test in apsis core (the default `None` registers
fine), but its `provenance()` line is missing from the rendered
block. Reviewers reading the supplementary material see the gap.

### Reproducibility envelope

`crate_version` + `commit_hash` together pin the implementation to a
specific source state. Two runs that report identical provenance
blocks ran the same Rust code, modulo platform-level f64 variance
(`apsis::contract` § *What this contract does NOT guarantee* names
this scope). The `commit_hash` is `Option` because not every build
comes from a git checkout (CI from tarball, vendored source); the
crate's `build.rs` decides population policy.

---

## Alternatives rejected

| Alternative | Reason rejected |
|---|---|
| Citation as a separate registry consulted by name | Drifts away from the actual operator instance. A custom researcher operator (no published crate) can't register in a global table. The citation must travel with the operator. |
| Structured BibTeX fields (separate author / title / journal / year) | Operators with multiple references would need an array of structured records. paper.md consumers want the raw BibTeX entry, not a structured form they have to re-serialise. The single `&'static str` field is what the consumers actually use. |
| Capture commit hash in apsis core's `build.rs`, not the operator's | Wrong source state. apsis core's commit reflects core's version, not the publishing crate's. A 1PN operator pinned to apsis-1pn `0.1.0` should report apsis-1pn's commit, not apsis core's. The macro pattern in `Citation::PROVENANCE_RECIPE` documents the right shape. |
| Mandatory citation (no `Option`) | Forces test fakes and internal tooling to fabricate citations. The default `None` lets the federation members take the obligation while keeping the trait usable for non-publishing operators. |
| Citation as a JSON / TOML manifest file alongside the crate | An out-of-band file drifts from the source. The runtime API and the file would need synchronisation, with no compile-time link between them. The trait method puts the citation in the same compilation unit as the operator. |
| Render the integrator's and kernel's citations alongside operators | Integrator and kernel live on different traits (`Integrator`, `Kernel`) that don't expose `citation()` yet. Future expansion will fold them in so the full reference list comes from one call. The current `provenance()` is documented as covering the operator stack only. |

---

## Consequences

**Good:**
- Operator stack is self-documenting. The simulation's full reference
  list reads off `System.provenance()` at runtime — no manual
  bookkeeping, no drift between code and bibliography.
- `crate_version` + `commit_hash` give a reproducibility envelope
  pinned to source state. Two provenance blocks that match
  guarantee bit-equal Rust code (modulo f64 platform variance).
- Federation members publishing without a citation are visible:
  the omission shows up in the rendered block as "this operator
  is contributing force without naming its derivation".
- Python bindings produce identical layout to Rust; cross-language
  diff works.

**Neutral:**
- Operators must add 4–10 lines (the `citation()` impl) and a
  `build.rs` to capture commit hash. Both are template-copyable
  (`apsis-1pn` and `apsis-radiation` are the references).
- `Citation` is `Copy` (all `&'static str`) so the per-operator
  cost is one pointer-pair lookup at aggregation time.

**Watch out:**
- An operator that hand-writes the `crate_name` / `crate_version`
  literals (instead of `env!`) drifts away from the actual crate
  state. The macro pattern in `Citation::PROVENANCE_RECIPE` and the
  template crates document the right shape; reviewers should
  enforce.
- `commit_hash` from an operator built outside a git checkout is
  `None`, so two builds of the same crate (one from git, one from
  a tarball) report different provenance lines. Not a bug — the
  two runs used different inputs.
- Integrator / kernel citations are not aggregated yet. Stated
  scope, not a bug; future work will extend the trait surface.
