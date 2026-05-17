# ADR-012 — Apsis Record v0.2: diagnostics, diff, resume, enrichment

**Status:** Draft
**Date:** 2026-05-17
**Extends:** ADR-011 (apsis record v0.1)

---

## Context

ADR-011 shipped apsis records as a reproducibility certificate: header
with provenance, frame stream with snapshots and events, BLAKE3 trailer
over the frames. The format is sufficient to reproduce a run from
`{record, Cargo.lock}` by re-running from initial state — but it
cannot resume a partial run, it ships no diagnostic emitter for the
reserved `Diagnostic` frame kind, and consumers comparing two records
fall back to byte diffs against a structured format.

Three gaps surfaced when comparing v0.1 against REBOUND's
SimulationArchive:

- **Mid-run resume.** REBOUND restores integrator scratch (IAS15
  b-coefficients, Kepler step state, Mercurius hybrid mode) from any
  snapshot index; apsis cannot. For long runs that crashed, "re-run
  from initial state" costs the wall-clock the failed run already
  consumed.
- **Diagnostic emitter.** Energy and angular-momentum drift over time
  is the first numerical-stability question any reader asks. Reserving
  a `Diagnostic` frame kind without shipping an emitter ships
  infrastructure without value.
- **Semantic diff.** Comparing two records means running an external
  byte diff on a structured format. The header is typed TOML;
  differences are categorisable (operator versions, integrator config,
  kernel variant, trajectory), and the right place for that logic is
  the library that owns the format.

Two smaller items track alongside these:

- The build-time `APSIS_GIT_COMMIT` capture in `crates/apsis/build.rs`
  re-runs on `.git/index` changes but not on workdir-only edits, so
  the `-dirty` suffix can be silently stale.
- Header provenance pins apsis version + `Cargo.lock` hash but not
  `rustc` version. f64 codegen varies between rustc releases; the
  reproducibility claim is incomplete without it.

## Decision

Ship apsis record format version 2 with the items below. v0.1 records
remain readable through the v0.1 code path; the reader dispatches on
`FORMAT_VER`.

### Mid-run snapshot resume

A `ResumeState` frame kind (`0x04`) carries per-integrator scratch
state. The writer emits one alongside every Snapshot frame when the
hook is constructed with `RecordHook::with_resume_capture(true)`.
Default omits ResumeState frames — they are several KB per snapshot
for IAS15, so opting in is explicit.

`Record::resume_from(snapshot_idx) -> Result<System, _>` reads the
`Snapshot` + `ResumeState` pair at the given index and rebuilds a
`System` whose subsequent `step()` calls produce a bit-equal
continuation of the original run.

Per-integrator state serialisation lives in each integrator's own
module via a new `Integrator::resume_state(&self) -> ResumeStateBytes`
+ `Integrator::restore_from(&mut self, ResumeStateBytes)` trait pair
(default no-op for stateless integrators).

### Diagnostic emission opt-in on `RecordHook`

`RecordHook` gains a `with_diagnostics(DiagnosticCadence)` builder.
When enabled, the hook emits a `Diagnostic` frame (kind `0x03`) at
the configured cadence carrying ΔE/E and ΔLz/Lz — the two drifts
named in §Context. Linear momentum is omitted: every integrator
apsis ships uses pairwise Newton-third-law forces, so `Σ m·Δv` is
conserved to roundoff and is a sanity gate rather than a drift
diagnostic; sims set up with `move_to_center_of_momentum()` also
have `|P₀| ≈ 0`, making any relative normalisation degenerate.

Diagnostic cadence is intentionally orthogonal to snapshot cadence
(`RecordPolicy`): a `BookendsAndEvents` policy with
`DiagnosticCadence::EveryNSteps(100)` is the expected long-run
configuration. Cadence options mirror `RecordPolicy`:
`Off | EveryNSteps(u32) | EveryTime(f64)`, default `Off`.

### Semantic `Record::diff` API

```rust
pub fn diff(a: &Record, b: &Record) -> RecordDiff;

pub struct RecordDiff {
    pub header: Vec<HeaderChange>,
    pub frames: FrameStreamDiff,
}

pub enum HeaderChange {
    OperatorAdded { name, version, crate_hash },
    OperatorRemoved { name, version },
    OperatorVersionChanged { name, before, after },
    IntegratorChanged { field, before, after },
    KernelChanged { before, after },
    SeedChanged { before, after },
    ...
}

pub struct FrameStreamDiff {
    pub event_count_delta: (usize, usize),
    pub trajectory_rms_at_final: f64,
    pub trailer_blake3_match: bool,
}
```

Categorises differences by source. A reviewer sees "operator X v1 →
v2 + trajectory rms 3.2e-7 at t_final" rather than "31 bytes differ
at offset 4287".

CLI surface piggybacks on the `apsis` binary entry (re-introduced
narrowly for this — see §"Out of scope" in ADR-011 about the prior
removal of the headless CLI): `apsis diff a.apsis b.apsis` renders
the `RecordDiff` as a human-readable report.

### Header enrichment

Three new fields in `[apsis]` section:

```toml
[apsis]
version = "0.1.0"           # existing
git_sha = "..."             # existing
created_utc = "..."         # existing
rustc_version = "1.95.0 (cd5fbb0 2026-04-12)"   # NEW
generated_by = "apsis 0.1.0"                     # NEW
[apsis.platform]                                 # NEW (optional)
os_kernel = "Linux 6.8.0-31-generic"             # NEW (omitted if detection unavailable)
```

`rustc_version` captured by extending the existing `build.rs`. The
field is mandatory in v0.2 readers; absent values cause
`UnknownRustcVersion` (warning, not error — old records can still be
opened, just without rustc provenance).

`generated_by` defaults to `format!("{} {}", CARGO_PKG_NAME,
CARGO_PKG_VERSION)`. Future external tooling (Python apsis-py
wrapper, custom Rust binaries) can override.

`os_kernel` is optional and per-platform: Linux reads from `uname`,
macOS likewise, Windows from `GetVersionEx` (or skipped if no
straightforward read). `hostname` deliberately omitted — privacy risk
in published research.

### Reader-side operator/kernel validation

`Record::open` already validates the format. v0.2 adds a second pass:
for each `[[operators]]` entry, run the same `KernelRequirements`
satisfaction check the registry runs at hook registration time, against
the record's `[kernel]` block. Mismatches surface as
`RecordError::IncompatibleKernel { operator, required, actual }`.

This catches "record claims apsis-1pn + Plummer kernel" malformations
at open time rather than at re-run time. No new types needed; the
existing `KernelRequirements::satisfies` logic is reused.

### `build.rs` dirty detection fix

Drop the `cargo:rerun-if-changed=../../.git/index` line and rerun the
build script every cargo invocation (via no `rerun-if-changed` for the
SHA-capture inputs). The recompile cost is two `git` calls per
incremental build, dominated by linker time; the win is that
workdir-only edits get the correct `-dirty` suffix.

## Consequences

**Format version bump:** records become `FORMAT_VER = 2`. v0.1 readers
on v0.2 files emit `UnsupportedFormatVersion(2)`. v0.2 readers on v0.1
files dispatch on the version field and read the older schema; the
reserved `0x02–0x0F` core-kind range absorbs the new `Diagnostic` and
`ResumeState` kinds.

**RecordHook gains a second responsibility:** writing snapshots +
events (provenance) and emitting periodic conservation drifts
(diagnostics). Both belong to "what this record captures from one
run", which is the cohesive scope of the certificate. A second hook
publishing Diagnostic frames was considered and rejected: two hooks
writing the same `.apsis` file would need `Arc<Mutex<RecordWriter>>`
coordination for a single-file output that has no other producer.

**Paper claim strengthens:** the v0.1 claim was
"`{record, Cargo.lock}` is the content-addressable closure of an apsis
run". v0.2 extends to "`{record, Cargo.lock}` reproduces the run from
initial state OR from any captured snapshot". The federation thesis
gains the parity with REBOUND that ADR-011 explicitly omitted.

**Out of scope (still):**

- Alternative formats (HDF5, Parquet, NetCDF)
- Compression of dense trajectories
- User-defined Event kinds (reserved range stays reserved; first
  in-tree producer triggers wiring)
- Cloud / remote record sinks
- Multi-file records

## References

- ADR-005 — federated perturbation operators
- ADR-008 — Kernel as System parameter
- ADR-011 — apsis record v0.1
- REBOUND `SimulationArchive` v5.0 — prior art for the resume capability
- BLAKE3 specification
