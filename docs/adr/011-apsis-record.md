# ADR-011 — Apsis Record: reproducibility-first binary certificate

**Status:** Accepted
**Date:** 2026-05-16
**Supersedes:** the `.grav` save format (`crates/apsis/src/io/snapshot.rs`)

---

## Context

The `.grav` save format, present from the first commit, was shaped for
the apsis-app save browser: `save_id` doubled as the filename timestamp,
the trail buffer was serialized in GPU std430 stride for direct upload,
schema versions tracked GUI feature additions. Post apsis-app extraction
(ADR-010), `.grav` has no in-tree consumer.

The federation thesis (ADR-005) claims `Cargo.lock` captures a
simulation's full physical model. Without a corresponding artifact for
the numerical output, the reproducibility claim is incomplete: a
reviewer cannot independently replay a published run.

## Decision

Introduce **Apsis Record** — a binary certificate of a simulation run.

A `.apsis` file contains a TOML header, a stream of binary frames, and a
BLAKE3 trailer. The header captures provenance (apsis git sha, BLAKE3
hash of `Cargo.lock`, per-operator crate name + version + checksum +
declared `KernelRequirements`, unit system, integrator config, kernel
variant + softening, seed) and body metadata. Frames record dynamic
state (Snapshot) and material physical events (collisions, escapes).
The default policy emits initial + final bookend snapshots and every
event; dense trajectory capture is opt-in.

The writer is a `RecordHook` implementing `SimHook` — zero new
extension surface, the existing hook system serves as the integration
point. The reader is `Record::open(path)` returning lazy iterators
over header, events, and optional dense trajectory.

### File format

```text
File:
  MAGIC      "APSR"           4 B
  FMT_VER    u16 LE = 1       2 B
  RSVD       u16 LE = 0       2 B
  HDR_LEN    u64 LE           8 B
  HEADER     UTF-8 TOML       N B
  FRAME 0    (Snapshot @ t=0)        ← initial bookend, required
  FRAME 1..N (Event | Snapshot per policy)
  FRAME N+1  (Snapshot @ t=Tfinal)   ← final bookend, required
  FRAME N+2  (Trailer)               ← required for valid record

Frame:
  KIND  u8      1 B    0=Snapshot · 1=Event · 0xFF=Trailer
                       0x02–0x0F reserved for core (bump FMT_VER on use)
                       0x10–0xFE reserved for operator-emitted events
                                 (opt-in via plugin protocol, post-v0.1)
  T     f64 LE  8 B    sim time
  LEN   u32 LE  4 B    payload bytes
  DATA  …       N B

Snapshot payload:
  n_bodies u32 LE
  per body i = 0..n_bodies (order matches Header.bodies.list):
    pos_x, pos_y, pos_z, vel_x, vel_y, vel_z   f64 LE × 6

Event payload:
  event_kind u8
    0 = Collision: body_a_idx u32, body_b_idx u32, distance f64
    1 = Escape:    body_idx u32,    radius f64

Trailer payload:
  step_count   u64 LE
  frame_count  u64 LE   number of data frames (Snapshot + Event), excludes the trailer
  blake3       32 B    BLAKE3 of the frame stream only (Snapshot + Event bytes)
```

The trailer's `blake3` covers the frame stream, not the header or the
trailer itself. The header is plaintext + re-parseable, and
`created_utc` is per-run wall-clock metadata; hashing it would couple
the content-addressable trailer to the timestamp and break the
"same `{seed, config}` → byte-equal frame stream + trailer"
contract that grounds the paper claim.

Header TOML schema:

```toml
[apsis]
version = "0.1.0"
git_sha = "abc123..."
created_utc = "2026-05-16T11:23:45Z"

[reproducibility]
cargo_lock_blake3 = "deadbeef..."   # 64-char hex
seed = 42

[unit_system]
g = 1.0
length = "AU"; mass = "M_sun"; time = "yr/2pi"

[integrator]
kind = "Ias15"
dt_mode = "Adaptive"
initial_dt = 0.01
[integrator.params]
epsilon = 1.0e-9        # kind-specific; nested so unknown integrators
                        # can carry arbitrary JSON without schema churn

[kernel]
variant = "Newton"           # or "Plummer"
softening = 0.0              # required when variant = "Plummer"; absent otherwise

[[operators]]
name = "apsis-1pn"
version = "0.1.0"
crate_hash = "..."           # 64-char hex from Cargo.lock checksum
[operators.requirements]      # what this operator declares it needs from the kernel
kernel_exactness = "exact"   # "exact" | "softened_ok"
kernel_continuity = "smooth" # "smooth" | "c0_ok"
# Future requirement enum variants are added here; reader treats unknown
# fields as a forward-compat warning, not an error.

[bodies]
count = 9
[[bodies.list]]
name = "sun"
mass = 1.0
density = 1.408
physical_radius = 4.65e-3
color = [255, 233, 100]
q_pr = 0.0
albedo = 0.5
class = "Star"
```

### API

```rust
// Writer
use apsis::records::{Header, RecordHook, RecordPolicy, provenance::header_from_system};
let mut sys = System::new(bodies, units);
let header = header_from_system(&sys, seed, /* lock_path */ None)?;
let hook = RecordHook::with_header("run.apsis", header, RecordPolicy::default())?;
sys.hooks_mut().register(0, Box::new(hook));

// Reader
use apsis::records::Record;
let rec = Record::open("run.apsis")?;
rec.header().reproducibility.seed;
for ev in rec.events()? { /* … */ }
let (initial, final_) = rec.bookends()?;
```

The `Header` argument is mandatory and explicit so the lockfile-hash
+ operator-provenance gathering is a separate, testable step.
`provenance::header_from_system` is the canonical builder; tooling that
constructs a `Header` from a manifest can bypass it.

```rust
pub enum RecordPolicy {
    BookendsAndEvents,   // default — initial + final + events only
    EveryNSteps(u32),    // bookends + Snapshot when steps() % N == 0
    EveryTime(f64),      // bookends + Snapshot when t crosses k·Δt
    Dense,               // every step — debug
}
```

Cadence edge cases:

- `EveryNSteps(N)`: Snapshot emitted after the post-step hook fires when
  `system.steps() % N == 0`. The initial bookend at `t=0` is separate
  from this counter.
- `EveryTime(Δt)`: tracks the last snapshot's `t_last`; emits a Snapshot
  when `t >= t_last + Δt` after a step completes. If a single step
  crosses multiple intervals (large adaptive `dt`), one Snapshot is
  emitted at the post-step time — not back-filled per interval.
- `Dense`: a Snapshot for every step; the initial bookend is the t=0
  snapshot, the final bookend is the last step's snapshot.

Python:

```python
import apsis

# Write — single keyword API on System. Default policy: bookends + events.
sys = apsis.System(bodies=..., units=..., integrator="ias15", dt=1e-3)
sys.attach_record("run.apsis", seed=42)
# Variants: sys.attach_record(path, every_steps=100)
#           sys.attach_record(path, every_time=0.1)
#           sys.attach_record(path, dense=True)

# Read — apsis.Record exposes header (TOML string), events(), snapshot_count()
rec = apsis.Record("run.apsis")
```

Writing in Python is mediated by the Rust hook — the lockfile-hash
provenance is computed at file creation time by the Rust binary. There
is no Python-side `RecordHook` class; `System.attach_record` is the
sole entry point.

### Deliberate format decisions

- **Trailer required for `Record::open` to succeed.** An incomplete
  record is not a certificate.
- **BLAKE3 trailer over the frame stream only.** Header is plaintext
  TOML carrying per-run metadata (`created_utc`); covering it would
  couple the content hash to wall-clock time and break the "same
  `{seed, config}` → byte-equal trailer" contract. The same primitive
  hashes `Cargo.lock` separately in the header (`cargo_lock_blake3`):
  one algorithm, two scopes, neither composed.
- **Body metadata in Header only, not Snapshot.** Snapshot frames carry
  pos + vel only. Static body inventory is the current contract;
  mass-dynamic perturbations would add a `MassUpdate` frame kind
  without breaking format.
- **Little-endian everywhere** (`f64::to_le_bytes`). Portable across x86
  and ARM.
- **`format_ver = 1`** with policy "bump on breaking change". No
  backward compat promise pre-apsis-1.0.

## Consequences

**Architectural wins:**

- `{record, Cargo.lock}` is the content-addressable closure of an apsis
  run: re-running with the same lockfile reproduces the byte-equal
  frame stream. The federation thesis extends from the input side
  (composition of crates) to the output side (the record).
- The hook system gains a flagship in-tree consumer, validating the
  extension surface for downstream uses (instrumentation, monitoring,
  custom recording).
- The bookend-by-default policy reflects the framing of a record as a
  certificate, not a recording: small, diff-friendly, citable.
- Frame kinds `0x10-0xFE` are reserved for operator-emitted events,
  giving per-operator event taxonomies a stable space without forcing
  a format renegotiation when they are introduced.

**Cleanup carried by the same PR:**

| path | LOC | reason |
|---|---|---|
| `crates/apsis/src/core/physics_thread.rs` | 1512 | real-time GUI orchestration |
| `crates/apsis/src/core/precision_run/` | — | cascade with physics_thread |
| `crates/apsis/src/core/trail/` | 27 | replaced by Snapshot frames |
| `crates/apsis/src/core/system/snapshot.rs` | 87 | `.grav` bridge |
| `crates/apsis/src/io/snapshot.rs` | ~1000 | `.grav` save format |
| `crates/apsis/src/domain/field/` | 21 | UI overlay trait orphan |

Total: ~2700 LOC of GUI-runtime infrastructure retired from the core.

**Migration cost:**

- `.grav` files become unreadable in v0.1. No published consumers
  exist; apsis-app vendored its own `.grav` reader at extraction time.
  External migration cost: zero.
- `io::headless` migrates its snapshot output from `.grav` to `.apsis`;
  the `run_config.snapshot_interval` field is reinterpreted as
  `RecordPolicy::EveryTime`. CSV export (`io::recorder`) is unchanged.

**Paper claim added (v0.1):**

New subsection §"Reproducibility certificate" under Implementation, plus
one sentence appended to §Summary:

> The full simulation, from physical model (`Cargo.lock`-pinned operator
> crates) to numerical output (an Apsis Record), is bit-exactly
> reproducible from a single hash-pinned configuration.

**Roadmap (next iteration):**

| item | shape |
|---|---|
| `Diagnostic` frame emitter | second built-in `SimHook` (`EnergyTrackerHook`) that consumes the reserved frame kind and writes ΔE, ΔL, ΔP at a configurable cadence |
| Mid-run snapshot resume | header schema extended with per-integrator scratch state (IAS15 b-coefficients, dt history, Mercurius hybrid mode) so a record can be opened and the simulation continued bit-exactly from any Snapshot |
| Semantic `Record::diff` API | structured comparison categorising differences by operator (name/version/hash), integrator config, kernel variant, and trajectory; richer than byte-diff because the header is typed |
| Header enrichment | `rustc_version`, `generated_by`, OS `kernel_version` for full build-environment provenance |
| `Record::open` operator/kernel compatibility validation | reader runs the same `KernelRequirements` check the registry runs, against the header's operator list |
| User-defined Event kinds | reserved range (`0x10–0xFE`) is in place; plugin protocol wiring follows when an in-tree producer needs it |

**Out of scope:**

| | rationale |
|---|---|
| Alternative formats (HDF5, Parquet, NetCDF) | one canonical binary; alternative formats are downstream tooling territory |
| Compression of dense trajectories | uncompressed is sufficient at the body counts the validation portfolio addresses |

## Naming

`io::recorder::Recorder` (CSV stream) and `records::Record` (binary
certificate) coexist. Different module, type, purpose; the CSV recorder
streams data points (verb), the Apsis Record is an artifact (noun).

## References

- ADR-005 — federated perturbation operators (thesis the record extends)
- ADR-008 — Kernel as System parameter (record captures kernel + softening)
- ADR-009 — consolidated Python distribution
- ADR-010 — apsis-app extraction (removed the historical `.grav` consumer)
- REBOUND `SimulationArchive` v5.0 (2026-05) — prior art for binary
  record format in the N-body context
- BLAKE3 — <https://github.com/BLAKE3-team/BLAKE3>
