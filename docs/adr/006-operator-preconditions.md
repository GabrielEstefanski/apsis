# ADR-006 — Operator Preconditions: Kernel, Regime, Units

**Status:** Accepted
**Date:** 2026-05-13
**PRs:** [#88](https://github.com/GabrielEstefanski/apsis/pull/88) (trait surface),
[#89](https://github.com/GabrielEstefanski/apsis/pull/89) (observable-constructor + UnitSystem binding),
[#90](https://github.com/GabrielEstefanski/apsis/pull/90) (regime-of-validity static + dynamic),
[#91](https://github.com/GabrielEstefanski/apsis/pull/91) (apsis-1pn-py getter fix),
[#93](https://github.com/GabrielEstefanski/apsis/pull/93) (typed `UnitSystemMismatch` Err)

---

## Context

Every perturbation has preconditions. 1PN assumes the unperturbed
gravity is bit-exact 1/r (softening invalidates the derivation). It
assumes `m_secondary / m_primary ≪ 1` (test-particle pairwise
expansion drops EIH cross-terms beyond ~1 %). Its `c` value depends
on which unit system the simulation runs in.

Pre-#88 those preconditions lived as comments in rustdoc and as
implicit assumptions in the operator's author's head. The simulator
had no machinery to surface a violation: register a 1PN operator
built for `UnitSystem::solar()` (IAU, G≈4π²) into a
`UnitSystem::solar_canonical()` system (G=1) and the integration
proceeded silently with a `c` value off by ~0.009 %. Register on a
softened kernel and the numerical apsidal precession from softening
alone exceeded the GR signal by ~2 × 10³ — energy and angular
momentum stayed conserved at machine precision while the trajectory
was physically wrong.

The locked memory entry [[feedback_fine_physics_traps]] catalogues
this pattern: precision physics derived around a clean baseline can
fail catastrophically and silently when the baseline assumption is
violated. Surfacing the violation is the contract; auto-fixing would
erase it.

Three precondition families surfaced once the trait split landed:

1. **Kernel preconditions** — what shape of gravitational baseline
   the operator's derivation requires (Exactness, Continuity).
2. **Regime-of-validity preconditions** — what physical bounds the
   body state must satisfy for the derivation to apply (mass ratio,
   eccentricity, v/c, periapse, …).
3. **Unit-system preconditions** — what `UnitSystem` the operator's
   dimensional parameters were computed for.

Each family needed structural enforcement at registration; the
regime family additionally needed dynamic re-checking because body
state evolves.

---

## Decision

Each operator declares its preconditions through `Operator` trait
methods, the `System` enforces them at registration and (for regime)
at a configurable cadence during integration, and violations surface
through one of two channels chosen by severity.

### Kernel preconditions — `Operator::kernel_requirements`

Returns a `KernelRequirements { required_exactness, min_continuity }`.
The `System::add_*_perturbation` methods compare against the active
kernel's `properties(...)` and emit one structured `warn_diag` per
violated invariant. Soft channel — the operator still registers and
integration proceeds; the warning is the contract.

* Exactness violation message includes the count of softened bodies
  and the maximum softening, so the user sees what to fix.
* Continuity violation names the discontinuity tier the kernel
  reports and what the operator required.
* `KernelRequirements::none()` is the default; operators must opt in
  to declare a requirement.

### Regime-of-validity — `Operator::check_regime` + `regime_check_cadence`

* `check_regime(bodies, t) -> Vec<RegimeViolation>` returns one entry
  per crossed bound at the current state. Empty vector when within
  the operator's envelope.
* `regime_check_cadence() -> usize` declares how many outer
  integration steps between dynamic checks (default 100; 1PN
  overrides to 15 000 because its only bound is mass-ratio, which
  is static).
* Each `(operator, bound)` pair fires *exactly one* `warn_diag` per
  `System` lifetime via `regime_warnings_emitted` dedup. A persistent
  violation does not respam the bus.
* `Severity` (`Approaching`, `Exceeded`, `Hard`) is `#[non_exhaustive]`
  so future tiers don't break consumers; one bound key per severity
  prevents an Approaching → Hard escalation from being silently
  suppressed by the dedup state.
* `System::reset_regime_warnings()` re-arms the bus when the caller
  deliberately changes scenario (loaded a snapshot, replaced bodies).

Soft channel — same as kernel violations. Integration proceeds with
the user's choice; the warning records that the derivation no longer
strictly applies.

### Unit-system binding — `Operator::declared_units` + typed `Err`

* `declared_units() -> Option<UnitSystem>` returns `Some(units)` for
  operators whose parameters are dimensional (1PN's `c`, radiation's
  per-body β cached prefactor, future J2's `J₂` coefficient). Default
  `None` for unit-agnostic operators (constant pushes, dimensionless
  couplings).
* `System::add_*_perturbation` returns `Result<(),
  Box<UnitSystemMismatch>>`. On mismatch the operator is **not**
  registered; no kernel-precondition warnings fire; no regime check
  runs; `*_perturbations.push` does not happen.
* Hard channel — caller owns the policy. Propagate with `?`, log and
  skip, swap operator, fall back, or `.expect(...)` for end-of-line
  scripts. The Python binding maps the boxed `Err` to
  `apsis.UnitSystemMismatchError` (subclass of `Exception`,
  preserving `operator`, `operator_units`, `system_units` attributes).

`UnitSystemMismatch` is `Box`ed in the `Err` variant per
`clippy::result_large_err` — `UnitSystem` is ~70 bytes per copy, two
of them blow past the 128-byte default.

### Why two channels — soft warn vs hard `Err`

`UnitSystemMismatch` makes integration silently *wrong*: a 1PN
operator built for one unit system applied in another produces
internally-consistent dynamics that are physically incorrect. The
caller has no way to recover by inspection — the trajectory will
look fine. Hard `Err`.

Kernel and regime violations make integration *less applicable to
the operator's derivation* but the trajectory is still a
self-consistent integration of the registered force model. The user
may legitimately want to proceed (running a softened-gravity
scenario knowing 1PN is now decorative; running a binary mass ratio
to see what happens past the test-particle bound). Soft warn.

The two-tier semantics is locked in `apsis::contract` § *Failure
model*.

### The observable-constructor convention

Operators with dimensional parameters expose three constructor
families (locked in `apsis::contract` § *Observable constructor
convention*):

* **Regime-based constructor** (`for_units(units)`,
  `for_<regime>(...)`): the most common path. Unit system in,
  derived parameters out. `PostNewtonian1PN::for_units(units)`
  derives `c` exactly from `c_SI · T_scale / L_scale`.
* **Observable-inversion constructor** (`from_<observable>(...)`):
  takes a measured observable and inverts to the operator's
  coupling. Reserved for future operators (central-force law from
  apsidal precession measurement, J2 from secular nodal regression).
* **Raw escape** (`from_raw_<param>(...)`): explicit "I computed
  this externally". Pinned to a `UnitSystem` so the registration
  check still validates frame consistency.

> **Naming note.** Earlier revisions of this ADR called these the
> "Pattern A" and "Pattern B" constructors. The descriptive names
> are preferred going forward to avoid the impression of provisional
> taxonomy and to keep the public API self-explanatory.

Every constructor family carries the `UnitSystem` so the operator
publishes `Some(units)` from `declared_units`. Unit confusion is
structurally impossible past the registration boundary.

### Bug discovered while building demo

PR #91 fixed a separate getter bug in `apsis-1pn-py`:
`call_method0("length_scale_si")` resolved to `getattr(...)()` and
failed with `'float' object is not callable` because the Python
`UnitSystem` exposes the scale fields as `#[getter]` properties, not
methods. The fix (`getattr` + `extract::<f64>`) preserves the
duck-typed extraction across the cdylib boundary. Discovered while
running the unit-binding demo end-to-end; included as evidence that
the cross-language plugin path needs the same precondition checking
the Rust path has.

---

## Alternatives rejected

| Alternative | Reason rejected |
|---|---|
| Kernel and regime violations as hard `panic!` | Domain validation panics break composability. A long pipeline with one bad operator should not crash the whole run. |
| Unit-system mismatch as soft `warn_diag` | Silent wrongness is worse than verbose failure here. The trajectory looks fine while being physically inconsistent — there is no recovery via inspection. Hard `Err` forces the caller to acknowledge. |
| Auto-correct softening to zero on registering an operator that requires Exact | Surfacing the violation is the contract; auto-fixing erases the user's choice and hides the bug. |
| Single `Result<(), PreconditionError>` covering all three families | Conflates orthogonal failure modes. The caller policy for "softening violates 1PN's exactness requirement" (acknowledge and proceed) differs from "operator built for IAU solar but System runs canonical solar" (must fix). One typed error per family keeps the policy choice clean. |
| Static-only regime checks | A regime can change dynamically (eccentricity growth, mass ratio change via planet–planetesimal merger). Static-only would catch bad initial conditions but miss in-run drift. |
| No dedup on regime warnings | A persistent violation across 10⁴ outer steps generates 10⁴ warnings — drowns the log. The dedup makes the warning a one-time signal: the *first* time the bound is crossed, the user sees it. |
| Regime-based constructor only, no raw escape | Raw escape exists for hypothetical experiments ("what if c were 5 % larger?") and for cases where `c` (or any param) is computed by neighbouring code so cross-checking with `units` would be redundant. The regime-based path alone forces a workaround. |

---

## Consequences

**Good:**
- Three orthogonal precondition families have one consistent
  surface (`Operator` trait method per family).
- Unit-system mismatch is structurally impossible to ignore;
  silent-wrongness class of bug eliminated.
- Soft channel keeps integration usable when the user knows what
  they're doing; hard channel forces acknowledgement when the
  alternative is silent wrongness.
- Regime warn-once dedup keeps the bus signal-rich.

**Neutral:**
- Operator authors must opt in to declare each family
  (`KernelRequirements::none()` is default; `declared_units()`
  default `None`; `check_regime` default empty). The federation
  thesis depends on authors taking this seriously, but the contract
  is permissive of unit-agnostic / regime-agnostic operators.
- Python callers must catch `apsis.UnitSystemMismatchError`
  explicitly; `try/except Exception` works.

**Watch out:**
- Adding a new precondition family (numerical stability bounds, FFT
  resolution requirements for spectral operators) requires extending
  the trait surface. The current shape reserves room
  (`#[non_exhaustive]` on `Severity`, `Potential`,
  `ConservationClass`).
- Log channel separation (`Source::Domain` vs
  `Source::Performance`) is open backlog: regime warnings currently
  share a bus with perf logs. Operators don't see the bleed; consumers
  do. Tracked for a future PR.
