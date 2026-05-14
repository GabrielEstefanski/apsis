# ADR-005 — Federated Perturbation Operators

**Status:** Accepted
**Date:** 2026-05-13
**PRs:** [#88](https://github.com/GabrielEstefanski/apsis/pull/88) (trait split),
[#94](https://github.com/GabrielEstefanski/apsis/pull/94) (first non-conservative publisher)
**Plugin protocol:** `crates/apsis-py-core` (rlib transport),
`crates/apsis-1pn-py` (template consumer)

---

## Context

Pre-#88 the simulator carried a single `PerturbationForce` trait that
attempted to express every kind of force-contributing operator: the
gravitational corrections derived from a Hamiltonian, the dissipative
forces that have no Hamiltonian, and the pure observers that only
read state. One trait, three different contracts mashed together.

That conflation broke as the perturbation library grew:

* **`System::total_energy`** wanted to sum each operator's potential
  contribution. Pure-force operators (test-particle 1PN with no
  closed-form V; future drag terms with no V at all) had nothing to
  return, but the trait demanded an answer.
* **Symplectic integrators** want to know at registration time
  whether a non-conservative operator was attached, so they can warn
  that conservation invariants no longer hold. The single trait
  carried no signal to distinguish.
* **The conservation report** classifies the system into
  `FullyConservative` / `HamiltonianForceOnly` / `Dissipative` based
  on what's registered. With a single trait, classification needed a
  per-operator boolean lookup — fragile and easy to forget.
* **Future plugins** (radiation pressure split into a Hamiltonian
  half + a dissipative half — see ADR-006 / `apsis-radiation`) need
  to register two distinct objects sharing nothing but their
  derivation paper. The single trait forced bundling.

Separately, the federation thesis ([[project_thesis_anchor]],
[[project_paper_positioning]]) argues that perturbations should be
first-class scientific artifacts — independently versioned,
independently validated, independently citable. That requires a
trait surface stable enough that an out-of-tree crate can target it
without depending on apsis internals.

---

## Decision

Replace the single `PerturbationForce` trait with three composable
traits with explicit responsibilities, and ship a plugin transport
that lets out-of-tree consumers (Rust *and* Python) register
operators against the same trait surface.

### Trait split

| Trait | Adds | Used for |
|---|---|---|
| `Operator` | base — `name`, `declared_units`, `kernel_requirements`, `check_regime`, `regime_check_cadence`, `citation`, `observe` | Pure observers; base for the other two |
| `HamiltonianOperator: Operator` | `accumulate_force` + `potential` (default `NotAvailable`) | Conservative force contributions; potential summed into `System::total_energy` when `Value(_)` |
| `NonConservativeOperator: Operator` | `accumulate_force` only | Dissipative forces (drag, PR drag, radiation reaction); register triggers a `warn_diag` on symplectic integrators |

Registration is split to match: `System::add_hamiltonian_perturbation`,
`System::add_non_conservative_perturbation`, `System::register_observer`.
Each method returns `Result<(), Box<UnitSystemMismatch>>` (see ADR-006
for the typed-error rationale).

### Conservation classification

`System::conservation_report()` reads the registered stack and emits
one of:

* `FullyConservative` — all Hamiltonian operators expose
  `Potential::Value`; energy conservation is the gate.
* `HamiltonianForceOnly` — at least one Hamiltonian operator returns
  `Potential::NotAvailable`; force-derivation conservation is the
  gate, energy reporting excludes the unavailable contribution.
* `Dissipative` — at least one `NonConservativeOperator` registered;
  energy is expected to drift.

The classification is a function of *what is registered*, not of
runtime behaviour. Misclassification cannot happen because the trait
choice at registration time is the classification.

### Plugin transport

Out-of-tree consumers (crate authors writing perturbations against
the public API) get two paths:

* **Rust side** — depend on `apsis = { ... }`, implement
  `HamiltonianOperator` or `NonConservativeOperator`, build with the
  observable-constructor convention, ship the crate (`apsis-1pn`,
  `apsis-radiation` are the templates).
* **Python side** — wrap the Rust operator in a typed `PyCapsule`
  via `apsis_py_core::box_into_capsule`, return a pure-Python
  `apsis.Perturbation` instance carrying the capsule and a label.
  `System.add_hamiltonian_perturbation(...)` unwraps the capsule
  back into Rust. The Python wrapper owns nothing physics; the
  capsule is single-shot ownership transfer (consumed at registration).

`apsis-1pn-py` is the reference Python binding; future bindings
follow the same shape (the README of `apsis-1pn-py` is the
template).

### Federation invariants the contract guarantees

The `apsis::contract` rustdoc module names the guarantees with one
test per claim:

* **Determinism (system-level).** `(bodies, perturbations, dt) →
  trajectory` is a pure function bit-for-bit.
* **Newtonian consistency under attach.** Attaching a no-op operator
  produces a trajectory bit-equal to the bare run.
* **Read-only access to base dynamics.** `accumulate_force` takes
  `&self` and `&[Body]`; mutation is structurally impossible without
  interior mutability (which is a contract violation, gated by
  `tests::invariant_perturbation_is_pure_function_of_state`).
* **Commutativity (bit-exact for N = 2).** Order of registration
  does not affect the trajectory at the IEEE-754 level for two
  perturbations.
* **Associativity within ULP for N ≥ 3.** Three perturbations
  iterated against the same buffer in two orders agree within the
  IEEE-754 summation envelope (not bit-equal — that would be
  physically wrong).
* **Additive composition.** Each operator contributes by `+=`;
  sentinel-checked in `tests::composition_perturbation_is_additive_via_sentinel`.

---

## Alternatives rejected

| Alternative | Reason rejected |
|---|---|
| Keep one `PerturbationForce` trait with `is_conservative()` boolean | Boolean opt-out is forgettable; symplectic-warning logic and conservation classification stay fragile. The trait split makes the choice structural. |
| Generic `Perturbation<C: Conservative>` with marker types | Compile-time clean but doesn't survive type erasure into `Box<dyn …>`, which is what registration needs. Runtime dispatch over heterogeneous operator stacks needs distinct trait objects. |
| Single trait + dispatch on `Option<Potential>` return | Same problem as the boolean: a non-conservative operator returning `None` is indistinguishable from a Hamiltonian one whose author forgot the impl. The split removes the ambiguity. |
| Plugin as Python-only ABI (no Rust `Box<dyn Operator>` exposure) | Forces every Python perturbation to be re-implemented in Rust before it can run; the federation thesis wants the Python boundary to be plumbing, not a re-derivation. |
| Plugin transport via `serde`-serialised force functions | Closures aren't serializable; capturing operator state across the FFI requires real ownership transfer. PyCapsule with a typed capsule name is the standard PyO3 pattern. |

---

## Consequences

**Good:**
- Out-of-tree perturbation crates target a stable trait surface.
  `apsis-1pn` and `apsis-radiation` compile against the public API
  alone — no `pub(crate)` access, no patches to core sources.
- Conservation classification and symplectic-warning logic are
  driven by which trait the operator implements, not by a runtime
  flag.
- The Hamiltonian-vs-non-conservative split inside a single physical
  regime works (radiation pressure + PR drag in `apsis-radiation`,
  registered as two distinct operators sharing one paper).
- Rust↔Python boundary has one well-typed transport
  (`apsis-py-core`); future plugin Python bindings copy
  `apsis-1pn-py` line for line.

**Neutral:**
- Three trait names instead of one. The naming makes the contract
  obvious at the call site (`add_hamiltonian_perturbation` reads
  what it does).
- Operator authors pick the right trait at construction time.
  Picking wrong (registering a dissipative force as
  `HamiltonianOperator`) compiles but breaks the conservation
  classification — caught by the conservation report tests, not by
  the integrator.

**Watch out:**
- A future trait that doesn't fit either Hamiltonian or
  non-conservative (e.g., a stochastic forcing with bounded
  variance) needs a fourth trait or a `Stochastic` variant in the
  classification enum. The current `#[non_exhaustive]` on
  `ConservationClass` reserves the room.
- Plugin authors who bypass the observable-constructor convention
  (ADR-006) and publish operators with un-bound units force every
  consumer to know the operator's unit system out-of-band. This is
  a style violation, not a contract violation; the registration
  check still catches mismatches at runtime.
