# ADR-008 — Kernel as System Parameter, Not Body Field

**Status:** Accepted
**Date:** 2026-05-14
**PR:** [#121](https://github.com/GabrielEstefanski/apsis/pull/121)

---

## Context

Pre-PR-#121 every `Body` carried a `softening: f64` field, computed at
construction by `default_softening(mass) = EPS_BASE · m^{1/3}` and
applied pairwise via `ε²_ij = (ε²_i + ε²_j) / 2` inside force kernels.
Two `Kernel` trait impls existed: `PlummerKernel` (default,
softening-aware) and `TruncatedPlummerKernel` (counter-test fixture).

The default was Plummer with non-zero per-body ε. Paper-grade work
(Mercury 1PN, apsis-radiation gates, federation tests) routinely called
`Body::unsoftened()` to zero the field — and `System::with_exact_gravity()`
existed as a "unsoften the whole system" convenience. Both were
workarounds for the wrong default.

Three concrete failures motivated revisiting the architecture:

1. The Lagrange equilateral and Euler collinear three-body presets
   (PR #120's review surfaced these) showed bodies escaping on
   hyperbolae rather than orbiting their analytical equilibrium —
   the default softening contaminated the fine force balance and
   the apparent fix (correct ω formula) didn't help while softening
   silently broke the contract.

2. The Mercury 1PN federation gate, the observable-inversion
   central-force round-trip gate, and most apsis-radiation tests had to opt out of
   the default explicitly. The most-used scenarios fought the default
   the most.

3. `feedback_fine_physics_traps` (memory) catalogued the pattern:
   "precision physics derived around a clean baseline fails
   catastrophically when the baseline assumption is violated".
   The default was the baseline-violating case.

## Decision

Drop `softening` from `Body` entirely. Make the gravity kernel a
single `NewtonKernel { epsilon: f64 }`-parameterised system-level
configuration. Default `epsilon = 0` is exact `1/r²` Newton; `epsilon > 0`
is the Plummer-softened regularisation. The two are mathematically the
same kernel (continuous in the `ε → 0` limit), so there is no separate
`PlummerKernel` impl.

```rust
// Before
struct Body { /* ... */ pub softening: f64 }

enum Kernel {
    Plummer { /* per-body softening read from bodies */ },
    TruncatedPlummer { ... },
}

// After
struct Body { /* ... */ /* no softening */ }

trait Kernel {
    fn potential(&self, r_squared: f64) -> f64;
    fn acceleration_factor(&self, r_squared: f64) -> f64;
    fn epsilon_squared(&self) -> f64;  // 0 for exact Newton
    /* ... */
}

struct NewtonKernel { pub epsilon: f64 }  // ε=0 exact, ε>0 Plummer
struct TruncatedPlummerKernel { /* counter-test only */ }
```

The default `System` constructs with `NewtonKernel::exact()` (ε=0).
Cluster / cosmological scenarios that need softening opt in via
`System::with_kernel(Arc::new(NewtonKernel::new(0.01)))`.

`System::with_exact_gravity()` is removed (default is exact). The
`Body::unsoftened()` builder is removed (no field to zero).
`default_softening()` and `pair_eps2()` are removed. The snapshot
codec drops the per-body softening byte and the system-level
`softening_scale` header field; schema bumps `v11 → v12` with reader
fallback that discards the bytes when reading older saves.

## Validation cross-check

REBOUND uses the same architecture:

- `reb_simulation` carries a single `double softening` field
  (default `0.0`)
- `reb_particle` (Body equivalent) has no softening field
- Force kernels read `r->softening` and use `r² + softening²`

Orekit (high-precision satellite mechanics) does not have a softening
concept at all — its gravity force models are exact `1/r²` exclusively.

## Consequences

**Architectural wins:**

- `Body` is leaner (8 bytes per body × N less)
- `BodyArrays` SoA snapshot drops one column (40 → 32 bytes per row)
- Force evaluators take a single `r²` argument (kernel reads ε from
  state); `pair_eps2` helper is gone
- Pre-existing tests with `.unsoftened()` calls become redundant —
  bulk removal across ~60 sites
- Newton/Plummer collapse into one impl: less branching, less
  pattern matching, less serialisation surface, no edge cases
- `ε → 0` limit is continuous in the math and in the type system

**Backward compatibility:**

- `.grav` snapshot files at schema v11 still load correctly; the
  reader discards the dropped bytes
- Templates that used to rely on default softening to mask
  close-encounter blowups (Pythagorean, Chaotic Ejection) now declare
  `suggested_integrator: Some(Mercurius/Ias15)` so the integrator
  layer handles the close encounters via adaptive substepping
- Python binding: `Body.<material>(...)` no longer accepts a
  `softening=` kwarg; `unsoftened()` method is gone; `softening`
  getter is gone

**Future kernels:**

The single `Kernel` trait stays open for genuinely different physics:
`YukawaKernel { range: f64 }`, `MONDKernel { a0: f64 }`, etc. Each is
a new impl with its own state. Newton/Plummer are a single
`NewtonKernel` because they are the same family parameterised by ε.

**Barnes–Hut acceptance criterion at ε > 0:**

The BH walk's `s/d < θ` test now uses pure geometric distance
`d = √(dx² + dy² + dz²)` rather than the prior softened
`d = √(dx² + dy² + dz² + body_eps²)`. For the new exact default
(ε = 0) this is a no-op. For callers that opt into
`NewtonKernel::new(ε > 0)` for cluster work, the BH acceptance
becomes uniformly slightly more aggressive at close range (smaller
`d` → harder to accept the multipole approximation), which is
conservative for accuracy. Existing BH parity tests run at ε = 0 so
this behaviour change is not regression-gated; cluster workloads at
ε > 0 should re-baseline opening-angle settings if they depended on
the prior softened acceptance.

## Anti-pattern this rules out

A future PR proposing "add `softening` back to `Body` for per-body
heterogeneous softening" should be rejected. The lesson learned:
softening is a property of the gravity computation, not of any
individual body. Per-body softening creates the silent-baseline-
violation failure mode this ADR removes; if some genuine use case
needs per-particle force-law variation, the right design is a
`Kernel::variable_epsilon(...)` impl that takes the per-body data
explicitly, not a field hidden on `Body`.

## References

- Plummer (1911). *Mon. Not. R. Astron. Soc.* 71, 460–470.
- Dehnen & Read (2011). *Eur. Phys. J. Plus* 126, 55. (softening review)
- REBOUND source: `rebound/src/rebound.h` (system-level `softening` field).
- `feedback_fine_physics_traps` (project memory entry on softening
  contamination of fine physics).
