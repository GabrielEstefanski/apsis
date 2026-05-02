# Material as physics component — design notes (deferred)

**Date:** 2026-04-29
**Status:** decision recorded; implementation deferred to post-Phase-0
(after 3D port and inclined Mercury smoke test land).

## Context

Today every `Body` constructor forces an astrophysical commitment:
`Body::rocky(m)`, `Body::star(m)`, `Body::asteroid(m)`. There is no
`Body::point_mass(m)` — a user running a Kepler 2-body parity test or
a theoretical figure-8 must pick a material that does not affect the
test. The forced choice is friction without payoff.

At the same time, `Material` is **not** UI cosmetics. It carries
real physics:

* density → physical radius (used by collision Q*, disruption models)
* `q_pr` (radiation pressure efficiency, consumed by `RadiationField`)
* luminosity model (`Star` → main-sequence, `WhiteDwarf` → cooling, ...)

Future physics that would extend this surface: black-hole specific
behaviour (event horizon, accretion), neutron-star equation of state,
icy-body sublimation under Poynting–Robertson, planetesimal
disruption thresholds. **Material as a first-class physics component
is the right vision.** The problem to fix is "forced choice when
nothing material applies", not "Material exists".

## Decision

Introduce a wrapper enum on `Body.material`:

```rust
pub enum MaterialKind {
    /// Abstract point mass — no material commitments. Used for
    /// theoretical scenarios, parity tests, scaffolding.
    PointMass,
    /// Body with explicit astrophysical material.
    Physical(Material),
}
```

And a new constructor:

```rust
impl Body {
    pub fn point_mass(mass: f64) -> Self {
        Self {
            mass,
            material: MaterialKind::PointMass,
            // density / physical_radius / color / luminosity defaulted
            // to values that mean "no extent / no photonic interaction"
            ..
        }
    }
}
```

Existing material constructors (`Body::rocky`, `Body::star`, …) stay
as convenience and produce `MaterialKind::Physical(_)`.

## Why not `Option<Material>`

`Option` says "this value may be absent". Here we have a *different
kind of body*, not a missing material. Encoding domain variation as
`Option` invites:

* defensive `if let Some(m) = ...` branches scattered across
  consumers,
* sentinel-value drift in fields that can no longer be filled in
  meaningfully (density, radius, color when material is `None`),
* breakage every time a future variant (`Composite { core, mantle }`,
  `BlackHole { ... }`, `NeutronStar { ... }`) needs to be added,
  because every consumer matches `Option<Material>` exhaustively.

A wrapper enum makes "this body is a point mass" a tagged state of
the domain, not the absence of one.

## The transversal-coupling rule

The companion rule that prevents the wrapper enum from becoming the
*next* anti-pattern: **never expose the discriminant where a method
will do.** Consumers ask the material a specific question; the
material answers `Option<T>` if the question makes no sense for some
variants:

```rust
impl MaterialKind {
    pub fn radiation_efficiency(&self) -> Option<f64> { ... }
    pub fn solar_luminosity(&self, mass_to_solar: f64, r_solar: f64) -> Option<f64> { ... }
    pub fn density(&self) -> Option<f64> { ... }
    pub fn physical_radius(&self) -> Option<f64> { ... }
}
```

Consumer code reads naturally:

```rust
if let Some(q_pr) = body.material.radiation_efficiency() { /* … */ }
```

The match-on-variant happens *inside* each method. Adding
`MaterialKind::Composite { … }` later updates one place per query;
no consumer outside `material.rs` needs to change. This is the
diff between "tagged-union exposed" (consumer drag) and "tagged-union
encapsulated" (free evolution).

## Evolution path

When five or more methods on `MaterialKind` return `Option<…>` and
all of them return `None` for `PointMass` simultaneously, the correlation
is the signal to decompose:

```rust
pub struct Body {
    pub mass: f64,
    pub physical: Option<PhysicalAttrs>,    // density, radius, structural
    pub photonic: Option<PhotonicAttrs>,    // q_pr, luminosity
    pub composition: Option<Material>,      // narrative astrophysical tag
}
```

Each perturbation declares the attribute it needs. A body with only
`physical` does not enter the radiation pipeline. A body with only
`composition` is a point mass with a descriptive label. This is the
fully-federated form: `apsis-blackhole` can declare a `BlackHoleAttrs`
component without touching core.

The `MaterialKind` step is the cheap intermediate that delivers
ergonomic relief now and signals the structural direction without
paying for the full decomposition.

## Cost / when

* ~1–2h of work
* Touches: `domain::body` (struct field, constructors), every site
  that pattern-matches on `Material` (radiation `q_pr` lookup,
  luminosity assignment in `update_luminosity`, snapshot serialisation
  if material participates), maybe 4–5 unit tests that assumed a
  default material.
* Schedule: **after** the 3D port (commit 7 lands), **before** the
  paper submission push. Out of scope for the current 3D port
  series.

## What this fixes for the federated-perturbation thesis

The thesis claims "perturbations as first-class artifacts via
federation". `apsis-1pn` is a separate crate, the contract works.
But `q_pr` (a per-body attribute consumed by `RadiationField`)
lives in a closed enum `Material` inside core. A reviewer reading
`paper.md` and then `domain/body/material.rs` can ask: *how does
a downstream perturbation declare a custom per-body attribute?*

This refactor is the first move toward an honest answer:

1. `MaterialKind` introduces "not every body declares all
   attributes" as a tagged state.
2. The methods-not-discriminant pattern shows the encapsulation
   discipline.
3. The future per-attribute decomposition (`PhysicalAttrs`,
   `PhotonicAttrs`, …) is the structural realisation.

Step 1 alone is a meaningful answer. Steps 2 and 3 follow when
demand arrives.
