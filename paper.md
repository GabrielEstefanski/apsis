---
title: 'APSIS: Verified Extension Contracts for N-Body Simulation in Rust'
tags:
  - Rust
  - N-body
  - gravitational dynamics
  - general relativity
  - simulation software
authors:
  - name: Gabriel Braga Estefanski
    affiliation: 1
affiliations:
 - name: Independent researcher
   index: 1
date: 24 April 2026
bibliography: paper.bib
---

# Summary

`apsis` is a Rust library providing verified extension contracts for
gravitational N-body simulation. Physical preconditions of perturbation
forces — for example, whether a correction assumes an unsoftened `1/r`
potential — are promoted from informal documentation to type-level
declarations, checked at extension registration and enforced through an
out-of-tree companion crate whose compilation runs as a
continuous-integration gate. The mechanism is demonstrated end-to-end
by `apsis-1pn`, which implements the first-post-Newtonian Schwarzschild
correction; as evidence that the integrator stack resolves 1PN-scale
effects at the accuracy the verification claim requires, the
demonstration reproduces Mercury's perihelion precession to within
4.4×10⁻⁶ of the general-relativistic prediction over 500 orbital
periods.

The solver provides four integration schemes — second-order Velocity
Verlet, fourth-order Yoshida composition, Wisdom–Holman mixed-variable,
and the adaptive Gauss–Radau IAS15 scheme [@ReinSpiegel2015] — alongside
stable public traits for user-registered force models and perturbations.
The library's scope is narrow by intent: the solver is two-dimensional
and targets small-to-medium body counts (N ≤ 10³). Large-N collisionless
dynamics, stellar evolution, and hybrid close-encounter regimes — the
domains of REBOUND [@ReinLiu2012], MERCURIUS [@ReinTamayoHernandezPapaloizou2019],
and NBODY6/7 [@Aarseth2003] — remain outside the library's claims.

# Statement of need

Extension mechanisms are a central design feature of gravitational
N-body codes. A base integrator is augmented with conservative
corrections — general-relativistic precession, J2 oblateness, tidal
dissipation — or with non-gravitational forces such as radiation
pressure and gas drag. Each extension carries implicit preconditions
about the base integrator: the softening model, the force-determinism
guarantee, the units of `G`, `c`, and `M`. When those preconditions
are violated, the integrator reports no error and continues to satisfy
conservation invariants to machine precision. The only signal that
something is wrong is a quantitative comparison against an analytic
reference — the step a researcher is most likely to skip when every
other indicator reports health.

These mechanisms are well-established in the N-body literature.
REBOUNDx [@Tamayo2020] is the canonical example, adding conservative
and dissipative forces to the symplectic integrations that REBOUND
[@ReinLiu2012] produces; similar extension patterns exist in MERCURIUS
[@ReinTamayoHernandezPapaloizou2019] and NBODY6/7 [@Aarseth2003]. A
concrete instance sharpens the failure mode: a first-post-Newtonian
correction implicitly assumes exact `1/r` gravity, so if the base
simulation applies Plummer softening — common for numerical stability
— the softening produces a numerical apsidal precession that for
Sun–Mercury parameters exceeds the relativistic effect by three orders
of magnitude with the wrong sign, while energy and angular-momentum
conservation remain satisfied to machine precision.

`apsis` promotes this class of precondition from prose to the type
level. Extension points declare their physical assumptions as methods
on the `PerturbationForce` trait; registering an extension whose
assumptions are not satisfied by the current system emits a
structured diagnostic event with per-body statistics identifying the
violating bodies. A second design commitment makes the extension
surface a *buildable* contract rather than a documented one:
extensions live in independent Cargo crates that depend only on the
library's published interface, and their compilation runs as a
continuous-integration gate. These two properties — type-expressed
preconditions and out-of-tree verified extensions — are not, to the
author's knowledge, combined in any existing N-body code.

The contribution is to the *methodology* of extending an N-body
simulator rather than to the inventory of simulators. A research
group already running REBOUND, MERCURIUS, or an equivalent production
code is not served by replacing it with `apsis`; the claim of `apsis`
is orthogonal to the claim those codes make. The narrow scope (2D,
N ≤ 10³) is a deliberate trade: ship a verification infrastructure
with a complete physical demonstration, rather than a wider
simulation platform with verification deferred to later work.

# Design and validation

The library rests on two design commitments. First, the physical
preconditions of a perturbation force are part of that force's *type*,
not of its prose documentation. The public trait `PerturbationForce`
carries an optional method `requires_exact_gravity()` whose return
value tells the library whether the Newtonian base must be a bit-exact
`1/r` potential for the perturbation to produce a meaningful result.
Second, the public API boundary is a *buildable* contract rather than
a documented one. The companion crate `apsis-1pn` lives beside the
library in the Cargo workspace and imports `apsis` only through its
published interface, with no access to `pub(crate)` internals. The
consequence is that any change to `apsis` that would break an external
consumer's compilation fails the continuous-integration build of
`apsis-1pn`, not a manual review.

The precondition mechanism surfaces through the library's structured
diagnostic channel. When `System::add_perturbation(force)` is invoked and
`force.requires_exact_gravity()` returns `true`, the library counts
registered bodies with non-zero Plummer softening and — if any are
found — emits a warning that names the largest softening length in
use. The caller dismisses the warning by invoking `Body::unsoftened()`
per body or `System::with_exact_gravity()` for the whole system; both
are no-ops when already satisfied and therefore safe to include
unconditionally in research scripts.

The mechanism is validated against a configuration chosen for having
both a closed-form relativistic prediction and a sharp counter-test:
the 1PN Schwarzschild correction applied to Sun–Mercury. The
`apsis-1pn` crate integrates this configuration for 500 orbital
periods under IAS15 [@ReinSpiegel2015] and compares the accumulated
longitude of periastron against the closed form `6πGM/(c²a(1-e²))`
per orbit [@Will1993]. With both bodies unsoftened — the type-level
precondition satisfied — the drift is 42.983 arcseconds per century
against the predicted 43.000, a relative error of 4.4×10⁻⁶, stable
over the integration window and monotonic in step count. With the
library's default Plummer softening left in place — the precondition
violated — the drift is −83 128 arcseconds per century, three orders
of magnitude larger than the relativistic effect and of the wrong
sign, while energy and angular momentum are conserved to machine
precision throughout. The first result establishes that the integrator
resolves 1PN-scale effects at the accuracy the mechanism requires; the
second establishes that the type-level precondition catches a real
and severe failure mode rather than a cosmetic check. Both are
asserted as continuous-integration gates — 1 % tolerance on the GR
agreement, a non-negotiable warning emission on the counter-test —
and the full suite completes in under ten seconds on commodity
hardware.

# Availability and reproducibility

`apsis` is available under the Apache License 2.0 at
<https://github.com/gabrielbragaestefanski/apsis>. The Mercury
validation described above reproduces on a clean checkout with a
single command,

```bash
cargo test --release -p apsis-1pn --tests -- --ignored
```

after installing a Rust 1.85+ toolchain. The continuous-integration
configuration additionally compiles every example crate, rejects
warnings under `cargo clippy --all-targets`, and verifies that the
library crate resolves no user-interface dependency. A pinned snapshot
of the source archive corresponding to this paper is deposited at
Zenodo (DOI forthcoming).

# Acknowledgements

I thank the authors of REBOUND, REBOUNDx, MERCURIUS, and NBODY6/7 for
setting the standards of rigour against which this library's narrower
claim is positioned.

# References
