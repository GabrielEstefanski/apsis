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

`apsis` is a Rust library for gravitational N-body simulation. The
library's public extension API promotes physical preconditions of
perturbation forces — for example, whether a correction assumes an
unsoftened `1/r` potential — from informal documentation to type-level
declarations checked at extension registration and enforced in continuous
integration. The mechanism is demonstrated by an out-of-tree companion
crate, `apsis-1pn`, which implements the first-post-Newtonian
Schwarzschild correction and reproduces Mercury's perihelion precession
within 4.4×10⁻⁶ of the general-relativistic prediction over 500 orbital
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

Mature N-body simulation software in planetary-scale astrophysics
exposes an extension mechanism for augmenting the Newtonian base physics
with additional effects. REBOUNDx [@Tamayo2020] is the canonical
example, adding conservative and dissipative forces to the symplectic
integrations produced by REBOUND; similar extension points exist in most
N-body codes used in production research.

The failure mode that `apsis` addresses is narrower than any of those
tools sets out to solve. In existing extension mechanisms the physical
preconditions a correction implicitly assumes about the base integrator
— softening model, force determinism, the particular units of `G`, `c`,
and `M` — are carried in prose documentation and enforced, when at all,
by runtime checks at registration. A misapplication is therefore silent:
the simulation runs, the integrator reports no error, and standard
conservation diagnostics (energy, angular momentum) remain satisfied.

A concrete case motivates the mechanism described here. An extension
computing a first-post-Newtonian correction implicitly assumes an exact
`1/r` gravitational potential. If the base simulation applies Plummer
softening — common for numerical stability with arbitrary initial
conditions — the softening introduces a numerical apsidal precession.
For Mercury-scale parameters, with a default softening length that is
physically reasonable elsewhere, that numerical precession exceeds the
relativistic effect by roughly three orders of magnitude and has the
wrong sign. The only upstream signal is a quantitative comparison
against an analytic prediction, which is precisely the step a user is
likely to skip when every other indicator reports health.

`apsis` promotes this class of precondition to the type level. Extension
points declare their physical assumptions as trait methods; registering
an extension whose assumptions are not satisfied by the current system
emits a structured diagnostic event through the library's log bus, with
per-body statistics identifying the specific bodies that violate the
contract. The companion crate `apsis-1pn` is built, tested, and CI-gated
as an independent Cargo crate that depends only on the library's
published interface; the public API boundary is therefore enforced by
the compiler and the continuous-integration build rather than by review.
These two properties — compile-visible preconditions and out-of-tree
verified extensions — are not, to the authors' knowledge, combined in
any existing N-body library.
