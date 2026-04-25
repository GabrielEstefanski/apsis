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
potential or a smooth Hamiltonian required by symplectic integration —
are promoted from informal documentation to type-level declarations,
checked at extension registration and enforced through an out-of-tree
companion crate whose compilation runs as a continuous-integration
gate. The mechanism is demonstrated end-to-end by `apsis-1pn`, which
implements the first-post-Newtonian Schwarzschild correction; as
evidence that the integrator stack resolves 1PN-scale effects at the
accuracy the verification claim requires, the demonstration reproduces
Mercury's perihelion precession to within $4.4 \times 10^{-6}$ of the
general-relativistic prediction over 500 orbital periods. The same
match mechanism catches a distinct invariant-class violation (a
discontinuous kernel) with an independent observable (impulsive
energy-error events in one-to-one correspondence with cutoff-radius
crossings), demonstrating that the contract is compositional rather
than specialised to softening.

The solver provides four integration schemes — second-order Velocity
Verlet, fourth-order Yoshida composition, Wisdom–Holman mixed-variable,
and the adaptive Gauss–Radau IAS15 scheme [@ReinSpiegel2015] — alongside
stable public traits for user-registered force models and perturbations.
The library's scope is narrow by intent: the solver is two-dimensional
and targets small-to-medium body counts ($N \le 10^3$). Large-N
collisionless dynamics, stellar evolution, and hybrid close-encounter
regimes — the domains of REBOUND [@ReinLiu2012], MERCURIUS
[@ReinTamayoHernandezPapaloizou2019], and NBODY6/7 [@Aarseth2003] —
remain outside the library's claims.

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
$N \le 10^3$) is a deliberate trade: ship a verification infrastructure
with a complete physical demonstration, rather than a wider
simulation platform with verification deferred to later work.

# Design and validation

The library rests on two design commitments. First, the physical
preconditions of an extension are part of that extension's *type*, not
of its prose documentation. Two extension points exercise this pattern:
a perturbation force declares, via `PerturbationForce::kernel_requirements`,
the invariants the gravitational kernel must satisfy for the perturbation's
derivation to be meaningful; a kernel implementation declares, via
`Kernel::properties`, the invariants it in fact satisfies for the current
bodies. Second, the public API boundary is a *buildable* contract rather
than a documented one. The companion crate `apsis-1pn` lives beside the
library in the Cargo workspace and imports `apsis` only through its
published interface, with no access to `pub(crate)` internals. The
consequence is that any change to `apsis` that would break an external
consumer's compilation fails the continuous-integration build of
`apsis-1pn`, not a manual review.

Let $K: \mathbb{R}_+ \to \mathbb{R}_+$ denote the scalar kernel determining
the pair potential $U_{ij} = -G \cdot m_i \cdot m_j \cdot K(r)$, where
$r = |x_i - x_j|$. The library encodes two formal invariants of $K$.
**Exactness**: a kernel is *Exact* if $K(r) = 1/r$, *Softened* if
$K(r) = 1/\sqrt{r^2 + \varepsilon^2}$ with non-trivial $\varepsilon$, and
*Modified* otherwise. **Continuity**: a kernel is in $C^n$ if the force
$-dK/dr$ belongs to $C^n(\mathbb{R}_+)$, and *Smooth* if $C^\infty$. A
perturbation declares the minimum invariants it requires (typed as
`KernelRequirements`); a kernel implementation declares the invariants
it provides for the current body configuration (typed as
`KernelProperties`); a mismatch on any field is identified at extension
registration.

These invariants are not ad-hoc labels. Exactness violation is a
statement about the derivation: the 1PN correction is obtained by
expanding the geodesic equation around the Newtonian Hamiltonian
$H_N = p^2/2m - GMm/r$, and substituting a softened potential invalidates
the expansion itself — the observed apsidal drift is the signature of
applying a Taylor series on top of a different unperturbed system, not
a numerical artifact. Continuity violation is a statement about
phase-space geometry: symplectic integration relies on the Hamiltonian
flow preserving phase-space volume, which requires a smooth $H$; force
discontinuities produce impulsive accelerations that cannot be
represented within any symplectic splitting scheme, independent of
integrator order or step control.

The mechanism surfaces through the library's structured diagnostic
channel. When `System::add_perturbation(force)` is invoked, the active
kernel's properties are computed from the current bodies and matched
field-by-field against `force.kernel_requirements()`; every invariant
violation emits a `warn_diag!` event naming the specific invariant,
the value required, and the value provided. A Plummer kernel with
every body `.unsoftened()` reports Exactness::Exact dynamically, so
a correctly configured run stays silent. `System::with_exact_gravity()`
and per-body `Body::unsoftened()` are idempotent helpers safe to
include unconditionally in research scripts.

Two counter-tests exercise the two invariants separately. The **Exactness**
counter-test is the Sun–Mercury configuration integrated for 500
orbital periods under the adaptive Gauss–Radau IAS15 scheme
[@ReinSpiegel2015]. With both bodies unsoftened — Exactness satisfied —
the accumulated longitude of periastron drifts by 42.983 arcseconds
per century against the closed-form general-relativistic prediction
$6\pi GM / (c^2 a (1 - e^2))$ = 43.000 arcseconds per century
[@Will1993], a relative error of $4.4 \times 10^{-6}$, stable over the
integration window and monotonic in step count. With the library's
default Plummer softening left in place — Exactness violated — the
drift is $-83\,128$ arcseconds per century: three orders of magnitude
larger than the relativistic effect and of the wrong sign, while energy
and angular momentum remain conserved to machine precision throughout.

The **Continuity** counter-test is a distinct configuration designed
to exercise the second invariant on a distinct observable. An
equal-mass two-body orbit ($a = 1$, $e = 0.5$) is integrated under a
truncated-Plummer kernel that matches the standard Plummer profile
inside a cutoff radius $R_c = 1$ (semi-major-axis units) and switches
to a scaled Plummer outside, with the outside scale $\alpha = 0.8$
chosen so that $K$ is continuous at $R_c$, the force $-dK/dr$ has a
finite jump of $(1 - \alpha) \cdot R_c / (R_c^2 + \varepsilon^2)^{3/2} = 0.2$
there, and the trajectory remains reliably bound (the orbit's apoapse
sits near $r \approx 2.06$, well inside the marginal-binding threshold
at $\alpha \approx 0.5$ for these parameters). Under fourth-order
Yoshida composition at fixed timestep $10^{-3} \cdot T$, the orbit
crosses $R_c$ eleven times over 60 simulation units and the integrator
produces impulsive energy-error events of magnitude $4.7 \times 10^{-6}$
to $2.0 \times 10^{-4}$ — in one-to-one correspondence with the
crossings, each event temporally matched to its crossing within
$10 \cdot dt$, and no events between crossings. A reference run with
the smooth PlummerKernel on the same bodies exhibits no events above
$2.7 \times 10^{-14}$ per step, separating the Continuity signature
from the symplectic round-off floor by roughly eight orders of
magnitude. The observed signature is a consequence of the continuity
violation itself, not of the specific `TruncatedPlummerKernel` used
to exhibit it: any kernel whose declared properties include
`Continuity::C0` and whose orbital configuration places the
discontinuity within the radial range of the trajectory produces the
same class of observable.

Both counter-tests are asserted as continuous-integration gates — 1 %
relative-error tolerance on the GR agreement, exact bijection between
crossing and spike events with $10 \cdot dt$ temporal matching on the
continuity measurement, and non-negotiable warning-emission on both
registrations. The full suite completes in under twenty seconds on
commodity hardware.

**Reproducibility.** All measurements correspond to: IAS15 with
initial timestep $10^{-4} \cdot T$ and adaptivity enabled for the
Exactness counter-test (Sun–Mercury standard orbital elements,
$\varepsilon = 0$ for the satisfied case, $\varepsilon \approx 0.02$ AU
for the violated case, 500-period integration); fourth-order Yoshida
at fixed $dt = 10^{-3} \cdot T$ for the Continuity counter-test
(equal-mass two-body $a = 1$, $e = 0.5$, both bodies unsoftened,
$R_c = 1$, $\alpha = 0.8$, 60 simulation-unit integration). Sources
at `crates/apsis-1pn/tests/mercury_precession_gate.rs` and
`crates/apsis-1pn/tests/kernel_continuity_counter_test.rs`; both
reproduce on a clean checkout per the §Availability command.

Two formally distinct invariants (Exactness, Continuity), when
violated, produce two formally distinct and quantitatively separable
observable signatures, each caught independently by the registration
check. This — not empirical superiority in any numerical regime — is
the claim the mechanism supports.

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
