---
title: 'APSIS: A Federated Model for Composable N-Body Force Artifacts'
tags:
  - Rust
  - N-body
  - gravitational dynamics
  - general relativity
  - simulation software
authors:
  - name: Gabriel Braga Estefanski
    orcid: 0009-0009-1041-2715
    affiliation: 1
affiliations:
 - name: Independent researcher
   index: 1
author: "Gabriel Braga Estefanski^[Independent researcher. ORCID: 0009-0009-1041-2715]"
date: 22 May 2026
bibliography: paper.bib
geometry:
  - margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{microtype}
  - \usepackage{float}
  - \floatplacement{figure}{!ht}
  - \sloppy
---

# Summary

`apsis` is a Rust library implementing a *federated perturbation
model* for gravitational N-body simulation: each force perturbation
is published as an independent Cargo crate, versioned, citable, and
composed into a simulation through the library's public extension
API. A simulation's full physical model is captured as a
`Cargo.lock` file — reproducible bit-for-bit at the force-composition
level. The library is the runtime for composing force artifacts; it
is not the artifact.

Every force crate carries a contract surface. A perturbation
declares its physical preconditions on the gravitational kernel —
exact `1/r` versus softened, smooth versus $C^0$ — as type-level
`KernelRequirements`; the library matches these against the active
kernel's properties at registration and emits a structured
diagnostic for every violated invariant. Two formally distinct
invariants (Exactness, Continuity), when violated, produce two
formally distinct and quantitatively separable observables, each
caught independently by the registration check; the contract is
compositional, not specialised to softening.

The mechanism is demonstrated end-to-end by `apsis-1pn`, an
out-of-tree crate implementing the first-post-Newtonian
Schwarzschild correction. With the contract enforced, `apsis-1pn`
reproduces Mercury's perihelion precession to within $2.8 \times
10^{-5}$ (28 ppm) of the closed-form general-relativistic
prediction over 500 orbital periods under the adaptive Gauss–Radau
IAS15 scheme [@ReinSpiegel2015], reproduced bit-identically across
Windows and Linux on x86_64. With the Plummer softening contract
violated, the same machinery measures a precession three orders of
magnitude larger and of the wrong sign — caught by the registration
warning, never as a numerical artifact. The federation architecture
is additionally exercised with validated radiation [@Burns1979] and
central-force [@Tamayo2020] modules; the contract surface is not
specialised to Hamiltonian gravitational perturbations.

The solver provides Velocity Verlet, Yoshida fourth-order,
Wisdom–Holman in democratic-heliocentric coordinates
[@WisdomHolman1991] (uncorrected leapfrog with Kepler drifts; the
order-17 symplectic corrector of [@Wisdom2006] is tracked as future
work), implicit midpoint, the MERCURIUS hybrid
[@ReinTamayoHernandezPapaloizou2019], and the adaptive Gauss–Radau
IAS15 scheme [@ReinSpiegel2015], alongside stable public traits for
user-registered force models and perturbations.
Scope is narrow by intent: $N \le 10^3$ in the current validated
regime. Large-N collisionless dynamics, stellar evolution, and
dense close-encounter regimes — the domains of REBOUND
[@ReinLiu2012] and NBODY6/7 [@Aarseth2003] — remain outside this
library's claims. `apsis` does not attempt to replace mature
integrators or optimize numerical performance; its contribution is
orthogonal: defining how physical models are structured, published,
and composed.

The full simulation, from physical model (`Cargo.lock`-pinned
operator crates) to numerical output (an *Apsis Record*), is
bit-exactly reproducible from a single hash-pinned configuration.
This claim is operationalised by an empirical cross-platform test
over the full v0.1 federation portfolio (see §Cross-platform
reproducibility).

# Statement of need

In current practice, a force perturbation in a published N-body
simulation lives in the methods section of a paper and, sometimes,
in a fork of an established framework. The fork is not a citable
artifact of its own, the prose drifts as it is restated by
subsequent work, and the next group reimplements the same effect
from scratch. The framework — REBOUND [@ReinLiu2012], REBOUNDx
[@Tamayo2020], MERCURIUS [@ReinTamayoHernandezPapaloizou2019],
NBODY6/7 [@Aarseth2003] — is mature, citable, and validated, but
it absorbs every extension into a single binary with one citation
covering everything.

Each extension carries implicit preconditions about the base
integrator: the softening model, the force-determinism guarantee,
the units of `G`, `c`, and `M`. When those preconditions are
violated, the integrator reports no error and continues to satisfy
conservation invariants to machine precision. The only signal that
something is wrong is a quantitative comparison against an analytic
reference — the step a researcher is most likely to skip when every
other indicator reports health. A concrete instance sharpens the
failure mode: a first-post-Newtonian correction implicitly assumes
exact `1/r` gravity, so if the base simulation applies Plummer
softening — common for numerical stability — the softening produces
a numerical apsidal precession that for Sun–Mercury parameters
exceeds the relativistic effect by three orders of magnitude with
the wrong sign, while energy and angular-momentum conservation
remain satisfied to machine precision.

`apsis` replaces this publication path with a *federated
perturbation model*. A force is a Cargo crate that declares its
physical preconditions on the gravitational kernel through the
`KernelRequirements` type — `apsis-1pn` declares
`exact_and_smooth()`; future crates declare a different combination
of exactness and continuity invariants depending on the physics.
The library matches the declared requirements against the active
kernel at `System::add_hamiltonian_perturbation` (or
`add_non_conservative_perturbation`) and emits a structured
diagnostic for each violated invariant. Forgetting a precondition
surfaces as a registration warning, not as a wrong number in a
paper.

A simulation's physical model is therefore not embedded in code,
but in its dependency graph: `Cargo.toml` declares the forces a
paper uses, `Cargo.lock` pins them bit-precisely. Composition at
the call site is a single registration per operator:

```rust
let mut sys = System::new(bodies, UnitSystem::solar_canonical())
    .with_integrator(IntegratorKind::IAS15);

let gr = PostNewtonian1PN::for_units(UnitSystem::solar_canonical());
sys.add_hamiltonian_perturbation(gr)?;            // apsis-1pn
sys.add_non_conservative_perturbation(drag)?;     // apsis-radiation
```

A follow-up paper extending the model adds one line plus one entry
in `Cargo.toml`. This is reproducibility at the force-composition
level, distinct from script-level reproducibility — the latter
captures the configuration but not the physics implementation.

The contribution is to the *methodology* of extending an N-body
simulator rather than to the inventory of simulators. A research
group already running REBOUND, MERCURIUS, or an equivalent
production code is not served by replacing it with `apsis`. The
narrow scope ($N \le 10^3$ in the validated regime) is a deliberate
trade: ship a verification infrastructure with a complete physical
demonstration,
rather than a wider simulation platform with verification deferred
to later work. These two properties — type-expressed preconditions
and out-of-tree verified federated extensions — are not, to the
author's knowledge, combined in any existing N-body code.

# Design and validation

The library rests on two design commitments. First, the physical
preconditions of an extension are part of that extension's *type*, not
of its prose documentation. Two extension points exercise this pattern:
an operator declares, via `Operator::kernel_requirements`, the
invariants the gravitational kernel must satisfy for the operator's
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
the pair potential

$$U_{ij} = -G \cdot m_i \cdot m_j \cdot K(r), \qquad r = |x_i - x_j|.$$

The library encodes two formal invariants of $K$. **Exactness**: a
kernel is *Exact* if $K(r) = 1/r$, *Softened* if

$$K(r) = \frac{1}{\sqrt{r^2 + \varepsilon^2}}$$

with non-trivial $\varepsilon$, and *Modified* otherwise. **Continuity**:
a kernel is in $C^n$ if the force $-dK/dr$ belongs to $C^n(\mathbb{R}_+)$,
and *Smooth* if $C^\infty$. A
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
channel. When `System::add_hamiltonian_perturbation(operator)` is
invoked (or its non-conservative counterpart), the active kernel's
properties are matched field-by-field against
`operator.kernel_requirements()`; every invariant violation emits a
`warn_diag!` event naming the specific invariant, the value required,
and the value provided. The default kernel is `NewtonKernel::exact()`
($\varepsilon = 0$), which reports `Exactness::Exact`, so a correctly
configured run stays silent. Cluster-scale work that opts into a
softened kernel via `System::with_kernel(Arc::new(NewtonKernel::new(eps)))`
with $\varepsilon > 0$ triggers the diagnostic at registration of any
Exactness-requiring operator.

Two counter-tests exercise the two invariants separately. The **Exactness**
counter-test is the Sun–Mercury configuration integrated for 500
orbital periods under the adaptive Gauss–Radau IAS15 scheme
[@ReinSpiegel2015]. Under the default `NewtonKernel::exact()`
($\varepsilon = 0$) — Exactness satisfied — the measured cumulative
perihelion advance is
51.7705 arcsec, matching the closed-form general-relativistic
prediction $\Delta\omega_{\text{orbit}} = 6\pi GM / (c^2 a (1 - e^2))$
summed over 500 orbits (51.7720 arcsec in canonical f64 evaluation)
within $2.8 \times 10^{-5}$ relative agreement. The per-century rate
is 42.991 arcsec; the historical 43 arcsec/century [@Will1993] is
matched to four significant figures. The measurement reproduces
bit-identically across Windows and Linux on x86_64; the
hardware-specific reproducibility detail is in §Cross-platform
reproducibility.

Extending the same scenario to 4153 orbits (figure below) shows
the cumulative perihelion advance tracking the GR prediction
linearly across a 1000-year horizon. The per-orbit precession rate
inherits the f64 round-off floor of the integrator; the residual at
the end of the run sits at $2.2 \times 10^{-4}$ relative to the
predicted total, consistent with the per-orbit precision reported in
§Cross-platform reproducibility scaled to the longer horizon.

![Mercury perihelion precession under `apsis` (IAS15 + apsis-1pn) versus the closed-form Schwarzschild GR prediction over 4153 orbits ($\sim 1000$ years). Top: cumulative $\Delta\omega$ as a function of orbit number; the measured trajectory (solid) is visually indistinguishable from the GR prediction (dashed) on this scale. Bottom: the residual measured $-$ GR, showing a linear secular drift consistent with the per-orbit precision floor reported in §Cross-platform reproducibility.](paper/figures/mercury_1pn_long_horizon.pdf){#fig:mercury-1pn-long-horizon width=85%}

With a Plummer-softened kernel (`NewtonKernel::new(eps)` with
$\varepsilon \approx 0.02$ AU, the cluster-scale softening for a
solar-mass body) opted in — Exactness
violated — the drift is $-83\,128$ arcseconds per century: three
orders of magnitude larger than the relativistic effect and of the
wrong sign, while energy and angular momentum remain conserved to
machine precision throughout.

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
continuity measurement, and warning emission required on both
registrations. The full suite completes in under twenty seconds on
a 2024-class x64 workstation.

**Federation evidence beyond gravity.** Two additional first-party
crates exercise the federation API on physics distinct from the
Hamiltonian gravitational regime of `apsis-1pn`. `apsis-radiation`
implements radiation-pressure and Poynting–Robertson forces per
[@Burns1979] as a `NonConservativeOperator`, validated against the
analytic $\beta$-parameter decay law for a Sun-orbiting dust grain
(`crates/apsis-radiation/tests/dust_decay_gate.rs`).
`apsis-central` implements the [@Tamayo2020] central-force
perturbation with both a regime-based forward constructor and an
observable-inversion constructor (`from_apsidal_rate`) that selects
the coupling amplitude such that the resulting precession reproduces
a target observable apsidal rate
(`crates/apsis-central/tests/round_trip_gate.rs`). Together with
`apsis-1pn`, the three operator crates exercise both operator
categories defined by the contract surface (`HamiltonianOperator`
for `apsis-1pn` and `apsis-central`, `NonConservativeOperator` for
`apsis-radiation`) and both supported constructor patterns
(regime-based for `apsis-1pn` and `apsis-radiation`, observable-
inversion via `from_apsidal_rate` for `apsis-central`). Each crate
is an independent Cargo artifact, pinned by version in the
simulation's `Cargo.lock` and matched against the active kernel at
registration through the same `KernelRequirements` machinery.

**Executable contract surface.** The kernel-precondition mechanism
demonstrated above is one of three guarantee classes the library
publishes to a perturbation author. The full surface is formalised in
`apsis::contract`: every guarantee is named in the module-level
documentation, every guarantee is gated by a continuous-integration
test whose name matches the guarantee, and every test is co-located
with the prose. Reading the module top-to-bottom reads the contract;
running `cargo test -p apsis --lib contract` verifies it. The library
distinguishes itself from comparable surfaces in REBOUND/REBOUNDx
[@ReinLiu2012; @Tamayo2020] not on test count — REBOUND has a wider
validation portfolio measured by problem count — but on **shape**:
that a reviewer can mechanically check the claims the contract makes.

The three classes are:

*Kernel invariants* — the simulation is deterministic at the system
level (the bare integrator and any registered perturbations together);
attaching a no-op perturbation produces a trajectory bit-equal to the
bare-Newton run; perturbation evaluation is a pure function of
`(bodies, scratch_acc)`. Four tests, including a negative test that
proves the determinism check observes trajectory state rather than
returning a fixed value.

*Composition rules* — registration is commutative at the IEEE-754
accumulator step for $N = 2$; associative within the IEEE-754 summation
envelope for $N \ge 3$; additive (perturbations contribute by `+=`,
never overwrite, verified by sentinel pre-population of the
accumulator); the system's effective `KernelRequirements` is the
set-union of the individual perturbations'. Four tests. The
trajectory-level corollary of associativity (registering $[A, B, C]$
versus $[C, B, A]$ produces equivalent science) is not asserted at the
contract level — adaptive integrators amplify ULP-level acceleration
differences through chaotic substep selection, so a trajectory-level
gate would measure integrator behaviour rather than the composition
operator. The associativity claim therefore lives at the per-call
accumulator level, where the IEEE-754 statement is well-defined.

*Failure model* — exactly one warning per violated invariant per
registration; repeated registration produces a faithful audit trail
rather than silent coalescing; emission is unconditional on subscriber
state, so a registration with no consumer attached completes normally
and a subsequent subscriber-attached registration still observes the
warning. Four tests. The two demonstrated counter-tests above are
specific instances of the first guarantee (one Exactness diagnostic on
softened-kernel violation; one Continuity diagnostic on truncated-Plummer
violation); the remaining tests pin the audit-trail and
no-silent-acceptance properties that the kernel-precondition
demonstration alone would not exhibit.

Twelve tests in total at `crates/apsis/src/contract.rs`. The same
file holds the prose statement of every guarantee, the rationale for
the invariants the contract does *not* extend to (cross-platform
bit-exactness, cross-thread determinism, build-flag invariance), and
the load-bearing iteration-order property of the perturbation storage
that a future refactor must not break.

**Run configuration.** All measurements correspond to: IAS15 with
initial timestep $10^{-4} \cdot T$ and adaptivity enabled for the
Exactness counter-test (Sun–Mercury standard orbital elements,
$\varepsilon = 0$ for the satisfied case, $\varepsilon \approx 0.02$ AU
for the violated case, 500-period integration); fourth-order Yoshida
at fixed $dt = 10^{-3} \cdot T$ for the Continuity counter-test
(equal-mass two-body $a = 1$, $e = 0.5$, default exact kernel,
$R_c = 1$, $\alpha = 0.8$, 60 simulation-unit integration). Sources
at `crates/apsis-1pn/tests/mercury_precession_gate.rs` and
`crates/apsis-1pn/tests/kernel_continuity_counter_test.rs`; both
reproduce on a clean checkout per the §Availability command. The
twelve composition-contract tests run under
`cargo test -p apsis --lib contract` and live in
`crates/apsis/src/contract.rs`.

Two formally distinct invariants (Exactness, Continuity), when
violated, produce two formally distinct and quantitatively separable
observable signatures, each caught independently by the registration
check. This — not empirical superiority in any numerical regime — is
the claim the mechanism supports.

# Cross-platform reproducibility

The `Cargo.lock`-as-experiment claim of §Summary is operationalised
by an empirical cross-platform test, distinct from the formal
contract surface of §Design (which is scoped to single-host
invariants per `crates/apsis/src/contract.rs`). With the lockfile,
the rustc toolchain version, and the source commit pinned, the full
v0.1 federation portfolio — `apsis-1pn` Mercury long-horizon under
IAS15, the MERCURIUS hybrid on a Sun + four outer planets
configuration with a Jupiter-crossing test particle, Wisdom–Holman
on the same outer-planets initial conditions, the `apsis-central`
round-trip gate, and the four-scenario REBOUND parity portfolio
(Kepler, figure-8, Pythagorean, retrograde) — produces byte-
identical trajectory output on heterogeneous x86_64 hosts (Windows
on AMD Zen 4 against Linux on Intel Ice Lake). Verification covers
per-column ULP agreement, file size, and SHA256 of the captured
trajectory CSVs.

The mechanism is the routing of every libc-bound transcendental on
an integration-critical path (`sin`, `cos`, `cbrt`, `pow`, `log`,
...) through the pure-Rust `libm` crate. Hardware-rounded arithmetic
($+$, $-$, $\times$, $\div$, $\sqrt{}$) is bit-identical across IEEE
754-conformant x86_64 implementations by mandate and does not
require replacement. The IAS15 step-size controller's `pow` call was
measured against an IEEE-754 correctly-rounded mpmath oracle on the
42,662 unique inputs generated during the Mercury 1PN run: Windows
UCRT and the `libm` crate match the oracle on 96.97 % and 95.29 %
of inputs respectively, both within IEEE-754 tolerance for
transcendentals (which permits but does not require correctly-
rounded results). The 1-ULP differences in the remaining cases
propagate through the controller's substep-cadence selection to a
0.002 arcsec/century absolute shift in cumulative Mercury $\Delta\omega$ over
500 orbits, separating the $4.4 \times 10^{-6}$ relative agreement
on Windows UCRT (scenario-specific accidental cancellation in
UCRT's rounding distribution) from the $2.8 \times 10^{-5}$ result
reproduced bit-identically by `libm` and glibc. Both values sit
several orders of magnitude below the current observational
precision of Mercury's perihelion advance [@VermaFienga2014] and
below the experimental constraint on the PPN $\gamma$ parameter
from Cassini Doppler tracking [@BertottiIessTortora2003]. The
choice prioritises trajectory determinism across deployment
platforms over the scenario-specific rounding alignment that the
Windows-only result reflected.

The methodology and per-implementation analysis are recorded in
`paper/notebooks/2026-05-20-cross-platform-determinism.md` (federation
portfolio bit-equality protocol) and
`paper/notebooks/2026-05-22-controller-pow-implementations.md`
(controller `pow` ULP-distribution analysis).

# Reproducibility certificate

Each apsis simulation emits an *Apsis Record* — a binary certificate
documenting the run's full provenance and physical state. The record
contains the apsis-core git commit, a BLAKE3 hash of the project's
`Cargo.lock`, the unit system, integrator configuration, kernel
variant (and softening parameter when applicable), seed, every
registered operator with its crate name + version + lockfile hash +
declared `KernelRequirements`, and bookend snapshots of body state.
Material physical events (collisions, escapes) are recorded inline.
The certificate is a single binary file with a human-readable TOML
header and a binary frame stream.

The federation thesis — that a simulation's physical model is
captured in `Cargo.lock` — is here extended to the run itself:
`{record, Cargo.lock}` is the content-addressable closure of the
experiment. The record's frame stream and the trailer's BLAKE3 are
bit-exactly reproducible across replays with the same configuration;
the only per-run difference is the header's wall-clock
`created_utc` field, which is metadata and excluded from the content
hash. A reviewer with both files reproduces the run. The default
policy emits initial + final bookend snapshots and every material
event, which keeps records small and diff-friendly; dense trajectory
capture is an explicit policy opt-in. The bit-equal trajectory
reproduction demonstrated in §Cross-platform reproducibility applies
to the data the Record wraps; the binary file therefore reproduces
across the same heterogeneous hosts modulo the `created_utc`
metadata field.

# Availability and reproducibility

`apsis` is available under the Apache License 2.0 at
<https://github.com/GabrielEstefanski/apsis>. The Mercury
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

The author used an AI assistant for documentation drafting, code
review, refactoring, and design-discussion support. All algorithmic
decisions, physics interpretations, validation methodology, numerical
implementations, and reported results are the author's responsibility.
The AI assistant was not used for novel physics, mathematical
derivations, or scientific decision-making.

# References
