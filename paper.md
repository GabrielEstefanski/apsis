---
title: 'APSIS: A Federated Model for Composable N-Body Force Artifacts'
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
abstract: |
  We introduce a federated perturbation model for N-body simulation
  in which force operators are independently published, versioned
  Cargo crates that declare their physical preconditions on the
  gravitational kernel as type-level requirements. The library
  matches declared requirements against the active kernel at
  registration and emits a structured diagnostic for every violated
  invariant; the same registration step is the load-bearing
  mechanism for the contract surface and for the per-run *Apsis
  Record* that pins the operator stack to a `Cargo.lock` for end-to-
  end reproducibility. We exercise the architecture with three
  first-party operator crates — `apsis-1pn` (first-post-Newtonian
  Schwarzschild correction), `apsis-radiation` (radiation pressure
  and Poynting–Robertson drag, validated to 0.7 % of the Burns 1979
  analytic law), and `apsis-central` (parametric central-force
  precession with an observable-inversion constructor, round-trip
  agreement to 2.7 %) — spanning Hamiltonian and non-conservative
  categories and both supported constructor patterns. Among these,
  `apsis-1pn` carries the deepest demonstration of the contract
  mechanism: it reproduces Mercury's perihelion precession to
  within $2.8 \times 10^{-5}$ of
  the closed-form general-relativistic prediction over 500 orbits
  under the adaptive Gauss–Radau IAS15 scheme, reproduced bit-
  identically across Windows and Linux on x86_64; violating the
  kernel-exactness precondition with Plummer softening surfaces a
  registration warning and yields a drift more than four orders of
  magnitude larger and of the wrong sign, never as a numerical
  artifact. The
  federated model and the reproducibility certificate — not
  empirical superiority in any numerical regime — are the
  contributions the architecture supports.
geometry:
  - margin=1in
fontsize: 11pt
header-includes:
  - \usepackage{microtype}
  - \usepackage{inconsolata}
  - \usepackage{hyphenat}
  - \usepackage{float}
  - \floatplacement{figure}{!ht}
  - \setlength{\emergencystretch}{3em}
  - \sloppy
---

\newpage

# Introduction

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
the units of `G`, `c`, and `M`, and the step-control assumptions.
When any of these is violated, the integrator reports no error and
continues to satisfy conservation invariants to machine precision.
The only signal that something is wrong is a quantitative
comparison against an analytic reference — the step a researcher
is most likely to skip when every other indicator reports health.

Softening is one such precondition. It is routine in cluster-scale
dynamics (where it prevents unphysical close-encounter heating) and
absent in two-body solar-system regimes (where exact $1/r$ is
required for derivation-level consistency with most post-Newtonian
and tidal corrections), yet both regimes share the same library,
the same default integrator, and the same prose documentation
describing the operator interface. Force-determinism, unit
convention, and step-control assumptions are others, all reachable
through the same mechanism. The worked example with Sun–Mercury
parameters appears in §Results, alongside two non-gravitational
operator crates that exercise the same machinery on physics
outside the Hamiltonian softening pair.

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
demonstration, rather than a wider simulation platform with
verification deferred to later work. These two properties — type-
expressed preconditions and out-of-tree verified federated
extensions — are not, to the author's knowledge, combined in any
existing N-body code.

The remainder of the paper is organised as follows. §Methods
introduces the kernel formalism, the diagnostic mechanism that
matches `KernelRequirements` against `KernelProperties` at
extension registration, the executable contract surface that
publishes the library's compositional guarantees, and the citable
operator stack that emits a BibTeX block for the registered force
composition. §Results opens with the federation evidence across
three operator crates spanning both operator categories
(Hamiltonian and non-conservative) and both supported constructor
patterns, then exercises the contract mechanism in depth on
`apsis-1pn` (Exactness counter-test on Mercury 1PN, Continuity
counter-test on a truncated Plummer kernel), checks external-
implementation parity against REBOUND IAS15 on four canonical
scenarios, and closes with the cross-platform reproducibility
test. §Discussion covers the scope
of the cross-platform claim (including ARM64 hardware not yet
verified), solver-family guidance, and a prioritised future-work
agenda. The per-run *Apsis Record* binary format specification is
in Appendix A. Source, validation harnesses, and a pinned snapshot
are available as detailed in §Data and code availability.

# Methods {#sec:methods}

The library rests on three design commitments. First, the
physical preconditions of an extension are part of that
extension's *type*, not of its prose documentation. Two extension
points exercise this pattern: an operator declares, via
`Operator::kernel_requirements`, the invariants the gravitational
kernel must satisfy for the operator's derivation to be
meaningful; a kernel implementation declares, via
`Kernel::properties`, the invariants it in fact satisfies for the
current bodies. Second, the public API boundary is a *buildable*
contract rather than a documented one. The companion crate
`apsis-1pn` lives beside the library in the Cargo workspace and
imports `apsis` only through its published interface, with no
access to `pub(crate)` internals. The consequence is that any
change to `apsis` that would break an external consumer's
compilation fails the continuous-integration build of `apsis-1pn`,
not a manual review. Third, the registered-operator list that
drives the contract surface is also the citation surface of the
simulation: every operator with its version, DOI, and lockfile
hash is emittable as a BibTeX block from the live `System` or
from the persisted *Apsis Record*, so the act of composing the
physical model and the act of citing it share the same source of
truth.

## Kernel formalism

The library tracks two invariants of the pair potential. They are
not numerical-tolerance labels: each names a condition on the
derivation of operators that depend on it, and a kernel that fails
either invariant invalidates the perturbative or symplectic ground
the dependent operator was built on.

Let $K: \mathbb{R}_+ \to \mathbb{R}_+$ denote the scalar kernel
determining the pair potential

$$U_{ij} = -G \cdot m_i \cdot m_j \cdot K(r), \qquad r = |x_i - x_j|.$$

The library encodes two formal invariants of $K$. **Exactness**: a
kernel is *Exact* if $K(r) = 1/r$, *Softened* if

$$K(r) = \frac{1}{\sqrt{r^2 + \varepsilon^2}}$$

with non-trivial $\varepsilon$, and *Modified* otherwise.
**Continuity**: a kernel is in $C^n$ if the force $-dK/dr$ belongs
to $C^n(\mathbb{R}_+)$, and *Smooth* if $C^\infty$. A perturbation
declares the minimum invariants it requires (typed as
`KernelRequirements`); a kernel implementation declares the
invariants it provides for the current body configuration (typed
as `KernelProperties`); a mismatch on any field is identified at
extension registration.

Exactness violation is a statement about the derivation: the 1PN
correction is obtained by expanding the geodesic equation around
the Newtonian Hamiltonian
$H_N = p^2/2m - GMm/r$, and substituting a softened potential
invalidates the expansion itself — the observed apsidal drift is
the signature of applying a Taylor series on top of a different
unperturbed system, not a numerical artifact. Continuity violation
is a statement about phase-space geometry: symplectic integration
relies on the Hamiltonian flow preserving phase-space volume, which
requires a smooth $H$; force discontinuities produce impulsive
accelerations that cannot be represented within any symplectic
splitting scheme, independent of integrator order or step control.

## Diagnostic mechanism

The mechanism surfaces through the library's structured diagnostic
channel. When `System::add_hamiltonian_perturbation(operator)` is
invoked (or its non-conservative counterpart), the active kernel's
properties are matched field-by-field against
`operator.kernel_requirements()`; every invariant violation emits a
`warn_diag!` event naming the specific invariant, the value
required, and the value provided. The default kernel is
`NewtonKernel::exact()` ($\varepsilon = 0$), which reports
`Exactness::Exact`, so a correctly configured run stays silent.
Cluster-scale work that opts into a softened kernel via
`System::with_kernel(Arc::new(NewtonKernel::new(eps)))` with
$\varepsilon > 0$ triggers the diagnostic at registration of any
Exactness-requiring operator.

## Executable contract surface {#sec:executable-contract}

The kernel-precondition mechanism above is one of three guarantee
classes the library publishes to a perturbation author. The full
surface is formalised in `apsis::contract`: every guarantee is
named in the module-level documentation, every guarantee is gated
by a continuous-integration test whose name matches the guarantee,
and every test is co-located with the prose. Reading the module
top-to-bottom reads the contract; running
`cargo test -p apsis --lib contract` verifies it. The library
distinguishes itself from comparable surfaces in REBOUND/REBOUNDx
[@ReinLiu2012; @Tamayo2020] not on test count — REBOUND has a
wider validation portfolio measured by problem count — but on
**shape**: that a reviewer can mechanically check the claims the
contract makes.

The three classes are:

*Kernel invariants* — the simulation is deterministic at the system
level (the bare integrator and any registered perturbations
together); attaching a no-op perturbation produces a trajectory
bit-equal to the bare-Newton run; perturbation evaluation is a pure
function of `(bodies, scratch_acc)`. Four tests, including a
negative test that proves the determinism check observes trajectory
state rather than returning a fixed value.

*Composition rules* — registration is commutative at the IEEE-754
accumulator step for $N = 2$; associative within the IEEE-754
summation envelope for $N \ge 3$; additive (perturbations
contribute by `+=`, never overwrite, verified by sentinel pre-
population of the accumulator); the system's effective
`KernelRequirements` is the set-union of the individual
perturbations'. Four tests. The trajectory-level corollary of
associativity (registering $[A, B, C]$ versus $[C, B, A]$ produces
equivalent science) is not asserted at the contract level —
adaptive integrators amplify ULP-level acceleration differences
through chaotic substep selection, so a trajectory-level gate
would measure integrator behaviour rather than the composition
operator. The associativity claim therefore lives at the per-call
accumulator level, where the IEEE-754 statement is well-defined.

*Failure model* — exactly one warning per violated invariant per
registration; repeated registration produces a faithful audit trail
rather than silent coalescing; emission is unconditional on
subscriber state, so a registration with no consumer attached
completes normally and a subsequent subscriber-attached
registration still observes the warning. Four tests. The two
counter-tests demonstrated in §Results are specific instances of
the first guarantee (one Exactness diagnostic on softened-kernel
violation; one Continuity diagnostic on truncated-Plummer
violation); the remaining tests pin the audit-trail and no-silent-
acceptance properties that the kernel-precondition demonstration
alone would not exhibit.

Twelve tests in total at `crates/apsis/src/contract.rs`. The same
file holds the prose statement of every guarantee, the rationale
for the invariants the contract does *not* extend to (cross-
platform bit-exactness, cross-thread determinism, build-flag
invariance), and the load-bearing iteration-order property of the
perturbation storage that a future refactor must not break.

## Citable operator stack

The registered-operator list that drives the contract surface
above is also the input to `System::cite()`, a Python accessor
that emits a BibTeX block keyed by every operator crate's name,
version, source-repository URL, build commit when available, and
the `Cargo.lock` hash that pins the dependency closure. For a
simulation registering `apsis-1pn` and `apsis-radiation`, the
emitted block carries one `@software` entry per crate. The
Hamiltonian operator from `apsis-1pn` produces:

```bibtex
@software{apsis-1pn_0.1.0,
  author  = {Estefanski, G. B.},
  title   = {apsis-1pn},
  version = {0.1.0},
  commit  = {fb7218a},
  url     = {https://github.com/GabrielEstefanski/apsis},
  note    = {First-post-Newtonian Schwarzschild correction.
             Cargo.lock blake3: e3f3...e94e;
             kernel_requirements: exact_and_smooth},
}
```

followed by `apsis-radiation`, whose entry carries the same
shape but a different `kernel_requirements` tag because Burns 1979
imposes no constraint on the gravitational kernel it composes with:

```bibtex
@software{apsis-radiation_0.1.0,
  author  = {Estefanski, G. B.},
  title   = {apsis-radiation},
  version = {0.1.0},
  commit  = {fb7218a},
  url     = {https://github.com/GabrielEstefanski/apsis},
  note    = {Radiation pressure and Poynting--Robertson drag after Burns 1979.
             Cargo.lock blake3: e3f3...e94e;
             kernel_requirements: unconstrained},
}
```

The output is a self-describing citation for the simulation's full
physical model: a follow-up paper using the same operator stack
reproduces the BibTeX entries by calling `System::cite()` on a
live `System` before integration, or by reading them back from the
*Apsis Record* (Appendix A) that captures the same header. Forces
published as versioned crates are citable artifacts at the
granularity at which they are published; the cite-generator is
the operational form of that claim.

# Results {#sec:results}

## Federation evidence across operator categories

The three first-party operator crates collectively exercise both
operator categories defined by the contract surface
(`HamiltonianOperator` for `apsis-1pn` and `apsis-central`,
`NonConservativeOperator` for `apsis-radiation`) and both supported
constructor patterns (regime-based for `apsis-1pn` and
`apsis-radiation`, observable-inversion via `from_apsidal_rate`
for `apsis-central`). Each crate is an independent Cargo artifact,
pinned by version in the simulation's `Cargo.lock` and matched
against the active kernel at registration through the same
`KernelRequirements` machinery.

`apsis-radiation` implements radiation pressure and Poynting–
Robertson drag per [@Burns1979] as a `NonConservativeOperator`.
The gate (`crates/apsis-radiation/tests/dust_decay_gate.rs`)
integrates a Sun-orbiting dust grain ($\beta = 0.5$, initial
circular orbit at $r_0 = 1$ AU, IAS15 at $dt = 10^{-3}$) under
direct gravity plus radiation forces for ten orbital periods,
measures the total specific-energy change $\Delta E$ between
endpoints, and compares against the constant-$r$ analytic
prediction $\Delta E_{\text{analytic}} = -\beta\, G M v^2 m\, T /
(r^2 c)$ derived from the tangential PR force at the initial
epoch; empirical agreement is 0.7 % at the gated 5 % tolerance,
with the gate width set to absorb the $\sim 2$ % bias of the
constant-$r$ approximation as the orbit drifts $\sim 0.5$ %
inward over the window. A counter-test on the same orbit without
the operator registered confirms IAS15 conserves energy to
$10^{-12}$, isolating the measured $\Delta E$ to the operator.

`apsis-central` implements the [@Tamayo2020] central-force
perturbation $\propto 1/r^\gamma$ with two constructors: a
regime-based one parameterised by coupling amplitude and exponent,
and `from_apsidal_rate`, an observable-inversion constructor that
selects the coupling amplitude so the resulting precession
reproduces a target observable apsidal rate. The round-trip gate
(`crates/apsis-central/tests/round_trip_gate.rs`) registers a
target rate $\dot\omega_{\text{in}} = 1.5\times 10^{-3}$ rad per
Gaussian time unit at $\gamma = -3$, $e = 0.1$, integrates for
fifty orbital periods sampled phase-locked at integer multiples
of $2\pi$, fits the linear secular drift of $\omega$ from the
unwrapped phase-locked series, and compares against the input:
agreement is 2.7 % at the same 5 % regression bound, with the
documented Tamayo-formula bias ($\sim 11$ % from using the
instantaneous separation rather than the secular semi-major axis,
partially cancelled by the $(1-e^2)$ correction) dominating the
residual. A Keplerian baseline counter-test holds the apsidal
drift below $10^{-7}$ when the operator is not registered.

The two gates use disjoint physics (radiation back-reaction and
apsidal precession from a $1/r^\gamma$ central force) and
disjoint test scaffolding from `apsis-1pn`'s Mercury gate,
supporting the claim that the federation contract does not bind
to a specific operator family.

The remaining subsections demonstrate the contract mechanism in
depth on a single operator (`apsis-1pn`): the Exactness counter-
test on Mercury's perihelion precession, the Continuity counter-
test on a truncated-Plummer kernel, and the cross-platform bit-
exactness of the entire portfolio.

## Exactness counter-test: Mercury 1PN

The Exactness counter-test is the Sun–Mercury configuration
integrated for 500 orbital periods under the adaptive Gauss–Radau
IAS15 scheme [@ReinSpiegel2015]. Under the default
`NewtonKernel::exact()` ($\varepsilon = 0$) — Exactness satisfied —
the measured cumulative perihelion advance is 51.7705 arcsec,
matching the closed-form general-relativistic prediction
$\Delta\omega_{\text{orbit}} = 6\pi GM / (c^2 a (1 - e^2))$ summed
over 500 orbits (51.7720 arcsec in canonical f64 evaluation) within
$2.8 \times 10^{-5}$ relative agreement. The per-century rate is
42.991 arcsec; the historical 43 arcsec/century [@Will1993] is
matched to four significant figures. The measurement reproduces bit-
identically across Windows and Linux on x86_64; the hardware-
specific reproducibility detail is in §Cross-platform reproducibility.

Extending the same scenario to 4153 orbits (figure below) shows the
cumulative perihelion advance tracking the GR prediction linearly
across a 1000-year horizon. The per-orbit precession rate inherits
the f64 round-off floor of the integrator; the residual at the end
of the run sits at $2.2 \times 10^{-4}$ relative to the predicted
total, consistent with the per-orbit precision reported in
§Cross-platform reproducibility scaled to the longer horizon.

![Mercury perihelion precession under `apsis` (IAS15 + apsis-1pn) versus the closed-form Schwarzschild GR prediction over 4153 orbits ($\sim 1000$ years). Top: cumulative $\Delta\omega$ as a function of orbit number; the measured trajectory (solid) is visually indistinguishable from the GR prediction (dashed) on this scale. Bottom: the residual measured $-$ GR, showing a linear secular drift consistent with the per-orbit precision floor reported in §Cross-platform reproducibility.](paper/figures/mercury_1pn_long_horizon.pdf){#fig:mercury-1pn-long-horizon width=85%}

With a Plummer-softened kernel (`NewtonKernel::new(eps)` with
$\varepsilon \approx 0.02$ AU, the cluster-scale softening for a
solar-mass body) opted in — Exactness violated — the measured per-orbit
cumulative drift agrees with the closed-form softened-Plummer
apsidal-precession prediction
$\Delta\varpi_\text{orbit} = -3\pi\varepsilon^2 / [a^2(1-e^2)^2]$
to 2.55 % (5 % gated), derived from the Plummer pair potential in the
companion lab notebook and cross-checked against an independent scipy
DOP853 integration at 3.2 % residual; the 5 % envelope absorbs both
the apsis and the scipy residuals with $\sim 2\times$ margin. Scaled to
arcseconds per Earth century, the prediction is
$\dot\varpi_\text{Plummer} \approx -2.35\times 10^6$ arcsec/century —
$\sim 5\times 10^4$ times the relativistic effect and of the wrong
sign, while energy and angular momentum remain conserved to machine
precision throughout. The contract enforcement returns this
quantitatively-bracketed signature at registration time, not as a
numerical artifact emerging only under post-hoc analytic comparison.

## Continuity counter-test: TruncatedPlummer

The Continuity counter-test is a distinct configuration designed
to exercise the second invariant on a distinct observable. An
equal-mass two-body orbit ($a = 1$, $e = 0.5$) is integrated under
a truncated-Plummer kernel that matches the standard Plummer
profile inside a cutoff radius $R_c = 1$ (semi-major-axis units)
and switches to a scaled Plummer outside, with the outside scale
$\alpha = 0.8$ chosen so that $K$ is continuous at $R_c$, the
force $-dK/dr$ has a finite jump of
$(1 - \alpha) \cdot R_c / (R_c^2 + \varepsilon^2)^{3/2} = 0.2$
there, and the trajectory remains reliably bound (the orbit's
apoapse sits near $r \approx 2.06$, well inside the marginal-
binding threshold at $\alpha \approx 0.5$ for these parameters).

Under fourth-order Yoshida composition at fixed timestep
$\mathrm{d}t = 10^{-3}$ in canonical time units, the orbit crosses
$R_c$ eleven times over 60 simulation units and the integrator
produces impulsive energy-error events of magnitude
$4.7 \times 10^{-6}$ to $2.0 \times 10^{-4}$ — in one-to-one
correspondence with the crossings, each event temporally matched
to its crossing within $10 \cdot \mathrm{d}t$, and no events between
crossings. Every spike falls within the a-priori envelope
$|\Delta E|/|E_0| \le \Delta F \cdot v_\text{cross} \cdot \mathrm{d}t / |E_0|
= 4.00\times 10^{-4}$ derived from the shadow-Hamiltonian breakdown
at the discontinuity, where $|E_0| = G(m_1+m_2)/(2a) = 0.5$ is the
specific energy of the relative motion (closed form in the companion
lab notebook); the measured peak sits at 50 % of this bound, with the remaining
$2\times$ margin and the magnitude spread both reflecting the
seven Yoshida-4 substep weights, which partially cancel the
wrong-side work whenever the crossing is not centred in the step.
A reference run with the smooth PlummerKernel on the same bodies
exhibits no events above $2.7 \times 10^{-14}$ per step, separating
the Continuity signature from the symplectic round-off floor by
roughly eight orders of magnitude.

The observed signature is a consequence of the continuity
violation itself, not of the specific `TruncatedPlummerKernel`
used to exhibit it: any kernel whose declared properties include
`Continuity::C0` and whose orbital configuration places the
discontinuity within the radial range of the trajectory produces
the same class of observable.

Both counter-tests are asserted as continuous-integration gates — 1
% relative-error tolerance on the GR agreement, exact bijection
between crossing and spike events with $10 \cdot dt$ temporal
matching on the continuity measurement, and warning emission
required on both registrations. The full suite completes in under
twenty seconds on a 2024-class x64 workstation.

**Run configuration.** All measurements correspond to: IAS15 with
initial timestep $10^{-4} \cdot T$ and adaptivity enabled for the
Exactness counter-test (Sun–Mercury standard orbital elements,
$\varepsilon = 0$ for the satisfied case, $\varepsilon \approx 0.02$
AU for the violated case, 500-period integration); fourth-order
Yoshida at fixed $\mathrm{d}t = 10^{-3}$ (canonical units) for the
Continuity counter-test (equal-mass two-body $a = 1$, $e = 0.5$,
default exact kernel, $R_c = 1$, $\alpha = 0.8$, 60 simulation-unit
integration). Sources at
`crates/apsis-1pn/tests/mercury_precession_gate.rs` and
`crates/apsis-1pn/tests/kernel_continuity_counter_test.rs`; both
reproduce on a clean checkout per the command in §Data and code
availability.

## REBOUND parity portfolio {#sec:rebound-parity}

apsis IAS15 is checked against REBOUND IAS15 [@ReinSpiegel2015] on
four canonical scenarios spanning regular (Kepler, figure-8
choreography), chaotic (Pythagorean three-body, Burrau 1913), and
long-horizon (retrograde Kepler over $10^4$ orbits) regimes. Both
implementations are exercised at matched IAS15 tolerance
($\epsilon_b = 10^{-9}$, the REBOUND default) with identical
initial conditions; trajectories are sampled at one position per
orbital period for the regular scenarios and at fixed $\Delta t =
0.1$ for the Pythagorean run, with total energy sampled at the
same cadence, and compared post-hoc. Per-scenario harnesses live
at `validation/rebound-parity/`.

The configuration-space comparison (Fig. \ref{fig:rebound-traj})
shows apsis and REBOUND trajectories overlapping at the ULP floor
for the two regular scenarios. In the Pythagorean three-body the
two implementations track each other at the ULP floor through the
regular regime; in the chaotic close-encounter cluster, both
diverge from initial conditions at the IAS15 f64 round-off rate,
yielding $|\Delta E|/|E_0| \sim 10^{-10}$ at $T=70$ — the
regime-limited precision rather than an apsis–REBOUND divergence.

![apsis IAS15 trajectories (filled / coloured) overlaid by REBOUND
IAS15 (dotted black) across three configuration-space scenarios.
Kepler $e=0.5$ samples are stroboscopic at periapsis; the
analytical ellipse is shown for reference and the maximum
apsis–REBOUND position residual over the run is annotated as
$|\Delta r|_{\max}$ in the panel callout. Figure-8 choreography
over 10 periods; the three apsis bodies trace the same closed
curve and are visually indistinguishable in the choreography
phase. Pythagorean three-body (Burrau 1913) integrated to $T=70$
through the close-encounter cluster; bodies 0, 1, 2 carry masses
3, 4, 5 in canonical units.](paper/figures/rebound_parity_trajectories.pdf){#fig:rebound-traj width=100%}

The fourth scenario — the retrograde Kepler over $10^4$ orbits —
is the long-horizon parity check and is shown separately in
Fig. \ref{fig:rebound-brouwer} because its log-log energy-error
structure does not share axes with the configuration-space panels
above. Both implementations track Brouwer's $\sqrt{N}$ random-walk
law [@Brouwer1937] and agree to $|\Delta E|/|E_0| = 2.6\times10^{-14}$
at the long-horizon checkpoint — consistent with the per-step
rounding propagation expected at this horizon.

![Relative energy drift $|E(t) - E_0|/|E_0|$ on the retrograde Kepler
scenario over $10^4$ orbits, apsis IAS15 (solid) and REBOUND IAS15
(dotted) on the same axes. The Brouwer $\sqrt{N}$ reference is
shown as a dashed line for visual anchoring. Cross-implementation
agreement at the long-horizon checkpoint is annotated.](paper/figures/rebound_parity_brouwer.pdf){#fig:rebound-brouwer}

## Cross-platform reproducibility {#sec:cross-platform}

The `Cargo.lock`-as-experiment claim is operationalised by an
empirical cross-platform test, distinct from the executable
contract surface of §Methods (which is scoped to single-host
invariants per `crates/apsis/src/contract.rs`). With
the lockfile, the rustc toolchain version, and the source commit
pinned, the full v0.1 federation portfolio — `apsis-1pn` Mercury
long-horizon under IAS15, the MERCURIUS hybrid on a Sun + four
outer planets configuration with a Jupiter-crossing test particle,
Wisdom–Holman on the same outer-planets initial conditions, the
`apsis-central` round-trip gate, and the four REBOUND parity
scenarios reported in §REBOUND parity portfolio — produces
byte-identical trajectory output on heterogeneous x86_64 hosts
(Windows on AMD Zen 4 against Linux on Intel Ice Lake).
Verification covers per-column ULP agreement, file size, and
SHA256 of the captured trajectory CSVs.

The mechanism is the routing of every libc-bound transcendental on
an integration-critical path (`sin`, `cos`, `cbrt`, `pow`, `log`,
...) through the pure-Rust `libm` crate. Hardware-rounded
arithmetic ($+$, $-$, $\times$, $\div$, $\sqrt{}$) is bit-identical
across IEEE 754-conformant x86_64 implementations by mandate and
does not require replacement. The IAS15 step-size controller's
`pow` call was measured against an IEEE-754 correctly-rounded
mpmath oracle on the 42,662 unique inputs generated during the
Mercury 1PN run: Windows UCRT and the `libm` crate match the
oracle on 96.97 % and 95.29 % of inputs respectively, both within
IEEE-754 tolerance for transcendentals (which permits but does not
require correctly-rounded results). The 1-ULP differences in the
remaining cases propagate through the controller's substep-cadence
selection to a 0.002 arcsec/century absolute shift in cumulative
Mercury $\Delta\omega$ over 500 orbits, separating the
$4.4 \times 10^{-6}$ relative agreement on Windows UCRT (scenario-
specific accidental cancellation in UCRT's rounding distribution)
from the $2.8 \times 10^{-5}$ result reproduced bit-identically by
`libm` and glibc. Both values sit several orders of magnitude below
the current observational precision of Mercury's perihelion
advance [@VermaFienga2014] and below the experimental constraint on
the PPN $\gamma$ parameter from Cassini Doppler tracking
[@BertottiIessTortora2003]. The choice prioritises trajectory
determinism across deployment platforms over the scenario-specific
rounding alignment that the Windows-only result reflected.

The methodology and per-implementation analysis are recorded in
`paper/notebooks/2026-05-20-cross-platform-determinism.md`
(federation portfolio bit-equality protocol) and
`paper/notebooks/2026-05-22-controller-pow-implementations.md`
(controller `pow` ULP-distribution analysis).

**Synthesis.** Across the five Results subsections, the library
demonstrates the federated perturbation model along five axes:
three operator crates (`apsis-1pn`, `apsis-radiation`,
`apsis-central`), two operator categories (Hamiltonian and non-
conservative), two constructor patterns (regime-based and
observable-inversion), two kernel invariants (Exactness and
Continuity, each caught independently by the registration check),
and bit-exact reproduction across deployment platforms. This — not
empirical superiority in any numerical regime — is the claim the
architecture supports.

# Discussion {#sec:discussion}

## Scope of the cross-platform claim

The bit-identical cross-platform reproduction reported in
§Cross-platform reproducibility is verified on x86_64 only — specifically
Windows on AMD Zen 4 against Linux on Intel Ice Lake. ARM64 hosts
(including Apple Silicon, increasingly common in the astrophysics
research community) are not in the v0.1 portfolio. The mechanism —
routing libc-bound transcendentals through the pure-Rust `libm`
crate while relying on IEEE 754 hardware arithmetic guarantees —
should extend to ARM64 in principle, but the empirical verification
has not been performed and is therefore not claimed.

The v0.1 hardware-rounded arithmetic guarantee covers FMA-free,
single-precision-free, default-rounding-mode operation. Code paths
using fused-multiply-add explicitly (none in the current
integration-critical surface) or non-default rounding modes would
not inherit the bit-equality.

## Choice of solver family

The library provides Velocity Verlet, Yoshida fourth-order, Wisdom–
Holman in democratic-heliocentric coordinates [@WisdomHolman1991]
(uncorrected leapfrog with Kepler drifts; the order-17 symplectic
corrector of [@Wisdom2006] is tracked as future work), implicit
midpoint, the MERCURIUS hybrid
[@ReinTamayoHernandezPapaloizou2019], and the adaptive Gauss–Radau
IAS15 scheme [@ReinSpiegel2015], alongside stable public traits
for user-registered force models and perturbations. For the
validation portfolio reported here, IAS15 was the appropriate
choice for the Exactness counter-test (adaptive resolution for the
1PN perturbation at periapsis); fourth-order Yoshida exercised the
Continuity counter-test (fixed-step symplectic with clean energy
bookkeeping for the impulsive-error signature). MERCURIUS
[@ReinTamayoHernandezPapaloizou2019] is the appropriate choice in
regimes with frequent close encounters between massive bodies,
where its hybrid switching to symplectic Kepler drifts during
quiescent phases recovers performance that IAS15's adaptive
substep selection cannot match once the controller is forced to
sustain small steps for extended windows. The uncorrected Wisdom–
Holman is retained for compatibility with the existing literature
on leapfrog comparisons; new work should prefer IAS15 or MERCURIUS
unless a specific historical reproduction motivates otherwise.

## Future work

The work below is prioritised: the first tier quantifies a current
floor of the satisfied 1PN result, the second generalises the
cross-platform analysis. (The theory-confirmed counter-test
prerequisite identified in earlier drafts has closed; the §3.2
Plummer and §3.3 continuity predictions are now reported as
closed-form expressions backed by independent scipy verification.)

*Error budget for Mercury 1PN agreement.* The 28 ppm relative
agreement reported in §Results is currently undisaggregated. A
decomposition into IAS15 truncation, `libm` transcendental
tolerance, and the unmodelled $v^4/c^4$ next-order post-Newtonian
correction (estimable from Mercury's orbital parameters and the 1PN
expansion structure) would identify which floor the v0.1
implementation actually sits on and which is the next obstacle to
tightening the claim.

*Cross-platform ULP-distribution analysis.* The summary statistics
in §Cross-platform reproducibility (UCRT 96.97 % oracle-exact, libm
95.29 %) admit an underlying ULP-distribution histogram that
visualises the rounding behaviour of each implementation against
the IEEE-754 correctly-rounded reference. A propagation-sensitivity
derivation linking per-call rounding to the IAS15 controller's
substep-cadence variation closes the chain from per-call ULP
deviation to per-trajectory arcsec/century drift.

# Data and code availability {#sec:availability}

`apsis` is available under the Apache License 2.0 at
<https://github.com/GabrielEstefanski/apsis>. The Mercury validation
described above reproduces on a clean checkout with a single
command,

```bash
cargo test --release -p apsis-1pn --tests -- --ignored
```

after installing a Rust 1.85+ toolchain. The continuous-integration
configuration additionally compiles every example crate, rejects
warnings under `cargo clippy --all-targets`, and verifies that the
library crate resolves no user-interface dependency. A pinned
snapshot of the source archive corresponding to this paper is
deposited at Zenodo (DOI forthcoming).

# Acknowledgements

I thank the authors of REBOUND, REBOUNDx, MERCURIUS, and NBODY6/7
for setting the standards of rigour against which this library's
narrower claim is positioned.

The author used an AI assistant for documentation drafting, code
review, refactoring, and design-discussion support. All algorithmic
decisions, physics interpretations, validation methodology,
numerical implementations, and reported results are the author's
responsibility. The AI assistant was not used for novel physics,
mathematical derivations, or scientific decision-making.

# References

::: {#refs}
:::

# Appendix A: Apsis Record binary format

Each apsis simulation emits an *Apsis Record* — a binary
certificate documenting the run's full provenance and physical
state. The certificate is a single binary file consisting of a
human-readable TOML header followed by a binary frame stream and a
BLAKE3 trailer. The header is emitted at `attach_record` time; the
following block is the verbatim header of a run with `apsis-1pn`
and `apsis-radiation` registered on Sun + Mercury under
`SOLAR_CANONICAL`:

```toml
[apsis]
version = "0.1.0"
git_sha = "9d6e1f50449d72f5499ee520daa049451d4d24cb-dirty"
created_utc = "2026-05-27T20:56:31Z"
rustc_version = "rustc 1.94.1 (e408947bf 2026-03-25)"
generated_by = "apsis 0.1.0"

[reproducibility]
cargo_lock_blake3 = "e3f3742765d9ade1ff9fddfa26bcb050a6f162043c4fc0b37dc560282856e94e"
seed = 42

[unit_system]
g = 1.0000000000000002
length = "AU"
mass = "Msun"
time = "T_G"
density = "Msun/AU3"

[integrator]
kind = "IAS15 (15th, adaptive)"
dt_mode = "Fixed"
initial_dt = 0.001

[kernel]
variant = "Newton"
exactness = "exact"
continuity = "smooth"

[[operators]]
name = "apsis-1pn"
version = "0.1.0"
crate_hash = "workspace:9d6e1f50449d72f5499ee520daa049451d4d24cb-dirty"

[operators.requirements]
kernel_exactness = "exact"
kernel_continuity = "smooth"

[[operators]]
name = "apsis-radiation"
version = "0.1.0"
crate_hash = "workspace:9d6e1f50449d72f5499ee520daa049451d4d24cb-dirty"

[operators.requirements]

[bodies]
count = 2

[[bodies.list]]
name = "Sun"
mass = 1.0
density = 2370030.08
physical_radius = 0.004652851346847559
color = [
    255,
    220,
    100,
]
q_pr = 0.0
albedo = 0.0
class = "Star"

[[bodies.list]]
name = "Mercury"
mass = 0.000000166
density = 6652806.314052459
physical_radius = 0.000018127511930821086
color = [
    139,
    90,
    43,
]
q_pr = 0.0
albedo = 0.3
class = "Planet"
```

Each `[[bodies.list]]` entry records the physical and rendering
metadata for one body; `q_pr` and `albedo` are the radiation-
coupling coefficients consumed by `apsis-radiation`. Numeric fields
are in the canonical units declared above (`length = AU`,
`mass = Msun`, `time = T_G`); the explicit `unit_system` block lets
a replay convert to SI or to a different canonical system without
ambiguity. The `g` field is computed from the canonical SI scales
rather than hardcoded, which records the 1-ULP residual that
`time_s = sqrt(AU³/(G_SI·MSUN))` produces under f64 — the
integrator runs with this `G_code` and the replay must see the same
value. The `dt_mode = "Fixed"` field indicates that `initial_dt`
was supplied explicitly by the caller rather than auto-derived;
IAS15's adaptive substep selection remains active per the
declared integrator kind, so the value records the seed step the
controller takes before the first adaptive adjustment, not a
fixed-step mode override. Material physical events (collisions,
escapes) are recorded inline in the binary frame stream that
follows.

The federation thesis — that a simulation's physical model is
captured in `Cargo.lock` — is here extended to the run itself:
`{record, Cargo.lock}` is the content-addressable closure of the
experiment. The record's frame stream and the trailer's BLAKE3 are
bit-exactly reproducible across replays with the same
configuration; the only per-run difference is the header's wall-
clock `created_utc` field, which is metadata and excluded from the
content hash. A reviewer with both files reproduces the run. The
default policy emits initial + final bookend snapshots and every
material event, which keeps records small and diff-friendly; dense
trajectory capture is an explicit policy opt-in. The bit-equal
trajectory reproduction demonstrated in §Cross-platform
reproducibility applies to the data the Record wraps; the binary
file therefore reproduces across the same heterogeneous hosts
modulo the `created_utc` metadata field.

The Record is the persistent counterpart of the live
`System::cite()` accessor described in §Methods: the same
operator-stack metadata (crate name, version, lockfile hash,
declared `KernelRequirements`) that drives runtime citation
generation is serialised into the Record header. A reviewer
holding a Record can confirm the dependency closure matches a
candidate build by comparing the lockfile hash, then re-run
`System::cite()` against an equivalent `System` to regenerate the
BibTeX block.
