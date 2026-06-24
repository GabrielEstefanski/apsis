# apsis-central

Out-of-tree perturbation crate for [`apsis`](../apsis). Generalized
central force `a = A · r^γ` per Tamayo, Rein, Shi & Hernandez
(2020, *MNRAS* 491, 2885).

**Observable-inversion exemplar.** This crate is the federation's
first implementation of the **observable-inversion constructor**
locked in ADR-005 / ADR-006: [`CentralForce::from_apsidal_rate`] takes a
measured (or desired) apsidal precession rate `ω̇` and inverts it
into the coupling `A` that produces it. The library is named after
the apsidal axis; the highest-leverage feature reproduces an apsidal
observable.

## Operator

| Field | Meaning | Convention |
|---|---|---|
| `source` | Body that sources the central force | Index into `&[Body]` |
| `A` | Coupling coefficient | `a = A · r^γ · r̂` (force per unit receiver mass) |
| `γ` | Radial power | `−2` is the gravity-like degenerate case (no precession) |

| `γ` | Force law | Notable use |
|---|---|---|
| `−3` | `A / r³` | Effective Schwarzschild precession (Nobili & Roxburgh 1986) |
| `−2` | `A / r²` | Degenerate — looks like gravity, no apsidal precession |
| `−1` | `A / r` | Logarithmic potential (galactic halo flat rotation) |
| `+1` | `A · r` | Hooke / harmonic oscillator |

Conservative — closed-form `V_central` published through
`HamiltonianOperator::potential`; `System::energy` accounts for the
radial contribution. Newton's third law applied to the source body
(recoil scaled by `−m_recv / m_src`); momentum conservation is the
test gate.

## Observable inversion

`CentralForce::from_apsidal_rate(source, target, ω̇, γ, &bodies, units)`
applies the Tamayo 2020 inversion:

```text
  A = G · M_source · ω̇ / [(1 + γ/2) · d^(γ+2) · n]
```

where `d` is the instantaneous separation and `n` the mean motion of
the target's current orbit. Errors:

| `Err` variant | Trigger |
|---|---|
| `DegenerateGamma` | `γ ≈ −2` — precession identically vanishes for `1/r²`, `A` diverges |
| `IndexOutOfRange` | source / target past the body vector end |
| `SourceEqualsTarget` | source == target |
| `UnboundOrbit` | target on a hyperbolic / parabolic flyby; mean motion undefined |

## Validation gate

`tests/round_trip_gate.rs` closes the loop end-to-end: register at
`ω̇ = 1.5 × 10⁻³ rad / Gaussian time` on an `e = 0.1` orbit, integrate
50 orbits with IAS15, fit `ω̇` from period-locked samples, assert
agreement within 5%. Empirical agreement when written: 2.38%. The
counter-test confirms IAS15 alone produces no apsidal drift on the
same baseline (`< 10⁻⁹`), attributing the measured drift to the
operator.

The 5% bound is set by the Tamayo inversion's near-circular
approximation (uses instantaneous `d` rather than `a`), not by IAS15
truncation or operator-implementation noise. At `e → 0` agreement
tightens to the IAS15 floor.

## Use

```rust
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_central::CentralForce;

let units = UnitSystem::solar_canonical();
let sun = Body::star(1.0);
let target = Body::rocky(1e-7).at(0.387, 0.0).with_velocity(0.0, 1.61);
let bodies = vec![sun, target];

// Observable inversion: pick a desired ω̇ and let the operator compute the coupling.
let force = CentralForce::from_apsidal_rate(
    0,                              // source
    1,                              // target
    5e-9,                            // rad / Gaussian time
    -3.0,                            // γ = −3 (Schwarzschild-effective)
    &bodies,
    units,
)?;

let mut sys = System::new(bodies, units)
    .with_integrator(IntegratorKind::Ias15)
    .with_dt(1e-3);
sys.add_hamiltonian_perturbation(Box::new(force))?;
```

## Reference

Tamayo, D., Rein, H., Shi, P., & Hernandez, D. M. (2020). REBOUNDx: a
library for adding conservative and dissipative forces to otherwise
symplectic N-body integrations. *MNRAS* 491, 2885–2901.
DOI: [10.1093/mnras/stz2870](https://doi.org/10.1093/mnras/stz2870).

[`CentralForce::from_apsidal_rate`]: src/lib.rs
