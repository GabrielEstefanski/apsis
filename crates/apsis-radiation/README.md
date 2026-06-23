# apsis-radiation

Out-of-tree perturbation crate for [`apsis`](../apsis). Radiation
pressure (Hamiltonian) and Poynting–Robertson drag (non-conservative)
per Burns, Lamy & Soter (1979).

Ships a
[`NonConservativeOperator`](../apsis/src/physics/integrator/operator.rs)
(Poynting–Robertson drag) alongside a conservative
[`HamiltonianOperator`](../apsis/src/physics/integrator/operator.rs)
(radiation pressure), both derived from Burns et al. (1979).

Compiles against the public `apsis` API alone — no `pub(crate)` access,
no patches to core sources.

## Operators

### `RadiationPressure` — Hamiltonian

Radial 1/r² force opposing gravity from the source body, scaled per
receiver by `β` (Burns et al. eq. 19, the dimensionless
radiation-to-gravity ratio):

```text
  F_rad(i ← source) = + β_i · G · M_source · m_i / r² · r̂
```

Conservative — derives from the closed-form potential `V_rad = − β·M·m/r`,
which the operator publishes through
[`HamiltonianOperator::potential`](../apsis/src/physics/integrator/operator.rs).
`System::energy()` therefore accounts for radiation-pressure energy on
the same footing as gravity.

### `PoyntingRobertsonDrag` — non-conservative

Relativistic angular-momentum loss from re-emitted radiation (Burns
et al. eq. 7):

```text
  F_PR(i ← source) = − (β_i · G · M_source / r²) · [ (2·v_r/c) · r̂ + v / c ]
```

where `v_r = v · r̂`. Symplectic integrators registered alongside this
operator emit a structured warning at registration time — the
energy-drift signal is the physical effect, not a numerical artifact.

For a circular orbit at semi-major axis `a` the analytic decay
timescale is

```text
  τ_PR = a² · c² / (4 · β · G · M_source)
```

## Validation signal

The CI gate in [`tests/dust_decay_gate.rs`](tests/dust_decay_gate.rs)
integrates a β = 0.5 dust grain on a circular orbit for ten orbits with
PR drag enabled. The measured energy drift agrees with the Burns 1979
prediction to within the gate's 5 % tolerance (typical agreement is
< 1 %; the bound covers ULP variance and the constant-r approximation
in the analytic formula). A counter-test confirms IAS15 alone preserves
energy to 10⁻¹² on the same baseline, so the measured drift is
attributable to the operator and not to integrator noise.

## Convention

- `bodies[source]` is the radiating body (typically `bodies[0]`); it
  must itself feel no radiation force, so `betas[source]` must be 0.
  The invariant is checked at registration via `Operator::check_regime`
  and surfaced as a `Hard` violation on breach.
- `betas[i] ≥ 0`. `β > 1` is the "blowout" regime where radiation
  overpowers gravity; the operator integrates the unbound trajectory
  without warning.
- `r̂` points from receiver to source — matches the apsis sign
  convention shared with `apsis-1pn`.

## Use

```rust
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_radiation::{PoyntingRobertsonDrag, RadiationPressure};

let units = UnitSystem::solar_canonical();
let sun = Body::star(1.0);
// β = 0.1 dust grain on a circular orbit at 1 AU.
let dust = Body::rocky(1e-15).at(1.0, 0.0).with_velocity(0.0, 1.0);
let mut sys = System::new(vec![sun, dust], units)
    .with_integrator(IntegratorKind::Ias15)
    .with_dt(1e-3);

sys.add_hamiltonian_perturbation(Box::new(
    RadiationPressure::from_raw_betas(0, vec![0.0, 0.1], units),
))?;
sys.add_non_conservative_perturbation(Box::new(
    PoyntingRobertsonDrag::from_raw_betas(0, vec![0.0, 0.1], units),
))?;
```

## Reference

Burns, J. A., Lamy, P. L., & Soter, S. (1979). Radiation forces on
small particles in the solar system. *Icarus* 40, 1–48.
DOI: [10.1016/0019-1035(79)90050-2](https://doi.org/10.1016/0019-1035(79)90050-2).
