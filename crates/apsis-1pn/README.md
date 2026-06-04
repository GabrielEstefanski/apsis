# apsis-1pn

Out-of-tree perturbation crate for [`apsis`](../apsis). Adds the 1PN
(test-particle Schwarzschild) correction as a `HamiltonianOperator`.

Compiles against the public `apsis` API alone — no `pub(crate)` access,
no patches to core sources.

## Extension contract

Perturbations registered through [`System::add_hamiltonian_perturbation`](../apsis/src/core/system/perturbations.rs) (Hamiltonian-class operators like 1PN) or `System::add_non_conservative_perturbation` (drag, radiation reaction) must:

- **operate on the exact Newtonian kernel** when their derivation requires it (declared via [`kernel_requirements()`](../apsis/src/physics/integrator/operator.rs));
- **be additive** — accumulate into the supplied buffer, never overwrite or modify the base Hamiltonian;
- **declare their physical preconditions at the type level** so the kernel-vs-operator contract can be checked at registration, not at publication time.

Preconditions are expressed at the type level and surfaced at runtime, ensuring that invalid physical configurations are detectable without coupling perturbations to the kernel. This crate is the reference implementation of that contract.

## ⚠️ Critical precondition

1PN is derived around the bit-exact Newtonian potential. Default `System::new(...)` uses `NewtonKernel::exact()` (ε = 0) and the registration is silent. Attaching 1PN on top of a softened kernel **invalidates the physical model**.

For Mercury-like orbits, the numerical apsidal precession induced by Plummer softening alone is **~2000× larger than the relativistic signal, and inverts its sign**. Energy and angular momentum stay conserved at machine precision while the trajectory is physically wrong.

**This is not a numerical error — it is a model violation.**

The contract is enforced once, in the core: opting into a softened kernel via `System::with_kernel(Arc::new(NewtonKernel::new(ε > 0)))` emits a structured warning at `add_hamiltonian_perturbation` time naming the failed invariant. The warning is the deliberate behaviour — apsis does not silently correct invalid physical configurations. Surfacing the violation is the contract; auto-fixing would erase it.

## Validation signal

With the contract enforced, this crate reproduces Mercury's textbook 43 arcsec/century rate to within **100 ppm** of the GR prediction, gated in CI and bit-identical across Windows and Linux on x86_64. Example output (the recommended `for_units` API path):

```text
$ cargo run --release -p apsis-1pn --example mercury_perihelion
...
── GR comparison over 500 orbits ──
  predicted Δω      = +2.509332e-04 rad  (+51.7587 arcsec)
  measured  Δω      = +2.509130e-04 rad  (+51.7545 arcsec)
  relative error    = -8.053e-05
  rate              = 42.978 arcsec/century  (GR expects 43)
```

The Mercury agreement is gated in CI at 100 ppm (`mercury-gate` job in [`.github/workflows/rust.yml`](../../.github/workflows/rust.yml)). The cross-platform deterministic floor is documented in [`paper/notebooks/2026-05-22-controller-pow-implementations.md`](../../paper/notebooks/2026-05-22-controller-pow-implementations.md).

## Why this matters

Splitting each perturbation into an independently citable crate enables:

- **reproducible scientific workflows** — a paper's perturbation set is its dependency set, versioned and citable;
- **independent validation** — each force is tested in isolation against an analytic reference;
- **composition without kernel modification** — additive perturbations stack via the trait; no privileged extension is hardcoded into the simulator core.

## Use

```rust
use apsis::core::system::System;
use apsis::domain::body::Body;
use apsis::physics::integrator::IntegratorKind;
use apsis::units::UnitSystem;
use apsis_1pn::PostNewtonian1PN;

let units = UnitSystem::solar_canonical();
let sun = Body::star(1.0);
let mercury = Body::rocky(1.66e-7).at(0.387, 0.0).with_velocity(0.0, 1.61);

let mut sys = System::new(vec![sun, mercury], units)
    .with_integrator(IntegratorKind::Ias15)
    .with_dt(1e-4);
sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::for_units(units)));
sys.integrate_for(100.0);
```

## Scope

- Test-particle Schwarzschild limit only. Full Einstein–Infeld–Hoffmann
  cross-terms (multi-body PN) are out of scope.
- Dependency set restricted to `apsis`.

## References

- Rein, H., & Spiegel, D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics. *MNRAS*, 446, 1424–1437.
- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics*. Cambridge.
- Damour, T., & Schäfer, G. (1985). General relativistic equations of motion. *General Relativity and Gravitation*, 17, 879.
