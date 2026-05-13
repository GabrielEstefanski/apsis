# apsis-1pn

Out-of-tree perturbation crate for [`apsis`](../apsis).

**This crate proves that the perturbation extension contract is buildable, not just documented.** It compiles against the public API alone — no `pub(crate)` access, no patches to core sources, no dependency other than `apsis` itself. A future change to that API that breaks this crate fails CI loudly rather than quietly.

Use it as the **template** when writing new perturbation crates (radiation pressure, J2, drag, …).

## Extension contract

Perturbations registered through [`System::add_hamiltonian_perturbation`](../apsis/src/core/system/perturbations.rs) (Hamiltonian-class operators like 1PN) or `System::add_non_conservative_perturbation` (drag, radiation reaction) must:

- **operate on the exact Newtonian kernel** when their derivation requires it (declared via [`kernel_requirements()`](../apsis/src/physics/integrator/operator.rs));
- **be additive** — accumulate into the supplied buffer, never overwrite or modify the base Hamiltonian;
- **declare their physical preconditions at the type level** so the kernel-vs-operator contract can be checked at registration, not at publication time.

Preconditions are expressed at the type level and surfaced at runtime, ensuring that invalid physical configurations are detectable without coupling perturbations to the kernel. This crate is the reference implementation of that contract.

## ⚠️ Critical precondition

Attaching 1PN to a softened-gravity system **invalidates the physical model**.

For Mercury-like orbits, the numerical apsidal precession induced by Plummer softening alone is **~2000× larger than the relativistic signal, and inverts its sign**. Energy and angular momentum stay conserved at machine precision while the trajectory is physically wrong.

**This is not a numerical error — it is a model violation.**

Call `Body::unsoftened()` on every body or `System::with_exact_gravity()` system-wide. The contract is enforced once, in the core: a violation emits a structured warning at `add_hamiltonian_perturbation` time naming the failed invariant. The warning is the deliberate behaviour — apsis does not silently correct invalid physical configurations. Surfacing the violation is the contract; auto-fixing would erase it.

## Validation signal

With the contract enforced, this crate reproduces Mercury's textbook 43 arcsec/century rate to **4.4 ppm**:

```text
$ cargo run --release -p apsis-1pn --example mercury_perihelion
...
── GR comparison over 500 orbits ──
  predicted Δω      = +2.509427e-04 rad  (+51.7606 arcsec)
  measured  Δω      = +2.509438e-04 rad  (+51.7609 arcsec)
  relative error    = +4.449e-06
  rate              = 42.983 arcsec/century  (GR expects 43)
```

The number is gated in CI — see the `mercury-gate` job in [`.github/workflows/rust.yml`](../../.github/workflows/rust.yml).

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

let sun = Body::star(1.0).unsoftened();
let mercury = Body::rocky(1.66e-7).at(0.387, 0.0).with_velocity(0.0, 1.61).unsoftened();

let mut sys = System::new(vec![sun, mercury], UnitSystem::canonical())
    .with_integrator(IntegratorKind::Ias15)
    .with_dt(1e-4);
sys.add_hamiltonian_perturbation(Box::new(PostNewtonian1PN::solar_units()));
sys.integrate_for(100.0);
```

## What this is NOT

- **Distribution is intentionally minimal.** The dependency set never grows beyond `apsis` — that constraint *is* the proof.
- **Not a production GR engine.** Full Einstein–Infeld–Hoffmann cross-terms are out of scope; only the test-particle Schwarzschild limit is implemented.
- **Different architectural axis from REBOUND/REBOUNDx.** Mature codes solve "integrate this Solar System with extra forces"; this project explores strict kernel/perturbation separation and citable force composition. Research beyond the test-particle regime, and the full breadth of REBOUNDx (gr_full, spin-orbit, GW emission), should use REBOUNDx.

## References

- Rein, H., & Spiegel, D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics. *MNRAS*, 446, 1424–1437.
- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics*. Cambridge.
- Damour, T., & Schäfer, G. (1985). General relativistic equations of motion. *General Relativity and Gravitation*, 17, 879.
