# apsis-1pn-py

Python bindings for [`apsis-1pn`](../apsis-1pn).

The 1PN operator constructed in `apsis-1pn` crosses the Rust → Python
boundary via a typed [`PyCapsule`](https://docs.python.org/3/c-api/capsule.html)
carrying a `Box<dyn HamiltonianOperator>`. The capsule transport lives
in [`apsis-py-core`](../apsis-py-core) and is consumed here.

## Extension contract

Perturbations registered through `apsis.System.add_hamiltonian_perturbation` must:

- **operate on the exact Newtonian kernel** when their derivation requires it (declared via `kernel_requirements()`);
- **be additive** — accumulate into the scratch buffer, never modify the base Hamiltonian;
- **be attachable at runtime** via the typed PyCapsule transport, with single-consume ownership semantics enforced at the FFI boundary.

Preconditions are expressed at the type level and surfaced at runtime, ensuring that invalid physical configurations are detectable without coupling perturbations to the kernel. This crate is the reference implementation of that contract at the Python boundary.

## ⚠️ Critical precondition

Attaching 1PN to a softened-gravity system **invalidates the physical model**.

For Mercury-like orbits, the numerical apsidal precession induced by Plummer softening alone is **~2000× larger than the relativistic signal, and inverts its sign**. Energy and angular momentum stay conserved at machine precision while the trajectory is physically wrong.

**This is not a numerical error — it is a model violation.**

Pass `exact_gravity=True` to `apsis.System(...)` or call `Body.<material>(...).unsoftened()` on every body. The kernel-vs-perturbation contract is enforced once, in the core: a violation emits a structured warning at `add_hamiltonian_perturbation` time naming the failed invariant. The warning is the deliberate behaviour — apsis does not silently correct invalid physical configurations. Surfacing the violation is the contract; auto-fixing would erase it.

## Why this matters

Splitting each perturbation into an independently citable crate enables:

- **reproducible scientific workflows** — a paper's perturbation set is its dependency set, versioned and citable;
- **independent validation** — each force is tested in isolation against an analytic reference (1PN: 4.4 ppm of GR for Mercury's perihelion);
- **composition without kernel modification** — additive perturbations stack via the trait; no privileged extension is hardcoded into the simulator core.

## Use

```python
import apsis
import apsis_1pn

sun = apsis.Body.star(mass=1.0).unsoftened()
mercury = (apsis.Body.rocky(mass=1.66e-7)
           .at((0.387, 0.0))
           .with_velocity((0.0, 1.61))
           .unsoftened())

sys = apsis.System(
    bodies=[sun, mercury],
    units=apsis.units.SOLAR_CANONICAL,
    integrator="ias15",
    dt=1e-3,
    exact_gravity=True,
)
# Same UnitSystem on both sides — registration check passes.
sys.add_hamiltonian_perturbation(
    apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
)
sys.integrate_for(100.0)
```

## Scope

- Plumbing only. The 1PN formula is implemented once in
  [`apsis-1pn`](../apsis-1pn); this crate carries it across the
  Rust → Python boundary.

## References

- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics*. Cambridge.
- Damour, T., & Schäfer, G. (1985). General relativistic equations of motion. *General Relativity and Gravitation*, 17, 879.
