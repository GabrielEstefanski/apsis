# apsis-1pn

Python binding for [`apsis-1pn`](../apsis-1pn): the first post-Newtonian gravitational correction as a pluggable [`apsis.Perturbation`](../apsis-py).

Independent package from `apsis` itself. The split mirrors the Rust workspace and reflects the project's architectural thesis: every additional force is a separately citable component, not a hardcoded option of the simulator core.

## Install

```bash
pip install apsis apsis-1pn
```

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
    units=apsis.units.SOLAR,
    integrator="ias15",
    dt=1e-3,
    exact_gravity=True,
)
sys.add_perturbation(apsis_1pn.PostNewtonian1PN.solar_units())
sys.integrate_for(100.0)
```

## Kernel preconditions

1PN is derived around the bit-exact Newtonian potential. Attaching it to a softened-gravity system (Plummer kernel with ε > 0) substitutes a different unperturbed Hamiltonian whose apsidal precession alone is ~2 × 10³ larger than the 1PN signal for a Mercury-like orbit, silently inverting the sign of the measured precession. Either pass `exact_gravity=True` to `apsis.System(...)` or call `Body.<material>(...).unsoftened()` on every body.

The wrapper does not check this — the precondition is surfaced as a structured warning at `add_perturbation` time so the kernel-vs-perturbation contract is enforced once, in the core.

## References

- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics*. Cambridge.
- Damour, T., & Schäfer, G. (1985). General relativistic equations of motion. *General Relativity and Gravitation*, 17, 879.
