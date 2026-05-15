"""Python binding for ``apsis-1pn``.

**This package proves that the apsis perturbation extension model is
preserved across both Rust and Python boundaries** — without duplicating
physics, breaking ownership semantics, or requiring kernel modification.
The 1PN formula is implemented exactly once, in the Rust ``apsis-1pn``
crate; this package is plumbing only. New Python perturbation crates
(radiation pressure, J2, drag, …) follow this package's shape.

⚠ Critical precondition
-----------------------

1PN is derived around the bit-exact Newtonian potential. Default
``apsis.System(...)`` uses an exact ``NewtonKernel`` (ε = 0) and the
registration is silent. Attaching 1PN on top of a softened kernel
**invalidates the physical model**: the numerical apsidal precession
from Plummer softening alone is ~2000× larger than the relativistic
signal at Mercury's orbit *and inverts its sign* — energy and angular
momentum stay conserved at machine precision while the trajectory is
physically wrong.

**This is not a numerical error — it is a model violation.**

The kernel-requirement check emits a structured warning at
``add_hamiltonian_perturbation`` time if a softened kernel is in
place (currently reachable only from the Rust side via
``System::with_kernel(NewtonKernel::new(ε > 0))``).

Quick start
-----------

.. code-block:: bash

    pip install apsis apsis-1pn

.. code-block:: python

    import apsis
    import apsis_1pn

    sun = apsis.Body.star(mass=1.0)
    mercury = (apsis.Body.rocky(mass=1.66e-7)
               .at((0.387, 0.0))
               .with_velocity((0.0, 1.61)))

    sys = apsis.System(
        bodies=[sun, mercury],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    sys.add_hamiltonian_perturbation(
        apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
    )
    sys.integrate_for(100.0)

``PostNewtonian1PN`` is the only export. Use :meth:`for_units` to derive
``c`` from any :class:`apsis.UnitSystem` (the recommended path), or
:meth:`from_raw_c` to pass an explicit ``c`` value pinned to a unit
system. The ``System`` registration check panics if the perturbation's
unit system disagrees with the ``System``'s own.

References
----------

- Will, C. M. (1993). *Theory and Experiment in Gravitational Physics*. Cambridge.
- Damour, T., & Schäfer, G. (1985). General relativistic equations of motion.
  *General Relativity and Gravitation*, 17, 879.
"""

from apsis_1pn._native import (
    C_SOLAR_UNITS,
    PostNewtonian1PN,
    __version__,
)

__all__ = [
    "C_SOLAR_UNITS",
    "PostNewtonian1PN",
    "__version__",
]
