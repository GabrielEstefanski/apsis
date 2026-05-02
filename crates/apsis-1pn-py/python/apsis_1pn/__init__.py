"""Python binding for ``apsis-1pn``.

**This package proves that the apsis perturbation extension model is
preserved across both Rust and Python boundaries** — without duplicating
physics, breaking ownership semantics, or requiring kernel modification.
The 1PN formula is implemented exactly once, in the Rust ``apsis-1pn``
crate; this package is plumbing only. New Python perturbation crates
(radiation pressure, J2, drag, …) follow this package's shape.

⚠ Critical precondition
-----------------------

Attaching 1PN to a softened-gravity system **invalidates the physical
model**. For Mercury-like orbits, the numerical apsidal precession
from Plummer softening alone is ~2000× larger than the relativistic
signal *and inverts its sign* — energy and angular momentum stay
conserved at machine precision while the trajectory is physically
wrong.

**This is not a numerical error — it is a model violation.**

Pass ``exact_gravity=True`` to ``apsis.System(...)`` or call
``Body.<material>(...).unsoftened()`` on every body. A violation emits
a structured warning at ``add_perturbation`` time naming the failed
invariant.

Quick start
-----------

.. code-block:: bash

    pip install apsis apsis-1pn

.. code-block:: python

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

``PostNewtonian1PN`` is the only export. Use :meth:`solar_units` for
the canonical solar-system unit system or :meth:`with_c` for any other
unit choice (geometric units, SI, custom).

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
