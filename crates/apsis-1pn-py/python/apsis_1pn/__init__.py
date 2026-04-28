"""First post-Newtonian gravitational correction as an apsis Perturbation plugin.

Independent package from :mod:`apsis` itself, mirroring the Rust workspace
split between the simulator core and its physics extensions. A researcher
who needs the relativistic correction installs both:

.. code-block:: bash

    pip install apsis apsis-1pn

and writes:

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

The ``PostNewtonian1PN`` class is the headline export. Use
:meth:`solar_units` for the canonical solar-system unit system or
:meth:`with_c` for any other unit choice (geometric units, SI, custom).

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
