"""Python bindings for the apsis N-body simulation library.

The :mod:`apsis` Python package is a thin façade over the Rust crate of
the same name. The bindings expose a researcher-first API for setting
up and integrating gravitational systems; the underlying numerical
work (integrators, force models, conservation diagnostics) is performed
in Rust and is never duplicated on the Python side.

Quickstart
----------

Build a custom system explicitly::

    >>> import apsis
    >>> sun = apsis.Body.star(mass=1.0).unsoftened()
    >>> mercury = (apsis.Body.rocky(mass=3e-6)
    ...            .at((0.307, 0.0))
    ...            .with_velocity((0.0, 1.98))
    ...            .unsoftened())
    >>> sys = apsis.System(
    ...     bodies=[sun, mercury],
    ...     integrator="ias15",
    ...     dt=1e-3,
    ... )
    >>> sys.integrate_for(100.0)
    >>> abs(sys.energy_delta) < 1e-12
    True

References
----------

For the algorithmic specification of the IAS15 integrator, see Rein &
Spiegel (2015). For the cross-implementation parity portfolio that
validates this binding's underlying Rust core against REBOUND on
canonical scenarios (Kepler, figure-8 choreography), see the
``docs/experiments/`` directory in the project repository.

The Rust crate ``apsis`` is the source of truth for behaviour; this
package follows it bit-for-bit.
"""

from apsis._native import (
    AdaptiveStats,
    Body,
    IntegratorKind,
    Stats,
    System,
    Trajectory,
    UnitSystem,
    __version__,
)
from apsis._native import units

__all__ = [
    "AdaptiveStats",
    "Body",
    "IntegratorKind",
    "Stats",
    "System",
    "Trajectory",
    "UnitSystem",
    "units",
    "__version__",
]
