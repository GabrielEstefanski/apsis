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


class Perturbation:
    """User-facing wrapper for a Hamiltonian-class perturbation plugin.

    Researchers never construct ``Perturbation`` directly. Each
    perturbation crate (``apsis_1pn``, future J2 / tidal packages)
    provides factory methods that return a fully formed instance:

    .. code-block:: python

        import apsis
        import apsis_1pn

        sys = apsis.System(bodies=[...], units=apsis.units.SOLAR, ...)
        sys.add_hamiltonian_perturbation(
            apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
        )

    Non-conservative operators (drag, radiation reaction) travel in a
    separate capsule type with its own registration entry point; that
    surface is not yet exposed across the FFI.

    Why a pure-Python class rather than a PyO3 ``#[pyclass]``:
    cross-extension type identity in PyO3 is unreliable (each
    ``cdylib`` registers its own Python class object even when the
    underlying Rust type is shared via an ``rlib``). Defining the
    user-facing class once, here, gives every perturbation crate the
    same ``apsis.Perturbation`` to pass into
    ``System.add_hamiltonian_perturbation``. The boxed Rust trait
    object travels in the ``_capsule`` attribute (an opaque
    ``PyCapsule``) which the
    ``System.add_hamiltonian_perturbation`` boundary unwraps via the
    shared helpers in the ``apsis-py-core`` Rust crate.
    """

    __slots__ = ("_capsule", "_label")

    def __init__(self, _capsule: object, _label: str) -> None:
        self._capsule = _capsule
        self._label = _label

    @property
    def label(self) -> str:
        """Human-readable label set by the constructing crate."""
        return self._label

    def __repr__(self) -> str:
        return f"Perturbation(label={self._label!r})"


__all__ = [
    "AdaptiveStats",
    "Body",
    "IntegratorKind",
    "Perturbation",
    "Stats",
    "System",
    "Trajectory",
    "UnitSystem",
    "units",
    "__version__",
]
