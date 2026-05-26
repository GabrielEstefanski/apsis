"""Python bindings for the apsis N-body simulation library.

Operator submodules live under ``apsis.<name>`` (e.g. ``apsis.gr``).

```python
import apsis
from apsis.gr import PostNewtonian1PN

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
    PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
)
sys.integrate_for(100.0)
```

The Rust crate ``apsis`` is the source of truth for behaviour.
"""

from apsis._native import (
    AdaptiveStats,
    Body,
    IntegratorKind,
    Record,
    Stats,
    System,
    Trajectory,
    UnitSystem,
    UnitSystemMismatchError,
    __version__,
    units,
)


class Perturbation:
    """User-facing wrapper for a Hamiltonian-class perturbation.

    Constructed by operator-submodule factories (``apsis.gr.PostNewtonian1PN.for_units``,
    external ``apsis-plugin-X`` packages); never instantiated directly
    by researchers. Carries the boxed Rust operator in the ``_capsule``
    attribute (opaque ``PyCapsule``) until
    ``System.add_hamiltonian_perturbation`` consumes it.
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
    "Record",
    "Stats",
    "System",
    "Trajectory",
    "UnitSystem",
    "UnitSystemMismatchError",
    "__version__",
    "units",
]
