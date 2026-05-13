"""Type stubs for the public ``apsis`` package surface.

The class-level signatures are owned by ``apsis/_native.pyi``, which
mirrors the runtime API exposed by the Rust extension module
``apsis._native``. This file simply re-exports the public symbols to
match what ``apsis/__init__.py`` does at runtime; it is the single
type-checker entry point a researcher's IDE consults when they write
``import apsis``.

Each subsequent PR that adds a class or free function to the Rust
side adds a matching declaration in ``_native.pyi`` and a re-export
line here in the same commit; type checking is not a follow-up task.
"""

from apsis._native import (
    AdaptiveStats as AdaptiveStats,
    Body as Body,
    IntegratorKind as IntegratorKind,
    Stats as Stats,
    System as System,
    Trajectory as Trajectory,
    UnitSystem as UnitSystem,
    __version__ as __version__,
)
from apsis._native import units as units


class Perturbation:
    """Pure-Python wrapper for a Hamiltonian-class perturbation plugin.

    Constructed only by perturbation crates (``apsis_1pn``, ...). Pass to
    :meth:`System.add_hamiltonian_perturbation` exactly once.
    """

    def __init__(self, _capsule: object, _label: str) -> None: ...
    @property
    def label(self) -> str: ...
    def __repr__(self) -> str: ...


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
