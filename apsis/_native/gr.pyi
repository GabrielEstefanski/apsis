"""Type stubs for ``apsis._native.gr``."""

from __future__ import annotations

from apsis import Perturbation, UnitSystem

C_SOLAR_UNITS: float

class PostNewtonian1PN:
    @staticmethod
    def for_units(*, units: UnitSystem) -> Perturbation: ...
    @staticmethod
    def from_raw_c(*, c: float, units: UnitSystem) -> Perturbation: ...
