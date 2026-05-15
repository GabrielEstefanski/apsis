"""Type stubs for ``apsis._native.central``."""

from __future__ import annotations

from collections.abc import Sequence

from apsis import Body, Perturbation, UnitSystem

class CentralForce:
    @staticmethod
    def from_raw(
        *,
        source: int,
        a_central: float,
        gamma: float,
        units: UnitSystem,
    ) -> Perturbation: ...
    @staticmethod
    def from_apsidal_rate(
        *,
        source: int,
        target: int,
        omega_dot: float,
        gamma: float,
        bodies: Sequence[Body],
        units: UnitSystem,
    ) -> Perturbation: ...
