"""Type stubs for the Rust extension module ``apsis_1pn._native``."""

from __future__ import annotations

from apsis import Perturbation

__version__: str

C_SOLAR_UNITS: float
"""Speed of light in the canonical solar-system unit system (AU per
``year/2π``). Compile-time-derived from CODATA SI constants."""

class PostNewtonian1PN:
    """First post-Newtonian gravitational correction (Schwarzschild,
    test-particle form, applied pairwise).

    Constructed via the named factories below. Each returns a fully-formed
    :class:`apsis.Perturbation` ready for ``System.add_perturbation(...)``.
    """

    @staticmethod
    def solar_units() -> Perturbation:
        """1PN calibrated for canonical solar-system units (G = 1, AU, M_sun)."""

    @staticmethod
    def with_c(*, c: float) -> Perturbation:
        """1PN with an explicit speed of light in the caller's unit system."""
