"""Type stubs for the Rust extension module ``apsis_1pn._native``."""

from __future__ import annotations

from apsis import Perturbation, UnitSystem

__version__: str

C_SOLAR_UNITS: float
"""Speed of light in the canonical solar-system unit system (AU per
``year/2π``). Compile-time-derived from CODATA SI constants."""

class PostNewtonian1PN:
    """First post-Newtonian gravitational correction (Schwarzschild,
    test-particle form, applied pairwise).

    Constructed via the named factories below. Each returns a fully-formed
    :class:`apsis.Perturbation` ready for
    ``System.add_hamiltonian_perturbation(...)``.

    The factories follow the observable constructor convention documented
    in ``apsis::contract``. Every constructor binds the operator to a
    :class:`apsis.UnitSystem`; the ``System`` registration check panics
    on unit-system mismatch, so you cannot silently mix an operator
    built for one unit system into a ``System`` integrating in another.

    - :meth:`for_units` is the recommended path: derive ``c`` from the
      same ``UnitSystem`` you passed to ``apsis.System(...)``.
    - :meth:`from_raw_c` is the explicit-``c`` escape: useful when ``c``
      is computed by neighbouring code or intentionally non-physical.
      ``units`` is still required so the registration check remains
      load-bearing.
    """

    @staticmethod
    def for_units(*, units: UnitSystem) -> Perturbation:
        """1PN with ``c`` derived from the supplied :class:`apsis.UnitSystem`."""

    @staticmethod
    def from_raw_c(*, c: float, units: UnitSystem) -> Perturbation:
        """1PN with an explicit ``c`` value, pinned to ``units`` for the
        registration-time unit-system check.
        """
