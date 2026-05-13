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
    in ``apsis::contract``:

    - Named-regime factories (:meth:`solar_units`, :meth:`for_units`)
      derive ``c`` from a known physical setup or a supplied
      :class:`apsis.UnitSystem`. No raw numeric input.
    - Raw-escape factories (:meth:`from_raw_c`,
      :meth:`from_raw_c_validated`) accept ``c`` directly. The
      ``_validated`` form cross-checks against a unit system and raises
      ``ValueError`` on mismatch.
    """

    @staticmethod
    def solar_units() -> Perturbation:
        """1PN calibrated for canonical solar-system units (G = 1, AU, M_sun)."""

    @staticmethod
    def for_units(*, units: UnitSystem) -> Perturbation:
        """1PN with ``c`` derived from the supplied :class:`apsis.UnitSystem`."""

    @staticmethod
    def from_raw_c(*, c: float) -> Perturbation:
        """1PN with an explicit speed of light, no validation."""

    @staticmethod
    def from_raw_c_validated(*, c: float, units: UnitSystem) -> Perturbation:
        """1PN with an explicit ``c``, cross-checked against ``units``.

        Raises ``ValueError`` when the relative error between the
        supplied ``c`` and the ``c`` derived from ``units`` exceeds
        ``1e-9``. Use when ``c`` originates from an external source
        and you want the simulator to confirm it before attaching.
        """
