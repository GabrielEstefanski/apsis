from apsis import Perturbation, UnitSystem

C_SOLAR_UNITS: float
"""Speed of light in the canonical solar unit system (AU per ``year/2π``)."""

class PostNewtonian1PN:
    """First post-Newtonian Schwarzschild correction (test-particle form, applied pairwise).

    Constructed via the named factories below; each returns an
    :class:`apsis.Perturbation` ready for
    ``System.add_hamiltonian_perturbation(...)``. The ``System``
    registration check rejects mismatched ``UnitSystem`` between
    operator and system.
    """

    @staticmethod
    def for_units(*, units: UnitSystem) -> Perturbation:
        """1PN with ``c`` derived from ``units``. Recommended path —
        pass the same ``UnitSystem`` used to build the ``System``."""

    @staticmethod
    def from_raw_c(*, c: float, units: UnitSystem) -> Perturbation:
        """1PN with an explicit ``c`` value. ``units`` is still
        required so the registration unit-system check remains
        load-bearing. Use when ``c`` is computed by neighbouring code
        or intentionally non-physical."""
