from collections.abc import Sequence

from apsis import Body, Perturbation, UnitSystem

class CentralForce:
    """Central-potential perturbation `∝ r^γ` from a central body
    (Tamayo et al. 2019)."""

    @staticmethod
    def from_raw(
        *,
        source: int,
        a_central: float,
        gamma: float,
        units: UnitSystem,
    ) -> Perturbation:
        """Construct from explicit ``a_central`` and ``γ``."""
    @staticmethod
    def from_apsidal_rate(
        *,
        source: int,
        target: int,
        omega_dot: float,
        gamma: float,
        bodies: Sequence[Body],
        units: UnitSystem,
    ) -> Perturbation:
        """Invert a desired apsidal rate ``ω̇`` for ``target`` orbiting
        ``source`` to produce ``a_central``. ``bodies`` is the list the
        ``System`` will be built from; ``target`` must be on a bound
        orbit around ``source``."""
