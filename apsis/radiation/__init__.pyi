from collections.abc import Sequence

from apsis import Perturbation, UnitSystem

class RadiationPressure:
    """Radiation pressure indexed by per-body `β = F_rad / F_grav` (Burns 1979)."""

    @staticmethod
    def from_raw_betas(
        *,
        source: int,
        betas: Sequence[float],
        units: UnitSystem,
    ) -> Perturbation:
        """Construct from an explicit β per body. ``source`` is the
        radiating body's index; ``betas[i]`` is the β applied to body
        ``i``. Set ``betas[source] = 0.0`` — non-zero values are flagged
        as a ``self_radiation`` regime violation by the System's
        ``check_regime`` at registration."""
