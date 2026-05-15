"""Type stubs for ``apsis._native.radiation``."""

from __future__ import annotations

from collections.abc import Sequence

from apsis import Perturbation, UnitSystem

class RadiationPressure:
    @staticmethod
    def from_raw_betas(
        *,
        source: int,
        betas: Sequence[float],
        units: UnitSystem,
    ) -> Perturbation:
        """Construct from per-body β. ``betas[i]`` is applied to body ``i``;
        ``betas[source]`` should be ``0.0``. Non-zero source β is not a
        construction error but is flagged as a ``self_radiation`` regime
        violation by the System's ``check_regime`` at registration time."""
