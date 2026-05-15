"""Radiation pressure on dust grains and small bodies (Burns et al. 1979).

```python
from apsis.radiation import RadiationPressure
sys.add_hamiltonian_perturbation(
    RadiationPressure.from_raw_betas(
        source=0,
        betas=[0.0, 0.1],   # body 0 = sun (no β), body 1 = dust grain (β = 0.1)
        units=apsis.units.SOLAR_CANONICAL,
    ),
)
```

Each body is parameterised by its dimensionless `β = F_rad / F_grav`
ratio. The operator reduces effective central gravity on each receiver
by `(1 - β)` per Burns 1979.

Poynting–Robertson drag (the velocity-dependent companion) is
implemented in the `apsis-radiation` Rust crate but not yet exposed
to Python — the non-conservative-operator capsule transport is
pending.
"""

from apsis._native.radiation import RadiationPressure

__all__ = ["RadiationPressure"]
