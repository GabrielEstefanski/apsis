"""First post-Newtonian Schwarzschild correction (Anderson et al. 1975).

```python
from apsis.gr import PostNewtonian1PN
sys.add_hamiltonian_perturbation(
    PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
)
```

`C_SOLAR_UNITS` is the speed of light in the canonical solar
unit system (AU per `year/2π`), exposed for callers that need the
constant directly.
"""

from apsis._native.gr import C_SOLAR_UNITS, PostNewtonian1PN

__all__ = ["C_SOLAR_UNITS", "PostNewtonian1PN"]
