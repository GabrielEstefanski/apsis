"""Central-potential perturbations (Tamayo et al. 2019, Pattern B).

Adds a power-law radial force `∝ r^γ` from a central body, on top of
Newtonian gravity. Two construction patterns:

```python
from apsis.central import CentralForce

# Pattern A — explicit (a_central, γ):
sys.add_hamiltonian_perturbation(
    CentralForce.from_raw(
        source=0, a_central=5e-9, gamma=-3.0,
        units=apsis.units.SOLAR_CANONICAL,
    ),
)

# Pattern B — invert a desired apsidal rate ω̇ for the target body:
sys.add_hamiltonian_perturbation(
    CentralForce.from_apsidal_rate(
        source=0, target=1, omega_dot=5e-9, gamma=-3.0,
        bodies=[sun, planet],
        units=apsis.units.SOLAR_CANONICAL,
    ),
)
```

`γ = -3` recovers the Schwarzschild-effective near-circular
precession; arbitrary `γ` covers the broader central-potential
family.
"""

from apsis._native.central import CentralForce

__all__ = ["CentralForce"]
