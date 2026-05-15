"""First post-Newtonian Schwarzschild correction (Anderson et al. 1975).

```python
from apsis.gr import PostNewtonian1PN
sys.add_hamiltonian_perturbation(
    PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
)
```

`C_SOLAR_UNITS` is the speed of light in the canonical solar
unit system (AU per `year/2π`).

⚠ Critical precondition
-----------------------

1PN is derived around the bit-exact Newtonian potential. Default
``apsis.System(...)`` uses ``NewtonKernel::exact()`` (ε = 0); the
registration is silent. Attaching 1PN on top of a softened kernel
(opt-in via ``System::with_kernel(NewtonKernel::new(ε > 0))`` on the
Rust side, typically for cluster work — currently not reachable from
Python) makes the numerical apsidal precession from softening alone
~2 × 10³ larger than the relativistic signal at Mercury's orbit,
with the wrong sign — energy and angular momentum remain conserved
at machine precision while the trajectory is physically wrong. The
kernel-requirement check at registration emits a structured warning
when this happens.
"""

from apsis._native.gr import C_SOLAR_UNITS, PostNewtonian1PN

__all__ = ["C_SOLAR_UNITS", "PostNewtonian1PN"]
