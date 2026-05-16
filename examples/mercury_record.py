"""Mercury 1PN run recorded as an Apsis Record.

Builds the Mercury + 1PN system from `examples/mercury_perihelion.py`,
attaches a `RecordHook` (default policy: bookends + events), and runs
for 100 Mercury orbits. The resulting `mercury_1pn.apsis` carries the
full provenance — apsis git sha, BLAKE3 hash of `Cargo.lock`, every
registered operator with crate name + version + checksum + declared
`KernelRequirements`, integrator config, unit system, seed — plus
the initial and final body state. Pair it with the lockfile to
reproduce this run bit-exactly.

Run::

    python examples/mercury_record.py
    python examples/inspect_record.py mercury_1pn.apsis
"""

from __future__ import annotations

import math

import apsis
from apsis.gr import PostNewtonian1PN

A_MERCURY = 0.387_098
E_MERCURY = 0.205_63
M_MERCURY = 1.660_114e-7
M_SUN = 1.0
MU = M_SUN + M_MERCURY
N_ORBITS = 100


def main() -> None:
    sun = apsis.Body.star(mass=M_SUN)
    mercury_x = A_MERCURY * (1.0 - E_MERCURY)
    mercury_v = math.sqrt(MU * (1.0 + E_MERCURY) / (A_MERCURY * (1.0 - E_MERCURY)))
    mercury = (
        apsis.Body.rocky(mass=M_MERCURY)
        .at(mercury_x, 0.0)
        .with_velocity(0.0, mercury_v)
    )

    period = 2.0 * math.pi * math.sqrt(A_MERCURY**3 / MU)
    dt = period / 1000.0

    sys = apsis.System([sun, mercury], units=apsis.units.SOLAR_CANONICAL)
    sys.set_integrator(apsis.IntegratorKind.IAS15)
    sys.set_dt(dt)
    sys.add_hamiltonian_perturbation(
        PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    )

    # Default policy: initial + final bookends + events only.
    sys.attach_record("mercury_1pn.apsis", seed=42)

    duration = N_ORBITS * period
    while sys.t < duration:
        sys.step()

    # Dropping `sys` flushes the RecordHook trailer.
    del sys

    print("Record written to mercury_1pn.apsis")


if __name__ == "__main__":
    main()
