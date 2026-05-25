"""Mercury 1PN run recorded as an Apsis Record.

Builds the Mercury + 1PN system from ``python-examples/mercury_perihelion.py``
and attaches a record writer via :py:meth:`apsis.System.attach_record`.
The resulting ``mercury_1pn.apsis`` carries the full provenance — apsis
git sha, BLAKE3 hash of ``Cargo.lock``, every registered operator with
crate name + version + checksum + declared ``KernelRequirements``,
integrator config, unit system, seed — the initial and final body
state, periodic snapshots, ``Diagnostic`` frames (ΔE/E + ΔLz/Lz over
time), and ``ResumeState`` frames enabling mid-run restore via the
Rust-side ``apsis::records::restore_into``. Pair it with the lockfile
to reproduce this run bit-exactly.

Run::

    python python-examples/mercury_record.py
    python python-examples/inspect_record.py mercury_1pn.apsis

Paper-grade reproducibility note: rebuild with a clean working tree
(no ``-dirty`` suffix on ``git_sha``) before running publication
artifacts. A ``-dirty`` sha records that the binary was compiled
from uncommitted source state — the verifier cannot recover that
state from the lockfile + commit pair.
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
N_ORBITS = 10


def main() -> None:
    sun = apsis.Body.star(mass=M_SUN)
    mercury_x = A_MERCURY * (1.0 - E_MERCURY)
    mercury_v = math.sqrt(MU * (1.0 + E_MERCURY) / (A_MERCURY * (1.0 - E_MERCURY)))
    mercury = (
        apsis.Body.rocky(mass=M_MERCURY)
        .at((mercury_x, 0.0))
        .with_velocity((0.0, mercury_v))
    )

    period = 2.0 * math.pi * math.sqrt(A_MERCURY**3 / MU)
    dt = period / 1000.0

    sys = apsis.System(
        bodies=[sun, mercury],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=dt,
    )
    sys.add_hamiltonian_perturbation(
        PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
    )

    # Snapshots every 50 steps + Diagnostic frames every 50 steps for
    # post-hoc verification of the precession claim. capture_resume
    # opts into ResumeState frames so the run can be branched mid-flight.
    sys.attach_record(
        "mercury_1pn.apsis",
        seed=42,
        every_steps=50,
        diagnostics_every_steps=50,
        capture_resume=True,
    )

    sys.integrate_for(N_ORBITS * period)
    sys.finish()  # flush trailer; idempotent, also fires on GC

    print("Record written to mercury_1pn.apsis")


if __name__ == "__main__":
    main()
