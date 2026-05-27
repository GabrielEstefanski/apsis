"""Produce the canonical Apsis Record TOML header embedded in
``paper.md`` Appendix A.

Run from the workspace root after ``maturin develop --release``:

    python paper/scripts/dump_record_header.py

The script constructs a Sun + Mercury two-body system in solar
canonical units, registers ``apsis-1pn`` and ``apsis-radiation``,
attaches an Apsis Record with bookend snapshots, advances one step
to flush the header, then reopens the record read-only and prints
the raw TOML preamble. The Appendix A example block is replaced
verbatim with this output whenever the header schema, operator
metadata, or unit-system labels change.
"""

from __future__ import annotations

import os
import tempfile

import apsis
from apsis import gr, radiation


def build_system_with_both_operators() -> apsis.System:
    sun = apsis.Body.star(mass=1.0)
    mercury = (
        apsis.Body.rocky(mass=1.66e-7)
        .at((0.387, 0.0))
        .with_velocity((0.0, 1.61))
    )
    sys = apsis.System(
        bodies=[sun, mercury],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    sys.add_hamiltonian_perturbation(
        gr.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    )
    sys.add_hamiltonian_perturbation(
        radiation.RadiationPressure.from_raw_betas(
            source=0,
            betas=[0.0, 0.05],
            units=apsis.units.SOLAR_CANONICAL,
        )
    )
    return sys


def main() -> None:
    sys = build_system_with_both_operators()
    fd, path = tempfile.mkstemp(suffix=".apsis")
    os.close(fd)
    try:
        sys.attach_record(path, seed=42)
        sys.step()
        sys.finish()
        record = apsis.Record(path)
        print(record.header, end="")
    finally:
        os.unlink(path)


if __name__ == "__main__":
    main()
