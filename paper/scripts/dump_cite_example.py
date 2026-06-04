"""Produce the canonical ``apsis.System.cite()`` output embedded in
``paper.md`` §Citable operator stack.

Run from the workspace root after ``maturin develop --release``:

    python paper/scripts/dump_cite_example.py

The script constructs a Sun + Mercury-like two-body system in solar
canonical units, registers ``apsis-1pn`` (the Hamiltonian
first-post-Newtonian correction) and ``apsis-radiation`` (radiation
pressure on the secondary), and writes the resulting BibTeX block to
stdout. The paper.md example block is replaced verbatim with this
output whenever the renderer format or operator metadata changes.
"""

from __future__ import annotations

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
    print(sys.cite(), end="")


if __name__ == "__main__":
    main()
