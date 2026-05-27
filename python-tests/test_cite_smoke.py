"""Smoke tests for ``apsis.System.cite()`` — the Python wrapper of
``System::cite`` that emits a BibTeX ``@software`` block citing every
registered operator crate.

These tests register real first-party operators (``apsis.gr`` and
``apsis.radiation``) so the output asserted here is the same shape a
paper.bib consumer would get end-to-end. The fine-grained renderer
tests live in the Rust crate
(``crates/apsis/src/physics/integrator/citation.rs``).
"""

from __future__ import annotations

import apsis
from apsis import gr, radiation


def _sun_only_system() -> apsis.System:
    return apsis.System(
        bodies=[apsis.Body.star(mass=1.0)],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )


def test_cite_empty_when_no_operators_registered() -> None:
    sys = _sun_only_system()
    assert sys.cite() == ""


def test_cite_includes_apsis_1pn_software_entry() -> None:
    sys = _sun_only_system()
    sys.add_hamiltonian_perturbation(
        gr.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    )
    block = sys.cite()
    assert "@software{apsis-1pn_" in block
    assert "title   = {apsis-1pn}" in block
    assert "First-post-Newtonian Schwarzschild correction" in block
    assert "kernel_requirements: exact_and_smooth" in block


def test_cite_emits_one_entry_per_unique_crate_for_paper_md_case() -> None:
    sys = _sun_only_system()
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
    block = sys.cite()
    assert block.count("@software{") == 2
    apsis_1pn_at = block.find("@software{apsis-1pn_")
    apsis_rad_at = block.find("@software{apsis-radiation_")
    assert apsis_1pn_at >= 0 and apsis_rad_at >= 0
    assert apsis_1pn_at < apsis_rad_at, "registration order must be preserved"
    assert "Cargo.lock blake3:" in block
