"""Smoke and integration tests for `apsis.central`."""

from __future__ import annotations

from typing import Any, cast

import pytest

import apsis
from apsis import central


def _sun_planet_system() -> tuple[apsis.System, list[apsis.Body]]:
    sun = apsis.Body.star(mass=1.0)
    planet = apsis.Body.rocky(mass=1e-7).at((0.387, 0.0)).with_velocity((0.0, 1.61))
    bodies = [sun, planet]
    sys = apsis.System(
        bodies=bodies,
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    return sys, bodies


def test_submodule_imports() -> None:
    assert hasattr(central, "CentralForce")


def test_from_raw_returns_apsis_perturbation() -> None:
    p = central.CentralForce.from_raw(
        source=0,
        a_central=5e-9,
        gamma=-3.0,
        units=apsis.units.SOLAR_CANONICAL,
    )
    assert isinstance(p, apsis.Perturbation)
    assert "CentralForce" in p.label


def test_from_apsidal_rate_returns_apsis_perturbation() -> None:
    _, bodies = _sun_planet_system()
    p = central.CentralForce.from_apsidal_rate(
        source=0,
        target=1,
        omega_dot=5e-9,
        gamma=-3.0,
        bodies=bodies,
        units=apsis.units.SOLAR_CANONICAL,
    )
    assert isinstance(p, apsis.Perturbation)
    assert "CentralForce" in p.label


def test_factories_require_keyword_only_arguments() -> None:
    with pytest.raises(TypeError):
        cast(Any, central.CentralForce.from_raw)(
            0, 5e-9, -3.0, apsis.units.SOLAR_CANONICAL,
        )


def test_attach_to_system_advances_time() -> None:
    sys, _ = _sun_planet_system()
    p = central.CentralForce.from_raw(
        source=0,
        a_central=5e-9,
        gamma=-3.0,
        units=apsis.units.SOLAR_CANONICAL,
    )
    sys.add_hamiltonian_perturbation(p)
    sys.integrate_for(0.5)
    assert sys.t >= 0.5
    assert sys.steps > 0


def test_double_attach_raises_a_clear_error() -> None:
    sys, _ = _sun_planet_system()
    p = central.CentralForce.from_raw(
        source=0,
        a_central=5e-9,
        gamma=-3.0,
        units=apsis.units.SOLAR_CANONICAL,
    )
    sys.add_hamiltonian_perturbation(p)
    with pytest.raises(ValueError, match="already been attached"):
        sys.add_hamiltonian_perturbation(p)
