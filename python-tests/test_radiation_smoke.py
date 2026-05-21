"""Smoke and integration tests for `apsis.radiation`."""

from __future__ import annotations

import pytest

import apsis
from apsis import radiation


def _sun_dust_system() -> tuple[apsis.System, apsis.Perturbation]:
    sun = apsis.Body.star(mass=1.0)
    dust = apsis.Body.rocky(mass=1e-15).at((1.0, 0.0)).with_velocity((0.0, 1.0))
    sys = apsis.System(
        bodies=[sun, dust],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    p = radiation.RadiationPressure.from_raw_betas(
        source=0,
        betas=[0.0, 0.1],
        units=apsis.units.SOLAR_CANONICAL,
    )
    return sys, p


def test_submodule_imports() -> None:
    assert hasattr(radiation, "RadiationPressure")


def test_from_raw_betas_returns_apsis_perturbation() -> None:
    p = radiation.RadiationPressure.from_raw_betas(
        source=0,
        betas=[0.0, 0.05],
        units=apsis.units.SOLAR_CANONICAL,
    )
    assert isinstance(p, apsis.Perturbation)
    assert "RadiationPressure" in p.label


def test_from_raw_betas_requires_keyword_only_arguments() -> None:
    with pytest.raises(TypeError):
        radiation.RadiationPressure.from_raw_betas(  # type: ignore[misc]
            0, [0.0, 0.05], apsis.units.SOLAR_CANONICAL,
        )


def test_attach_to_system_advances_time() -> None:
    sys, p = _sun_dust_system()
    sys.add_hamiltonian_perturbation(p)
    sys.integrate_for(0.5)
    assert sys.t >= 0.5
    assert sys.steps > 0


def test_double_attach_raises_a_clear_error() -> None:
    sys, p = _sun_dust_system()
    sys.add_hamiltonian_perturbation(p)
    with pytest.raises(ValueError, match="already been attached"):
        sys.add_hamiltonian_perturbation(p)
