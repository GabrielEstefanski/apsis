"""Smoke and integration tests for `apsis.gr` — first post-Newtonian
correction, consolidated into the apsis distribution from the legacy
apsis-1pn-py package.

Boundary checks (fast): the submodule exposes ``PostNewtonian1PN``,
``C_SOLAR_UNITS``, and the factories return objects compatible with
``apsis.System.add_hamiltonian_perturbation``. Behaviour checks (also
fast): a fresh ``Perturbation`` consumed by ``add_hamiltonian_perturbation``
cannot be reused; the integration with 1PN attached actually advances
time without raising.

The Mercury perihelion precession-rate is *not* re-validated here —
that's the parent crate's job (``crates/apsis-1pn/tests/mercury_precession_gate.rs``)
and it requires a release-mode run.
"""

from __future__ import annotations

import math

import pytest

import apsis
from apsis import gr


def _mercury_two_body_with_1pn() -> tuple[apsis.System, apsis.Perturbation]:
    """Sun + Mercury under SOLAR_CANONICAL with a 1PN perturbation."""
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
    p = gr.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    return sys, p


# ── Module surface ────────────────────────────────────────────────────────────


def test_submodule_imports() -> None:
    assert gr is not None
    assert hasattr(gr, "PostNewtonian1PN")
    assert hasattr(gr, "C_SOLAR_UNITS")


def test_c_solar_units_constant_exposed() -> None:
    """``C_SOLAR_UNITS`` is the speed of light in canonical solar units
    (~10065 AU per year/(2π))."""
    assert isinstance(gr.C_SOLAR_UNITS, float)
    assert math.isclose(gr.C_SOLAR_UNITS, 10065.13, rel_tol=1e-3)


# ── PostNewtonian1PN factories ───────────────────────────────────────────────


def test_for_units_returns_apsis_perturbation() -> None:
    p = gr.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    assert isinstance(p, apsis.Perturbation)
    assert "PostNewtonian1PN" in p.label
    assert "for_units" in p.label


def test_from_raw_c_accepts_explicit_speed_of_light() -> None:
    p = gr.PostNewtonian1PN.from_raw_c(c=10000.0, units=apsis.units.SOLAR_CANONICAL)
    assert isinstance(p, apsis.Perturbation)
    assert "c=10000" in p.label


def test_from_raw_c_rejects_invalid_speed_of_light() -> None:
    for bad in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError, match="c"):
            gr.PostNewtonian1PN.from_raw_c(c=bad, units=apsis.units.SOLAR_CANONICAL)


def test_from_raw_c_requires_keyword_only_arguments() -> None:
    with pytest.raises(TypeError):
        gr.PostNewtonian1PN.from_raw_c(10000.0, apsis.units.SOLAR_CANONICAL)  # type: ignore[misc]


# ── Integration with apsis.System ────────────────────────────────────────────


def test_attach_to_system_advances_time() -> None:
    sys, p = _mercury_two_body_with_1pn()
    sys.add_hamiltonian_perturbation(p)
    sys.integrate_for(1.0)
    assert sys.t >= 1.0
    assert sys.steps > 0


def test_double_attach_raises_a_clear_error() -> None:
    sys, p = _mercury_two_body_with_1pn()
    sys.add_hamiltonian_perturbation(p)
    with pytest.raises(ValueError, match="already been attached"):
        sys.add_hamiltonian_perturbation(p)


def test_two_systems_need_two_perturbations() -> None:
    sys_a, p_a = _mercury_two_body_with_1pn()
    sys_b, p_b = _mercury_two_body_with_1pn()
    sys_a.add_hamiltonian_perturbation(p_a)
    sys_b.add_hamiltonian_perturbation(p_b)
    sys_a.integrate_for(0.5)
    sys_b.integrate_for(0.5)
    assert sys_a.t >= 0.5
    assert sys_b.t >= 0.5
