"""Smoke and integration tests for the apsis-1pn Python binding.

Two layers:

- Boundary checks (fast): the binding exposes ``PostNewtonian1PN``,
  ``C_SOLAR_UNITS``, and the factories return objects compatible with
  ``apsis.System.add_hamiltonian_perturbation``.
- Behaviour checks (also fast): a fresh ``Perturbation`` consumed by
  ``add_hamiltonian_perturbation`` cannot be reused; the integration
  with 1PN attached actually advances time without raising.

The Mercury perihelion precession-rate is *not* re-validated here —
that's the parent crate's job (``crates/apsis-1pn/tests/mercury_precession_gate.rs``)
and it requires an hours-long release-mode run. This file proves only
that the binding correctly plumbs the Rust impl through the Python API.
"""

from __future__ import annotations

import math

import pytest

import apsis
import apsis_1pn


def _mercury_two_body_with_1pn() -> tuple[apsis.System, apsis.Perturbation]:
    """Build the Sun + Mercury system used by every behaviour test below.

    Both the System and the perturbation use ``apsis.units.SOLAR_CANONICAL``
    (G = 1 + AU + yr/(2π) + M☉) — the apsis-1pn validation portfolio
    convention. The registration check in ``add_hamiltonian_perturbation``
    requires the two unit systems to match exactly.
    """
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
    p = apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    return sys, p


# ── Module surface ────────────────────────────────────────────────────────────


def test_module_imports() -> None:
    assert apsis_1pn is not None
    assert apsis_1pn.__version__


def test_c_solar_units_constant_exposed() -> None:
    """``C_SOLAR_UNITS`` is the compile-time-derived speed of light in the
    canonical solar unit system. ~10065 AU / (yr/2π)."""
    assert isinstance(apsis_1pn.C_SOLAR_UNITS, float)
    assert math.isclose(apsis_1pn.C_SOLAR_UNITS, 10065.13, rel_tol=1e-3)


# ── PostNewtonian1PN factories ───────────────────────────────────────────────


def test_for_units_constructor_returns_apsis_perturbation() -> None:
    """The factory result is an ``apsis.Perturbation`` — type identity is
    shared across the apsis and apsis_1pn extensions because the class
    is pure-Python in ``apsis/__init__.py``."""
    p = apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL)
    assert isinstance(p, apsis.Perturbation)
    assert "PostNewtonian1PN" in p.label
    assert "for_units" in p.label


def test_from_raw_c_factory_accepts_explicit_speed_of_light() -> None:
    """``from_raw_c`` accepts an explicit ``c`` and pins the operator to
    the supplied UnitSystem for the registration-time check."""
    p = apsis_1pn.PostNewtonian1PN.from_raw_c(c=10000.0, units=apsis.units.SOLAR_CANONICAL)
    assert isinstance(p, apsis.Perturbation)
    assert "c=10000" in p.label


def test_from_raw_c_rejects_invalid_speed_of_light() -> None:
    """Zero, negative, infinite, and NaN values for ``c`` are rejected."""
    for bad in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError, match="c"):
            apsis_1pn.PostNewtonian1PN.from_raw_c(c=bad, units=apsis.units.SOLAR_CANONICAL)


def test_from_raw_c_requires_keyword_only_arguments() -> None:
    """``c`` and ``units`` are kwarg-only — positional args are a
    contract violation."""
    with pytest.raises(TypeError):
        apsis_1pn.PostNewtonian1PN.from_raw_c(10000.0, apsis.units.SOLAR_CANONICAL)  # type: ignore[misc]


# ── Integration with apsis.System ────────────────────────────────────────────


def test_attach_to_system_advances_time() -> None:
    """Attaching 1PN and integrating runs without raising; the system time
    advances. Behavioural correctness of the precession is the parent
    crate's gate, not this binding's."""
    sys, p = _mercury_two_body_with_1pn()

    sys.add_hamiltonian_perturbation(p)
    sys.integrate_for(1.0)

    assert sys.t >= 1.0
    assert sys.steps > 0


def test_double_attach_raises_a_clear_error() -> None:
    """The single-consume contract surfaces a precise error message
    rather than a use-after-free or silent double-attach."""
    sys, p = _mercury_two_body_with_1pn()
    sys.add_hamiltonian_perturbation(p)

    with pytest.raises(ValueError, match="already been attached"):
        sys.add_hamiltonian_perturbation(p)


def test_two_systems_need_two_perturbations() -> None:
    """The recommended pattern: build a fresh perturbation per system."""
    sys_a, p_a = _mercury_two_body_with_1pn()
    sys_b, p_b = _mercury_two_body_with_1pn()

    sys_a.add_hamiltonian_perturbation(p_a)
    sys_b.add_hamiltonian_perturbation(p_b)

    sys_a.integrate_for(0.5)
    sys_b.integrate_for(0.5)

    assert sys_a.t >= 0.5
    assert sys_b.t >= 0.5


# Kernel-requirement violation coverage lives in the parent crate
# (`crates/apsis-1pn/tests/mercury_precession_gate.rs::plummer_kernel_under_1pn_triggers_diagnostic`).
# The Python binding does not yet expose a softened-kernel constructor
# (`System::with_kernel`), so the warning path is currently unreachable
# from Python; the assertion would test absence-of-warning rather than
# emission, which is meaningless under the new exact-default model.
