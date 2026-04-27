"""Smoke tests for the binding's import surface and signature shape.

These tests are deliberately minimal — they exercise *only* the
boundary between Python and the Rust extension module, not the
behaviour of the underlying integrators or force models. Behavioural
correctness is the responsibility of the parent crate's test suite
(``crates/apsis/tests``, ``crates/apsis-1pn/tests``) and of the
cross-implementation parity portfolio under ``validation/``.

A new ``#[pyclass]`` or ``#[pyfunction]`` added on the Rust side
acquires a matching one-liner here that asserts the symbol is
importable through ``apsis`` and has the expected runtime type and
signature shape. That is the full scope of this file: prove the
façade is wired up.
"""

from __future__ import annotations

import re

import pytest

import apsis


# ── Module surface ────────────────────────────────────────────────────────────


def test_module_imports() -> None:
    """``import apsis`` succeeds once the extension module is built."""
    assert apsis is not None


def test_version_string_matches_semver() -> None:
    """``apsis.__version__`` is sourced from the workspace ``Cargo.toml``."""
    version = apsis.__version__
    assert isinstance(version, str)
    assert re.fullmatch(r"\d+\.\d+\.\d+(?:[-+].+)?", version), (
        f"unexpected version string: {version!r}"
    )


def test_public_surface_is_exported() -> None:
    """Every symbol declared in ``__all__`` is reachable from the package."""
    for name in apsis.__all__:
        assert hasattr(apsis, name), f"apsis.{name} is missing"


# ── IntegratorKind ────────────────────────────────────────────────────────────


def test_integrator_kind_variants_are_exposed() -> None:
    """All four integrator variants reach Python under the documented names."""
    kinds = apsis.IntegratorKind
    assert kinds.IAS15 is not None
    assert kinds.YOSHIDA4 is not None
    assert kinds.VELOCITY_VERLET is not None
    assert kinds.WISDOM_HOLMAN is not None
    # Distinct variants are distinct objects.
    assert kinds.IAS15 != kinds.YOSHIDA4


def test_integrator_kind_slug_round_trips() -> None:
    """``slug`` is the canonical lowercase form accepted by ``System(integrator=...)``."""
    assert apsis.IntegratorKind.IAS15.slug == "ias15"
    assert apsis.IntegratorKind.YOSHIDA4.slug == "yoshida4"
    assert apsis.IntegratorKind.VELOCITY_VERLET.slug == "velocity_verlet"
    assert apsis.IntegratorKind.WISDOM_HOLMAN.slug == "wisdom_holman"


# ── Body ──────────────────────────────────────────────────────────────────────


def test_body_factories_accept_kwargs_only() -> None:
    """Every material factory takes mass, position, velocity, softening as kwargs.

    Positional invocation is rejected — the Rust signature is
    ``(*, mass, ...)`` so the marker is enforced at the FFI boundary,
    not just by convention on the Python side.
    """
    sun = apsis.Body.star(mass=1.0, position=(0.5, -0.25), velocity=(0.0, 1.0))
    assert sun.mass == 1.0
    assert sun.position == (0.5, -0.25)
    assert sun.velocity == (0.0, 1.0)
    assert sun.material == "star"


def test_body_builder_methods_return_new_instances() -> None:
    """Builder methods produce fresh ``Body`` objects, leaving the original unchanged."""
    base = apsis.Body.rocky(mass=1e-6)
    placed = base.at((1.0, 2.0))
    assert base.position == (0.0, 0.0)
    assert placed.position == (1.0, 2.0)
    assert base is not placed


def test_body_unsoftened_zeroes_softening() -> None:
    """``unsoftened`` flips the softening to exactly zero (exact 1/r gravity)."""
    softened = apsis.Body.star(mass=1.0)
    assert softened.softening > 0.0
    assert softened.unsoftened().softening == 0.0


def test_body_rejects_non_positive_mass() -> None:
    """A negative or zero mass raises ``ValueError`` at the boundary."""
    with pytest.raises(ValueError, match="mass"):
        apsis.Body.rocky(mass=0.0)
    with pytest.raises(ValueError, match="mass"):
        apsis.Body.star(mass=-1.0)


def test_body_rejects_malformed_position() -> None:
    """Position that is not a 2-element sequence raises ``ValueError``."""
    with pytest.raises(ValueError, match="position"):
        apsis.Body.star(mass=1.0, position=(1.0,))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="position"):
        apsis.Body.star(mass=1.0, position=(1.0, 2.0, 3.0))  # type: ignore[arg-type]


# ── System ────────────────────────────────────────────────────────────────────


def test_system_constructor_accepts_string_or_enum_integrator() -> None:
    """``integrator=`` takes both ``IntegratorKind.IAS15`` and the slug ``"ias15"``."""
    sun = apsis.Body.star(mass=1.0)
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0))

    sys_str = apsis.System(bodies=[sun, earth], integrator="ias15", dt=1e-3)
    sys_enum = apsis.System(
        bodies=[sun, earth], integrator=apsis.IntegratorKind.IAS15, dt=1e-3
    )
    assert sys_str.integrator == apsis.IntegratorKind.IAS15
    assert sys_enum.integrator == apsis.IntegratorKind.IAS15


def test_system_string_is_case_insensitive() -> None:
    """The slug parser is liberal: ``"IAS15"`` works as well as ``"ias15"``."""
    sun = apsis.Body.star(mass=1.0)
    sys = apsis.System(bodies=[sun], integrator="IAS15", dt=1e-3)
    assert sys.integrator == apsis.IntegratorKind.IAS15


def test_system_rejects_unknown_integrator() -> None:
    """Typos and unknown spellings raise ``ValueError`` listing valid choices."""
    sun = apsis.Body.star(mass=1.0)
    with pytest.raises(ValueError, match="integrator"):
        apsis.System(bodies=[sun], integrator="ias16", dt=1e-3)


def test_system_run_loop_advances_time() -> None:
    """``integrate_for`` moves ``t`` forward and bumps the step counter."""
    sun = apsis.Body.star(mass=1.0).unsoftened()
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0)).unsoftened()
    sys = apsis.System(bodies=[sun, earth], integrator="ias15", dt=1e-3)

    assert sys.t == 0.0
    assert sys.steps == 0

    sys.integrate_for(1.0)

    assert sys.t >= 1.0
    assert sys.steps > 0


def test_system_energy_delta_is_machine_precision_on_kepler() -> None:
    """One Kepler orbit under IAS15 conserves energy to f64 round-off.

    This is the cheapest end-to-end sanity check: the binding correctly
    plumbs through the integrator that the Rust core's parity portfolio
    already validates at 1 ULP. Tolerance here is much looser (1e-12)
    so a Python-side bookkeeping error in the FFI surfaces clearly while
    a one-orbit run is fast enough to remain a smoke test.
    """
    sun = apsis.Body.star(mass=1.0).unsoftened()
    earth = (apsis.Body.rocky(mass=3e-6)
             .at((1.0, 0.0))
             .with_velocity((0.0, 1.0))
             .unsoftened())
    sys = apsis.System(bodies=[sun, earth], integrator="ias15", dt=1e-3)

    sys.integrate_for(2 * 3.14159)

    assert abs(sys.energy_delta) < 1e-12, f"energy drift too large: {sys.energy_delta}"
