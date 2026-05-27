"""Tests for `apsis.Body.with_name()` and the `name` property."""

from __future__ import annotations

import apsis


def test_name_is_none_on_fresh_body() -> None:
    sun = apsis.Body.star(mass=1.0)
    assert sun.name is None


def test_with_name_sets_name() -> None:
    sun = apsis.Body.star(mass=1.0).with_name("Sun")
    assert sun.name == "Sun"


def test_with_name_chains_with_other_builders() -> None:
    mercury = (
        apsis.Body.rocky(mass=1.66e-7)
        .at((0.387, 0.0))
        .with_velocity((0.0, 1.61))
        .with_name("Mercury")
    )
    assert mercury.name == "Mercury"
    assert mercury.position == (0.387, 0.0, 0.0)


def test_name_survives_system_round_trip() -> None:
    sys = apsis.System(
        bodies=[
            apsis.Body.star(mass=1.0).with_name("Sun"),
            apsis.Body.rocky(mass=1.66e-7).at((0.387, 0.0)).with_name("Mercury"),
        ],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    names = [b.name for b in sys.bodies]
    assert names == ["Sun", "Mercury"]


def test_unnamed_body_gets_auto_filled_after_registration() -> None:
    sys = apsis.System(
        bodies=[apsis.Body.star(mass=1.0)],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    assert sys.bodies[0].name == "Body 1"
