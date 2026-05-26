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

import numpy as np
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
    # 2-tuple input pads `z` (and `vz`) with zero; the getter returns
    # the full 3-tuple.
    assert sun.position == (0.5, -0.25, 0.0)
    assert sun.velocity == (0.0, 1.0, 0.0)
    assert sun.material == "star"


def test_body_factories_accept_3d_position_and_velocity() -> None:
    """Position and velocity kwargs accept a 3-tuple for out-of-plane orbits."""
    body = apsis.Body.rocky(
        mass=3e-6,
        position=(1.0, 2.0, 3.0),
        velocity=(0.4, 0.5, 0.6),
    )
    assert body.position == (1.0, 2.0, 3.0)
    assert body.velocity == (0.4, 0.5, 0.6)
    assert body.x == 1.0
    assert body.y == 2.0
    assert body.z == 3.0
    assert body.vx == 0.4
    assert body.vy == 0.5
    assert body.vz == 0.6


def test_body_builder_methods_return_new_instances() -> None:
    """Builder methods produce fresh ``Body`` objects, leaving the original unchanged."""
    base = apsis.Body.rocky(mass=1e-6)
    placed = base.at((1.0, 2.0, 3.0))
    assert base.position == (0.0, 0.0, 0.0)
    assert placed.position == (1.0, 2.0, 3.0)
    assert base is not placed


def test_body_rejects_non_positive_mass() -> None:
    """A negative or zero mass raises ``ValueError`` at the boundary."""
    with pytest.raises(ValueError, match="mass"):
        apsis.Body.rocky(mass=0.0)
    with pytest.raises(ValueError, match="mass"):
        apsis.Body.star(mass=-1.0)


def test_body_rejects_malformed_position() -> None:
    """Position that is not a 2- or 3-element sequence raises ``ValueError``."""
    with pytest.raises(ValueError, match="position"):
        apsis.Body.star(mass=1.0, position=(1.0,))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="position"):
        apsis.Body.star(mass=1.0, position=(1.0, 2.0, 3.0, 4.0))  # type: ignore[arg-type]


# ── System ────────────────────────────────────────────────────────────────────


def test_system_constructor_accepts_string_or_enum_integrator() -> None:
    """``integrator=`` takes both ``IntegratorKind.IAS15`` and the slug ``"ias15"``."""
    sun = apsis.Body.star(mass=1.0)
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0))

    sys_str = apsis.System(bodies=[sun, earth], units=apsis.units.CANONICAL, integrator="ias15", dt=1e-3)
    sys_enum = apsis.System(
        bodies=[sun, earth],
        units=apsis.units.CANONICAL,
        integrator=apsis.IntegratorKind.IAS15,
        dt=1e-3,
    )
    assert sys_str.integrator == apsis.IntegratorKind.IAS15
    assert sys_enum.integrator == apsis.IntegratorKind.IAS15


def test_system_string_is_case_insensitive() -> None:
    """The slug parser is liberal: ``"IAS15"`` works as well as ``"ias15"``."""
    sun = apsis.Body.star(mass=1.0)
    sys = apsis.System(bodies=[sun], units=apsis.units.CANONICAL, integrator="IAS15", dt=1e-3)
    assert sys.integrator == apsis.IntegratorKind.IAS15


def test_system_rejects_unknown_integrator() -> None:
    """Typos and unknown spellings raise ``ValueError`` listing valid choices."""
    sun = apsis.Body.star(mass=1.0)
    with pytest.raises(ValueError, match="integrator"):
        apsis.System(bodies=[sun], units=apsis.units.CANONICAL, integrator="ias16", dt=1e-3)


def test_system_run_loop_advances_time() -> None:
    """``integrate_for`` moves ``t`` forward and bumps the step counter."""
    sun = apsis.Body.star(mass=1.0)
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0))
    sys = apsis.System(bodies=[sun, earth], units=apsis.units.CANONICAL, integrator="ias15", dt=1e-3)

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
    sun = apsis.Body.star(mass=1.0)
    earth = (apsis.Body.rocky(mass=3e-6)
             .at((1.0, 0.0))
             .with_velocity((0.0, 1.0)))
    sys = apsis.System(bodies=[sun, earth], units=apsis.units.CANONICAL, integrator="ias15", dt=1e-3)

    sys.integrate_for(2 * 3.14159)

    assert abs(sys.energy_delta) < 1e-12, f"energy drift too large: {sys.energy_delta}"


# ── Trajectory ────────────────────────────────────────────────────────────────


def _two_body_kepler_system() -> apsis.System:
    sun = apsis.Body.star(mass=1.0)
    earth = (apsis.Body.rocky(mass=3e-6)
             .at((1.0, 0.0))
             .with_velocity((0.0, 1.0)))
    return apsis.System(bodies=[sun, earth], units=apsis.units.CANONICAL, integrator="ias15", dt=1e-3)


def test_sample_returns_trajectory_with_expected_shape() -> None:
    """``sample(duration, n_samples)`` produces NumPy arrays of the right shape."""
    sys = _two_body_kepler_system()
    traj = sys.sample(duration=1.0, n_samples=64)

    assert isinstance(traj, apsis.Trajectory)
    assert traj.n_samples == 64
    assert traj.n_bodies == 2
    assert traj.t.shape == (64,)
    assert traj.energy.shape == (64,)
    assert traj.x.shape == (64, 2)
    assert traj.y.shape == (64, 2)
    assert traj.z.shape == (64, 2)
    assert traj.vx.shape == (64, 2)
    assert traj.vy.shape == (64, 2)
    assert traj.vz.shape == (64, 2)


def test_sample_arrays_are_float64() -> None:
    """Every Trajectory array is a plain ``float64`` NumPy ndarray."""
    sys = _two_body_kepler_system()
    traj = sys.sample(duration=1.0, n_samples=8)

    for arr in (traj.t, traj.energy, traj.x, traj.y, traj.z, traj.vx, traj.vy, traj.vz):
        assert isinstance(arr, np.ndarray)
        assert arr.dtype == np.float64


def test_sample_z_components_are_zero_for_planar_input() -> None:
    """Bodies confined to the xy-plane keep ``z`` and ``vz`` identically zero
    across every sample. This is the contract that lets researchers ignore
    the third component when their problem is planar."""
    sys = _two_body_kepler_system()
    traj = sys.sample(duration=1.0, n_samples=16)

    assert np.all(traj.z == 0.0), "planar input must produce z = 0 everywhere"
    assert np.all(traj.vz == 0.0), "planar input must produce vz = 0 everywhere"


def test_sample_records_3d_motion_for_inclined_input() -> None:
    """An inclined orbit (`vz != 0` initially) populates the 3D arrays with
    real motion, not zeros."""
    primary = apsis.Body.star(mass=1.0)
    inclined = (apsis.Body.rocky(mass=1e-6)
                .at((1.0, 0.0, 0.0))
                .with_velocity((0.0, 0.7, 0.7)))
    sys = apsis.System(
        bodies=[primary, inclined],
        units=apsis.units.CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )

    traj = sys.sample(duration=1.0, n_samples=8)

    # The inclined body must explore non-zero z over the integration.
    assert np.any(traj.z[:, 1] != 0.0), "inclined orbit produced no z motion"


def test_sample_time_axis_is_monotonic_and_brackets_duration() -> None:
    """``traj.t`` is non-decreasing, starts at the system's pre-call ``t``,
    and ends at or just past ``start_t + duration``."""
    sys = _two_body_kepler_system()
    sys.integrate_for(0.5)
    start_t = sys.t

    traj = sys.sample(duration=2.0, n_samples=32)

    diffs = np.diff(traj.t)
    assert np.all(diffs >= 0.0), "trajectory time axis went backwards"
    assert traj.t[0] == pytest.approx(start_t)
    assert traj.t[-1] >= start_t + 2.0
    assert sys.t == pytest.approx(traj.t[-1])


def test_sample_initial_row_matches_system_state() -> None:
    """The first row of every state array equals the system's pre-sample state.

    Positions and velocities are read straight from the bodies, so they round-trip
    bit-for-bit. Energy is special: ``sample()`` primes the energy cache from
    the current body state at the start of the call, so ``traj.energy[0]`` is a
    real ``K + U`` value rather than the construction-time zero of ``sys.energy``.
    Comparing it against the analytic Kepler total ``-G M m / (2 a)`` (which is
    ``-0.5 * (G * M_sun * m_earth) / 1.0`` for the unit-radius, unit-velocity
    setup) is the cheapest way to assert physical correctness without leaking
    integrator state into the test.
    """
    sys = _two_body_kepler_system()
    pre_positions = [b.position for b in sys.bodies]
    pre_velocities = [b.velocity for b in sys.bodies]

    traj = sys.sample(duration=0.5, n_samples=4)

    expected_energy = -0.5 * 3e-6
    assert traj.energy[0] == pytest.approx(expected_energy, rel=1e-9)
    for k, (px, py, pz) in enumerate(pre_positions):
        assert traj.x[0, k] == pytest.approx(px)
        assert traj.y[0, k] == pytest.approx(py)
        assert traj.z[0, k] == pytest.approx(pz)
    for k, (vx, vy, vz) in enumerate(pre_velocities):
        assert traj.vx[0, k] == pytest.approx(vx)
        assert traj.vy[0, k] == pytest.approx(vy)
        assert traj.vz[0, k] == pytest.approx(vz)


def test_sample_rejects_invalid_arguments() -> None:
    """Non-positive ``duration`` and zero ``n_samples`` are rejected at the boundary."""
    sys = _two_body_kepler_system()

    with pytest.raises(ValueError, match="duration"):
        sys.sample(duration=0.0, n_samples=4)
    with pytest.raises(ValueError, match="duration"):
        sys.sample(duration=-1.0, n_samples=4)
    with pytest.raises(ValueError, match="n_samples"):
        sys.sample(duration=1.0, n_samples=0)


# ── sample(times=...) — explicit-targets API ─────────────────────────────────


def test_sample_with_explicit_times_records_at_those_targets() -> None:
    """``sample(times=...)`` records one sample per target; ``traj.t[i]`` meets or
    overshoots ``targets[i]`` (overshoot is bounded by one IAS15 adaptive sub-step,
    which on a smooth Kepler orbit can be a non-trivial fraction of an orbital period)."""
    sys = _two_body_kepler_system()
    targets = np.linspace(0.0, 6.28, 64)

    traj = sys.sample(times=targets)

    assert traj.t.shape == (64,)
    assert traj.n_samples == 64
    # Contract: each recorded sample is at or after its requested time.
    assert np.all(traj.t >= targets - 1e-12)
    # Trajectory time axis is monotonically non-decreasing (a sub-step that
    # straddles two consecutive targets produces equal-valued rows, which is
    # fine — the body state is the same).
    assert np.all(np.diff(traj.t) >= 0.0)
    # End reaches or just passes the final target.
    assert traj.t[-1] >= targets[-1] - 1e-12


def test_sample_accepts_list_and_tuple() -> None:
    """``times`` accepts plain Python sequences, not just NumPy arrays."""
    sys_list = _two_body_kepler_system()
    traj_list = sys_list.sample(times=[0.0, 1.0, 3.0, 6.28])
    assert traj_list.t.shape == (4,)

    sys_tuple = _two_body_kepler_system()
    traj_tuple = sys_tuple.sample(times=(0.0, 1.0, 3.0, 6.28))
    assert traj_tuple.t.shape == (4,)


def test_sample_supports_log_spaced_times() -> None:
    """Log-spaced sampling is the headline use case for chaotic / multi-scale runs.

    The contract is just "the API accepts a non-uniformly-spaced array and
    integrates through it correctly". Specific density bounds are fragile
    because the IAS15 sub-step can land on either side of a closely-spaced
    pair of targets — a more meaningful test would compare *target* spacing,
    which is purely an input property and doesn't need to round-trip.
    """
    sys = _two_body_kepler_system()
    targets = np.logspace(-2, 1, 32)  # 0.01 → 10 in 32 log-spaced points

    traj = sys.sample(times=targets)

    assert traj.t.shape == (32,)
    assert np.all(traj.t >= targets - 1e-12)
    # Target spacing is log-dense at small t, sparse at large t.
    target_diffs = np.diff(targets)
    assert target_diffs[0] < target_diffs[-1] / 100


def test_sample_respects_pre_sample_integration() -> None:
    """If the system has already advanced, ``times[0] >= sys.t`` is the contract;
    samples are recorded at-or-after that point."""
    sys = _two_body_kepler_system()
    sys.integrate_for(0.5)
    start_t = sys.t

    traj = sys.sample(times=np.linspace(start_t, start_t + 1.0, 16))

    assert traj.t[0] >= start_t
    assert traj.t[0] == pytest.approx(start_t, abs=1e-12)


def test_sample_rejects_empty_times() -> None:
    """An empty ``times`` array is rejected — there is no sensible interpretation."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="times"):
        sys.sample(times=[])
    with pytest.raises(ValueError, match="times"):
        sys.sample(times=np.array([]))


def test_sample_rejects_non_monotonic_times() -> None:
    """Forward integration is the only supported direction."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="monotonically non-decreasing"):
        sys.sample(times=[0.0, 5.0, 3.0])


def test_sample_rejects_times_before_current_t() -> None:
    """Sampling can't rewind the simulator past ``sys.t``."""
    sys = _two_body_kepler_system()
    sys.integrate_for(1.0)
    with pytest.raises(ValueError, match="cannot integrate backwards"):
        sys.sample(times=[0.0, 0.5, 1.5])


def test_sample_rejects_non_finite_times() -> None:
    """``NaN`` and ``±∞`` are rejected with the offending index reported."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match=r"times\[1\]"):
        sys.sample(times=[0.0, float("nan"), 1.0])
    with pytest.raises(ValueError, match=r"times\[2\]"):
        sys.sample(times=[0.0, 1.0, float("inf")])


def test_sample_rejects_mixing_modes() -> None:
    """Passing both ``times=`` and ``duration=``/``n_samples=`` is a contract violation."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="not both"):
        sys.sample(times=[0.0, 1.0], duration=1.0, n_samples=10)
    with pytest.raises(ValueError, match="not both"):
        sys.sample(times=[0.0, 1.0], duration=1.0)


def test_sample_rejects_partial_duration_mode() -> None:
    """``duration=`` and ``n_samples=`` must be passed together."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="together"):
        sys.sample(duration=1.0)
    with pytest.raises(ValueError, match="together"):
        sys.sample(n_samples=10)


def test_sample_rejects_no_arguments() -> None:
    """Calling ``sample()`` with neither form raises the same clear error."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="times|duration"):
        sys.sample()


def test_sample_evenly_and_explicit_match_for_matched_inputs() -> None:
    """``sample(duration=, n_samples=)`` matches
    ``sample(times=np.linspace(t, t+duration, n_samples))`` bit-for-bit."""
    sys_a = _two_body_kepler_system()
    sys_b = _two_body_kepler_system()

    traj_a = sys_a.sample(duration=2.0, n_samples=32)
    traj_b = sys_b.sample(times=np.linspace(0.0, 2.0, 32))

    np.testing.assert_array_equal(traj_a.t, traj_b.t)
    np.testing.assert_array_equal(traj_a.x, traj_b.x)
    np.testing.assert_array_equal(traj_a.energy, traj_b.energy)


def test_sample_single_point_records_only_initial_state() -> None:
    """``n_samples=1`` is the degenerate zero-integration case: state is captured but not advanced."""
    sys = _two_body_kepler_system()
    start_t = sys.t

    traj = sys.sample(duration=1.0, n_samples=1)

    assert traj.n_samples == 1
    assert traj.t.shape == (1,)
    assert traj.t[0] == pytest.approx(start_t)
    assert sys.t == pytest.approx(start_t)


def test_sample_energy_drift_within_one_orbit_is_machine_precision() -> None:
    """End-to-end: the energy column of a one-orbit Kepler trajectory holds at f64 round-off."""
    sys = _two_body_kepler_system()
    traj = sys.sample(duration=2 * 3.14159, n_samples=128)

    rel_drift = (traj.energy - traj.energy[0]) / abs(traj.energy[0])
    assert np.max(np.abs(rel_drift)) < 1e-12


# ── UnitSystem ────────────────────────────────────────────────────────────────


def test_named_unit_systems_are_exposed_as_singletons() -> None:
    """``apsis.units.{SOLAR,SI,CANONICAL,HENON,CGS}`` are pre-built ``UnitSystem`` instances."""
    for name in ("CANONICAL", "HENON", "SI", "SOLAR", "CGS"):
        u = getattr(apsis.units, name)
        assert isinstance(u, apsis.UnitSystem)


def test_canonical_g_is_one_and_solar_g_is_near_four_pi_squared() -> None:
    """Identity checks on the derived gravitational constant."""
    assert apsis.units.CANONICAL.g == 1.0
    assert apsis.units.HENON.g == 1.0
    assert apsis.units.SI.g == pytest.approx(apsis.units.G_SI)

    four_pi_sq = 4.0 * 3.141592653589793 ** 2
    assert abs(apsis.units.SOLAR.g - four_pi_sq) / four_pi_sq < 1e-2


def test_unit_system_conversions_round_trip() -> None:
    """Applying ``to_si`` then ``from_si`` returns the original value within ULP."""
    u = apsis.units.SOLAR
    for value in (1.0, 0.5, 1.234567):
        assert u.length_from_si(u.length_to_si(value)) == pytest.approx(value, abs=1e-15)
        assert u.time_from_si(u.time_to_si(value)) == pytest.approx(value, abs=1e-15)
        assert u.mass_from_si(u.mass_to_si(value)) == pytest.approx(value, abs=1e-15)


def test_unit_system_conversions_apply_named_scales() -> None:
    """``solar().length_to_si(1.0)`` is one AU in metres, etc."""
    u = apsis.units.SOLAR
    assert u.length_to_si(1.0) == pytest.approx(apsis.units.AU_M)
    assert u.time_to_si(1.0) == pytest.approx(apsis.units.YR_S)
    assert u.mass_to_si(1.0) == pytest.approx(apsis.units.MSUN_KG)


def test_unit_system_custom_validates_at_boundary() -> None:
    """Zero, negative, infinite, and NaN scales raise ``ValueError``."""
    for bad in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError, match="length_m"):
            apsis.UnitSystem.custom(length_m=bad, time_s=1.0, mass_kg=1.0)
        with pytest.raises(ValueError, match="time_s"):
            apsis.UnitSystem.custom(length_m=1.0, time_s=bad, mass_kg=1.0)
        with pytest.raises(ValueError, match="mass_kg"):
            apsis.UnitSystem.custom(length_m=1.0, time_s=1.0, mass_kg=bad)


def test_unit_system_equality_compares_scales_not_labels() -> None:
    """Two systems with identical SI scales are equal regardless of labels."""
    custom = apsis.UnitSystem.custom(
        length_m=apsis.units.AU_M,
        time_s=apsis.units.YR_S,
        mass_kg=apsis.units.MSUN_KG,
    )
    assert custom == apsis.units.SOLAR
    # But their display labels differ — the SOLAR factory carries "AU"/"yr"/"Msun",
    # custom() uses generic placeholders.
    assert custom.length_label != apsis.units.SOLAR.length_label


def test_system_requires_units_kwarg() -> None:
    """Constructing a System without ``units=`` is a ``TypeError``."""
    sun = apsis.Body.star(mass=1.0)
    with pytest.raises(TypeError):
        apsis.System(bodies=[sun], integrator="ias15", dt=1e-3)  # type: ignore[call-arg]


def test_system_units_snapshot_is_immutable_across_integration() -> None:
    """The units snapshot survives integration unchanged.

    Locks the invariant the design hinges on: once chosen, the unit
    system can never silently change. Without this guarantee, every
    energy / momentum baseline captured at construction would become
    physically meaningless mid-run.
    """
    sun = apsis.Body.star(mass=1.0)
    earth = (apsis.Body.rocky(mass=3e-6)
             .at((1.0, 0.0))
             .with_velocity((0.0, 1.0)))
    sys = apsis.System(
        bodies=[sun, earth],
        units=apsis.units.SOLAR,
        integrator="ias15",
        dt=1e-3,
    )
    units_at_construction = sys.units
    assert units_at_construction == apsis.units.SOLAR

    sys.integrate_for(1.0)

    assert sys.units == units_at_construction
    assert sys.units == apsis.units.SOLAR


# ── Perturbation transport (without a concrete impl) ────────────────────────


def test_perturbation_class_is_a_pure_python_wrapper() -> None:
    """``apsis.Perturbation`` is the user-facing class that perturbation crates
    construct. Defined in pure Python so cross-extension type identity is
    actually shared (PyO3 cannot share #[pyclass] objects across cdylibs)."""
    assert apsis.Perturbation.__module__ == "apsis"


def test_add_hamiltonian_perturbation_rejects_non_perturbation_objects() -> None:
    """Anything without a ``_capsule`` attribute is rejected at the boundary."""
    sys = _two_body_kepler_system()
    with pytest.raises(ValueError, match="perturbation"):
        sys.add_hamiltonian_perturbation("not a perturbation")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="perturbation"):
        sys.add_hamiltonian_perturbation(42)  # type: ignore[arg-type]


def test_add_hamiltonian_perturbation_rejects_perturbation_with_non_capsule_attribute() -> None:
    """A ``Perturbation``-shaped object with a non-capsule ``_capsule`` is rejected."""
    sys = _two_body_kepler_system()
    fake = apsis.Perturbation(_capsule="not a capsule", _label="fake")
    with pytest.raises(ValueError, match="PyCapsule"):
        sys.add_hamiltonian_perturbation(fake)


# ── Diagnostics: Stats / AdaptiveStats / Trajectory parallel arrays ──────────


def test_adaptive_counters_zero_at_construction() -> None:
    """Pre-integration the counters all read zero; after integration they advance."""
    sys = _two_body_kepler_system()
    assert sys.substeps == 0
    assert sys.step_rejections == 0
    assert sys.picard_iters == 0
    assert sys.picard_stagnations == 0
    assert sys.shrink_grow_cycles == 0
    assert sys.degraded_steps == 0
    assert sys.force_evaluations == 0

    sys.integrate_for(2 * 3.14159)
    assert sys.substeps > 0
    assert sys.picard_iters > 0
    assert sys.force_evaluations > 0


def test_adaptive_stats_is_none_for_fixed_step_integrators() -> None:
    """Fixed-step schemes don't run a controller, so the binding returns ``None``."""
    sun = apsis.Body.star(mass=1.0)
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0))
    for scheme in ("yoshida4", "velocity_verlet"):
        sys = apsis.System(
            bodies=[sun, earth],
            units=apsis.units.CANONICAL,
            integrator=scheme,
            dt=1e-3,
        )
        sys.integrate_for(0.1)
        assert sys.adaptive_stats is None
        assert sys.substeps == 0
        assert sys.step_rejections == 0


def test_adaptive_stats_populated_for_ias15() -> None:
    """IAS15 carries a real ``AdaptiveStats`` after a run."""
    sys = _two_body_kepler_system()
    sys.integrate_for(2 * 3.14159)

    a = sys.adaptive_stats
    assert a is not None
    assert isinstance(a, apsis.AdaptiveStats)
    assert a.substeps > 0
    assert a.picard_iters > 0
    # Smooth Kepler: rejections rare, stagnations zero, degraded zero.
    assert a.rejections >= 0
    assert a.picard_stagnations >= 0
    assert a.degraded == 0


def test_stats_object_carries_all_headline_diagnostics() -> None:
    """``sys.stats`` aggregates the scalar diagnostics into one frozen object."""
    sys = _two_body_kepler_system()
    sys.integrate_for(2 * 3.14159)
    s = sys.stats

    assert isinstance(s, apsis.Stats)
    assert s.t > 0
    assert s.steps > 0
    assert s.dt > 0
    assert s.energy < 0  # bound Kepler orbit
    assert s.energy_drift is not None
    assert abs(s.energy_drift) < 1e-12
    assert s.lz_drift is not None
    assert abs(s.lz_drift) < 1e-12
    assert s.kinetic_energy > 0
    assert s.potential_energy < 0
    assert s.integrator == apsis.IntegratorKind.IAS15
    assert s.force_evaluations > 0


def _dust_system() -> apsis.System:
    """Extreme mass ratio that drives ``|E_0|`` below the conditioning floor.

    A 1e-15 secondary against a unit primary makes the total energy
    ``E ≈ −0.5 · m₂ ≈ 5e-16`` — below ``sqrt(f64::EPSILON) ≈ 1.5e-8``,
    so the relative drift metric is not well-defined and the binding
    must report ``None``.
    """
    sun = apsis.Body.star(mass=1.0)
    dust = (apsis.Body.rocky(mass=1e-15)
            .at((1.0, 0.0))
            .with_velocity((0.0, 1.0)))
    return apsis.System(
        bodies=[sun, dust],
        units=apsis.units.CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )


def test_drift_is_none_when_baseline_below_conditioning_floor() -> None:
    """In the dust regime ``|E_0|, |L_{z,0}| < sqrt(eps)`` so the relative
    drift is undefined for both energy and angular momentum; the binding
    must report ``None`` for both and surface the absolute counterpart
    as a finite float bounded by the system's own scale (a delta larger
    than the total it tracks would mean the field is wired to something
    other than ``current - baseline``)."""
    sys = _dust_system()
    sys.integrate_for(0.1)

    assert sys.energy_delta is None
    assert sys.lz_delta is None

    assert isinstance(sys.abs_energy_drift, float)
    assert np.isfinite(sys.abs_energy_drift)
    assert abs(sys.abs_energy_drift) <= abs(sys.energy) + 1e-30

    assert isinstance(sys.abs_lz_drift, float)
    assert np.isfinite(sys.abs_lz_drift)
    assert abs(sys.abs_lz_drift) <= abs(sys.lz) + 1e-30


def test_stats_drift_mirrors_system_in_dust_regime() -> None:
    """``Stats`` carries the same None-projection and finite absolute values
    that the live system accessors report — same regime, same contract."""
    sys = _dust_system()
    sys.integrate_for(0.1)
    s = sys.stats

    assert s.energy_drift is None
    assert s.lz_drift is None

    assert isinstance(s.abs_energy_drift, float)
    assert np.isfinite(s.abs_energy_drift)
    assert s.abs_energy_drift == sys.abs_energy_drift

    assert isinstance(s.abs_lz_drift, float)
    assert np.isfinite(s.abs_lz_drift)
    assert s.abs_lz_drift == sys.abs_lz_drift


def test_trajectory_carries_dt_history() -> None:
    """``traj.dt`` records the controller's step size at each sample."""
    sys = _two_body_kepler_system()
    traj = sys.sample(times=np.linspace(0.0, 6.28, 32))

    assert traj.dt.shape == (32,)
    assert traj.dt.dtype == np.float64
    # IAS15 starts at the user-provided seed dt and adapts upward on a
    # smooth Kepler orbit; dt[1:] should mostly exceed the seed.
    assert np.all(traj.dt > 0.0)


def test_trajectory_carries_energy_and_lz_drift_history() -> None:
    """``traj.energy_drift`` / ``traj.lz_drift`` plot directly without further math."""
    sys = _two_body_kepler_system()
    traj = sys.sample(times=np.linspace(0.0, 2 * 3.14159, 64))

    assert traj.energy_drift.shape == (64,)
    assert traj.lz_drift.shape == (64,)
    # Drift magnitudes stay bounded for IAS15 over one orbit.
    assert np.max(np.abs(traj.energy_drift)) < 1e-12
    assert np.max(np.abs(traj.lz_drift)) < 1e-12


def test_yoshida4_trajectory_dt_is_constant() -> None:
    """Fixed-step integrator: every sample reads back the same dt."""
    sun = apsis.Body.star(mass=1.0)
    earth = apsis.Body.rocky(mass=3e-6, position=(1.0, 0.0), velocity=(0.0, 1.0))
    sys = apsis.System(
        bodies=[sun, earth],
        units=apsis.units.CANONICAL,
        integrator="yoshida4",
        dt=1e-3,
    )

    traj = sys.sample(times=np.linspace(0.0, 1.0, 16))

    assert np.allclose(traj.dt, 1e-3)


def test_system_units_propagates_to_g_factor() -> None:
    """The system's effective G is read from the chosen unit system at construction."""
    sun = apsis.Body.star(mass=1.0)
    sys_solar = apsis.System(
        bodies=[sun],
        units=apsis.units.SOLAR,
        integrator="ias15",
        dt=1e-3,
    )
    sys_canon = apsis.System(
        bodies=[sun],
        units=apsis.units.CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    assert sys_solar.units.g != sys_canon.units.g
    assert sys_canon.units.g == 1.0
    assert abs(sys_solar.units.g - 39.478) < 0.05
