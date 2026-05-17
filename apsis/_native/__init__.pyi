"""Type stubs for the Rust extension module ``apsis._native``.

Mirrors the runtime API exposed by the PyO3 ``#[pymodule]`` in
``crates/apsis-python/src/lib.rs``. The ``apsis`` user-facing
package re-exports from here. Submodule stubs (``gr``, future
``radiation``/``central``) live as sibling files in this directory.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import numpy.typing as npt

__version__: str

# A 1-D ``float64`` NumPy array, used for trajectory time and energy axes.
_F64Array1D = npt.NDArray[np.float64]
# A 2-D ``float64`` NumPy array, used for per-body trajectory state arrays.
_F64Array2D = npt.NDArray[np.float64]

# ── IntegratorKind ────────────────────────────────────────────────────────────

class IntegratorKind:
    """Numerical integration scheme applied to the simulation's body state.

    See ``docs/integrator.md`` in the project repository for the
    per-integrator contract (execution profile, force-model
    determinism requirement, selection rubric).
    """

    IAS15: IntegratorKind
    YOSHIDA4: IntegratorKind
    VELOCITY_VERLET: IntegratorKind
    WISDOM_HOLMAN: IntegratorKind

    @property
    def name(self) -> str:
        """Canonical Python enum name (``"IAS15"``, ``"YOSHIDA4"``, ...)."""

    @property
    def slug(self) -> str:
        """Lower-case canonical slug (``"ias15"``, ``"yoshida4"``, ...).

        The same string accepted by every ``integrator=`` kwarg, by the
        core ``IntegratorKind::FromStr`` impl, and by the project's
        ``run.toml`` config files.
        """

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

# ── Body ──────────────────────────────────────────────────────────────────────

class Body:
    """Point-mass body with kinematics, mass, softening, and a binding-
    layer material slug.

    Bodies are constructed through one of the nine preset factories
    (``Body.star``, ``Body.rocky``, ...). All factories share a
    kwargs-only signature; position and velocity default to zero.

    The underlying core body holds no material taxonomy field —
    physical defaults from the preset are applied once at construction
    and never re-read. The Python wrapper still surfaces a
    ``body.material`` slug for ergonomic introspection on freshly-built
    bodies; bodies retrieved via ``System.bodies`` lose that tag
    (they report ``"body"``) because the round-trip drops the
    binding-layer state.

    Builder methods (:meth:`at`, :meth:`with_velocity`,
    :meth:`with_density`) return a new ``Body`` — bodies are
    value-typed on the Python side.
    """

    @staticmethod
    def star(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Main-sequence luminous body."""

    @staticmethod
    def brown_dwarf(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Brown dwarf — sub-stellar, deuterium-burning regime."""

    @staticmethod
    def white_dwarf(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """White dwarf — compact stellar remnant."""

    @staticmethod
    def gas_giant(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Gas giant — Jupiter-class hydrogen/helium envelope."""

    @staticmethod
    def ice_giant(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Ice giant — Neptune-class water/methane envelope."""

    @staticmethod
    def rocky(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Rocky body — terrestrial planet or large rocky satellite."""

    @staticmethod
    def icy(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Icy body — water-dominated composition (outer satellites, KBOs)."""

    @staticmethod
    def asteroid(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Asteroid — rocky minor body."""

    @staticmethod
    def comet(
        *,
        mass: float,
        position: tuple[float, float] | None = None,
        velocity: tuple[float, float] | None = None,
        softening: float | None = None,
    ) -> Body:
        """Comet — volatile-rich minor body."""

    # ── Builder ──
    def at(self, position: tuple[float, float]) -> Body:
        """Place the body at ``position = (x, y)``. Returns a new ``Body``."""

    def with_velocity(self, velocity: tuple[float, float]) -> Body:
        """Set the body's velocity to ``(vx, vy)``. Returns a new ``Body``."""

    def with_density(self, density: float) -> Body:
        """Override the material-default density. Recomputes physical radius."""

    # ── Properties ──
    @property
    def mass(self) -> float: ...
    @property
    def position(self) -> tuple[float, float]: ...
    @property
    def velocity(self) -> tuple[float, float]: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def vx(self) -> float: ...
    @property
    def vy(self) -> float: ...
    @property
    def softening(self) -> float: ...
    @property
    def density(self) -> float: ...
    @property
    def radius(self) -> float: ...
    @property
    def material(self) -> str:
        """Construction-time preset slug (e.g. ``"star"``, ``"rocky"``,
        ``"gas_giant"``). Round-tripped bodies (those returned by
        ``System.bodies``) report ``"body"`` because the slug is
        binding-layer state and the core body holds no taxonomy
        reference."""
    @property
    def luminosity(self) -> float: ...
    @property
    def q_pr(self) -> float:
        """Radiation-pressure receiver coefficient ``Q_pr``. Positive on
        radiation receivers (asteroids, comets, icy grains); zero on
        emitters and large planets."""

    def __repr__(self) -> str: ...

# ── Errors ────────────────────────────────────────────────────────────────────

class UnitSystemMismatchError(Exception):
    """Raised by ``System.add_*_perturbation`` when the operator's
    declared :class:`UnitSystem` disagrees with the System's own.

    Attributes (set by the Rust binding for programmatic recovery):
      - ``operator``: ``str`` — operator identifier (e.g.
        ``"apsis_1pn::PostNewtonian1PN"``).
      - ``operator_units``: ``str`` — Display of the operator's
        :class:`UnitSystem`.
      - ``system_units``: ``str`` — Display of the System's
        :class:`UnitSystem`.

    The string message itself includes both unit systems and a hint;
    catching by class is the recommended path:

    .. code-block:: python

        try:
            sys.add_hamiltonian_perturbation(p)
        except apsis.UnitSystemMismatchError as e:
            log(f"{e.operator}: {e.operator_units} vs {e.system_units}")
    """

    operator: str
    operator_units: str
    system_units: str


# ── System ────────────────────────────────────────────────────────────────────

class System:
    """Orchestrator for the simulation: bodies, integrator, run loop.

    Construction is kwargs-only — every dial a researcher is likely to
    set is named at the call site, so the meaning of any code that
    builds a ``System`` is clear without consulting the docs.
    """

    def __init__(
        self,
        *,
        bodies: Sequence[Body],
        units: UnitSystem,
        integrator: IntegratorKind | str,
        dt: float,
        epsilon: float | None = None,
    ) -> None: ...

    # ── Run loop ──
    def step(self) -> None:
        """Advance by exactly one integrator step."""

    def integrate_for(self, duration: float) -> int:
        """Advance for ``duration`` time units. Returns step count."""

    def integrate_until(self, t_end: float) -> int:
        """Advance until ``t >= t_end``. Returns step count."""

    def finish(self) -> None:
        """Close attached records and fire each hook's lifecycle-end callback.
        Idempotent; also fires automatically on garbage collection."""

    def attach_record(
        self,
        path: str,
        *,
        seed: int | None = None,
        every_steps: int | None = None,
        every_time: float | None = None,
        dense: bool = False,
        diagnostics_every_steps: int | None = None,
        diagnostics_every_time: float | None = None,
    ) -> None:
        """Attach an Apsis Record writer to this system. Subsequent steps
        write to ``path``; the file is closed (with a trailer) on
        :meth:`finish` or garbage collection.

        Snapshot cadence — at most one of ``every_steps`` /
        ``every_time`` / ``dense`` may be set; default is bookend
        snapshots + events only.

        Diagnostic cadence — at most one of ``diagnostics_every_steps``
        / ``diagnostics_every_time`` may be set. When enabled, the
        record carries ``Diagnostic`` frames (``ΔE/E``, ``ΔLz/Lz``)
        readable via :meth:`Record.diagnostics`."""

    def sample(
        self,
        *,
        times: Sequence[float] | _F64Array1D | None = None,
        duration: float | None = None,
        n_samples: int | None = None,
    ) -> Trajectory:
        """Record the state at a set of target times, returning a
        :class:`Trajectory` of NumPy arrays.

        Two invocation forms — exactly one must be used:

        Explicit times (primary)::

            traj = sys.sample(times=np.linspace(0.0, 100.0, 1024))
            traj = sys.sample(times=np.logspace(-3, 2, 200))
            traj = sys.sample(times=[0.0, 1.0, 10.0, 100.0])

        Evenly spaced (convenience)::

            traj = sys.sample(duration=10.0, n_samples=128)

        ``times`` must be non-empty, finite, monotonically non-decreasing,
        and ``times[0] >= sys.t``. There is no interpolation: each
        recorded sample is the integrator's actual output at (or just
        past) the requested time. Sampling advances ``sys.t`` to
        ``traj.t[-1]``.
        """

    # ── Mutators ──
    def recenter_com(self) -> None:
        """Translate every body so the centre of mass is at the origin."""

    # ── Read-only properties (O(1)) ──
    @property
    def t(self) -> float: ...
    @property
    def steps(self) -> int: ...
    @property
    def bodies(self) -> list[Body]: ...
    @property
    def energy(self) -> float: ...
    @property
    def energy_delta(self) -> float:
        """Relative energy drift, ``(E - E_0) / |E_0|``."""
    @property
    def kinetic_energy(self) -> float: ...
    @property
    def potential_energy(self) -> float: ...
    @property
    def lz(self) -> float:
        """Total z-component of angular momentum."""
    @property
    def lz_delta(self) -> float:
        """Relative angular-momentum drift."""
    @property
    def dt(self) -> float:
        """Current controller time step (mutates with IAS15; constant for fixed-step schemes)."""
    @property
    def integrator(self) -> IntegratorKind: ...

    @property
    def units(self) -> UnitSystem:
        """The unit system this system was constructed against. Frozen."""

    # ── Adaptive controller counters (zero for fixed-step integrators) ──
    @property
    def substeps(self) -> int: ...
    @property
    def step_rejections(self) -> int: ...
    @property
    def picard_stagnations(self) -> int: ...
    @property
    def shrink_grow_cycles(self) -> int: ...
    @property
    def picard_iters(self) -> int: ...
    @property
    def degraded_steps(self) -> int: ...
    @property
    def force_evaluations(self) -> int:
        """Estimated total force evaluations (`steps × force_evals_per_step`)."""

    @property
    def stats(self) -> Stats:
        """Frozen snapshot of cumulative scalar diagnostics."""

    @property
    def adaptive_stats(self) -> AdaptiveStats | None:
        """Controller counters for adaptive integrators; ``None`` for fixed-step."""

    def add_hamiltonian_perturbation(self, perturbation: object) -> None:
        """Attach a Hamiltonian-class perturbation built by a downstream crate.

        Accepts an :class:`apsis.Perturbation` (pure-Python wrapper carrying
        an opaque capsule from a perturbation crate). The perturbation is
        consumed by the call; build a fresh one for each system.
        Kernel-requirement violations (e.g. attaching a 1PN correction to a
        softened-gravity system) emit structured warnings.
        """

    def __repr__(self) -> str: ...

# ── Trajectory ────────────────────────────────────────────────────────────────

class Trajectory:
    """Dense recording of a simulation interval, returned by :meth:`System.sample`.

    All arrays are ``float64`` NumPy ``ndarray`` views materialised once at
    construction time and handed out as zero-copy reads thereafter. The 1-D
    arrays (:attr:`t`, :attr:`energy`) have shape ``(n_samples,)``; the 2-D
    arrays (:attr:`x`, :attr:`y`, :attr:`vx`, :attr:`vy`) have shape
    ``(n_samples, n_bodies)`` with the body index on the second axis.

    A typical plot of body ``k``'s orbit is
    ``plt.plot(traj.x[:, k], traj.y[:, k])``; the energy-conservation
    diagnostic is
    ``plt.plot(traj.t, (traj.energy - traj.energy[0]) / abs(traj.energy[0]))``.
    """

    @property
    def n_samples(self) -> int:
        """Number of recorded samples (length of the time axis)."""

    @property
    def n_bodies(self) -> int:
        """Number of bodies tracked in this trajectory."""

    @property
    def t(self) -> _F64Array1D:
        """Sample times, shape ``(n_samples,)``."""

    @property
    def x(self) -> _F64Array2D:
        """Body x-coordinates, shape ``(n_samples, n_bodies)``."""

    @property
    def y(self) -> _F64Array2D:
        """Body y-coordinates, shape ``(n_samples, n_bodies)``."""

    @property
    def vx(self) -> _F64Array2D:
        """Body x-velocities, shape ``(n_samples, n_bodies)``."""

    @property
    def vy(self) -> _F64Array2D:
        """Body y-velocities, shape ``(n_samples, n_bodies)``."""

    @property
    def energy(self) -> _F64Array1D:
        """Total mechanical energy at each sample, shape ``(n_samples,)``."""

    @property
    def dt(self) -> _F64Array1D:
        """Controller dt at each sample. Constant for fixed-step integrators;
        traces the adaptive controller for IAS15."""

    @property
    def energy_drift(self) -> _F64Array1D:
        """Relative energy drift ``δE/E₀`` at each sample, shape ``(n_samples,)``."""

    @property
    def lz_drift(self) -> _F64Array1D:
        """Relative angular-momentum drift ``δLz/Lz₀`` at each sample."""

    def __repr__(self) -> str: ...

# ── Stats / AdaptiveStats ─────────────────────────────────────────────────────

class Stats:
    """Frozen snapshot of cumulative scalar diagnostics; returned by ``System.stats``."""

    @property
    def t(self) -> float: ...
    @property
    def steps(self) -> int: ...
    @property
    def dt(self) -> float: ...
    @property
    def energy(self) -> float: ...
    @property
    def energy_drift(self) -> float: ...
    @property
    def kinetic_energy(self) -> float: ...
    @property
    def potential_energy(self) -> float: ...
    @property
    def lz(self) -> float: ...
    @property
    def lz_drift(self) -> float: ...
    @property
    def integrator(self) -> IntegratorKind: ...
    @property
    def force_evaluations(self) -> int: ...
    def __repr__(self) -> str: ...

class AdaptiveStats:
    """Frozen snapshot of an adaptive integrator's controller counters.
    Returned by ``System.adaptive_stats`` for IAS15; ``None`` for fixed-step."""

    @property
    def substeps(self) -> int: ...
    @property
    def rejections(self) -> int: ...
    @property
    def rejections_picard(self) -> int: ...
    @property
    def rejections_truncation(self) -> int: ...
    @property
    def picard_iters(self) -> int: ...
    @property
    def picard_stagnations(self) -> int: ...
    @property
    def shrink_grow_cycles(self) -> int: ...
    @property
    def degraded(self) -> int: ...
    def __repr__(self) -> str: ...

# ── UnitSystem ────────────────────────────────────────────────────────────────

class UnitSystem:
    """Closed system of units for length, time, and mass.

    Construct via one of the named factories (:meth:`canonical`, :meth:`henon`,
    :meth:`si`, :meth:`solar`, :meth:`cgs`) or :meth:`custom`. Once chosen,
    a ``UnitSystem`` is immutable — there is no setter for any of its fields.

    All physics inputs (positions, velocities, masses, dt) passed to a
    :class:`System` are interpreted in the canonical units of the supplied
    ``UnitSystem``. The wrapper performs no dimensional checking — passing a
    value in the wrong unit is a silent physical error, matching REBOUND's
    convention.
    """

    @staticmethod
    def canonical() -> UnitSystem:
        """Hénon-style canonical N-body units: ``G = 1`` by construction."""

    @staticmethod
    def henon() -> UnitSystem:
        """Alias for :meth:`canonical` using the literature name."""

    @staticmethod
    def si() -> UnitSystem:
        """SI units: metre, second, kilogram."""

    @staticmethod
    def solar() -> UnitSystem:
        """Solar-system canonical units: AU, year, solar mass. ``G ≈ 39.478``."""

    @staticmethod
    def cgs() -> UnitSystem:
        """CGS units: centimetre, second, gram."""

    @staticmethod
    def custom(*, length_m: float, time_s: float, mass_kg: float) -> UnitSystem:
        """Build a custom unit system from explicit SI scales.

        All three scales must be strictly positive and finite; zero, negative,
        infinite, or NaN values raise ``ValueError`` at the boundary.
        """

    @property
    def length_scale_si(self) -> float:
        """SI metres per code-unit length."""
    @property
    def time_scale_si(self) -> float:
        """SI seconds per code-unit time."""
    @property
    def mass_scale_si(self) -> float:
        """SI kilograms per code-unit mass."""
    @property
    def g(self) -> float:
        """Newtonian gravitational constant in this system's canonical units."""
    @property
    def length_label(self) -> str: ...
    @property
    def time_label(self) -> str: ...
    @property
    def mass_label(self) -> str: ...

    def length_to_si(self, x: float) -> float: ...
    def length_from_si(self, x: float) -> float: ...
    def time_to_si(self, x: float) -> float: ...
    def time_from_si(self, x: float) -> float: ...
    def mass_to_si(self, x: float) -> float: ...
    def mass_from_si(self, x: float) -> float: ...

    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

# ── Record (Apsis Record reader) ──────────────────────────────────────────────

# Tuple shapes for events surfaced by :meth:`Record.events`:
#  - ("collision", t, body_a, body_b, distance)
#  - ("escape",    t, body, radius)
_CollisionTuple = tuple[str, float, int, int, float]
_EscapeTuple = tuple[str, float, int, float]
_EventTuple = _CollisionTuple | _EscapeTuple

class Record:
    """Read-only view of an Apsis Record (``.apsis``) file.

    Writing is mediated by :meth:`System.attach_record` (which builds
    a header from the live System state and registers a Rust-side
    writer hook). The reader surface here is for downstream analysis
    and post-hoc inspection.
    """

    def __init__(self, path: str) -> None:
        """Open a record. Raises ``IOError`` on missing magic, version
        mismatch, truncation, missing trailer, or missing bookend."""

    @property
    def header(self) -> str:
        """Raw TOML header as a string. Parse with ``tomllib`` (3.11+)
        for structured access."""

    def events(self) -> list[_EventTuple]:
        """All event tuples in time order."""

    def snapshot_count(self) -> int:
        """Number of dense Snapshot frames (initial bookend + per-policy
        snapshots + final bookend)."""

    def diagnostics(self) -> list[tuple[float, float, float]]:
        """Diagnostic frames as ``(t, d_energy_rel, d_lz_rel)`` tuples in
        time order. Empty when the record was written without a
        diagnostic cadence."""


# ── units submodule ───────────────────────────────────────────────────────────

class _UnitsModule:
    """Stub for ``apsis.units`` — exposes the named-factory singletons and
    the SI constants. Mirrors the ``register`` function in
    ``crates/apsis-py/src/units.rs``."""

    CANONICAL: UnitSystem
    HENON: UnitSystem
    SI: UnitSystem
    SOLAR: UnitSystem
    SOLAR_CANONICAL: UnitSystem
    CGS: UnitSystem

    G_SI: float
    AU_M: float
    YR_S: float
    MSUN_KG: float
    CM_M: float
    G_KG: float

    UnitSystem: type[UnitSystem]

units: _UnitsModule
