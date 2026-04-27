"""Type stubs for the Rust extension module ``apsis._native``.

Mirrors the runtime API exposed by the PyO3 ``#[pymodule]`` definition
in ``crates/apsis-py/src/lib.rs``. This file is the source of truth
for static type checking; ``apsis/__init__.pyi`` re-exports from here.

The ``_native`` namespace itself is private — researchers should
``import apsis`` and access the same classes through the package
façade rather than reaching into ``apsis._native``.

Each declaration here is checked against the runtime by the smoke
tests in ``tests/test_basic.py`` (signature matching) and against
the user-facing API by ``mypy --strict`` over the ``examples/``
directory in CI.
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
    """Point-mass body with kinematics, mass, softening, and material class.

    Bodies are constructed through one of the nine material factories
    (``Body.star``, ``Body.rocky``, ...). All factories share a
    kwargs-only signature; position and velocity default to zero.

    Builder methods (:meth:`at`, :meth:`with_velocity`,
    :meth:`with_density`, :meth:`unsoftened`) return a new ``Body`` —
    bodies are value-typed on the Python side.
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

    def unsoftened(self) -> Body:
        """Drop softening to zero, restoring exact 1/r gravity for this body."""

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
        """Material class slug (e.g. ``"star"``, ``"rocky"``, ``"gas_giant"``)."""
    @property
    def luminosity(self) -> float: ...

    def __repr__(self) -> str: ...

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
        integrator: IntegratorKind | str,
        dt: float,
        epsilon: float | None = None,
        exact_gravity: bool = False,
    ) -> None: ...

    # ── Run loop ──
    def step(self) -> None:
        """Advance by exactly one integrator step."""

    def integrate_for(self, duration: float) -> int:
        """Advance for ``duration`` time units. Returns step count."""

    def integrate_until(self, t_end: float) -> int:
        """Advance until ``t >= t_end``. Returns step count."""

    def sample(self, duration: float, n_samples: int) -> Trajectory:
        """Integrate for ``duration`` time units while recording ``n_samples``
        evenly spaced snapshots, returning a :class:`Trajectory` of NumPy
        arrays. The first sample captures the pre-integration state; the
        last sample is taken after advancing to ``current_t + duration``
        (within one adaptive sub-step under IAS15). Sampling advances the
        system state in place.
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

    def __repr__(self) -> str: ...
