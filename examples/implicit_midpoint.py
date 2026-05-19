"""ImplicitMidpoint integrator on Mercury + 1PN — symplectic + A-stable.

Same physical scenario as ``mercury_perihelion.py`` but with the
single-stage Gauss-Legendre implicit symplectic method
(``integrator="implicit_midpoint"``) instead of the adaptive
Gauss-Radau IAS15.

Demonstrates:

- **Energy bounded:** symplectic A-stable methods preserve a shifted
  Hamiltonian; ``|dE/E|`` stays bounded over the full integration
  window. The asymptotic floor is dominated by the implicit solver's
  per-step convergence tolerance, not by ``ε_machine``. This is the
  defining property of the method and is asserted.
- **Phase-drift trade-off:** at 2nd-order accuracy and this ``dt``,
  the orbit phase picks up spurious precession from integrator
  truncation that adds to the GR signal. The example prints the gap
  rather than asserting against it. Tighter ``dt`` or a higher-order
  integrator (IAS15, Yoshida4, WHFast) brings the precession
  measurement down to GR.

Run::

    python examples/implicit_midpoint.py
"""

from __future__ import annotations

import math

import apsis
from apsis.gr import C_SOLAR_UNITS, PostNewtonian1PN

A_MERCURY = 0.387_098
E_MERCURY = 0.205_63
M_MERCURY = 1.660_114e-7
M_SUN = 1.0

N_ORBITS = 100
DT = 1e-4   # ImplicitMidpoint is 2nd-order, so smaller dt than IAS15
MU = M_SUN + M_MERCURY


def perihelion_longitude(rx: float, ry: float, vx: float, vy: float, mu: float) -> float:
    h = rx * vy - ry * vx
    r = math.hypot(rx, ry)
    ex = (vy * h) / mu - rx / r
    ey = (-vx * h) / mu - ry / r
    return math.atan2(ey, ex)


def unwrap_radians(d: float) -> float:
    while d > math.pi:
        d -= 2.0 * math.pi
    while d <= -math.pi:
        d += 2.0 * math.pi
    return d


def main() -> None:
    sun = apsis.Body.star(mass=M_SUN)
    r_peri = A_MERCURY * (1.0 - E_MERCURY)
    v_peri = math.sqrt(MU * (2.0 / r_peri - 1.0 / A_MERCURY))
    mercury = apsis.Body.rocky(mass=M_MERCURY).at((r_peri, 0.0)).with_velocity((0.0, v_peri))

    sys = apsis.System(
        bodies=[sun, mercury],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="implicit_midpoint",
        dt=DT,
    )
    sys.add_hamiltonian_perturbation(
        PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
    )

    s, m = sys.bodies
    omega_0 = perihelion_longitude(m.x - s.x, m.y - s.y, m.vx - s.vx, m.vy - s.vy, MU)

    period = 2.0 * math.pi * math.sqrt(A_MERCURY**3 / MU)
    sys.integrate_for(period * N_ORBITS)

    s, m = sys.bodies
    omega_1 = perihelion_longitude(m.x - s.x, m.y - s.y, m.vx - s.vx, m.vy - s.vy, MU)
    measured = unwrap_radians(omega_1 - omega_0)

    c = C_SOLAR_UNITS
    predicted = 6.0 * math.pi * M_SUN / (c * c * A_MERCURY * (1.0 - E_MERCURY**2)) * N_ORBITS

    arcsec = (180.0 / math.pi) * 3600.0

    print(f"Mercury 1PN — {N_ORBITS} orbits under ImplicitMidpoint (dt = {DT})")
    print(f"  predicted Δω (GR)                 = {predicted * arcsec:+.4f} arcsec")
    print(f"  measured  Δω (GR + truncation)    = {measured  * arcsec:+.4f} arcsec")
    print(f"  excess from 2nd-order truncation  = {(measured - predicted) * arcsec:+.4f} arcsec")
    de_rel = sys.energy_delta
    assert de_rel is not None, "Mercury 1PN is well-conditioned; energy_delta must be Some"
    print(f"  |dE/E| = {abs(de_rel):.3e}  (symplectic: bounded)")

    # Symplectic A-stable: bounded |dE/E|, floor set by the implicit
    # solver's convergence tolerance (~10⁻¹³ measured).
    assert abs(de_rel) < 1e-10, f"|dE/E| = {abs(de_rel):.3e} not bounded"


if __name__ == "__main__":
    main()
