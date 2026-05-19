"""Mercury perihelion precession under 1PN — canonical GR test.

Sun + Mercury integrated under IAS15 with the apsis.gr first
post-Newtonian Schwarzschild correction attached. Measures the
accumulated longitude of periastron over 100 Mercury orbits and
compares against the closed-form Schwarzschild prediction
``Δω = 6π G M / (c² a (1 - e²))`` per orbit.

Run::

    python examples/mercury_perihelion.py
"""

from __future__ import annotations

import math

import apsis
from apsis.gr import C_SOLAR_UNITS, PostNewtonian1PN

# Mercury orbital elements in canonical solar units (G = 1, AU, year/2π, M_sun).
A_MERCURY = 0.387_098
E_MERCURY = 0.205_63
M_MERCURY = 1.660_114e-7
M_SUN = 1.0

N_ORBITS = 100
MU = M_SUN + M_MERCURY


def perihelion_longitude(rx: float, ry: float, vx: float, vy: float, mu: float) -> float:
    """Argument of periastron from a 2D state vector (planar orbit, μ = G·M_total)."""
    h = rx * vy - ry * vx                  # specific angular momentum (z)
    r = math.hypot(rx, ry)
    # Eccentricity vector e = (v × h)/μ - r̂ (in 2D, h is a scalar in z)
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
        integrator="ias15",
        dt=1e-3,
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

    arcsec_per_rad = (180.0 / math.pi) * 3600.0
    rel_err = (measured - predicted) / predicted

    print(f"Mercury 1PN — {N_ORBITS} orbits under IAS15 + apsis.gr.PostNewtonian1PN")
    print(f"  predicted Δω = {predicted * arcsec_per_rad:+.4f} arcsec  ({predicted:+.6e} rad)")
    print(f"  measured  Δω = {measured  * arcsec_per_rad:+.4f} arcsec  ({measured:+.6e} rad)")
    print(f"  relative error = {rel_err:+.3e}")
    de_rel = sys.energy_delta
    if de_rel is None:
        print(f"  |dE|   = {sys.abs_energy_drift:.3e}  (|E0| below conditioning floor)")
    else:
        print(f"  |dE/E| = {abs(de_rel):.3e}")

    # Demo-level GR recovery at IAS15 + dt = 1e-3.
    assert abs(rel_err) < 1e-2, f"Mercury 1PN precession off by {rel_err:.3e}"


if __name__ == "__main__":
    main()
