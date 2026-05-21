"""Central-potential precession — Tamayo 2019, Patterns A and B.

A test particle around a central body, with an additional radial
force ``∝ r^γ`` that drives apsidal precession. Demonstrates both
construction patterns of ``apsis.central.CentralForce``:

- **Pattern A (`from_raw`)**: pass the coupling ``a_central`` and
  exponent ``γ`` directly.
- **Pattern B (`from_apsidal_rate`)**: pass the desired apsidal rate
  ``ω̇`` and the operator inverts it to compute ``a_central``.

Run::

    python python-examples/central_precession.py
"""

from __future__ import annotations

import math

import apsis
from apsis.central import CentralForce

A_PLANET = 0.387
# Near-circular regime — `CentralForce::from_apsidal_rate` inverts an
# apsidal rate using the near-circular approximation; higher eccentricity
# picks up an O(e²) correction the constructor does not apply.
E_PLANET = 0.05
M_PLANET = 1e-7
M_SUN = 1.0
MU = M_SUN + M_PLANET

GAMMA = -3.0
A_CENTRAL = 5e-9
N_ORBITS = 50


def make_system() -> tuple[apsis.System, list[apsis.Body]]:
    sun = apsis.Body.star(mass=M_SUN)
    r_peri = A_PLANET * (1.0 - E_PLANET)
    v_peri = math.sqrt(MU * (2.0 / r_peri - 1.0 / A_PLANET))
    planet = apsis.Body.rocky(mass=M_PLANET).at((r_peri, 0.0)).with_velocity((0.0, v_peri))
    bodies = [sun, planet]
    sys = apsis.System(
        bodies=bodies,
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    return sys, bodies


def perihelion_longitude(rx: float, ry: float, vx: float, vy: float, mu: float) -> float:
    h = rx * vy - ry * vx
    r = math.hypot(rx, ry)
    ex = (vy * h) / mu - rx / r
    ey = (-vx * h) / mu - ry / r
    return math.atan2(ey, ex)


def measure_precession(sys: apsis.System, n_orbits: int) -> tuple[float, float]:
    period = 2.0 * math.pi * math.sqrt(A_PLANET**3 / MU)

    s, p = sys.bodies
    omega_0 = perihelion_longitude(p.x - s.x, p.y - s.y, p.vx - s.vx, p.vy - s.vy, MU)
    sys.integrate_for(period * n_orbits)
    s, p = sys.bodies
    omega_1 = perihelion_longitude(p.x - s.x, p.y - s.y, p.vx - s.vx, p.vy - s.vy, MU)

    d = omega_1 - omega_0
    while d > math.pi:
        d -= 2.0 * math.pi
    while d <= -math.pi:
        d += 2.0 * math.pi
    return d, period


def pattern_a() -> None:
    sys, _ = make_system()
    sys.add_hamiltonian_perturbation(
        CentralForce.from_raw(
            source=0,
            a_central=A_CENTRAL,
            gamma=GAMMA,
            units=apsis.units.SOLAR_CANONICAL,
        ),
    )
    measured, _ = measure_precession(sys, N_ORBITS)
    arcsec = measured * (180.0 / math.pi) * 3600.0
    print(f"  Pattern A — from_raw(a_central={A_CENTRAL:.1e}, γ={GAMMA})")
    print(f"    Δω over {N_ORBITS} orbits = {arcsec:+.4f} arcsec")


def pattern_b() -> None:
    sys, bodies = make_system()
    target_omega_dot = 5e-8  # rad / canonical time unit
    sys.add_hamiltonian_perturbation(
        CentralForce.from_apsidal_rate(
            source=0,
            target=1,
            omega_dot=target_omega_dot,
            gamma=GAMMA,
            bodies=bodies,
            units=apsis.units.SOLAR_CANONICAL,
        ),
    )
    measured, period = measure_precession(sys, N_ORBITS)
    expected = target_omega_dot * period * N_ORBITS
    rel_err = (measured - expected) / expected
    arcsec_meas = measured * (180.0 / math.pi) * 3600.0
    arcsec_exp = expected * (180.0 / math.pi) * 3600.0
    print(f"  Pattern B — from_apsidal_rate(ω̇={target_omega_dot:.1e}, γ={GAMMA})")
    print(f"    expected Δω = {arcsec_exp:+.4f} arcsec over {N_ORBITS} orbits")
    print(f"    measured Δω = {arcsec_meas:+.4f} arcsec  (rel err {rel_err:+.3e})")
    # Pattern B inversion is near-circular; ~4 % drift at e = 0.05, γ = -3
    # over 50 orbits sums O(e²), rate variation along the orbit, and
    # perturbation back-reaction.
    assert abs(rel_err) < 0.05, f"Pattern B round-trip off by {rel_err:.3e}"


def main() -> None:
    print("apsis.central.CentralForce — apsidal precession demo")
    pattern_a()
    pattern_b()


if __name__ == "__main__":
    main()
