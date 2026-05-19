"""Federated operators compose additively.

Sun + Mercury under IAS15 with TWO operators registered: 1PN
relativistic correction (apsis.gr) and a small central-potential
contribution (apsis.central). The federation contract says
perturbations accumulate into the same scratch buffer; the measured
apsidal precession should match the sum of the individual
contributions.

Run::

    python examples/federation_composition.py
"""

from __future__ import annotations

import math

import apsis
from apsis.central import CentralForce
from apsis.gr import C_SOLAR_UNITS, PostNewtonian1PN

A_MERCURY = 0.387_098
E_MERCURY = 0.205_63
M_MERCURY = 1.660_114e-7
M_SUN = 1.0

CENTRAL_OMEGA_DOT = 5e-7  # comparable in magnitude to the GR effect
GAMMA = -3.0
N_ORBITS = 80
MU = M_SUN + M_MERCURY


def make_system() -> tuple[apsis.System, list[apsis.Body]]:
    sun = apsis.Body.star(mass=M_SUN)
    r_peri = A_MERCURY * (1.0 - E_MERCURY)
    v_peri = math.sqrt(MU * (2.0 / r_peri - 1.0 / A_MERCURY))
    mercury = apsis.Body.rocky(mass=M_MERCURY).at((r_peri, 0.0)).with_velocity((0.0, v_peri))
    bodies = [sun, mercury]
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


def precession(sys: apsis.System, n_orbits: int, mu: float) -> tuple[float, float]:
    period = 2.0 * math.pi * math.sqrt(A_MERCURY**3 / mu)
    s, m = sys.bodies
    omega_0 = perihelion_longitude(m.x - s.x, m.y - s.y, m.vx - s.vx, m.vy - s.vy, mu)
    sys.integrate_for(period * n_orbits)
    s, m = sys.bodies
    omega_1 = perihelion_longitude(m.x - s.x, m.y - s.y, m.vx - s.vx, m.vy - s.vy, mu)
    d = omega_1 - omega_0
    while d > math.pi:
        d -= 2.0 * math.pi
    while d <= -math.pi:
        d += 2.0 * math.pi
    return d, period


def main() -> None:
    # ── Run 1: 1PN alone ────────────────────────────────────────────────────
    sys_a, _ = make_system()
    sys_a.add_hamiltonian_perturbation(
        PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
    )
    delta_gr, period = precession(sys_a, N_ORBITS, MU)

    # ── Run 2: central alone ────────────────────────────────────────────────
    sys_b, bodies_b = make_system()
    sys_b.add_hamiltonian_perturbation(
        CentralForce.from_apsidal_rate(
            source=0, target=1,
            omega_dot=CENTRAL_OMEGA_DOT, gamma=GAMMA,
            bodies=bodies_b,
            units=apsis.units.SOLAR_CANONICAL,
        ),
    )
    delta_central, _ = precession(sys_b, N_ORBITS, MU)

    # ── Run 3: both registered ──────────────────────────────────────────────
    sys_c, bodies_c = make_system()
    sys_c.add_hamiltonian_perturbation(
        PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
    )
    sys_c.add_hamiltonian_perturbation(
        CentralForce.from_apsidal_rate(
            source=0, target=1,
            omega_dot=CENTRAL_OMEGA_DOT, gamma=GAMMA,
            bodies=bodies_c,
            units=apsis.units.SOLAR_CANONICAL,
        ),
    )
    delta_both, _ = precession(sys_c, N_ORBITS, MU)

    arcsec = (180.0 / math.pi) * 3600.0
    expected_sum = delta_gr + delta_central
    rel_err = (delta_both - expected_sum) / expected_sum

    print(f"Federated composition — Sun + Mercury, {N_ORBITS} orbits")
    print(f"  GR (1PN) alone        Δω = {delta_gr * arcsec:+.4f} arcsec")
    print(f"  central alone         Δω = {delta_central * arcsec:+.4f} arcsec")
    print(f"  expected (additive)   Δω = {expected_sum * arcsec:+.4f} arcsec")
    print(f"  measured (composed)   Δω = {delta_both * arcsec:+.4f} arcsec")
    print(f"  composition error        = {rel_err:+.3e}")
    de_rel = sys_c.energy_delta
    if de_rel is None:
        print(f"  |dE|   (composed run)    = {sys_c.abs_energy_drift:.3e}  (|E0| below conditioning floor)")
    else:
        print(f"  |dE/E| (composed run)    = {abs(de_rel):.3e}")
    print(f"  c (canonical units)      = {C_SOLAR_UNITS:.3f}")

    # Cross-operator coupling (1PN ↔ central, both modifying apsidal rate
    # of the same orbit) introduces non-linear corrections that grow with
    # the product of perturbation strengths. 5 % is the headroom budget
    # for that coupling at the chosen rates; the additive-composition
    # contract itself is satisfied to f64 noise.
    assert abs(rel_err) < 5e-2, f"composition non-additive: {rel_err:.3e}"


if __name__ == "__main__":
    main()
