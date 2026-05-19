"""Radiation pressure on a dust grain — Burns 1979.

A dust grain on an initially circular Keplerian orbit at 1 AU, with
the Sun's radiation pressure reducing effective gravity by a factor
``(1 - β)`` where ``β = F_rad / F_grav``. Pure radiation pressure
(no Poynting–Robertson drag) is conservative: the grain settles
onto a new bound orbit determined by ``μ_eff = 1 - β`` rather than
``μ = 1``, with a slightly larger semi-major axis and modest
eccentricity inherited from the initial-state mismatch with the
new effective potential. The asserted invariant is that this new
semi-major axis is preserved over many orbits.

Run::

    python examples/radiation_dust.py
"""

from __future__ import annotations

import math

import apsis
from apsis.radiation import RadiationPressure

BETA = 0.1
N_ORBITS = 50


def kepler_period(a: float, mu: float) -> float:
    return 2.0 * math.pi * math.sqrt(a**3 / mu)


def main() -> None:
    sun = apsis.Body.star(mass=1.0)
    dust = apsis.Body.rocky(mass=1e-15).at((1.0, 0.0)).with_velocity((0.0, 1.0))

    sys = apsis.System(
        bodies=[sun, dust],
        units=apsis.units.SOLAR_CANONICAL,
        integrator="ias15",
        dt=1e-3,
    )
    sys.add_hamiltonian_perturbation(
        RadiationPressure.from_raw_betas(
            source=0,
            betas=[0.0, BETA],
            units=apsis.units.SOLAR_CANONICAL,
        ),
    )

    mu_eff = 1.0 - BETA  # central μ reduced by (1-β)
    period_eff = kepler_period(1.0, mu_eff)
    period_newton = kepler_period(1.0, 1.0)

    print(f"Dust grain at 1 AU, β = {BETA}")
    print(f"  Newtonian period      = {period_newton:.4f} (canonical units)")
    print(f"  effective period      = {period_eff:.4f}  (μ_eff = 1 - β = {mu_eff})")
    print(f"  expected period ratio = {period_eff / period_newton:.4f}")

    sys.integrate_for(period_newton * N_ORBITS)

    s, d = sys.bodies
    r_final = math.hypot(d.x - s.x, d.y - s.y)
    v_final = math.hypot(d.vx - s.vx, d.vy - s.vy)
    # Specific orbital energy ε = v²/2 - μ/r ; semi-major a = -μ/(2ε)
    # Initial state (r=1, v=1) was Newtonian-circular; with μ_eff = 1-β the
    # orbit is mildly eccentric with a_eff = -μ_eff/(2·(v²/2 - μ_eff/r)).
    # For β = 0.1: a_eff = 0.9 / (2·0.4) = 1.125.
    a_initial = -mu_eff / (2.0 * (0.5 * 1.0**2 - mu_eff / 1.0))
    eps = 0.5 * v_final**2 - mu_eff / r_final
    a_final = -mu_eff / (2.0 * eps)

    print(f"  after {N_ORBITS} Newtonian periods:")
    print(f"    r = {r_final:.4f}  v = {v_final:.4f}")
    print(f"    a (against μ_eff): initial = {a_initial:.4f}, final = {a_final:.4f}")
    de_rel = sys.energy_delta
    if de_rel is None:
        print(f"  |dE|   = {sys.abs_energy_drift:.3e}  (|E0| below conditioning floor)")
    else:
        print(f"  |dE/E| = {abs(de_rel):.3e}")

    # Smoke gate: catches unbinding and NaN integration.
    assert abs(a_final - a_initial) < 0.02, f"semi-major drift {a_final - a_initial:.3e}"


if __name__ == "__main__":
    main()
