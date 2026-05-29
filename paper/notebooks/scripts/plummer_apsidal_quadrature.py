"""Exact apsidal-angle quadrature for the Plummer-softened two-body orbit.

Method-independent oracle for the §3.2 counter-test. Computes the true
precession per radial period

    Δϖ = 2 ∫_{r_min}^{r_max} (L/r²) / sqrt(2(E − Φ(r)) − L²/r²) dr − 2π

for the FULL Plummer potential Φ(r) = −GM/sqrt(r²+ε²) (no ε-expansion),
via the central-potential apsidal-angle integral (Goldstein §3.6), plus
the radial period T_r = 2 ∫ dr / sqrt(...) for the Kepler-vs-radial
"orbit" bookkeeping.

Independent of the LRL-vector ω definition, of Kepler-vs-radial sampling,
and of any time-integrator — so it decomposes the ~3% gap between the
leading-order closed form and the time-integrated rate (scipy DOP853 /
apsis IAS15) into (1) closed-form next-order truncation and (2) the
measurement-definition difference. numpy-only (Gauss–Legendre +
bisection); no scipy.

Output is ASCII-only so it does not depend on the console encoding
(Windows cp1252, Linux UTF-8, redirected pipes/CI logs all work);
non-ASCII lives only in this docstring and comments, which Python reads
as UTF-8 source on every platform.

Sanity check: at eps=0 the orbit is pure Kepler and the precession is 0.
"""

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]

GM = 1.0
A = 0.387098
E_ECC = 0.20563

RAD_TO_ARCSEC = 180.0 / np.pi * 3600.0
T_MERCURY_DAYS = 87.969
DAYS_PER_CENTURY = 36_525.0
ORBITS_PER_CENTURY = DAYS_PER_CENTURY / T_MERCURY_DAYS

# Gauss–Legendre nodes/weights on [-1, 1], mapped to [0, pi] per call.
_GL = np.polynomial.legendre.leggauss(400)
_GL_NODES: FloatArray = _GL[0]
_GL_WEIGHTS: FloatArray = _GL[1]


def phi(r: float, eps: float) -> float:
    """Specific Plummer potential (GM-scaled)."""
    return -GM / np.sqrt(r * r + eps * eps)


def orbit_constants(eps: float) -> tuple[float, float, float]:
    """E, L, r_peri from the §3.2 IC: at periapsis r0 = a(1-e), tangential
    vis-viva speed using the Kepler mu (exactly as the harness sets it)."""
    r0 = A * (1.0 - E_ECC)
    v0 = np.sqrt(GM * (2.0 / r0 - 1.0 / A))
    energy = 0.5 * v0 * v0 + phi(r0, eps)
    angmom = r0 * v0
    return energy, angmom, r0


def _radial_kinetic(r: float, energy: float, angmom: float, eps: float) -> float:
    """Half rdot^2 = E - Phi(r) - L^2/2r^2. Zero at the turning points."""
    return energy - phi(r, eps) - angmom * angmom / (2.0 * r * r)


def turning_points(eps: float) -> tuple[float, float, float, float]:
    energy, angmom, r0 = orbit_constants(eps)
    r_min = r0  # IC has rdot = 0 by construction

    lo = A  # interior point, half-rdot^2 > 0
    hi = A * (1.0 + E_ECC)
    while _radial_kinetic(hi, energy, angmom, eps) > 0.0:
        hi *= 1.5
    for _ in range(200):  # bisection: kinetic(lo) > 0, kinetic(hi) < 0
        mid = 0.5 * (lo + hi)
        if _radial_kinetic(mid, energy, angmom, eps) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-15 * mid:
            break
    r_max = 0.5 * (lo + hi)
    return energy, angmom, r_min, r_max


def _grid(eps: float) -> tuple[float, FloatArray, FloatArray, FloatArray]:
    """Build the Gauss–Legendre grid on [r_min, r_max] via r = c - d*cos(u),
    u in [0, pi], which cancels the inverse-sqrt turning-point singularity
    (integrand then smooth -> GL is spectral). Returns (L, r, half-rdot^2,
    quadrature weight including the dr/du Jacobian)."""
    energy, angmom, r_min, r_max = turning_points(eps)
    c = 0.5 * (r_min + r_max)
    d = 0.5 * (r_max - r_min)
    u: FloatArray = 0.5 * np.pi * (_GL_NODES + 1.0)
    r: FloatArray = c - d * np.cos(u)
    pot: FloatArray = -GM / np.sqrt(r * r + eps * eps)
    rk2: FloatArray = np.clip(2.0 * (energy - pot - angmom * angmom / (2.0 * r * r)), 1e-300, None)
    weight: FloatArray = 0.5 * np.pi * _GL_WEIGHTS * (d * np.sin(u))
    return angmom, r, rk2, weight


def apsidal_sweep(eps: float) -> float:
    """Angle swept periapse -> apoapse."""
    angmom, r, rk2, weight = _grid(eps)
    return float(np.dot(weight, (angmom / (r * r)) / np.sqrt(rk2)))


def radial_period(eps: float) -> float:
    """Full radial period (periapse -> apoapse -> periapse)."""
    _, _, rk2, weight = _grid(eps)
    return float(2.0 * np.dot(weight, 1.0 / np.sqrt(rk2)))


def precession_per_orbit_rad(eps: float) -> float:
    return 2.0 * apsidal_sweep(eps) - 2.0 * np.pi


def closed_form_per_orbit_rad(eps: float) -> float:
    """Leading-order secular apsidal precession per radial period."""
    return -3.0 * np.pi * eps**2 / (A**2 * (1.0 - E_ECC**2) ** 2)


def next_order_a_rad(eps: float) -> float:
    """First-order contribution of the eps^4 potential term V4 = -3 mu eps^4 / 8r^5,
    via the Landau-Lifshitz precession formula dvarpi = d/dL (closed-loop integral
    of dU dt). One of the two O(eps^4) pieces; the other is the second-order
    contribution of V2, left to the exact quadrature here."""
    return 15.0 * np.pi * (4.0 + 3.0 * E_ECC**2) * eps**4 / (8.0 * A**4 * (1.0 - E_ECC**2) ** 4)


def kepler_period() -> float:
    return 2.0 * np.pi * A**1.5  # GM = 1


def to_arcsec_per_century(per_orbit_rad: float) -> float:
    return per_orbit_rad * RAD_TO_ARCSEC * ORBITS_PER_CENTURY


if __name__ == "__main__":
    print("=== Self-check: eps=0 must give precession = 0 (Kepler closes) ===")
    print(f"  precession(eps=0) = {precession_per_orbit_rad(0.0):.3e} rad/orbit  (should be ~0)")
    print()

    eps = 0.02
    t_rad = radial_period(eps)
    t_kep = kepler_period()
    opc = ORBITS_PER_CENTURY
    # Convention-free observable: the precession RATE (rad per unit time).
    # Each source is divided by ITS OWN period, so the per-radial-vs-per-
    # Kepler "which orbit" ambiguity cancels. The closed form's natural
    # output is a rate (secular dvarpi/dt x 2pi/n); the exact quadrature
    # gives the geometric precession per radial period; the integrators
    # report osculating-omega drift per Kepler period.
    rate_exact = precession_per_orbit_rad(eps) / t_rad
    rate_closed = closed_form_per_orbit_rad(eps) / t_kep
    rate_apsis = (-2.289e6 / (RAD_TO_ARCSEC * opc)) / t_kep
    rate_scipy = (-2.275e6 / (RAD_TO_ARCSEC * opc)) / t_kep
    print(f"=== eps = {eps}: apsidal precession RATE (rad/time) -- convention-free ===")
    print(f"  exact (full Plummer quadrature) = {rate_exact:.6e}   <- ground truth")
    print(f"  closed form (leading order)     = {rate_closed:.6e}  "
          f"({(rate_closed / rate_exact - 1.0) * 100:+.2f}%  next-order truncation)")
    print(f"  apsis IAS15 measurement         = {rate_apsis:.6e}  "
          f"({(rate_apsis / rate_exact - 1.0) * 100:+.2f}%)")
    print(f"  scipy DOP853 measurement        = {rate_scipy:.6e}  "
          f"({(rate_scipy / rate_exact - 1.0) * 100:+.2f}%)")
    print(f"  T_radial = {t_rad:.6f}   T_kepler = {t_kep:.6f}   "
          f"(differ by {(t_rad / t_kep - 1.0) * 100:+.2f}%)")
    # Gate reference: the integrators (and the Rust gate) report the
    # osculating-omega drift PER KEPLER PERIOD; the exact apsidal-precession
    # rate times T_kepler is the value to pin in
    # softened_plummer_precession_validation.rs.
    drift_per_kepler = rate_exact * t_kep
    print(f"  -> GATE REFERENCE (exact drift per Kepler period) = {drift_per_kepler:.8e} rad")
    print(f"     = {drift_per_kepler * RAD_TO_ARCSEC * opc:.6e} arcsec/century")
    print()

    print("=== eps-sweep: closed-form RATE error vs exact (true next-order eps^2) ===")
    print(f"  {'eps':>7} {'rate_closed/rate_exact':>24} {'(ratio-1)/eps^2':>16}")
    for e in (0.002, 0.005, 0.01, 0.02, 0.04, 0.08):
        rr = (closed_form_per_orbit_rad(e) / t_kep) / (precession_per_orbit_rad(e) / radial_period(e))
        print(f"  {e:>7.3f} {rr:>24.6f} {(rr - 1.0) / e**2:>16.3f}")
    print()

    print("=== Purist: analytic next-order piece (a) vs exact, per radial period ===")
    print("  residual after subtracting leading + piece(a) should fall as eps^4 -> the")
    print("  remainder is the second-order V2 piece (b). All per radial period.")
    print(f"  {'eps':>7} {'exact':>14} {'LO only gap':>13} {'LO+(a) gap':>13} {'(b) rem/eps^4':>14}")
    for e in (0.002, 0.005, 0.01, 0.02, 0.04):
        ex = precession_per_orbit_rad(e)
        lo = closed_form_per_orbit_rad(e)
        lo_a = lo + next_order_a_rad(e)
        gap_lo = (ex - lo) / lo  # relative
        gap_lo_a = (ex - lo_a) / lo
        rem = (ex - lo_a) / e**4
        print(f"  {e:>7.3f} {ex:>14.6e} {gap_lo * 100:>12.3f}% {gap_lo_a * 100:>12.3f}% {rem:>14.3f}")
