"""
error_budget_endpoint_numerical.py
==================================

Extended-precision discriminator for the endpoint structure of the
osculating omega derived in
`error_budget_endpoint_symbolic.py`, in the gate's measurement convention
(fixed-time endpoint, Newtonian e-vector, mu = 1).

Two questions, answered with no f64 integrator in the loop:

  1. **Exact endpoint** (t = N * P0, P0 = osculating Kepler period of the
     IC): does the truncated-1PN physics produce any O(eps) angle term?
     The f64 ensembles show a constant -0.56 eps; if that were physics it
     would survive here. Gate GN1: the fitted O(eps) coefficient alpha in
     residual = alpha*eps + beta*eps^2 must satisfy |alpha| < 0.05 —
     ten-fold below the f64 observation.

  2. **Overshoot grid** (t = N * P0 + dt): does the derived Q(nu) =
     eps [3 nu - (3/e - e) sin nu - (5/2) sin 2nu] predict the angle
     shift, parameter-free, using the *measured* endpoint true anomaly?
     Gate GN2: |Delta_omega_meas - (Q(nu2) - Q(nu1))| < 10 * 400 * eps^2
     at every grid point (second-order headroom; signal/bound ~ 20-80x).

Integrator and force: identical machinery to error_budget_k_numerical.py
(mpmath.odefun Taylor IVP, degree 20, dps 40; rhs mirrors
crates/apsis-1pn accumulate_force in the test-particle limit).

Run:  python paper/notebooks/scripts/error_budget_endpoint_numerical.py
"""

import sys
import time
from dataclasses import dataclass
from typing import Any

import mpmath
from mpmath import atan2, mp, mpf, pi, sin, sqrt

mp.dps = 40

MU = mpf("1")
A_ORB = mpf("1")
TAYLOR_DEGREE = 20


# ── Physics (mirrors error_budget_k_numerical.py) ─────────────────────────────

def c_from_eps(eps, a_orb, e):
    return sqrt(MU / (eps * a_orb * (1 - e**2)))


def periapsis_ic(a_orb, e):
    r0 = a_orb * (1 - e)
    v0 = sqrt(MU * (1 + e) / (a_orb * (1 - e)))
    return r0, mpf("0"), mpf("0"), v0


def make_rhs(mu, c):
    c2 = c * c

    def rhs(t, state):
        x, y, vx, vy = state
        r2 = x * x + y * y
        r = sqrt(r2)
        rinv = 1 / r
        nx, ny = x * rinv, y * rinv
        v2 = vx * vx + vy * vy
        ndv = nx * vx + ny * vy
        pref = mu / (c2 * r2)
        sc_n = 4 * mu * rinv - v2
        sc_v = 4 * ndv
        ax = -mu / r2 * nx + pref * (sc_n * nx + sc_v * vx)
        ay = -mu / r2 * ny + pref * (sc_n * ny + sc_v * vy)
        return [vx, vy, ax, ay]

    return rhs


# ── Osculating measurement (the gate's convention: Newtonian e-vector) ────────

def wrap_pi(x):
    while x > pi:
        x -= 2 * pi
    while x <= -pi:
        x += 2 * pi
    return x


def osculating_omega_nu(state):
    """omega = atan2(e_y, e_x); nu = position angle - omega (wrapped)."""
    x, y, vx, vy = state
    r = sqrt(x * x + y * y)
    h = x * vy - y * vx
    ex = (vy * h) / MU - x / r
    ey = -(vx * h) / MU - y / r
    omega = atan2(ey, ex)
    nu = wrap_pi(atan2(y, x) - omega)
    return omega, nu


def Q_of_nu(eps, e, nu):
    """Derived endpoint-offset function (error_budget_endpoint_symbolic.py)."""
    return eps * (3 * nu - (3 / e - e) * sin(nu) - mpf(5) / 2 * sin(2 * nu))


# ── One scenario: integrate N orbits, measure at fixed times ──────────────────

@dataclass
class GridPoint:
    """One displaced-endpoint measurement (all values mpmath mpf)."""

    frac: Any
    nu: Any
    meas: Any
    pred: Any
    diff: Any


@dataclass
class Scenario:
    """One (e, eps, N) integration with its endpoint measurements."""

    eps: Any
    n: int
    nu_end: Any
    residual_exact: Any
    q_end: Any
    grid: list[GridPoint]


def run_scenario(e, eps, n_orbits, dt_grid_frac) -> Scenario:
    """
    Integrate the 1PN EOM through t_end = N*P0 (+ max grid offset),
    measure (omega, nu) at t_end and at each grid offset.
    """
    c = c_from_eps(eps, A_ORB, e)
    x0, y0, vx0, vy0 = periapsis_ic(A_ORB, e)
    rhs = make_rhs(MU, c)
    P0 = 2 * pi * sqrt(A_ORB**3 / MU)
    t_end = P0 * n_orbits

    tol = mpf(10) ** (-(mp.dps - 6))
    t0 = time.time()
    sol = mpmath.odefun(rhs, 0, [x0, y0, vx0, vy0], tol=tol, degree=TAYLOR_DEGREE)

    # Ordered cache build out to the largest needed time.
    t_max = t_end + max(abs(f) for f in [*dt_grid_frac, mpf(0)]) * P0
    n_coarse = int(25 * n_orbits) + 2
    for k in range(1, n_coarse + 1):
        sol(t_max * mpf(k) / n_coarse)
    t1 = time.time()

    omega_end, nu_end = osculating_omega_nu(sol(t_end))
    predicted_secular = 6 * pi * eps * n_orbits
    residual_exact = omega_end - predicted_secular  # omega(0) = 0 by IC

    grid: list[GridPoint] = []
    for frac in dt_grid_frac:
        t_g = t_end + frac * P0
        om_g, nu_g = osculating_omega_nu(sol(t_g))
        meas_shift = om_g - omega_end
        pred_shift = Q_of_nu(eps, e, nu_g) - Q_of_nu(eps, e, nu_end)
        grid.append(
            GridPoint(frac=frac, nu=nu_g, meas=meas_shift, pred=pred_shift,
                      diff=meas_shift - pred_shift)
        )
    t2 = time.time()
    print(
        f"    [integration {t1 - t0:.1f}s, measurement {t2 - t1:.1f}s]  "
        f"nu_end={mpmath.nstr(nu_end, 8)}  residual={mpmath.nstr(residual_exact, 8)}",
        flush=True,
    )
    return Scenario(
        eps=eps, n=n_orbits, nu_end=nu_end,
        residual_exact=residual_exact, q_end=Q_of_nu(eps, e, nu_end), grid=grid,
    )


# ── Fit residual = alpha*eps + beta*eps^2 (A1's fit machinery) ────────────────

def fit_alpha_beta(eps_list, y_list):
    n = len(eps_list)
    A = mpmath.matrix(n, 2)
    b = mpmath.matrix(n, 1)
    for i, (ep, y) in enumerate(zip(eps_list, y_list)):
        A[i, 0] = ep
        A[i, 1] = ep**2
        b[i, 0] = y
    x = mpmath.lu_solve(A.T * A, A.T * b)
    return x[0], x[1]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("error_budget_endpoint_numerical.py")
    print(f"mp.dps={mp.dps}  mu={MU}  a_orb={A_ORB}  taylor_degree={TAYLOR_DEGREE}")

    e_mercury = mpf("0.20563")
    n_ladder = 3
    eps_ladder = [mpf("1e-5"), mpf("3e-6"), mpf("1e-6")]
    grid_frac = [mpf("-0.005"), mpf("-0.002"), mpf("0.001"), mpf("0.002"), mpf("0.005")]

    # ── GN1: exact-endpoint ladder at e = e_Mercury, N = 3 ───────────────────
    print(f"\nGN1 ladder: e={float(e_mercury)}  N={n_ladder}  exact endpoint")
    results = []
    for eps in eps_ladder:
        print(f"  eps={float(eps):.1e}", flush=True)
        results.append(run_scenario(e_mercury, eps, n_ladder, grid_frac))

    alpha, beta = fit_alpha_beta(
        [r.eps for r in results], [r.residual_exact for r in results]
    )
    gn1_pass = abs(alpha) < mpf("0.05")
    print("\n  [GN1] residual = alpha*eps + beta*eps^2:")
    print(f"        alpha = {mpmath.nstr(alpha, 8)}   (f64 ensembles: -0.56)")
    print(f"        beta  = {mpmath.nstr(beta, 8)}")
    print(f"  [GN1] |alpha| < 0.05: {'PASS' if gn1_pass else 'FAIL'}")
    if not gn1_pass:
        print("  BLOCKED: O(eps) endpoint term present in exact-endpoint physics —")
        print("  the f64 constant offset would NOT be an integrator artefact.")
        sys.exit(1)

    # ── GN2: overshoot grid — Q(nu) parameter-free, all rungs ────────────────
    print("\nGN2 overshoot grid: Delta_omega vs Q(nu_2) - Q(nu_1), per rung")
    gn2_pass = True
    for r in results:
        bound = 10 * 400 * r.eps ** 2
        print(f"  eps={float(r.eps):.1e}  bound={float(bound):.3e}")
        for g in r.grid:
            ok = abs(g.diff) < bound
            gn2_pass = gn2_pass and ok
            ratio = g.meas / g.pred if g.pred != 0 else mpf("nan")
            print(
                f"    dt/P0={float(g.frac):+.3f}  nu={mpmath.nstr(g.nu, 6)}  "
                f"meas={mpmath.nstr(g.meas, 8)}  pred={mpmath.nstr(g.pred, 8)}  "
                f"meas/pred={mpmath.nstr(ratio, 8)}  "
                f"diff={mpmath.nstr(g.diff, 4)}  {'PASS' if ok else 'FAIL'}"
            )
    print(f"  [GN2] all grid points within 4000*eps^2: {'PASS' if gn2_pass else 'FAIL'}")
    if not gn2_pass:
        print("  BLOCKED: derived Q(nu) does not reproduce the overshoot response.")
        sys.exit(1)

    # ── GN3: N-constancy probe at eps = 1e-5 (N = 10 vs N = 3) ───────────────
    print("\nGN3 probe: N=10 at eps=1e-5 (exact endpoint)")
    r10 = run_scenario(e_mercury, eps_ladder[0], 10, [])
    print(f"  residual(N=3)  = {mpmath.nstr(results[0].residual_exact, 8)}")
    print(f"  residual(N=10) = {mpmath.nstr(r10.residual_exact, 8)}")
    print("  (physics endpoint residual scales with N through nu_end and the")
    print("   eps^2 secular term — a constant-in-N O(eps) term would not.)")

    # ── GN4: eccentricity check of Q at e = 0.4, eps = 1e-5, N = 3 ───────────
    print("\nGN4: e=0.4 rung (validates the e-dependence of Q)")
    r_e4 = run_scenario(mpf("0.4"), eps_ladder[0], n_ladder, grid_frac)
    gn4_pass = True
    bound = 10 * 400 * eps_ladder[0] ** 2
    for g in r_e4.grid:
        ok = abs(g.diff) < bound
        gn4_pass = gn4_pass and ok
        ratio = g.meas / g.pred if g.pred != 0 else mpf("nan")
        print(
            f"    dt/P0={float(g.frac):+.3f}  meas/pred={mpmath.nstr(ratio, 8)}  "
            f"diff={mpmath.nstr(g.diff, 4)}  {'PASS' if ok else 'FAIL'}"
        )
    print(f"  [GN4] {'PASS' if gn4_pass else 'FAIL'}")
    if not gn4_pass:
        sys.exit(1)

    t_total = time.time() - t_start
    print(f"\nTotal wall-time: {t_total / 60:.1f} min")
    print()
    print("=" * 72)
    print("MACHINE-READABLE SUMMARY")
    print("=" * 72)
    print(f"GN1 PASS  alpha = {mpmath.nstr(alpha, 8)}  (|alpha| < 0.05; f64 shows -0.56)")
    print(f"GN1       beta  = {mpmath.nstr(beta, 8)}")
    print("GN2 PASS  Q(nu) reproduces the overshoot response, parameter-free")
    print(f"GN3       residual(N=3) = {mpmath.nstr(results[0].residual_exact, 6)}  "
          f"residual(N=10) = {mpmath.nstr(r10.residual_exact, 6)}")
    print("GN4 PASS  e-dependence of Q validated at e=0.4")
    print("=" * 72)


if __name__ == "__main__":
    main()
