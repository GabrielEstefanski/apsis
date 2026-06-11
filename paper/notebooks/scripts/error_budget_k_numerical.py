"""
error_budget_k_numerical.py
===========================

Numerically measures the second-order coefficient k in the 1PN apsidal advance:

    Delta_omega(eps) = 6*pi*eps * (1 + k*eps + O(eps^2))

where eps = mu / (c^2 * a_orb * (1 - e^2)).

Force from crates/apsis-1pn/src/lib.rs::accumulate_force (test-particle limit):

    a = -mu/r^2 * n_hat + (mu/(c^2*r^2)) * [(4*mu/r - v^2)*n_hat + 4*(n_hat.v)*v]

n_hat = (x, y)/r (source at origin to particle).

Integrator: mpmath.odefun (Taylor-series IVP solver, rigorous at working
precision), degree=20. Timing benchmark (this machine, dps=40, tol=1e-34):
degree=8 -> 28.8 s/eval; degree=14 -> 0.83 s/eval; degree=20 -> 0.30 s/eval.
The large Taylor order permits ~0.02-length internal steps, cutting internal
step count ~100x vs degree=8.

Periapsis detection: coarse ordered scan (builds odefun cache) + bisection
on rdot = x*vx + y*vy (cache lookups, fast). Angle: atan2 unwrapped along
the trajectory.
"""

import sys
import time
import mpmath
from mpmath import mp, mpf, sqrt, pi, atan2, log10

# ─────────────────────────── Precision ────────────────────────────────────────
mp.dps = 40

MU    = mpf('1')
A_ORB = mpf('1')

TAYLOR_DEGREE = 20


# ─────────────────────────── Physics ──────────────────────────────────────────

def c_from_eps(eps, a_orb, e):
    return sqrt(MU / (eps * a_orb * (1 - e**2)))

def periapsis_ic(a_orb, e):
    r0 = a_orb * (1 - e)
    v0 = sqrt(MU * (1 + e) / (a_orb * (1 - e)))
    return r0, mpf('0'), mpf('0'), v0   # x, y, vx, vy

def make_rhs(mu, c):
    """EXACT equations from accumulate_force (test-particle limit)."""
    c2 = c * c
    def rhs(t, state):
        x, y, vx, vy = state
        r2   = x*x + y*y
        r    = sqrt(r2)
        rinv = 1 / r
        nx, ny = x*rinv, y*rinv
        v2   = vx*vx + vy*vy
        ndv  = nx*vx + ny*vy
        pref = mu / (c2 * r2)
        sc_n = 4 * mu * rinv - v2
        sc_v = 4 * ndv
        ax   = -mu / r2 * nx + pref * (sc_n * nx + sc_v * vx)
        ay   = -mu / r2 * ny + pref * (sc_n * ny + sc_v * vy)
        return [vx, vy, ax, ay]
    return rhs

def rdot_unnorm(state):
    """x*vx + y*vy = r * rdot (same sign as rdot)."""
    return state[0]*state[2] + state[1]*state[3]


# ─────────────────────────── Core measurement ─────────────────────────────────

def measure_periapsis_advance(mu, c, a_orb, e, dps_override=None,
                               n_coarse=25, n_bisect=160, n_angle=250):
    """
    Integrate one radial period; return Delta_omega = swept_angle - 2*pi.
    """
    dps = dps_override if dps_override is not None else mp.dps
    saved = mp.dps
    mp.dps = dps
    try:
        x0, y0, vx0, vy0 = periapsis_ic(a_orb, e)
        rhs   = make_rhs(mu, c)
        T_kep = 2 * pi * sqrt(a_orb**3 / mu)
        T_end = T_kep * mpf('1.1')
        tol   = mpf(10) ** (-(dps - 6))

        t0 = time.time()
        sol = mpmath.odefun(rhs, 0, [x0, y0, vx0, vy0],
                            tol=tol, degree=TAYLOR_DEGREE)

        # ── Coarse scan (ordered evaluations build the cache) ─────────────────
        dt_c = T_end / mpf(n_coarse)
        past_apoapsis = False
        t_lo_b = None
        t_hi_b = None

        t_prev  = mpf('0')
        rd_prev = rdot_unnorm(sol(t_prev))   # ~0 at periapsis IC

        for k in range(1, n_coarse + 1):
            t_cur  = mpf(k) * dt_c
            rd_cur = rdot_unnorm(sol(t_cur))

            if not past_apoapsis:
                if rd_cur < 0:
                    past_apoapsis = True
            else:
                if rd_prev < 0 and rd_cur >= 0:
                    t_lo_b = t_prev
                    t_hi_b = t_cur
                    break

            t_prev  = t_cur
            rd_prev = rd_cur

        if t_lo_b is None:
            raise RuntimeError(
                f"Failed to bracket periapsis in [0, {float(T_end):.3f}] "
                f"(n_coarse={n_coarse}). past_apoapsis={past_apoapsis}"
            )

        t1 = time.time()
        print(f"scan {t1-t0:.1f}s [{float(t_lo_b):.4f},{float(t_hi_b):.4f}]",
              end='  ', flush=True)

        # ── Bisection (cached range → fast polynomial evaluation) ─────────────
        for _ in range(n_bisect):
            t_mid  = (t_lo_b + t_hi_b) / 2
            rd_mid = rdot_unnorm(sol(t_mid))
            if rd_mid < 0:
                t_lo_b = t_mid
            else:
                t_hi_b = t_mid
            if t_hi_b - t_lo_b < mpf(10)**(-36):
                break

        t_peri2 = (t_lo_b + t_hi_b) / 2
        t2 = time.time()
        print(f"bisect {t2-t1:.1f}s t_p2={mpmath.nstr(t_peri2, 10)}",
              end='  ', flush=True)

        # ── Angle measurement: unwrap along trajectory ────────────────────────
        dt_a      = t_peri2 / mpf(n_angle)
        phi_prev  = atan2(y0, x0)   # = 0 (periapsis on +x axis)
        phi_total = mpf('0')

        for k in range(1, n_angle + 1):
            t = mpf(k) * dt_a
            s = sol(t)
            phi_cur = atan2(s[1], s[0])
            dphi    = phi_cur - phi_prev
            while dphi >   pi: dphi -= 2*pi
            while dphi <= -pi: dphi += 2*pi
            phi_total += dphi
            phi_prev  = phi_cur

        Delta_omega = phi_total - 2 * pi
        t3 = time.time()
        print(f"angle {t3-t2:.1f}s dw={mpmath.nstr(Delta_omega, 10)}", flush=True)

        return Delta_omega

    finally:
        mp.dps = saved


# ─────────────────────────── Fitting ──────────────────────────────────────────

def fit_k(eps_list, y_list):
    """Least-squares fit y = k*eps + m*eps^2. Returns (k, m, jackknife spread)."""
    n = len(eps_list)
    A = mpmath.matrix(n, 2)
    b = mpmath.matrix(n, 1)
    for i, (eps, y) in enumerate(zip(eps_list, y_list)):
        A[i, 0] = eps
        A[i, 1] = eps**2
        b[i, 0] = y
    x = mpmath.lu_solve(A.T * A, A.T * b)
    k_fit, m_fit = x[0], x[1]

    if n > 2:
        k_jk = []
        for skip in range(n):
            ei = [eps_list[i] for i in range(n) if i != skip]
            yi = [y_list[i]   for i in range(n) if i != skip]
            nj = n - 1
            Aj = mpmath.matrix(nj, 2)
            bj = mpmath.matrix(nj, 1)
            for i2, (e2, y2) in enumerate(zip(ei, yi)):
                Aj[i2, 0] = e2
                Aj[i2, 1] = e2**2
                bj[i2, 0] = y2
            xj = mpmath.lu_solve(Aj.T * Aj, Aj.T * bj)
            k_jk.append(xj[0])
        k_spread = max(abs(kv - k_fit) for kv in k_jk)
    else:
        k_spread = mpf('0')

    return k_fit, m_fit, k_spread


# ─────────────────────────── Per-eccentricity run ─────────────────────────────

def run_ladder(e, eps_list, label):
    print(f"\n{'='*72}")
    print(f"e = {float(e):.4f}  ({label})")
    print(f"{'='*72}")

    dws = []
    ys  = []

    for eps in eps_list:
        c = c_from_eps(eps, A_ORB, e)
        print(f"  eps={float(eps):.2e}  c={float(c):.6e}  ", end='', flush=True)
        dw = measure_periapsis_advance(MU, c, A_ORB, e)
        dws.append(dw)
        y = dw / (6 * pi * eps) - 1
        ys.append(y)
        print(f"    Delta_omega = {mpmath.nstr(dw, 35)}")
        print(f"    y           = {mpmath.nstr(y, 20)}")

    # ── Gate G1: first-order reproduced at smallest eps ───────────────────────
    eps_min = min(eps_list)
    i_min   = eps_list.index(eps_min)
    g1_val  = abs(ys[i_min])
    g1_pass = g1_val < mpf('1e-3')
    print(f"\n  [G1] |y(eps={float(eps_min):.1e})| = {mpmath.nstr(g1_val, 10)} < 1e-3: "
          f"{'PASS' if g1_pass else 'FAIL'}")
    if not g1_pass:
        print("  BLOCKED: Gate G1 failed — first-order 1PN advance not reproduced")
        sys.exit(1)

    # ── Gate G2: dps=40 vs dps=55 agreement >= 25 significant digits ─────────
    c_min = c_from_eps(eps_min, A_ORB, e)
    print(f"  [G2] Re-running eps_min at dps=55 ... ", end='', flush=True)
    dw_55 = measure_periapsis_advance(MU, c_min, A_ORB, e, dps_override=55)
    dw_40 = dws[i_min]
    rel = abs(dw_55 - dw_40) / abs(dw_40) if dw_40 != 0 else mpf('0')
    sig = int(-log10(rel)) if rel > 0 else 55
    g2_pass = sig >= 25
    print(f"  [G2] rel={mpmath.nstr(rel, 6)}  sig_digits={sig}: {'PASS' if g2_pass else 'FAIL'}")
    if not g2_pass:
        print(f"  BLOCKED: Gate G2 failed — {sig} digits agreement (need >=25)")
        sys.exit(1)

    # ── Gate G3: near-Newtonian limit (c=1e30) ────────────────────────────────
    print(f"  [G3] Near-Newtonian c=1e30 ... ", end='', flush=True)
    dw_newt = measure_periapsis_advance(MU, mpf('1e30'), A_ORB, e)
    g3_val  = abs(dw_newt)
    g3_pass = g3_val < mpf('1e-30')
    print(f"  [G3] |dw_Newt|={mpmath.nstr(g3_val, 6)} < 1e-30: {'PASS' if g3_pass else 'FAIL'}")
    if not g3_pass:
        print(f"  BLOCKED: Gate G3 failed — {float(g3_val):.3e} >= 1e-30")
        sys.exit(1)

    # ── Fit ───────────────────────────────────────────────────────────────────
    k_fit, m_fit, k_unc = fit_k(eps_list, ys)
    print(f"\n  k_fit           = {mpmath.nstr(k_fit, 20)}")
    print(f"  m_fit           = {mpmath.nstr(m_fit, 20)}")
    print(f"  k_jackknife_unc = {mpmath.nstr(k_unc, 10)}")

    return dws, ys, k_fit, m_fit, k_unc


# ─────────────────────────── Main ─────────────────────────────────────────────

def main():
    t_start = time.time()
    print("error_budget_k_numerical.py")
    print(f"mp.dps={mp.dps}  mu={MU}  a_orb={A_ORB}  taylor_degree={TAYLOR_DEGREE}")

    e1   = mpf('0.2')
    eps1 = [mpf('1e-4'), mpf('5e-5'), mpf('2e-5'), mpf('1e-5'),
            mpf('5e-6'), mpf('2e-6'), mpf('1e-6')]
    dws1, ys1, k1, m1, ku1 = run_ladder(e1, eps1, "primary  e=0.2")

    e2   = mpf('0.4')
    eps2 = [mpf('1e-4'), mpf('1e-5'), mpf('1e-6')]
    dws2, ys2, k2, m2, ku2 = run_ladder(e2, eps2, "secondary  e=0.4")

    # k = k0 + k1_coef * e^2 (two-eccentricity solve)
    de2     = e1**2 - e2**2
    k1_coef = (k1 - k2) / de2
    k0_coef = k1 - k1_coef * e1**2

    t_total = time.time() - t_start
    print(f"\nTotal wall-time: {t_total/60:.1f} min")

    print()
    print("=" * 72)
    print("MACHINE-READABLE SUMMARY")
    print("=" * 72)
    print(f"mp.dps={mp.dps}  mu={MU}  a_orb={A_ORB}  taylor_degree={TAYLOR_DEGREE}")
    print()
    print("# e = 0.2 ladder")
    for eps, dw, y in zip(eps1, dws1, ys1):
        print(f"  eps={mpmath.nstr(eps,6)}  "
              f"Delta_omega={mpmath.nstr(dw,30)}  "
              f"y={mpmath.nstr(y,20)}")
    print()
    print("# e = 0.4 ladder")
    for eps, dw, y in zip(eps2, dws2, ys2):
        print(f"  eps={mpmath.nstr(eps,6)}  "
              f"Delta_omega={mpmath.nstr(dw,30)}  "
              f"y={mpmath.nstr(y,20)}")
    print()
    print(f"k(e=0.2)           = {mpmath.nstr(k1, 20)}")
    print(f"k(e=0.2) jk_unc    = {mpmath.nstr(ku1, 10)}")
    print(f"k(e=0.4)           = {mpmath.nstr(k2, 20)}")
    print(f"k(e=0.4) jk_unc    = {mpmath.nstr(ku2, 10)}")
    print()
    print("k = k0 + k1_coef * e^2  (two-eccentricity solve):")
    print(f"  k0      = {mpmath.nstr(k0_coef, 20)}")
    print(f"  k1_coef = {mpmath.nstr(k1_coef, 20)}")
    print()
    print("Gates (all must PASS):")
    print("  G1 PASS  |y(eps_min)| < 1e-3 at both eccentricities")
    print("  G2 PASS  dps=40 vs dps=55 agreement >= 25 significant digits")
    print("  G3 PASS  near-Newtonian (c=1e30) advance < 1e-30 rad")
    print("=" * 72)


if __name__ == "__main__":
    main()
