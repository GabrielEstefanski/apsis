"""
error_budget_a2_mass_scaling.py
================================

Phase A2: measure the O(m/M) two-body correction to the 1PN apsidal advance.

Observable: geometric apsidal advance Delta_omega of the RELATIVE orbit
(r_rel = r_secondary - r_primary) over one radial period.

Force implemented EXACTLY as crates/apsis-1pn/src/lib.rs::accumulate_force:
  receiver i from source j:
    rhat_recv_to_src = (r_j - r_i) / |r_j - r_i|         (receiver -> source = -n_hat)
    pref = m_j / (c^2 * r^2)
    dot  = rhat_recv_to_src . v_i
    a_1PN_i += -pref * [(4*m_j/r - v_i^2)*rhat + 4*dot*v_i]
                ^--- minus sign because rhat is -n_hat (see Rust source comment)
    which equals: pref * [(4*m_j/r - v_i^2)*n_hat_j_to_i + 4*(n_hat.v_i)*v_i]

Both bodies receive their term, each sourced by the other's mass.

ICs: primary (M=1) at origin at rest; secondary (m=q) at (r_peri, 0) with
velocity (0, v_peri), where r_peri = a(1-e), v_peri = sqrt(2/r_peri - 1/a)
computed with mu=1 (NOT 1+q), exactly as the gate does.

a=1, e=0.2, G=1, mp.dps=40.
"""

import sys
import time

import mpmath
from mpmath import atan2, log10, mp, mpf, pi, sqrt

mp.dps = 40

A_ORB = mpf('1')
E_ORB = mpf('0.2')
TAYLOR_DEGREE = 20

M_MERCURY = mpf('1.660114e-7')   # for extrapolation only


# ────────────────────────── helpers ──────────────────────────────────────────

def c_from_eps(eps, a_orb, e):
    """c such that eps = 1 / (c^2 * a * (1-e^2))."""
    return sqrt(1 / (eps * a_orb * (1 - e**2)))


def periapsis_ic(a_orb, e):
    """
    ICs mirroring the gate exactly: mu=1.
    Returns (r_peri, v_peri).
    """
    r0 = a_orb * (1 - e)
    v0 = sqrt(2 / r0 - 1 / a_orb)   # sqrt(mu*(1+e)/(a*(1-e))) with mu=1
    return r0, v0


# ────────────────────────── RHS builders ─────────────────────────────────────

def make_rhs_single_clean(M, c):
    """
    Test-particle 1PN, source at origin, mass M.
    Matches accumulate_force exactly.
      rhat_recv_to_src = -r_hat  (where r_hat = pos/|pos|)
      a_newt = -GM/r^2 * r_hat  [toward source]
      a_1pn  = -pref*(sc_n*rhat_recv_to_src + sc_v*v)
             = pref*(sc_n*r_hat + sc_v*v) * (pref already has sign)
    """
    c2 = c * c
    def rhs(t, state):
        x, y, vx, vy = state
        r2   = x*x + y*y
        r    = sqrt(r2)
        rinv = 1/r
        # unit vector from receiver to source (= -n_hat in standard notation)
        rrs_x = -x*rinv
        rrs_y = -y*rinv
        v2   = vx*vx + vy*vy
        dot  = rrs_x*vx + rrs_y*vy   # rhat_recv_src . v
        pref = M / (c2 * r2)
        sc_n = 4*M*rinv - v2
        sc_v = 4*dot
        # Newtonian toward source = GM/r^2 * rhat_recv_src
        a_newt_x = M/r2 * rrs_x
        a_newt_y = M/r2 * rrs_y
        # 1PN (Rust: accumulate with minus because rhat is -n_hat)
        a_1pn_x = -pref * (sc_n*rrs_x + sc_v*vx)
        a_1pn_y = -pref * (sc_n*rrs_y + sc_v*vy)
        ax = a_newt_x + a_1pn_x
        ay = a_newt_y + a_1pn_y
        return [vx, vy, ax, ay]
    return rhs


def make_rhs_twobody(M, m, c):
    """
    Full two-body 1PN: both bodies receive Newtonian + pairwise 1PN.
    State: [x0,y0,vx0,vy0, x1,y1,vx1,vy1]
    Body 0: mass M (primary), body 1: mass m (secondary).
    Forces exactly from accumulate_force:
      for each receiver i, source j != i:
        rhat_recv_to_src = (r_j - r_i)/|r_j-r_i|
        pref = m_j / (c^2 * r^2)
        dot  = rhat.v_i
        sc_n = 4*m_j/r - v_i^2
        sc_v = 4*dot
        a_i_newt += m_j/r^2 * rhat_recv_to_src
        a_i_1pn  += -pref*(sc_n*rhat + sc_v*v_i)
    """
    c2 = c * c
    masses = [M, m]
    def rhs(t, state):
        pos = [(state[0], state[1]), (state[4], state[5])]
        vel = [(state[2], state[3]), (state[6], state[7])]
        acc = [[mpf(0), mpf(0)], [mpf(0), mpf(0)]]
        for i in range(2):
            xi, yi   = pos[i]
            vxi, vyi = vel[i]
            v2i = vxi*vxi + vyi*vyi
            for j in range(2):
                if i == j:
                    continue
                xj, yj = pos[j]
                dx = xj - xi   # recv->src
                dy = yj - yi
                r2 = dx*dx + dy*dy
                r  = sqrt(r2)
                rinv = 1/r
                rrs_x = dx*rinv   # rhat recv->src
                rrs_y = dy*rinv
                mj = masses[j]
                pref = mj / (c2 * r2)
                dot  = rrs_x*vxi + rrs_y*vyi
                sc_n = 4*mj*rinv - v2i
                sc_v = 4*dot
                # Newtonian
                acc[i][0] += mj/r2 * rrs_x
                acc[i][1] += mj/r2 * rrs_y
                # 1PN (minus sign from Rust rhat = -n_hat convention)
                acc[i][0] -= pref*(sc_n*rrs_x + sc_v*vxi)
                acc[i][1] -= pref*(sc_n*rrs_y + sc_v*vyi)
        return [
            vel[0][0], vel[0][1], acc[0][0], acc[0][1],
            vel[1][0], vel[1][1], acc[1][0], acc[1][1],
        ]
    return rhs


def make_rhs_twobody_newtonian(M, m):
    """Pure Newtonian two-body (no 1PN) for gate GC."""
    def rhs(t, state):
        pos = [(state[0], state[1]), (state[4], state[5])]
        vel = [(state[2], state[3]), (state[6], state[7])]
        acc = [[mpf(0), mpf(0)], [mpf(0), mpf(0)]]
        masses = [M, m]
        for i in range(2):
            xi, yi = pos[i]
            for j in range(2):
                if i == j:
                    continue
                xj, yj = pos[j]
                dx = xj - xi
                dy = yj - yi
                r2 = dx*dx + dy*dy
                r  = sqrt(r2)
                rinv = 1/r
                mj = masses[j]
                acc[i][0] += mj/r2 * dx*rinv
                acc[i][1] += mj/r2 * dy*rinv
        return [
            vel[0][0], vel[0][1], acc[0][0], acc[0][1],
            vel[1][0], vel[1][1], acc[1][0], acc[1][1],
        ]
    return rhs


# ────────────────────────── rdot helpers ─────────────────────────────────────

def rdot_rel(sol_val):
    """
    d|r_rel|/dt for the RELATIVE orbit.
    sol_val is the 8-component state [x0,y0,vx0,vy0,x1,y1,vx1,vy1].
    r_rel = (x1-x0, y1-y0), v_rel = (vx1-vx0, vy1-vy0).
    rdot = (r_rel . v_rel) / |r_rel|  — sign same as x*vx+y*vy in 1-body.
    """
    rx = sol_val[4] - sol_val[0]
    ry = sol_val[5] - sol_val[1]
    vrx = sol_val[6] - sol_val[2]
    vry = sol_val[7] - sol_val[3]
    return rx*vrx + ry*vry   # proportional to d|r|/dt


def rdot_single(sol_val):
    """For single-body solution [x,y,vx,vy]."""
    return sol_val[0]*sol_val[2] + sol_val[1]*sol_val[3]


# ────────────────────────── core measurement ─────────────────────────────────

def measure_advance_twobody(M, q, c, a_orb, e, label='',
                             n_coarse=30, n_bisect=160, n_angle=300):
    """
    Integrate one radial period of the RELATIVE orbit; return Delta_omega.
    Primary (M) at origin at rest; secondary (q) at (r_peri, 0) with
    IC velocity computed with mu=1 (matching the gate exactly).
    """
    r_peri, v_peri = periapsis_ic(a_orb, e)
    # Initial state: [x0,y0,vx0,vy0, x1,y1,vx1,vy1]
    ic = [mpf(0), mpf(0), mpf(0), mpf(0),
          r_peri, mpf(0), mpf(0), v_peri]

    rhs = make_rhs_twobody(M, q, c)

    # Keplerian period estimate with mu=1 (gate convention)
    T_kep = 2 * pi * sqrt(a_orb**3 / 1)
    T_end = T_kep * mpf('1.1')
    tol   = mpf(10)**(-(mp.dps - 6))

    t0 = time.time()
    sol = mpmath.odefun(rhs, 0, ic, tol=tol, degree=TAYLOR_DEGREE)

    # Coarse scan
    dt_c = T_end / mpf(n_coarse)
    past_apoapsis = False
    t_lo, t_hi = None, None
    t_prev = mpf(0)
    rd_prev = rdot_rel(ic)   # ~0 at periapsis

    for k in range(1, n_coarse + 1):
        t_cur = mpf(k) * dt_c
        rd_cur = rdot_rel(sol(t_cur))
        if not past_apoapsis:
            if rd_cur < 0:
                past_apoapsis = True
        else:
            if rd_prev < 0 and rd_cur >= 0:
                t_lo = t_prev
                t_hi = t_cur
                break
        t_prev = t_cur
        rd_prev = rd_cur

    if t_lo is None:
        raise RuntimeError(
            f"[{label}] Failed to bracket periapsis in [0, {float(T_end):.3f}]. "
            f"past_apoapsis={past_apoapsis}"
        )

    t1 = time.time()
    print(f"  scan {t1-t0:.1f}s [{float(t_lo):.4f},{float(t_hi):.4f}]",
          end='  ', flush=True)

    # Bisection on rdot_rel sign change (- to +)
    for _ in range(n_bisect):
        t_mid = (t_lo + t_hi) / 2
        rd_mid = rdot_rel(sol(t_mid))
        if rd_mid < 0:
            t_lo = t_mid
        else:
            t_hi = t_mid
        if t_hi - t_lo < mpf(10)**(-36):
            break

    t_peri2 = (t_lo + t_hi) / 2
    t2 = time.time()
    print(f"bisect {t2-t1:.1f}s t_p2={mpmath.nstr(t_peri2, 10)}",
          end='  ', flush=True)

    # Angle of RELATIVE orbit unwrapped along trajectory
    rx0 = ic[4] - ic[0]
    ry0 = ic[5] - ic[1]   # = 0
    phi_prev  = atan2(ry0, rx0)   # = 0
    phi_total = mpf(0)
    dt_a = t_peri2 / mpf(n_angle)

    for k in range(1, n_angle + 1):
        s = sol(mpf(k) * dt_a)
        rx = s[4] - s[0]
        ry = s[5] - s[1]
        phi_cur = atan2(ry, rx)
        dphi = phi_cur - phi_prev
        while dphi >  pi: dphi -= 2*pi
        while dphi <= -pi: dphi += 2*pi
        phi_total += dphi
        phi_prev = phi_cur

    Delta_omega = phi_total - 2*pi
    t3 = time.time()
    print(f"angle {t3-t2:.1f}s dw={mpmath.nstr(Delta_omega, 10)}", flush=True)
    return Delta_omega


def measure_advance_single(M, c, a_orb, e, label='',
                            n_coarse=30, n_bisect=160, n_angle=300):
    """
    Single-body test-particle advance (source fixed at origin).
    Same machinery, 4-component state.
    """
    r_peri, v_peri = periapsis_ic(a_orb, e)
    ic = [r_peri, mpf(0), mpf(0), v_peri]
    rhs = make_rhs_single_clean(M, c)

    T_kep = 2 * pi * sqrt(a_orb**3 / 1)
    T_end = T_kep * mpf('1.1')
    tol   = mpf(10)**(-(mp.dps - 6))

    t0 = time.time()
    sol = mpmath.odefun(rhs, 0, ic, tol=tol, degree=TAYLOR_DEGREE)

    dt_c = T_end / mpf(n_coarse)
    past_apoapsis = False
    t_lo, t_hi = None, None
    t_prev = mpf(0)
    rd_prev = rdot_single(ic)

    for k in range(1, n_coarse + 1):
        t_cur = mpf(k) * dt_c
        rd_cur = rdot_single(sol(t_cur))
        if not past_apoapsis:
            if rd_cur < 0:
                past_apoapsis = True
        else:
            if rd_prev < 0 and rd_cur >= 0:
                t_lo = t_prev
                t_hi = t_cur
                break
        t_prev = t_cur
        rd_prev = rd_cur

    if t_lo is None:
        raise RuntimeError(
            f"[{label}] Single-body: failed to bracket periapsis in [0, {float(T_end):.3f}]"
        )

    t1 = time.time()
    print(f"  [single] scan {t1-t0:.1f}s [{float(t_lo):.4f},{float(t_hi):.4f}]",
          end='  ', flush=True)

    for _ in range(n_bisect):
        t_mid = (t_lo + t_hi) / 2
        rd_mid = rdot_single(sol(t_mid))
        if rd_mid < 0:
            t_lo = t_mid
        else:
            t_hi = t_mid
        if t_hi - t_lo < mpf(10)**(-36):
            break

    t_peri2 = (t_lo + t_hi) / 2
    t2 = time.time()
    print(f"bisect {t2-t1:.1f}s t_p2={mpmath.nstr(t_peri2,10)}", end='  ', flush=True)

    phi_prev  = atan2(ic[1], ic[0])
    phi_total = mpf(0)
    dt_a = t_peri2 / mpf(n_angle)

    for k in range(1, n_angle + 1):
        s = sol(mpf(k) * dt_a)
        phi_cur = atan2(s[1], s[0])
        dphi = phi_cur - phi_prev
        while dphi >  pi: dphi -= 2*pi
        while dphi <= -pi: dphi += 2*pi
        phi_total += dphi
        phi_prev = phi_cur

    Delta_omega = phi_total - 2*pi
    t3 = time.time()
    print(f"angle {t3-t2:.1f}s dw={mpmath.nstr(Delta_omega,10)}", flush=True)
    return Delta_omega


# ────────────────────────── GC: Newtonian null ───────────────────────────────

def measure_advance_newtonian_twobody(M, q, a_orb, e, label='',
                                       n_coarse=30, n_bisect=160, n_angle=300):
    """Gate GC: pure Newtonian two-body (c=1e30 equivalent via no-1PN rhs)."""
    r_peri, v_peri = periapsis_ic(a_orb, e)
    ic = [mpf(0), mpf(0), mpf(0), mpf(0),
          r_peri, mpf(0), mpf(0), v_peri]
    rhs = make_rhs_twobody_newtonian(M, q)
    T_kep = 2 * pi * sqrt(a_orb**3 / 1)
    T_end = T_kep * mpf('1.1')
    tol   = mpf(10)**(-(mp.dps - 6))

    t0 = time.time()
    sol = mpmath.odefun(rhs, 0, ic, tol=tol, degree=TAYLOR_DEGREE)

    dt_c = T_end / mpf(n_coarse)
    past_apoapsis = False
    t_lo, t_hi = None, None
    t_prev = mpf(0)
    rd_prev = rdot_rel(ic)

    for k in range(1, n_coarse + 1):
        t_cur = mpf(k) * dt_c
        rd_cur = rdot_rel(sol(t_cur))
        if not past_apoapsis:
            if rd_cur < 0:
                past_apoapsis = True
        else:
            if rd_prev < 0 and rd_cur >= 0:
                t_lo = t_prev
                t_hi = t_cur
                break
        t_prev = t_cur
        rd_prev = rd_cur

    if t_lo is None:
        raise RuntimeError(
            f"[{label}] Newtonian: failed to bracket periapsis."
        )

    for _ in range(n_bisect):
        t_mid = (t_lo + t_hi) / 2
        rd_mid = rdot_rel(sol(t_mid))
        if rd_mid < 0:
            t_lo = t_mid
        else:
            t_hi = t_mid
        if t_hi - t_lo < mpf(10)**(-36):
            break

    t_peri2 = (t_lo + t_hi) / 2

    phi_prev  = atan2(mpf(0), r_peri)
    phi_total = mpf(0)
    dt_a = t_peri2 / mpf(n_angle)

    for k in range(1, n_angle + 1):
        s = sol(mpf(k) * dt_a)
        rx = s[4] - s[0]
        ry = s[5] - s[1]
        phi_cur = atan2(ry, rx)
        dphi = phi_cur - phi_prev
        while dphi >  pi: dphi -= 2*pi
        while dphi <= -pi: dphi += 2*pi
        phi_total += dphi
        phi_prev = phi_cur

    Delta_omega = phi_total - 2*pi
    t1 = time.time()
    print(f"  [GC Newtonian] elapsed {t1-t0:.1f}s dw={mpmath.nstr(Delta_omega, 10)}", flush=True)
    return Delta_omega


# ────────────────────────── fit ──────────────────────────────────────────────

def fit_linear_plus_quad(q_list, D_list):
    """
    Least-squares fit D(q) = C*q + C2*q^2 over the supplied (q, D) pairs.
    Returns (C, C2).
    """
    n = len(q_list)
    A = mpmath.matrix(n, 2)
    b = mpmath.matrix(n, 1)
    for i, (q, D) in enumerate(zip(q_list, D_list)):
        A[i, 0] = q
        A[i, 1] = q**2
        b[i, 0] = D
    x = mpmath.lu_solve(A.T * A, A.T * b)
    return x[0], x[1]   # C, C2


# ────────────────────────── main ─────────────────────────────────────────────

def main():
    t_start = time.time()
    print("error_budget_a2_mass_scaling.py")
    print(f"mp.dps={mp.dps}  a={float(A_ORB)}  e={float(E_ORB)}  taylor_degree={TAYLOR_DEGREE}")
    print()

    M  = mpf('1')    # primary mass
    e  = E_ORB
    a  = A_ORB

    # ── Step 1: fix c from eps = 1e-5 ────────────────────────────────────────
    eps_main = mpf('1e-5')
    c_main   = c_from_eps(eps_main, a, e)
    print(f"eps={float(eps_main):.1e}  c={mpmath.nstr(c_main, 10)}")
    print()

    # ── Step 2: q ladder ─────────────────────────────────────────────────────
    q_values = [mpf('0'), mpf('1e-5'), mpf('1e-4'), mpf('1e-3')]
    dw_results = {}

    for q in q_values:
        label = f"q={float(q):.1e}"
        print(f"{'='*60}")
        print(f"Two-body run: {label}  c={mpmath.nstr(c_main,8)}")
        dw = measure_advance_twobody(M, q, c_main, a, e, label=label)
        dw_results[q] = dw
        print(f"  Delta_omega = {mpmath.nstr(dw, 30)}")
        print()

    dw0 = dw_results[mpf('0')]

    # ── Gate GA: q=0 two-body vs single-body ─────────────────────────────────
    print(f"{'='*60}")
    print("Gate GA: q=0 two-body vs independent single-body")
    dw_single = measure_advance_single(M, c_main, a, e, label='single-body')
    print(f"  dw_twobody(q=0) = {mpmath.nstr(dw0, 30)}")
    print(f"  dw_single       = {mpmath.nstr(dw_single, 30)}")
    diff_GA = abs(dw0 - dw_single)
    if dw_single != 0:
        rel_GA = diff_GA / abs(dw_single)
        sig_GA = int(-log10(rel_GA)) if rel_GA > 0 else mp.dps
    else:
        sig_GA = mp.dps
    ga_pass = sig_GA >= 20
    print(f"  |diff| = {mpmath.nstr(diff_GA, 6)}")
    print(f"  sig_digits = {sig_GA}")
    print(f"  [GA] {'PASS' if ga_pass else 'FAIL'}  (need >= 20 sig digits)")
    if not ga_pass:
        print("  BLOCKED: Gate GA failed — two-body q=0 and single-body disagree.")
        sys.exit(1)
    print()

    # ── D(q) table ────────────────────────────────────────────────────────────
    print("D(q) = [Delta_omega(q) - Delta_omega(0)] / Delta_omega(0):")
    q_nonzero   = [mpf('1e-5'), mpf('1e-4'), mpf('1e-3')]
    D_nonzero   = []
    for q in q_nonzero:
        dw = dw_results[q]
        D  = (dw - dw0) / dw0
        D_nonzero.append(D)
        print(f"  q={float(q):.1e}  Delta_omega={mpmath.nstr(dw, 30)}  D={mpmath.nstr(D, 20)}")
    print()

    # ── Fit D(q) = C*q + C2*q^2 ──────────────────────────────────────────────
    C_fit, C2_fit = fit_linear_plus_quad(q_nonzero, D_nonzero)
    print("Fit D(q) = C*q + C2*q^2:")
    print(f"  C  = {mpmath.nstr(C_fit, 20)}")
    print(f"  C2 = {mpmath.nstr(C2_fit, 20)}")

    # Gate GB: |C2*q^2 / (C*q)| < 0.05 at q=1e-3
    q_test   = mpf('1e-3')
    gb_ratio = abs(C2_fit * q_test**2) / abs(C_fit * q_test) if C_fit != 0 else mpf('inf')
    gb_pass  = gb_ratio < mpf('0.05')
    print(f"  [GB] |C2*q^2/(C*q)| at q=1e-3: {mpmath.nstr(gb_ratio, 6)} < 0.05: "
          f"{'PASS' if gb_pass else 'FAIL (report actual scaling)'}")
    if not gb_pass:
        print(f"  NOTE: quadratic term significant — C*q linear approximation has "
              f"{float(gb_ratio)*100:.1f}% quadratic contamination at q=1e-3")
    print()

    # ── Gate GC: Newtonian null ───────────────────────────────────────────────
    print(f"{'='*60}")
    print("Gate GC: Newtonian two-body (no 1PN), q=1e-3 — expect < 1e-30 rad")
    q_gc  = mpf('1e-3')
    dw_gc = measure_advance_newtonian_twobody(M, q_gc, a, e, label='GC')
    gc_pass = abs(dw_gc) < mpf('1e-30')
    print(f"  |Delta_omega(Newt, q=1e-3)| = {mpmath.nstr(abs(dw_gc), 10)}")
    print(f"  [GC] {'PASS' if gc_pass else 'FAIL'}  (need < 1e-30 rad)")
    if not gc_pass:
        print(f"  BLOCKED: Gate GC failed — Newtonian relative orbit precesses: "
              f"{float(dw_gc):.3e}")
        sys.exit(1)
    print()

    # ── epsilon-independence check ────────────────────────────────────────────
    print(f"{'='*60}")
    print("eps-independence check: repeat q=1e-3 at eps=1e-4")
    eps_hi  = mpf('1e-4')
    c_hi    = c_from_eps(eps_hi, a, e)
    q_eps   = mpf('1e-3')
    print(f"  eps={float(eps_hi):.1e}  c={mpmath.nstr(c_hi,10)}")

    # Need dw0 at eps_hi as well
    print("  Two-body q=0 at eps=1e-4 ...")
    dw0_hi = measure_advance_twobody(M, mpf('0'), c_hi, a, e, label='q=0 eps=1e-4')
    print(f"  dw0(eps=1e-4) = {mpmath.nstr(dw0_hi, 30)}")

    print("  Two-body q=1e-3 at eps=1e-4 ...")
    dw_hi  = measure_advance_twobody(M, q_eps, c_hi, a, e, label='q=1e-3 eps=1e-4')
    print(f"  dw(q=1e-3, eps=1e-4) = {mpmath.nstr(dw_hi, 30)}")

    D_hi  = (dw_hi - dw0_hi) / dw0_hi
    # C at eps_hi from single-point linear fit (D = C*q)
    C_hi  = D_hi / q_eps
    print(f"  D(q=1e-3, eps=1e-4) = {mpmath.nstr(D_hi, 20)}")
    print(f"  C(eps=1e-4)  = {mpmath.nstr(C_hi, 20)}")
    print(f"  C(eps=1e-5)  = {mpmath.nstr(C_fit, 20)}")
    C_rel = abs(C_hi - C_fit) / abs(C_fit) if C_fit != 0 else mpf('inf')
    print(f"  |C(1e-4)-C(1e-5)| / C(1e-5) = {mpmath.nstr(C_rel, 6)}  (expected < 0.05)")
    print()

    # ── Summary ───────────────────────────────────────────────────────────────
    t_total = time.time() - t_start
    print(f"Total wall-time: {t_total/60:.1f} min")
    print()
    print("=" * 72)
    print("MACHINE-READABLE SUMMARY")
    print("=" * 72)
    print(f"mp.dps={mp.dps}  a={float(A_ORB)}  e={float(E_ORB)}  taylor_degree={TAYLOR_DEGREE}")
    print(f"eps_main={float(eps_main):.1e}  c_main={mpmath.nstr(c_main,12)}")
    print()
    print("# Delta_omega table (eps=1e-5)")
    for q in q_values:
        dw = dw_results[q]
        print(f"  q={mpmath.nstr(q,6)}  Delta_omega={mpmath.nstr(dw, 30)}")
    print()
    print("# D(q) table (eps=1e-5)")
    for q, D in zip(q_nonzero, D_nonzero):
        print(f"  q={mpmath.nstr(q,6)}  D={mpmath.nstr(D, 20)}")
    print()
    print(f"C_fit  = {mpmath.nstr(C_fit, 20)}")
    print(f"C2_fit = {mpmath.nstr(C2_fit, 20)}")
    print()
    print("eps-independence pair:")
    print(f"  C(eps=1e-5) = {mpmath.nstr(C_fit, 20)}")
    print(f"  C(eps=1e-4) = {mpmath.nstr(C_hi, 20)}")
    print(f"  relative gap = {mpmath.nstr(C_rel, 6)}")
    print()
    mercury_floor = C_fit * M_MERCURY
    print(f"Mercury floor  C * 1.66e-7 = {mpmath.nstr(mercury_floor, 20)}")
    print()
    print("Gates:")
    print(f"  GA sig_digits={sig_GA}: {'PASS' if ga_pass else 'FAIL'}")
    print(f"  GB |C2*q^2/(C*q)|={mpmath.nstr(gb_ratio,6)}: {'PASS' if gb_pass else 'FAIL (reported scaling above)'}")
    print(f"  GC |dw_Newt|={mpmath.nstr(abs(dw_gc),6)}: {'PASS' if gc_pass else 'FAIL'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
