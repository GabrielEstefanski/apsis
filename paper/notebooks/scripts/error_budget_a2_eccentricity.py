"""
error_budget_a2_eccentricity.py
===============================

A2 follow-up: the eccentricity dependence of the two-body coefficient C.

The closed Phase-A2 result measured C = 8.291 at e = 0.2; its e-dependence
was "not probed, bounds the coefficient at the percent level"
(error_budget_a2_mass_scaling.py). The Q-corrected error budget leaves a
constant -3.6e-8 relative deficit, and the only floor not evaluated at
Mercury's e = 0.20563 is A2 (A1's k_osc(e) is closed-form, already at
Mercury's e). This script measures C(e) on a grid — reusing the A2
machinery verbatim — to test whether C(0.20563) accounts for the deficit.

The measurement is blind to the pre-registered target C(0.20563) ~ 8.074
(validation/audit/ledger.md, 2026-06-12): it integrates the grid and
reports. C(e=0.2) reproducing the committed 8.291 is gate EA.

Run:  python paper/notebooks/scripts/error_budget_a2_eccentricity.py
"""

import time

from error_budget_a2_mass_scaling import (
    M_MERCURY,
    c_from_eps,
    fit_linear_plus_quad,
    measure_advance_twobody,
)
from mpmath import mp, mpf, nstr

mp.dps = 40

A_ORB = mpf("1")
EPS = mpf("1e-5")              # amplified, matching the closed A2 run
Q_LADDER = [mpf("1e-4"), mpf("1e-3")]
E_GRID = [mpf("0.15"), mpf("0.2"), mpf("0.20563"), mpf("0.25"), mpf("0.3")]
E_MERCURY = mpf("0.20563")

# Closed A2 result at e=0.2 (cross-check target for gate EA).
C_AT_E020 = mpf("8.291")

# Empirical target: the Q-corrected floor the ensemble produces, relative,
# at N=500 (validation/mercury-1pn-error-budget, b1). A1+A2 should sum here.
MEASURED_CORRECTED_FLOOR = mpf("1.2461e-6")


def measure_C(e):
    """C(e) from D(q) = C q + C2 q^2 over the q-ladder at fixed eps."""
    c = c_from_eps(EPS, A_ORB, e)
    dw0 = measure_advance_twobody(mpf("1"), mpf("0"), c, A_ORB, e, label=f"e={float(e)} q=0")
    D = []
    for q in Q_LADDER:
        lbl = f"e={float(e)} q={float(q):.0e}"
        dw = measure_advance_twobody(mpf("1"), q, c, A_ORB, e, label=lbl)
        D.append((dw - dw0) / dw0)
    C, C2 = fit_linear_plus_quad(Q_LADDER, D)
    return C, C2


def k_osc(e):
    """Osculating-convention second-order coefficient (A1, closed form)."""
    return mpf("-11") / 6 - 8 * e + e**2 / 3


def main():
    t0 = time.time()
    print("error_budget_a2_eccentricity.py")
    print(f"mp.dps={mp.dps}  eps={float(EPS):.0e}  q_ladder={[float(q) for q in Q_LADDER]}")
    print()

    results = {}
    for e in E_GRID:
        print(f"{'=' * 60}\ne = {float(e)}")
        C, C2 = measure_C(e)
        results[e] = C
        print(f"  C(e={float(e)}) = {nstr(C, 12)}   C2 = {nstr(C2, 8)}")
        print()

    # ── Gate EA: reproduce the committed C(e=0.2) = 8.291 ────────────────────
    c_020 = results[mpf("0.2")]
    ea_rel = abs(c_020 - C_AT_E020) / C_AT_E020
    ea_pass = ea_rel < mpf("1e-3")
    print(f"[EA] C(e=0.2) = {nstr(c_020, 10)} vs committed 8.291: "
          f"rel={nstr(ea_rel, 4)} < 1e-3: {'PASS' if ea_pass else 'FAIL'}")
    print()

    # ── C(e) trend and local slope at Mercury's e ────────────────────────────
    print("C(e) grid:")
    for e in E_GRID:
        print(f"  e={float(e):.5f}  C={nstr(results[e], 12)}")
    e_lo, e_hi = mpf("0.2"), mpf("0.25")
    dCde = (results[e_hi] - results[e_lo]) / (e_hi - e_lo)
    print(f"\n  local slope dC/de near 0.2 = {nstr(dCde, 8)}")
    print()

    # ── Budget closure at Mercury's e (direct measurement, no interpolation) ─
    c_merc = results[E_MERCURY]
    q_merc = M_MERCURY
    a2_old = C_AT_E020 * q_merc
    a2_new = c_merc * q_merc
    eps_merc = 1 / (mpf("10065.3201686") ** 2 * mpf("0.387098") * (1 - E_MERCURY**2))
    a1 = k_osc(E_MERCURY) * eps_merc
    pred_old = a1 + a2_old
    pred_new = a1 + a2_new
    deficit_old = MEASURED_CORRECTED_FLOOR - pred_old
    deficit_new = MEASURED_CORRECTED_FLOOR - pred_new

    print(f"{'=' * 60}")
    print("MACHINE-READABLE SUMMARY")
    print(f"{'=' * 60}")
    ea_tag = "PASS" if ea_pass else "FAIL"
    print(f"C(e=0.20563)        = {nstr(c_merc, 12)}   (pre-registered target ~8.074)")
    print(f"C(e=0.2)            = {nstr(c_020, 12)}   (committed 8.291; EA {ea_tag})")
    print(f"dC/de near 0.2      = {nstr(dCde, 8)}")
    print()
    print(f"eps(Mercury,raw_c)  = {nstr(eps_merc, 8)}")
    print(f"A1 floor k_osc*eps  = {nstr(a1, 8)}")
    print(f"A2 floor C(0.2)*q   = {nstr(a2_old, 8)}   (e=0.2 C)")
    print(f"A2 floor C(0.20563)*q = {nstr(a2_new, 8)}   (Mercury-e C)")
    print()
    print(f"measured corrected floor = {nstr(MEASURED_CORRECTED_FLOOR, 8)}")
    print(f"predicted (e=0.2 C)      = {nstr(pred_old, 8)}   deficit = {nstr(deficit_old, 6)}")
    print(f"predicted (Mercury-e C)  = {nstr(pred_new, 8)}   deficit = {nstr(deficit_new, 6)}")
    print()
    shrink = abs(deficit_new) / abs(deficit_old) if deficit_old != 0 else mpf("inf")
    print(f"deficit |new|/|old| = {nstr(shrink, 6)}")
    closed = abs(deficit_new) < abs(deficit_old) / 3
    print(f"VERDICT: e-dependence of C {'ACCOUNTS FOR' if closed else 'DOES NOT ACCOUNT FOR'} "
          f"the deficit (>3x shrink: {closed})")
    print(f"\nTotal wall-time: {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
