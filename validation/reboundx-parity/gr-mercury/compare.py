"""Comparator for REBOUNDx parity — Sun-Mercury 1PN (apsis-1pn vs REBOUNDx gr).

Gates the *parity* (cross-implementation agreement); the analytic 43"/century
is reported as the physics anchor but NOT gated here — accuracy is owned by
`crates/apsis-1pn/tests/mercury_precession_gate.rs`. Following the
rebound-parity lesson, |Δr| at fixed times conflates orbital phase drift with
geometric drift under adaptive IAS15, so the gated metrics are orbital
invariants (a, e, h), cross-implementation energy, and — for 1PN — the secular
precession rate (the gauge-invariant observable). |Δr| is informational.

Gate bounds (a priori; see paper/notebooks/2026-05-29-reboundx-parity-gr.md):
  invariants/energy cross-impl   ≤ 1e-13  (ULP floor, as the rebound-parity scenarios)
  precession apsis-vs-reboundx   ≤ 2e-5   (gr mode; empirical formulation floor ~2e-6, 10x)

Exit: 0 all gated metrics pass; 1 input error; 2 a gated metric failed.

Run:
    python compare.py --apsis out/apsis.csv --rebound out/rebound.csv --mode gr
    python compare.py --apsis out/apsis_kepler.csv --rebound out/rebound_kepler.csv --mode control
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

MU: float = 1.0 + 1.660114e-7
A_MERCURY: float = 0.387098
E_MERCURY: float = 0.20563
C_LIGHT: float = 1.006513002441656681e4  # C_SOLAR_UNITS (matches the apsis side)

ORBITS_PER_CENTURY: float = 36525.0 / 87.969
RAD_TO_ARCSEC: float = 180.0 * 3600.0 / math.pi

TOL_INVARIANT_CROSS: float = 1.0e-13
TOL_ENERGY_CROSS: float = 1.0e-13
TOL_PRECESSION_PARITY: float = 2.0e-5  # gr mode only; relative apsis-vs-reboundx


def load(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found at {path}")
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        return [{k: float(v) for k, v in row.items()} for row in reader]


def elements(s: dict[str, float]) -> tuple[float, float, float, float]:
    """Osculating (a, e, ω, h) of Mercury relative to the Sun."""
    rx, ry = s["x1"] - s["x0"], s["y1"] - s["y0"]
    vrx, vry = s["vx1"] - s["vx0"], s["vy1"] - s["vy0"]
    r = math.hypot(rx, ry)
    v2 = vrx * vrx + vry * vry
    a = -MU / (2.0 * (0.5 * v2 - MU / r))
    h = rx * vry - ry * vrx
    e = math.sqrt(max(0.0, 1.0 - h * h / (MU * a)))
    rdv = rx * vrx + ry * vry
    ev_x = ((v2 - MU / r) * rx - rdv * vrx) / MU
    ev_y = ((v2 - MU / r) * ry - rdv * vry) / MU
    return a, e, math.atan2(ev_y, ev_x), h


def precession_arcsec_per_century(samples: list[dict[str, float]]) -> float:
    omegas = [elements(s)[2] for s in samples]
    total = 0.0
    for prev, cur in zip(omegas, omegas[1:]):
        step = cur - prev
        while step > math.pi:
            step -= 2.0 * math.pi
        while step <= -math.pi:
            step += 2.0 * math.pi
        total += step
    return (total / (len(samples) - 1)) * ORBITS_PER_CENTURY * RAD_TO_ARCSEC


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        # TextIO stub lacks reconfigure (TextIOWrapper-only) — known typeshed gap.
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # pyright: ignore[reportAttributeAccessIssue, reportUnknownMemberType]
    args = parse_args()
    try:
        apsis = load(Path(args.apsis).resolve())
        rebound = load(Path(args.rebound).resolve())
    except FileNotFoundError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1
    if len(apsis) != len(rebound):
        print(f"ERROR: sample mismatch apsis={len(apsis)} rebound={len(rebound)}", file=sys.stderr)
        return 1

    ea = [elements(s) for s in apsis]
    er = [elements(s) for s in rebound]
    a0, h0 = ea[0][0], ea[0][3]
    max_da = max(abs(x[0] - y[0]) / abs(a0) for x, y in zip(ea, er))
    max_de = max(abs(x[1] - y[1]) for x, y in zip(ea, er))
    max_dh = max(abs(x[3] - y[3]) / abs(h0) for x, y in zip(ea, er))

    e0a, e0r = apsis[0]["e_total"], rebound[0]["e_total"]
    drift_a = max(abs(s["e_total"] - e0a) / abs(e0a) for s in apsis)
    drift_r = max(abs(s["e_total"] - e0r) / abs(e0r) for s in rebound)
    cross_e = max(abs(a["e_total"] - r["e_total"]) / abs(e0a) for a, r in zip(apsis, rebound))
    max_dr = max(math.hypot(a["x1"] - r["x1"], a["y1"] - r["y1"]) for a, r in zip(apsis, rebound))

    prec_a = precession_arcsec_per_century(apsis)
    prec_r = precession_arcsec_per_century(rebound)
    prec_analytic = (6.0 * math.pi / (C_LIGHT**2 * A_MERCURY * (1.0 - E_MERCURY**2))) \
        * ORBITS_PER_CENTURY * RAD_TO_ARCSEC

    # ── Gated metrics ───────────────────────────────────────────────────────
    gated: list[tuple[str, float, float]] = [
        ("|Δa|/a cross", max_da, TOL_INVARIANT_CROSS),
        ("|Δe| cross", max_de, TOL_INVARIANT_CROSS),
        ("|Δh|/h cross", max_dh, TOL_INVARIANT_CROSS),
        ("|ΔE|/E0 cross", cross_e, TOL_ENERGY_CROSS),
    ]
    if args.mode == "gr":
        prec_parity = abs(prec_a / prec_r - 1.0)
        gated.append(("precession apsis-vs-reboundx", prec_parity, TOL_PRECESSION_PARITY))

    results = [(name, obs, tol, obs <= tol) for name, obs, tol in gated]
    all_passed = all(p for _, _, _, p in results)

    print(f"\n=== REBOUNDx parity — Sun-Mercury — {args.mode} ({len(apsis)} samples) ===")
    print(f"  {'gated metric':<32} {'observed':>12} {'tol':>12}  verdict")
    print(f"  {'-' * 32} {'-' * 12} {'-' * 12}  -------")
    for name, obs, tol, passed in results:
        print(f"  {name:<32} {obs:>12.3e} {tol:>12.3e}  {'pass' if passed else 'FAIL'}")
    print("  reported (not gated):")
    print(f"    |Δr| (phase-contaminated)   = {max_dr:.3e}")
    print(f"    |ΔE/E0| per side: apsis={drift_a:.3e} reboundx={drift_r:.3e}")
    print(f"    precession arcsec/cy: apsis={prec_a:+.4f} reboundx={prec_r:+.4f} "
          f"analytic={prec_analytic:+.4f}")
    if args.mode == "gr":
        print(f"    apsis vs analytic={prec_a / prec_analytic - 1.0:+.2e}  "
              f"reboundx vs analytic={prec_r / prec_analytic - 1.0:+.2e}")
    print(f"  -> {'ALL GATED METRICS PASS' if all_passed else 'GATE FAILURE'}\n")

    report = {
        "mode": args.mode, "n_samples": len(apsis), "all_passed": all_passed,
        "gated": [{"name": n, "observed": o, "tolerance": t, "passed": p} for n, o, t, p in results],
        "reported": {
            "max_dr": max_dr, "drift_apsis": drift_a, "drift_rebound": drift_r,
            "precession_apsis": prec_a, "precession_rebound": prec_r,
            "precession_analytic": prec_analytic,
        },
    }
    Path(args.apsis).resolve().parent.joinpath(f"comparison_{args.mode}.json").write_text(
        json.dumps(report, indent=2)
    )
    return 0 if all_passed else 2


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="REBOUNDx parity comparator (Sun-Mercury 1PN).")
    p.add_argument("--apsis", default="out/apsis.csv")
    p.add_argument("--rebound", default="out/rebound.csv")
    p.add_argument("--mode", choices=["gr", "control"], default="gr")
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(main())
