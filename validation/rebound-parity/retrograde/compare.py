"""Comparator for REBOUND parity — Retrograde Kepler e=0.5.

Loads `out/apsis.csv` and `out/rebound.csv`, computes the metrics defined
in the protocol's §Hypothesis (Tier 1 magnitude invariants + Tier 2 sign-
consistency binary checks + Tier 3 informational |Δr|) over the 10^4-orbit
long horizon, and exits 0 iff every gated metric is within tolerance.
Prints the §Decision rules outcome interpretation at the end.

## Tier structure (from the protocol notebook §Hypothesis)

Tier 1 — magnitude invariants (gated, 7 metrics, identical to Kepler-prograde):
    - |Δa|/a max over both sides
    - |Δe| max over both sides
    - |Δω| max over both sides (rad)
    - ||h| - |h_0|| / |h_0| max over both sides   (magnitude only)
    - |ΔE/E_0| apsis
    - |ΔE/E_0| rebound
    - cross-impl |ΔE|/|E_0|

Tier 2 — orientation invariants (gated, 3 binary checks):
    - apsis sign(h) consistency: sign(h(t)) = sign(h_0) ∀t  AND  |h(t)| > eps_floor ∀t
    - rebound sign(h) consistency: same on REBOUND side
    - cross-impl sign agreement: sign(h_apsis(t)) = sign(h_rebound(t)) ∀t

Tier 3 — geometric coherence (informational, NOT gated):
    - |Δr| max over horizon

## Two-horizon evaluation

The CSV holds 10001 samples spanning 10^4 orbits; the comparator gates the
full sample set. (A 100-orbit checkpoint was removed — it was bit-identical
to the prograde Kepler scenario, so it validated nothing new.)

## Exit codes

- 0 — all gated metrics within tolerance.
- 1 — input file error (missing file, sample count mismatch).
- 2 — at least one gated metric exceeded tolerance at either horizon.

## Run

    python compare.py
    python compare.py --apsis-csv path/to/apsis.csv --rebound-csv path/to/rebound.csv

## Protocol notebook

    docs/experiments/2026-05-01-rebound-parity-retrograde.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

# ── Protocol constants ──────────────────────────────────────────────────── #

M_PRIMARY: float = 1.0
M_SECONDARY: float = 1.0e-6
MU: float = M_PRIMARY + M_SECONDARY

N_ORBITS_LONG: int = 10_000
N_ORBITS_CHECKPOINT: int = 100  # baseline horizon for the √N round-off floor

# ── Tolerances ─────────────────────────────────────────────────────────── #
# Single long-horizon gate. The IAS15 cross-impl round-off floor grows as √N
# (Brouwer's law; Rein & Spiegel 2015): F = 13·EPS·√(N_long/N_check) ≈ 130·EPS.
# Gate 5×F — kept tight to detect departure from √N (a ∝N drift); ω adds the
# atan2 1/e condition factor.
EPS: float = 2.220446049250313e-16
E_ECC: float = 0.5  # IC eccentricity (mirrors the Rust example)
_F_LONG: float = 13.0 * EPS * math.sqrt(N_ORBITS_LONG / N_ORBITS_CHECKPOINT)

TOL_RELATIVE_SEMIAXIS: float = 5.0 * _F_LONG
TOL_ECCENTRICITY: float = 5.0 * _F_LONG
TOL_PERIAPSIS_OMEGA: float = 5.0 * _F_LONG / E_ECC
TOL_RELATIVE_ANGULAR_MOMENTUM: float = 5.0 * _F_LONG
TOL_RELATIVE_ENERGY_DRIFT: float = 5.0 * _F_LONG
TOL_CROSS_IMPL_ENERGY: float = 5.0 * _F_LONG

# Near-zero sign guard: h flips sign if retrograde handling breaks. ~10^10
# below |h_0| ≈ 0.866 — defensive, not a routine threshold.
EPS_FLOOR_H: float = 1.0e-10


# ── Data records ───────────────────────────────────────────────────────── #


@dataclass
class Sample:
    orbit: int
    t: float
    x0: float
    y0: float
    vx0: float
    vy0: float
    x1: float
    y1: float
    vx1: float
    vy1: float
    e_total: float


@dataclass
class OrbitalState:
    """`(a, e, ω, h)` of the secondary's orbit relative to the primary."""

    a: float
    e: float
    omega: float
    h: float  # signed; sign-checks are performed against the signed value


@dataclass
class MetricResult:
    name: str
    observed: float
    tolerance: float
    passed: bool
    detail: dict = field(default_factory=dict)


@dataclass
class BinaryCheckResult:
    name: str
    passed: bool
    detail: dict = field(default_factory=dict)


@dataclass
class HorizonReport:
    label: str
    n_samples: int
    t_final: float
    tier1: list[MetricResult]
    tier2: list[BinaryCheckResult]
    tier3_info: MetricResult
    all_gated_passed: bool


# ── Orbital element extraction ─────────────────────────────────────────── #


def relative_elements(s: Sample, mu: float) -> OrbitalState:
    """Compute orbital elements of the secondary relative to the primary.

    2D Kepler reduction; identical to the Kepler-prograde comparator. The
    `h` returned here is the signed z-component; Tier 1 metrics consume
    `|h|`, Tier 2 metrics consume the signed value.
    """
    rx = s.x1 - s.x0
    ry = s.y1 - s.y0
    vrx = s.vx1 - s.vx0
    vry = s.vy1 - s.vy0

    r = math.sqrt(rx * rx + ry * ry)
    v_sq = vrx * vrx + vry * vry

    eps = 0.5 * v_sq - mu / r
    a = -mu / (2.0 * eps)
    h = rx * vry - ry * vrx
    e_sq = max(0.0, 1.0 - h * h / (mu * a))
    e = math.sqrt(e_sq)
    r_dot_v = rx * vrx + ry * vry
    ev_x = ((v_sq - mu / r) * rx - r_dot_v * vrx) / mu
    ev_y = ((v_sq - mu / r) * ry - r_dot_v * vry) / mu
    omega = math.atan2(ev_y, ev_x)

    return OrbitalState(a=a, e=e, omega=omega, h=h)


def angle_diff(a: float, b: float) -> float:
    """Smallest signed difference `a − b` wrapped into [−π, π]."""
    d = math.fmod(a - b + math.pi, 2.0 * math.pi)
    if d < 0.0:
        d += 2.0 * math.pi
    return d - math.pi


# ── Per-horizon evaluation ─────────────────────────────────────────────── #


def evaluate_horizon(
    label: str,
    apsis: list[Sample],
    rebound: list[Sample],
    elem_apsis: list[OrbitalState],
    elem_rebound: list[OrbitalState],
) -> HorizonReport:
    """Evaluate Tier 1 + Tier 2 + Tier 3 over a sample range."""
    a0 = elem_apsis[0].a
    h0_signed = elem_apsis[0].h
    h0_mag = abs(h0_signed)

    # ── Tier 1: magnitude invariants ────────────────────────────────────── #

    # Δa per side, then take max-of-sides for cross-impl reporting style
    # consistent with Kepler-prograde.
    max_da_rel = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_da_rel = max(max_da_rel, abs(ea.a - er.a) / abs(a0))
    m_a = MetricResult(
        name="|Δa|/a (semi-major axis)",
        observed=max_da_rel,
        tolerance=TOL_RELATIVE_SEMIAXIS,
        passed=max_da_rel <= TOL_RELATIVE_SEMIAXIS,
    )

    max_de = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_de = max(max_de, abs(ea.e - er.e))
    m_e = MetricResult(
        name="|Δe| (eccentricity)",
        observed=max_de,
        tolerance=TOL_ECCENTRICITY,
        passed=max_de <= TOL_ECCENTRICITY,
    )

    max_domega = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_domega = max(max_domega, abs(angle_diff(ea.omega, er.omega)))
    m_om = MetricResult(
        name="|Δω| (periapsis orient.)",
        observed=max_domega,
        tolerance=TOL_PERIAPSIS_OMEGA,
        passed=max_domega <= TOL_PERIAPSIS_OMEGA,
    )

    # |h| magnitude drift — uses absolute values on both sides, separately
    # from the sign tier. The cross-impl form below mirrors the prograde
    # comparator's |Δh|/h structure but with absolute values everywhere.
    max_dh_mag_rel = 0.0
    for ea, er in zip(elem_apsis, elem_rebound):
        max_dh_mag_rel = max(
            max_dh_mag_rel, abs(abs(ea.h) - abs(er.h)) / h0_mag
        )
    m_h_mag = MetricResult(
        name="||h|−|h_0||/|h_0| (cross-impl)",
        observed=max_dh_mag_rel,
        tolerance=TOL_RELATIVE_ANGULAR_MOMENTUM,
        passed=max_dh_mag_rel <= TOL_RELATIVE_ANGULAR_MOMENTUM,
    )

    e0_apsis = apsis[0].e_total
    e0_rebound = rebound[0].e_total
    max_drift_apsis = max(abs(s.e_total - e0_apsis) / abs(e0_apsis) for s in apsis)
    max_drift_rebound = max(abs(s.e_total - e0_rebound) / abs(e0_rebound) for s in rebound)
    max_cross_e = max(
        abs(a.e_total - r.e_total) / abs(e0_apsis)
        for a, r in zip(apsis, rebound)
    )
    m_ea = MetricResult(
        name="|ΔE/E_0| apsis",
        observed=max_drift_apsis,
        tolerance=TOL_RELATIVE_ENERGY_DRIFT,
        passed=max_drift_apsis <= TOL_RELATIVE_ENERGY_DRIFT,
    )
    m_er = MetricResult(
        name="|ΔE/E_0| rebound",
        observed=max_drift_rebound,
        tolerance=TOL_RELATIVE_ENERGY_DRIFT,
        passed=max_drift_rebound <= TOL_RELATIVE_ENERGY_DRIFT,
    )
    m_cross_e = MetricResult(
        name="cross-impl |ΔE|/|E_0|",
        observed=max_cross_e,
        tolerance=TOL_CROSS_IMPL_ENERGY,
        passed=max_cross_e <= TOL_CROSS_IMPL_ENERGY,
    )

    tier1 = [m_a, m_e, m_om, m_h_mag, m_ea, m_er, m_cross_e]

    # ── Tier 2: sign consistency (binary checks) ────────────────────────── #

    sign_h0 = 1 if h0_signed > 0 else (-1 if h0_signed < 0 else 0)

    apsis_sign_ok = True
    apsis_floor_ok = True
    apsis_first_failure_orbit = None
    for ea in elem_apsis:
        if (math.copysign(1.0, ea.h) > 0) != (sign_h0 > 0):
            apsis_sign_ok = False
            if apsis_first_failure_orbit is None:
                apsis_first_failure_orbit = elem_apsis.index(ea)
                break
        if abs(ea.h) <= EPS_FLOOR_H:
            apsis_floor_ok = False
            if apsis_first_failure_orbit is None:
                apsis_first_failure_orbit = elem_apsis.index(ea)
                break
    apsis_sign_consistency = apsis_sign_ok and apsis_floor_ok
    t2_apsis = BinaryCheckResult(
        name="apsis sign(h) consistency",
        passed=apsis_sign_consistency,
        detail={
            "sign_consistent": apsis_sign_ok,
            "floor_respected": apsis_floor_ok,
            "first_failure_orbit": apsis_first_failure_orbit,
            "h_0_sign": sign_h0,
            "eps_floor": EPS_FLOOR_H,
        },
    )

    rebound_sign_ok = True
    rebound_floor_ok = True
    rebound_first_failure_orbit = None
    for er in elem_rebound:
        if (math.copysign(1.0, er.h) > 0) != (sign_h0 > 0):
            rebound_sign_ok = False
            if rebound_first_failure_orbit is None:
                rebound_first_failure_orbit = elem_rebound.index(er)
                break
        if abs(er.h) <= EPS_FLOOR_H:
            rebound_floor_ok = False
            if rebound_first_failure_orbit is None:
                rebound_first_failure_orbit = elem_rebound.index(er)
                break
    rebound_sign_consistency = rebound_sign_ok and rebound_floor_ok
    t2_rebound = BinaryCheckResult(
        name="rebound sign(h) consistency",
        passed=rebound_sign_consistency,
        detail={
            "sign_consistent": rebound_sign_ok,
            "floor_respected": rebound_floor_ok,
            "first_failure_orbit": rebound_first_failure_orbit,
            "h_0_sign": sign_h0,
            "eps_floor": EPS_FLOOR_H,
        },
    )

    cross_sign_ok = True
    cross_first_disagreement = None
    for idx, (ea, er) in enumerate(zip(elem_apsis, elem_rebound)):
        if (math.copysign(1.0, ea.h) > 0) != (math.copysign(1.0, er.h) > 0):
            cross_sign_ok = False
            cross_first_disagreement = idx
            break
    t2_cross = BinaryCheckResult(
        name="cross-impl sign(h) agreement",
        passed=cross_sign_ok,
        detail={"first_disagreement_orbit": cross_first_disagreement},
    )

    tier2 = [t2_apsis, t2_rebound, t2_cross]

    # ── Tier 3: informational |Δr| ──────────────────────────────────────── #

    max_dr = 0.0
    max_dr_orbit = 0
    for a, r in zip(apsis, rebound):
        dx = a.x1 - r.x1
        dy = a.y1 - r.y1
        dr = (dx * dx + dy * dy) ** 0.5
        if dr > max_dr:
            max_dr = dr
            max_dr_orbit = a.orbit
    info = MetricResult(
        name="|Δr| (secondary, NOT gated)",
        observed=max_dr,
        tolerance=float("inf"),
        passed=True,
        detail={
            "argmax_orbit": max_dr_orbit,
            "note": "phase-drift contaminated; not invariant across adaptive integrators",
        },
    )

    all_passed = all(m.passed for m in tier1) and all(b.passed for b in tier2)

    return HorizonReport(
        label=label,
        n_samples=len(apsis),
        t_final=apsis[-1].t,
        tier1=tier1,
        tier2=tier2,
        tier3_info=info,
        all_gated_passed=all_passed,
    )


# ── Decision-rules interpretation ──────────────────────────────────────── #


def decision_rule(report: HorizonReport) -> tuple[str, str]:
    """Map a horizon's outcome to (label, action) per §Decision rules."""
    tier1_pass = all(m.passed for m in report.tier1)
    tier2_pass = all(b.passed for b in report.tier2)
    if tier1_pass and tier2_pass:
        return ("PASS", "Tier 1 + Tier 2 both pass — integrator + sign convention OK")
    if not tier1_pass and tier2_pass:
        # Detect Brouwer-law saturation case: only energy/h drift fails, by
        # being just above 1e-13.
        e_metrics = [m for m in report.tier1 if "|ΔE" in m.name or "|h|" in m.name]
        all_e_close = all(
            (not m.passed) and m.observed < 5.0 * m.tolerance for m in e_metrics if not m.passed
        )
        if all_e_close and any(not m.passed for m in e_metrics):
            return (
                "BROUWER-SATURATION",
                "Tier 1 magnitude bound approached/exceeded by < 5×; consistent with "
                "Brouwer-law saturation envelope. Document honestly; do NOT widen bound.",
            )
        return (
            "TIER1-FAIL",
            "Magnitude-drift bug (energy or radial bookkeeping). Halt, localise to "
            "inner force / integration loop. Re-run prograde at same horizon to "
            "discriminate retrograde-specific from regime-driven cause.",
        )
    if tier1_pass and not tier2_pass:
        return (
            "TIER2-FAIL",
            "Sign-convention bug. Halt. Inspect cross-product order, eccentricity-vector "
            "composition, atan2 argument ordering, controller sign assumptions, "
            "underflow/overflow paths. First failing sample localises bug expression time.",
        )
    return (
        "DEEP-FAIL",
        "Tier 1 + Tier 2 both fail. Likely IC handling or state representation defect. "
        "Halt. Verify IC bit-identicality at t=0 and COM-shift preservation of sign(h).",
    )


# ── Main ───────────────────────────────────────────────────────────────── #


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    args = parse_args()
    apsis_path = Path(args.apsis_csv).resolve()
    rebound_path = Path(args.rebound_csv).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        apsis = load_csv(apsis_path)
        rebound = load_csv(rebound_path)
    except FileNotFoundError as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1

    if len(apsis) != len(rebound):
        print(
            f"ERROR: sample count mismatch — apsis has {len(apsis)}, "
            f"rebound has {len(rebound)}",
            file=sys.stderr,
        )
        return 1
    if len(apsis) != N_ORBITS_LONG + 1:
        print(
            f"ERROR: expected {N_ORBITS_LONG + 1} samples, got {len(apsis)}",
            file=sys.stderr,
        )
        return 1

    # ── Compute orbital elements per side, full set ─────────────────────── #
    elem_apsis = [relative_elements(s, MU) for s in apsis]
    elem_rebound = [relative_elements(s, MU) for s in rebound]

    # ── Long-horizon evaluation (10^4 orbits) ───────────────────────────── #
    # The 100-orbit checkpoint was removed: it was bit-identical to the prograde
    # Kepler scenario (round-off is sign-agnostic), validating nothing new. The
    # long-horizon carries the unique value (Brouwer √N + sign convention).
    rep_long = evaluate_horizon(
        f"long-horizon gate ({N_ORBITS_LONG} orbits)",
        apsis,
        rebound,
        elem_apsis,
        elem_rebound,
    )

    print_report(rep_long, apsis[0].e_total, rebound[0].e_total)
    write_json_report(output_dir, rep_long)

    return 0 if rep_long.all_gated_passed else 2


# ── I/O ─────────────────────────────────────────────────────────────────── #


def load_csv(path: Path) -> list[Sample]:
    if not path.exists():
        raise FileNotFoundError(f"CSV not found at {path}")
    samples: list[Sample] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            samples.append(
                Sample(
                    orbit=int(row["orbit"]),
                    t=float(row["t"]),
                    x0=float(row["x0"]),
                    y0=float(row["y0"]),
                    vx0=float(row["vx0"]),
                    vy0=float(row["vy0"]),
                    x1=float(row["x1"]),
                    y1=float(row["y1"]),
                    vx1=float(row["vx1"]),
                    vy1=float(row["vy1"]),
                    e_total=float(row["e_total"]),
                )
            )
    return samples


def print_horizon(rep: HorizonReport) -> None:
    print()
    print(f"  ── {rep.label} ──")
    print(f"  samples: {rep.n_samples}    t_final: {rep.t_final:.6e}")
    print()
    print(f"  Tier 1 (magnitude, gated)")
    print(f"  {'metric':<32} {'observed':>14} {'tolerance':>14}  verdict")
    print(f"  {'-' * 32} {'-' * 14} {'-' * 14}  {'-' * 7}")
    for m in rep.tier1:
        verdict = "pass" if m.passed else "FAIL"
        print(f"  {m.name:<32} {m.observed:>14.3e} {m.tolerance:>14.3e}  {verdict}")
    print()
    print(f"  Tier 2 (sign(h), gated, binary)")
    print(f"  {'check':<32} {'verdict':>14}  detail")
    print(f"  {'-' * 32} {'-' * 14}  {'-' * 7}")
    for b in rep.tier2:
        verdict = "pass" if b.passed else "FAIL"
        first_fail = b.detail.get("first_failure_orbit") or b.detail.get("first_disagreement_orbit")
        detail = "ok" if b.passed else f"first fail at orbit {first_fail}"
        print(f"  {b.name:<32} {verdict:>14}  {detail}")
    print()
    print(f"  Tier 3 (informational, NOT gated)")
    print(f"  {'metric':<32} {'observed':>14}")
    print(f"  {'-' * 32} {'-' * 14}")
    info = rep.tier3_info
    print(f"  {info.name:<32} {info.observed:>14.3e}")
    argmax = info.detail.get("argmax_orbit", "n/a")
    print(f"    └── peak at orbit {argmax}; phase-drift contaminated")


def print_report(
    rep: HorizonReport,
    e0_apsis: float,
    e0_rebound: float,
) -> None:
    print()
    print("REBOUND parity — Retrograde Kepler e=0.5 — comparison report")
    print(f"  E_0 apsis   : {e0_apsis:+.18e}")
    print(f"  E_0 rebound : {e0_rebound:+.18e}")

    print_horizon(rep)

    print()
    print("  ── Decision rules outcome (per §Decision rules in protocol) ──")
    verdict, action = decision_rule(rep)
    print(f"  [{rep.label}] {verdict}: {action}")
    print()


def write_json_report(
    output_dir: Path,
    rep: HorizonReport,
) -> None:
    report = {
        "all_passed": rep.all_gated_passed,
        "horizons": {
            rep.label: horizon_to_dict(rep),
        },
    }
    (output_dir / "comparison.json").write_text(json.dumps(report, indent=2))


def horizon_to_dict(rep: HorizonReport) -> dict:
    return {
        "n_samples": rep.n_samples,
        "t_final": rep.t_final,
        "all_gated_passed": rep.all_gated_passed,
        "tier1_magnitude": [asdict(m) for m in rep.tier1],
        "tier2_sign": [asdict(b) for b in rep.tier2],
        "tier3_info": asdict(rep.tier3_info),
        "decision_rule": dict(zip(("verdict", "action"), decision_rule(rep))),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Comparator for REBOUND vs apsis parity — Retrograde Kepler e=0.5.",
    )
    parser.add_argument(
        "--apsis-csv",
        default="out/apsis.csv",
        help="Path to the apsis-side CSV. Default: out/apsis.csv.",
    )
    parser.add_argument(
        "--rebound-csv",
        default="out/rebound.csv",
        help="Path to the REBOUND-side CSV. Default: out/rebound.csv.",
    )
    parser.add_argument(
        "--output-dir",
        default="out",
        help="Directory for the JSON report. Default: out/.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
