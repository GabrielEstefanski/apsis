"""Comparator for REBOUND parity -- Plummer cluster (softened kernel, 3D).

Computes every invariant itself from the long-format snapshot CSVs -- the
softened Hamiltonian's energy (PE with K(r) = 1/sqrt(r^2+eps^2)), 3D angular
momentum, linear momentum, and centre of mass -- plus the Tier-3 statistical
observables (virial ratio, Lagrangian radii, per-body cross-impl drift).
Neither side's internal energy bookkeeping enters any metric (REBOUND's
sim.energy() omits the softening term; verified v4.6.0).

Phase 0 (--informational): every metric and BOTH candidate L/P floor models
are reported; the exit code never signals metric failure. The gated phase-1
verdict is enabled only after the floor model is frozen in the protocol
notebook.

Run:
    python compare.py --selftest
    python compare.py --ics ics_n256.csv --informational

Protocol notebook:
    paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

EPS64 = 2.220446049250313e-16


def read_ics(path: Path) -> tuple[np.ndarray, float, int]:
    eps = None
    masses: list[float] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("# eps="):
                    eps = float(line.split("=", 1)[1])
                continue
            if line.startswith("body,"):
                continue
            masses.append(float(line.split(",")[1]))
    if eps is None:
        raise ValueError(f"no '# eps=' header in {path}")
    return np.array(masses), eps, len(masses)


def load_snapshots(path: Path, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (times[S], pos[S,N,3], vel[S,N,3]) from the long-format CSV.

    Long CSV columns: sample,t,body,x,y,z,vx,vy,vz  (indices 0..8)
    Positions are columns 3:6; velocities are columns 6:9.
    """
    rows = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("sample,"):
                continue
            rows.append([float(v) for v in line.rstrip("\n").split(",")])
    arr = np.array(rows)
    n_samples = arr.shape[0] // n
    assert arr.shape[0] == n_samples * n, (
        f"row count {arr.shape[0]} not divisible by n={n}"
    )
    arr = arr.reshape(n_samples, n, 9)
    times = arr[:, 0, 1].copy()
    pos = arr[:, :, 3:6].copy()
    vel = arr[:, :, 6:9].copy()
    return times, pos, vel


def softened_pe(masses: np.ndarray, pos: np.ndarray, eps: float) -> float:
    n = len(masses)
    pe = 0.0
    for i in range(n):
        d = pos[i + 1 :] - pos[i]
        r2 = np.einsum("ij,ij->i", d, d)
        pe -= float(np.sum(masses[i] * masses[i + 1 :] / np.sqrt(r2 + eps * eps)))
    return pe


def kinetic(masses: np.ndarray, vel: np.ndarray) -> float:
    return 0.5 * float(np.sum(masses * np.einsum("ij,ij->i", vel, vel)))


def ang_mom(masses: np.ndarray, pos: np.ndarray, vel: np.ndarray) -> np.ndarray:
    return np.sum(masses[:, None] * np.cross(pos, vel), axis=0)


def lin_mom(masses: np.ndarray, vel: np.ndarray) -> np.ndarray:
    return np.sum(masses[:, None] * vel, axis=0)


def com(masses: np.ndarray, pos: np.ndarray) -> np.ndarray:
    return np.average(pos, axis=0, weights=masses)


def scale_p(masses: np.ndarray, vel: np.ndarray) -> float:
    """Wilkinson cancellation scale of P = sum m_i v_i (0 by IC)."""
    return float(np.linalg.norm(np.sum(masses[:, None] * np.abs(vel), axis=0)))


def scale_l(masses: np.ndarray, pos: np.ndarray, vel: np.ndarray) -> float:
    """Wilkinson cancellation scale of L = sum m_i (r_i x v_i)."""
    return float(
        np.sum(
            masses[:, None]
            * (
                np.abs(pos[:, [1, 2, 0]] * vel[:, [2, 0, 1]])
                + np.abs(pos[:, [2, 0, 1]] * vel[:, [1, 2, 0]])
            )
        )
    )


def lagrangian_radii(
    masses: np.ndarray, pos: np.ndarray, fracs=(0.1, 0.5, 0.9)
) -> list[float]:
    c = com(masses, pos)
    r = np.sort(np.sqrt(np.einsum("ij,ij->i", pos - c, pos - c)))
    n = len(r)
    return [float(r[max(0, math.ceil(f * n) - 1)]) for f in fracs]


def _selftest() -> None:
    # softened two-body with exact closed form: PE = -(0.5*0.5)/sqrt(0.6^2+0.8^2) = -0.25
    m2 = np.array([0.5, 0.5])
    p2 = np.array([[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]])
    assert softened_pe(m2, p2, 0.8) == -0.25, softened_pe(m2, p2, 0.8)

    # 3D two-body, integer-exact L, P, COM
    m = np.array([2.0, 3.0])
    p = np.array([[1.0, 2.0, 3.0], [-1.0, 0.0, 2.0]])
    v = np.array([[4.0, 5.0, 6.0], [0.0, -1.0, 1.0]])
    assert np.array_equal(lin_mom(m, v), np.array([8.0, 7.0, 15.0])), lin_mom(m, v)
    assert np.allclose(com(m, p), np.array([-0.2, 0.8, 2.4]), rtol=0, atol=1e-16), com(m, p)
    assert np.array_equal(ang_mom(m, p, v), np.array([0.0, 15.0, -3.0])), ang_mom(m, p, v)

    # KE: 0.5*2*(16+25+36) + 0.5*3*(0+1+1) = 77 + 3 = 80
    assert kinetic(m, v) == 80.0, kinetic(m, v)
    print("selftest: ok")


def main() -> int:
    args = parse_args()
    _selftest()
    if args.selftest:
        return 0

    masses, eps, n = read_ics(Path(args.ics))
    t_a, pos_a, vel_a = load_snapshots(Path(args.apsis_csv), n)
    t_r, pos_r, vel_r = load_snapshots(Path(args.rebound_csv), n)

    if len(t_a) != len(t_r):
        print(f"ERROR: sample count mismatch {len(t_a)} vs {len(t_r)}", file=sys.stderr)
        return 1
    if float(np.max(np.abs(t_a - t_r))) > 1e-12:
        print("ERROR: sample times disagree -- exact_finish_time failed", file=sys.stderr)
        return 1

    n_steps_a = read_steps(Path(args.apsis_stats), "substeps_total")
    n_steps_r = read_steps(Path(args.rebound_stats), "steps_done")
    n_pair = n * (n - 1) // 2

    n_s = len(t_a)
    # compute (ke, pe) once per side; derive energy and virial ratio from the cached pairs
    ke_pe_a = [(kinetic(masses, vel_a[s]), softened_pe(masses, pos_a[s], eps)) for s in range(n_s)]
    ke_pe_r = [(kinetic(masses, vel_r[s]), softened_pe(masses, pos_r[s], eps)) for s in range(n_s)]
    e_a = np.array([ke + pe for ke, pe in ke_pe_a])
    e_r = np.array([ke + pe for ke, pe in ke_pe_r])
    q_a = [-ke / pe for ke, pe in ke_pe_a]
    q_r = [-ke / pe for ke, pe in ke_pe_r]
    L_a = np.array([ang_mom(masses, pos_a[s], vel_a[s]) for s in range(n_s)])
    L_r = np.array([ang_mom(masses, pos_r[s], vel_r[s]) for s in range(n_s)])
    P_a = np.array([lin_mom(masses, vel_a[s]) for s in range(n_s)])
    P_r = np.array([lin_mom(masses, vel_r[s]) for s in range(n_s)])
    C_a = np.array([com(masses, pos_a[s]) for s in range(n_s)])
    C_r = np.array([com(masses, pos_r[s]) for s in range(n_s)])

    # Energy gate model: round-off walk (13 eps sqrt(N_steps)) + comparator summation
    # floor (eps sqrt(n_pair)), headroom 10. Protocol notebook section Hypothesis.
    def tol_e(n_steps: int) -> float | None:
        if n_steps <= 0:
            return None
        walk = 13.0 * EPS64 * math.sqrt(n_steps)
        return 10.0 * (walk + EPS64 * math.sqrt(n_pair))

    de_a = float(np.max(np.abs(e_a - e_a[0]) / abs(e_a[0])))
    de_r = float(np.max(np.abs(e_r - e_r[0]) / abs(e_r[0])))
    de_x = float(np.max(np.abs(e_a - e_r) / abs(e_a[0])))
    _tol_a = tol_e(n_steps_a)
    _tol_r = tol_e(n_steps_r)
    _tol_cross = (
        None if (_tol_a is None or _tol_r is None)
        else math.sqrt(2) * max(_tol_a, _tol_r)
    )

    sl_a = max(scale_l(masses, pos_a[s], vel_a[s]) for s in range(n_s))
    sl_r = max(scale_l(masses, pos_r[s], vel_r[s]) for s in range(n_s))
    sp_a = max(scale_p(masses, vel_a[s]) for s in range(n_s))
    sp_r = max(scale_p(masses, vel_r[s]) for s in range(n_s))

    dL_a = float(np.max(np.linalg.norm(L_a - L_a[0], axis=1)))
    dL_r = float(np.max(np.linalg.norm(L_r - L_r[0], axis=1)))
    dL_x = float(np.max(np.linalg.norm(L_a - L_r, axis=1)))
    dP_a = float(np.max(np.linalg.norm(P_a - P_a[0], axis=1)))
    dP_r = float(np.max(np.linalg.norm(P_r - P_r[0], axis=1)))
    dP_x = float(np.max(np.linalg.norm(P_a - P_r, axis=1)))
    dC_a = float(np.max(np.linalg.norm(C_a - C_a[0], axis=1)))
    dC_r = float(np.max(np.linalg.norm(C_r - C_r[0], axis=1)))
    dC_x = float(np.max(np.linalg.norm(C_a - C_r, axis=1)))

    t_final = float(t_a[-1])
    # floors fall back to sqrt(1) when stats are missing; tolerances go null instead
    sqrt_a = math.sqrt(max(n_steps_a, 1))
    sqrt_r = math.sqrt(max(n_steps_r, 1))
    floors = {
        "L": {
            "scale_only": {
                "apsis": 10 * EPS64 * sl_a,
                "rebound": 10 * EPS64 * sl_r,
                "cross": 10 * math.sqrt(2) * EPS64 * max(sl_a, sl_r),
            },
            "scale_sqrt_steps": {
                "apsis": 10 * EPS64 * sl_a * sqrt_a,
                "rebound": 10 * EPS64 * sl_r * sqrt_r,
                "cross": 10 * math.sqrt(2) * EPS64 * max(sl_a * sqrt_a, sl_r * sqrt_r),
            },
        },
        "P": {
            "scale_only": {
                "apsis": 10 * EPS64 * sp_a,
                "rebound": 10 * EPS64 * sp_r,
                "cross": 10 * math.sqrt(2) * EPS64 * max(sp_a, sp_r),
            },
            "scale_sqrt_steps": {
                "apsis": 10 * EPS64 * sp_a * sqrt_a,
                "rebound": 10 * EPS64 * sp_r * sqrt_r,
                "cross": 10 * math.sqrt(2) * EPS64 * max(sp_a * sqrt_a, sp_r * sqrt_r),
            },
        },
        "COM": {
            "drift_model": {
                "apsis": 1.5 * EPS64 * sp_a * t_final,
                "rebound": 1.5 * EPS64 * sp_r * t_final,
                "cross": 3.0 * EPS64 * max(sp_a, sp_r) * t_final,
            },
        },
    }

    half = n_s // 2
    dr = np.linalg.norm(pos_a - pos_r, axis=2)
    max_dr = float(np.max(dr))

    report = {
        "informational": args.informational,
        "n": n,
        "eps": eps,
        "n_samples": n_s,
        "t_final": t_final,
        "n_steps": {"apsis": n_steps_a, "rebound": n_steps_r},
        "tier1_energy": {
            "de_apsis": de_a,
            "tol_apsis": _tol_a,
            "de_rebound": de_r,
            "tol_rebound": _tol_r,
            "de_cross": de_x,
            "tol_cross": _tol_cross,
        },
        "tier1_L": {"apsis": dL_a, "rebound": dL_r, "cross": dL_x},
        "tier2_P": {"apsis": dP_a, "rebound": dP_r, "cross": dP_x},
        "tier2_COM": {"apsis": dC_a, "rebound": dC_r, "cross": dC_x},
        "candidate_floors": floors,
        "tier3": {
            "virial_q": {
                "apsis_t0": q_a[0],
                "apsis_mean_2nd_half": float(np.mean(q_a[half:])),
                "apsis_std_2nd_half": float(np.std(q_a[half:])),
                "rebound_mean_2nd_half": float(np.mean(q_r[half:])),
                "rebound_std_2nd_half": float(np.std(q_r[half:])),
            },
            "lagrangian_radii": {
                "apsis_t0": lagrangian_radii(masses, pos_a[0]),
                "apsis_tfinal": lagrangian_radii(masses, pos_a[-1]),
                "rebound_t0": lagrangian_radii(masses, pos_r[0]),
                "rebound_tfinal": lagrangian_radii(masses, pos_r[-1]),
            },
            "max_dr_cross_impl": max_dr,
        },
    }

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "comparison.json").write_text(
        json.dumps(report, indent=2, allow_nan=False), encoding="utf-8"
    )
    print(json.dumps(report, indent=2, allow_nan=False))

    if args.informational:
        print("\nphase 0 -- informational run: metrics reported, nothing gated")
        return 0
    print(
        "\nERROR: gated mode is blocked until the L/P floor model is frozen "
        "in the protocol notebook (phase 0 close-out)",
        file=sys.stderr,
    )
    return 2


def read_steps(path: Path, key: str) -> int:
    if not path.exists():
        return -1
    return int(json.loads(path.read_text(encoding="utf-8")).get(key, -1))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plummer cluster parity comparator (3D, softened PE)."
    )
    parser.add_argument("--selftest", action="store_true", help="run embedded self-tests and exit")
    parser.add_argument("--ics", default="ics_n256.csv", help="committed IC CSV (masses + eps)")
    parser.add_argument("--apsis-csv", default="out/apsis.csv", help="apsis-side snapshots")
    parser.add_argument("--rebound-csv", default="out/rebound.csv", help="REBOUND-side snapshots")
    parser.add_argument(
        "--apsis-stats", default="out/apsis_stats.json", help="apsis step counts"
    )
    parser.add_argument(
        "--rebound-stats", default="out/rebound_stats.json", help="REBOUND step counts"
    )
    parser.add_argument("--output-dir", default="out", help="directory for comparison.json")
    parser.add_argument(
        "--informational",
        action="store_true",
        help="phase-0 mode: report, never gate",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
