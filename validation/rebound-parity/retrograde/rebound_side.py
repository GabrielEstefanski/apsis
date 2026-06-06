"""REBOUND parity — Retrograde Kepler e=0.5, REBOUND IAS15 side.

Mirror of `crates/apsis/examples/rebound_parity_retrograde.rs`. Runs a
canonical retrograde Kepler two-body orbit at eccentricity 0.5 under
REBOUND's IAS15 for 10000 orbital periods, sampling state and total
energy at the same target times the apsis side actually landed at.

Single IC difference vs the prograde Kepler harness: the secondary's
tangential periapsis velocity is `-v_peri` instead of `+v_peri`. Every
other component, mass, length, and physical scale is held identical.

For parity, REBOUND lands at apsis's *actual* sample times rather than
the nominal `k * 2π` targets — apsis's adaptive controller may overshoot
the nominal target by one substep, and comparing at apsis's actual `t`
removes "two implementations sampled at slightly different physical
times" as a spurious source of |Δr|. REBOUND's `exact_finish_time=1`
forces its controller to land exactly at the requested time.

Run:
    python rebound_side.py
    python rebound_side.py --apsis-csv ./out/apsis.csv --output ./out/rebound.csv

Protocol notebook:
    paper/notebooks/2026-05-01-rebound-parity-retrograde.md

Constants below mirror those in the Rust example. A change here is a
protocol change — update the notebook in lockstep.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import rebound

# ── Protocol constants ──────────────────────────────────────────────────── #

A: float = 1.0
E: float = 0.5
M_PRIMARY: float = 1.0
M_SECONDARY: float = 1.0e-6
N_ORBITS: int = 10_000
DT_FRACTION_OF_PERIOD: float = 1.0e-3


def main() -> int:
    args = parse_args()

    apsis_csv = Path(args.apsis_csv).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Discover apsis's actual sample times ────────────────────────────── #
    sample_times = read_apsis_sample_times(apsis_csv)
    if len(sample_times) != N_ORBITS + 1:
        print(
            f"ERROR: apsis CSV has {len(sample_times)} samples; "
            f"expected {N_ORBITS + 1}",
            file=sys.stderr,
        )
        return 1
    if sample_times[0] != 0.0:
        print(
            f"ERROR: first apsis sample time should be 0; got {sample_times[0]}",
            file=sys.stderr,
        )
        return 1

    # ── Initial conditions (mirror Rust side bit-for-bit) ───────────────── #
    #
    # The same floating-point expression evaluated in Python's f64 should
    # produce the same bits as Rust's f64 on the same hardware. The single
    # sign-flip (`-v_peri` instead of `+v_peri`) is an exact bit-level
    # operation in IEEE-754 (toggles the sign bit), so it does not introduce
    # any new IC-rounding source vs the prograde test.
    r_peri = A * (1.0 - E)
    v_peri = math.sqrt((1.0 + E) / (A * (1.0 - E)))

    m_total = M_PRIMARY + M_SECONDARY
    primary_x = -(M_SECONDARY / m_total) * r_peri
    primary_vy = (M_SECONDARY / m_total) * v_peri  # sign flipped vs prograde
    secondary_x = (M_PRIMARY / m_total) * r_peri
    secondary_vy = -(M_PRIMARY / m_total) * v_peri  # sign flipped vs prograde

    # ── Build REBOUND simulation ────────────────────────────────────────── #
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    period = 2.0 * math.pi
    sim.dt = period * DT_FRACTION_OF_PERIOD
    sim.exact_finish_time = 1

    sim.add(
        m=M_PRIMARY,
        x=primary_x, y=0.0, z=0.0,
        vx=0.0, vy=primary_vy, vz=0.0,
    )
    sim.add(
        m=M_SECONDARY,
        x=secondary_x, y=0.0, z=0.0,
        vx=0.0, vy=secondary_vy, vz=0.0,
    )

    # ── Integrate to each target time, recording state ──────────────────── #
    rows: list[tuple] = []
    for orbit, t_target in enumerate(sample_times):
        if t_target > 0.0:
            sim.integrate(t_target)
        # else: orbit 0, just record initial state.

        b0 = sim.particles[0]
        b1 = sim.particles[1]
        try:
            e_total = sim.energy()  # REBOUND ≥ 4
        except AttributeError:
            e_total = sim.calculate_energy()  # REBOUND 3.x fallback
        rows.append((
            orbit, sim.t,
            b0.x, b0.y, b0.vx, b0.vy,
            b1.x, b1.y, b1.vx, b1.vy,
            e_total,
        ))

        if orbit > 0 and orbit % 1000 == 0:
            print(f"  rebound: {orbit}/{N_ORBITS} orbits completed (t={sim.t:.3e})", flush=True)

    # ── Write CSV with the apsis-side schema ────────────────────────────── #
    with output_path.open("w", newline="") as f:
        f.write("# REBOUND parity — Retrograde Kepler e=0.5 — REBOUND IAS15 side\n")
        f.write("# protocol: paper/notebooks/2026-05-01-rebound-parity-retrograde.md\n")
        f.write("# integrator: IAS15 (REBOUND)\n")
        f.write(f"# rebound_version: {rebound.__version__}\n")
        f.write("# units: canonical (G = 1)\n")
        f.write(
            f"# a={A}, e={E}, m_primary={M_PRIMARY}, m_secondary={M_SECONDARY:e}\n"
        )
        f.write("# orbit_sense: retrograde (secondary_vy = -v_peri)\n")
        f.write(
            f"# period={period:.18e}, dt0={sim.dt:.18e}, n_orbits={N_ORBITS}\n"
        )
        f.write("orbit,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,e_total\n")
        for row in rows:
            orbit_idx = row[0]
            floats = row[1:]
            f.write(str(orbit_idx))
            for v in floats:
                f.write(f",{v:.18e}")
            f.write("\n")

    print(f"wrote {len(rows)} samples to {output_path}", flush=True)
    return 0


def read_apsis_sample_times(path: Path) -> list[float]:
    """Read the `t` column of the apsis CSV, skipping comment lines."""
    if not path.exists():
        raise FileNotFoundError(
            f"apsis CSV not found at {path}. Run the apsis side first via "
            f"`cargo run --release --example rebound_parity_retrograde -p apsis`."
        )
    times: list[float] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            times.append(float(row["t"]))
    return times


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="REBOUND IAS15 side of the retrograde Kepler e=0.5 parity comparison.",
    )
    parser.add_argument(
        "--apsis-csv",
        default="out/apsis.csv",
        help="Path to the apsis-side CSV (used to discover sample times). "
        "Default: out/apsis.csv (relative to cwd).",
    )
    parser.add_argument(
        "--output",
        default="out/rebound.csv",
        help="Path for the REBOUND-side CSV output. Default: out/rebound.csv.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
