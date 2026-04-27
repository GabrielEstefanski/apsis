"""REBOUND parity — Figure-8 choreography, REBOUND IAS15 side.

Mirror of `crates/apsis/examples/rebound_parity_figure8.rs`. Runs the
canonical Chenciner–Montgomery figure-8 three-body orbit under REBOUND's
IAS15 for the same number of periods and at the same dense sampling
cadence as the apsis side, evaluating REBOUND's state at the *actual*
sample times the apsis side landed at.

For parity, REBOUND lands at apsis's actual sample times rather than the
nominal `n · T / 200` targets — apsis's adaptive controller may
overshoot the nominal target by one substep, and comparing at apsis's
actual `t` removes "two implementations sampled at slightly different
physical times" as a spurious source of disagreement. REBOUND's
``exact_finish_time = 1`` is used to force its own controller to land
exactly at the requested time.

Run:
    python rebound_side.py
    python rebound_side.py --apsis-csv ./out/apsis.csv --output ./out/rebound.csv

Protocol notebook:
    docs/experiments/2026-04-26-rebound-parity-figure8.md

Constants below mirror those in the Rust example. A change here is a
protocol change — update the notebook in lockstep.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import rebound

# ── Protocol constants ──────────────────────────────────────────────────── #

MASS: float = 1.0

# Chenciner & Montgomery (2000) ICs, 8-digit literature form. Same string
# literals on both sides → identical f64 bit pattern on the same hardware.
R1: tuple[float, float] = (-0.97000436, 0.24308753)
R2: tuple[float, float] = (0.97000436, -0.24308753)
R3: tuple[float, float] = (0.0, 0.0)
V1: tuple[float, float] = (0.4662036850, 0.4323657300)
V2: tuple[float, float] = (0.4662036850, 0.4323657300)
V3: tuple[float, float] = (-0.93240737, -0.86473146)

PERIOD: float = 6.3259139870
N_PERIODS_DEFAULT: int = 10
SAMPLES_PER_PERIOD: int = 200
DT_FRACTION_OF_PERIOD: float = 1.0e-3


def main() -> int:
    args = parse_args()

    apsis_csv = Path(args.apsis_csv).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Discover apsis's actual sample times ────────────────────────────── #
    sample_times = read_apsis_sample_times(apsis_csv)
    if not sample_times:
        print(f"ERROR: apsis CSV at {apsis_csv} is empty", file=sys.stderr)
        return 1
    if sample_times[0] != 0.0:
        print(
            f"ERROR: first apsis sample time should be 0; got {sample_times[0]}",
            file=sys.stderr,
        )
        return 1

    # ── Build REBOUND simulation ────────────────────────────────────────── #
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.dt = PERIOD * DT_FRACTION_OF_PERIOD
    # Force REBOUND to land exactly at the requested t rather than overshooting.
    sim.exact_finish_time = 1

    sim.add(m=MASS, x=R1[0], y=R1[1], z=0.0, vx=V1[0], vy=V1[1], vz=0.0)
    sim.add(m=MASS, x=R2[0], y=R2[1], z=0.0, vx=V2[0], vy=V2[1], vz=0.0)
    sim.add(m=MASS, x=R3[0], y=R3[1], z=0.0, vx=V3[0], vy=V3[1], vz=0.0)

    # ── Integrate to each target time, recording state ──────────────────── #
    rows: list[tuple] = []
    for sample_idx, t_target in enumerate(sample_times):
        if t_target > 0.0:
            sim.integrate(t_target)
        # else: sample 0, just record initial state.

        b0 = sim.particles[0]
        b1 = sim.particles[1]
        b2 = sim.particles[2]
        try:
            e_total = sim.energy()  # REBOUND ≥ 4
        except AttributeError:
            e_total = sim.calculate_energy()  # REBOUND 3.x fallback
        rows.append((
            sample_idx, sim.t,
            b0.x, b0.y, b0.vx, b0.vy,
            b1.x, b1.y, b1.vx, b1.vy,
            b2.x, b2.y, b2.vx, b2.vy,
            e_total,
        ))

    # ── Write CSV with the apsis-side schema ────────────────────────────── #
    n_periods = (len(sample_times) - 1) // SAMPLES_PER_PERIOD
    with output_path.open("w", newline="") as f:
        f.write("# REBOUND parity — Figure-8 choreography — REBOUND IAS15 side\n")
        f.write("# protocol: docs/experiments/2026-04-26-rebound-parity-figure8.md\n")
        f.write("# integrator: IAS15 (REBOUND)\n")
        f.write(f"# rebound_version: {rebound.__version__}\n")
        f.write("# units: canonical (G = 1)\n")
        f.write(f"# mass={MASS}, period={PERIOD:.18e}\n")
        f.write(
            f"# n_periods={n_periods}, samples_per_period={SAMPLES_PER_PERIOD}, "
            f"dt0={sim.dt:.18e}\n"
        )
        f.write("sample,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,x2,y2,vx2,vy2,e_total\n")
        for row in rows:
            sample_idx = row[0]
            floats = row[1:]
            f.write(str(sample_idx))
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
            f"`cargo run --release --example rebound_parity_figure8 -p apsis`."
        )
    times: list[float] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            times.append(float(row["t"]))
    return times


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="REBOUND IAS15 side of the Figure-8 parity comparison.",
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
