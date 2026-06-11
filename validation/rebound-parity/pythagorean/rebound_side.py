"""REBOUND parity — Pythagorean three-body (Burrau 1913), REBOUND IAS15 side.

Mirror of `crates/apsis/examples/rebound_parity_pythagorean.rs`. Runs the
canonical Burrau Pythagorean problem under REBOUND's IAS15 to the same
horizon and at the same dense sampling cadence as the apsis side,
evaluating REBOUND's state at the *actual* sample times the apsis side
landed at.

For parity, REBOUND lands at apsis's actual sample times rather than the
nominal `n / 30` targets — apsis's adaptive controller may overshoot the
nominal target by one substep, and comparing at apsis's actual `t`
removes "two implementations sampled at slightly different physical
times" as a spurious source of disagreement. REBOUND's
``exact_finish_time = 1`` is used to force its own controller to land
exactly at the requested time.

Run:
    python rebound_side.py
    python rebound_side.py --apsis-csv ./out/apsis.csv --output ./out/rebound.csv

Protocol notebook:
    paper/notebooks/2026-04-30-rebound-parity-pythagorean.md

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

# Burrau (1913) opposite-side convention: side opposite mass mᵢ has length mᵢ.
MASSES: tuple[float, float, float] = (3.0, 4.0, 5.0)

# Burrau (1913) ICs. Integer-valued floats → identical f64 bit pattern on
# both sides without any string-to-double rounding choice.
R1: tuple[float, float] = (1.0, 3.0)
R2: tuple[float, float] = (-2.0, -1.0)
R3: tuple[float, float] = (1.0, -1.0)
V1: tuple[float, float] = (0.0, 0.0)
V2: tuple[float, float] = (0.0, 0.0)
V3: tuple[float, float] = (0.0, 0.0)

HORIZON_DEFAULT: float = 70.0
SAMPLES_PER_TIME_UNIT: int = 30
DT_INITIAL: float = 1.0e-3


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
    sim.dt = DT_INITIAL
    # Force REBOUND to land exactly at the requested t rather than overshooting.
    sim.exact_finish_time = 1

    sim.add(m=MASSES[0], x=R1[0], y=R1[1], z=0.0, vx=V1[0], vy=V1[1], vz=0.0)
    sim.add(m=MASSES[1], x=R2[0], y=R2[1], z=0.0, vx=V2[0], vy=V2[1], vz=0.0)
    sim.add(m=MASSES[2], x=R3[0], y=R3[1], z=0.0, vx=V3[0], vy=V3[1], vz=0.0)

    # ── Integrate to each target time, recording state ──────────────────── #
    #
    # Track the smallest accepted IAS15 dt across the integration so we can
    # report it alongside the substep count for cross-implementation
    # comparison with apsis. Both diagnostics go to stderr at end of run.
    rows: list[tuple] = []
    min_dt_observed = float("inf")
    for sample_idx, t_target in enumerate(sample_times):
        if t_target > 0.0:
            sim.integrate(t_target)
        # else: sample 0, just record initial state.

        # `sim.dt` after `integrate()` reflects the last accepted timestep.
        # During quiescent intervals this is large; during close encounters
        # the IAS15 controller drives it down. Tracking min_dt over the
        # full run gives the floor REBOUND's controller actually reached.
        if sample_idx > 0 and sim.dt > 0.0:
            min_dt_observed = min(min_dt_observed, abs(sim.dt))

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
    horizon = sample_times[-1]
    with output_path.open("w", newline="") as f:
        f.write("# REBOUND parity — Pythagorean three-body (Burrau 1913) — REBOUND IAS15 side\n")
        f.write("# protocol: paper/notebooks/2026-04-30-rebound-parity-pythagorean.md\n")
        f.write("# integrator: IAS15 (REBOUND)\n")
        f.write(f"# rebound_version: {rebound.__version__}\n")
        f.write("# units: canonical (G = 1)\n")
        f.write(f"# masses=({MASSES[0]},{MASSES[1]},{MASSES[2]}), horizon={horizon:.18e}\n")
        f.write(
            f"# samples_per_t_unit={SAMPLES_PER_TIME_UNIT}, "
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

    # Diagnostic counters surfaced for cross-implementation comparison
    # with the apsis side. `steps_done` counts accepted IAS15 substeps;
    # `min_dt_observed` is the smallest accepted timestep across the
    # full integration. Stderr-only, debug telemetry — not part of the
    # CSV schema or the gated comparison.
    try:
        steps_done = sim.steps_done
    except AttributeError:
        steps_done = -1  # older REBOUND versions lacked this field
    print(
        f"[diag] rebound substeps total: {steps_done}",
        f"\n[diag] rebound min dt observed: {min_dt_observed:.6e}",
        file=sys.stderr,
    )
    return 0


def read_apsis_sample_times(path: Path) -> list[float]:
    """Read the `t` column of the apsis CSV, skipping comment lines."""
    if not path.exists():
        raise FileNotFoundError(
            f"apsis CSV not found at {path}. Run the apsis side first via "
            f"`cargo run --release --example rebound_parity_pythagorean -p apsis`."
        )
    times: list[float] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            times.append(float(row["t"]))
    return times


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="REBOUND IAS15 side of the Pythagorean parity comparison.",
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
