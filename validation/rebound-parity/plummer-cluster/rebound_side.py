"""REBOUND parity -- Plummer cluster, REBOUND IAS15 side.

Mirror of `crates/apsis/examples/rebound_parity_plummer_cluster.rs`: same
committed ICs, same softening (parsed from the IC header), sampled at the
*actual* times the apsis side landed at (`exact_finish_time = 1`), emitted
in the same long format. Step count goes to a JSON for the comparator's
energy-gate model.

Run:
    python rebound_side.py --ics ics_n256.csv

Protocol notebook:
    paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import rebound

DT_INITIAL = 1.0e-3


def read_ics(path: Path) -> tuple[list[tuple[float, ...]], float]:
    eps = None
    rows: list[tuple[float, ...]] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("# eps="):
                eps = float(line.split("=", 1)[1])
                continue
            if line.startswith("#") or line.startswith("body,") or not line.strip():
                continue
            parts = line.rstrip().split(",")
            rows.append(tuple(float(v) for v in parts[1:]))
    if eps is None:
        raise ValueError(f"no '# eps=' header in {path}")
    return rows, eps


def read_apsis_sample_times(path: Path) -> list[float]:
    if not path.exists():
        raise FileNotFoundError(
            f"apsis CSV not found at {path}. Run the apsis side first via "
            f"`cargo run --release --example rebound_parity_plummer_cluster -p apsis`."
        )
    times: list[float] = []
    seen: set[int] = set()
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("sample,"):
                continue
            parts = line.split(",", 2)
            s = int(parts[0])
            if s not in seen:
                seen.add(s)
                times.append(float(parts[1]))
    return times


def main() -> int:
    args = parse_args()
    ics, eps = read_ics(Path(args.ics))
    times = read_apsis_sample_times(Path(args.apsis_csv))
    if not times or times[0] != 0.0:
        print("ERROR: apsis CSV missing or first sample time != 0", file=sys.stderr)
        return 1

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.dt = DT_INITIAL
    sim.softening = eps
    sim.exact_finish_time = 1
    for m, x, y, z, vx, vy, vz in ics:
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    min_dt_observed = float("inf")
    with out.open("w", newline="", encoding="utf-8") as f:
        f.write("# REBOUND parity -- Plummer cluster -- REBOUND IAS15 side\n")
        f.write("# protocol: paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md\n")
        f.write(f"# rebound_version: {rebound.__version__}\n")
        f.write(f"# n={len(ics)}, eps={eps:.18e}, dt0={DT_INITIAL:.18e}\n")
        f.write("sample,t,body,x,y,z,vx,vy,vz\n")
        for s, t_target in enumerate(times):
            if t_target > 0.0:
                sim.integrate(t_target)
            if s > 0 and sim.dt > 0.0:
                min_dt_observed = min(min_dt_observed, abs(sim.dt))
            for i, p in enumerate(sim.particles):
                f.write(
                    f"{s},{sim.t:.18e},{i},{p.x:.18e},{p.y:.18e},{p.z:.18e},"
                    f"{p.vx:.18e},{p.vy:.18e},{p.vz:.18e}\n"
                )

    steps = int(getattr(sim, "steps_done", -1))
    Path(args.stats_output).write_text(json.dumps({"steps_done": steps}), encoding="utf-8")
    print(f"wrote {len(times)} samples to {out}", flush=True)
    print(f"[diag] rebound steps_done: {steps}", file=sys.stderr)
    print(f"[diag] rebound min dt observed: {min_dt_observed:.6e}", file=sys.stderr)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="REBOUND IAS15 side of the Plummer cluster parity run."
    )
    parser.add_argument("--ics", default="ics_n256.csv", help="committed IC CSV")
    parser.add_argument("--apsis-csv", default="out/apsis.csv", help="apsis CSV for sample times")
    parser.add_argument("--output", default="out/rebound.csv", help="REBOUND-side CSV output")
    parser.add_argument("--stats-output", default="out/rebound_stats.json", help="step-count JSON")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
