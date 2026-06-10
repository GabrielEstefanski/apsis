"""Single-pair softening-convention check (protocol section Softening).

Both sides integrate one softened pair (unit source mass, probe at r = 1)
from rest for t = 1e-6 and report a_x = dvx/t; each is gated at 1e-9
relative against the closed form -1/(1+eps^2)^{3/2}. A convention mismatch
(eps placed differently) shifts the value by O(eps^2) ~ 1e-2 -- seven
decades above the gate. Measured agreement at drafting time: ~3e-13.

Run:
    python smoke_pair.py --ics ics_n256.csv

Protocol notebook:
    paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import rebound

GATE = 1.0e-9
WORKSPACE_ROOT = Path(__file__).resolve().parents[3]


def read_eps(ics: Path) -> float:
    with ics.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("# eps="):
                return float(line.split("=", 1)[1])
    raise ValueError(f"no '# eps=' header in {ics}")


def rebound_side(eps: float) -> float:
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.softening = eps
    sim.integrator = "ias15"
    sim.dt = 1.0e-7
    sim.add(m=1.0, x=0.0, y=0.0, z=0.0)
    sim.add(m=1.0e-12, x=1.0, y=0.0, z=0.0)
    sim.integrate(1.0e-6)
    return sim.particles[1].vx / sim.t


def apsis_side(eps: float) -> float:
    try:
        out = subprocess.run(
            [
                "cargo", "run", "--release", "--example", "rebound_parity_plummer_cluster",
                "-p", "apsis", "--", "--smoke", "--eps", repr(eps),
            ],
            cwd=WORKSPACE_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(e.stderr, file=sys.stderr)
        raise
    return float(out.stdout.strip().splitlines()[-1])


def main() -> int:
    args = parse_args()
    eps = read_eps(Path(args.ics))
    expected = -1.0 / (1.0 + eps * eps) ** 1.5

    ok = True
    for name, fn in (("apsis", apsis_side), ("rebound", rebound_side)):
        a = fn(eps)
        rel = abs(a / expected - 1.0)
        verdict = "pass" if rel < GATE else "FAIL"
        if rel >= GATE:
            ok = False
        print(f"smoke {name:8s} a_x={a:+.17e}  expected={expected:+.17e}  rel={rel:.3e}  {verdict}")
    return 0 if ok else 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Single-pair softening smoke test.")
    parser.add_argument("--ics", default="ics_n256.csv", help="IC CSV providing the protocol eps")
    return parser.parse_args()


if __name__ == "__main__":
    sys.exit(main())
