"""REBOUNDx parity — Sun-Mercury 1PN, REBOUND + REBOUNDx `gr` side.

Mirror of `crates/apsis-1pn/examples/reboundx_parity_gr.rs`. Runs Sun-Mercury
under REBOUND IAS15 with the REBOUNDx `gr` effect (single-dominant-mass 1PN),
sampling state and total (Newtonian) energy at the *actual* sample times the
apsis side landed at — same protocol as the rebound-parity scenarios, so
"two implementations sampled at slightly different physical times" is not a
spurious source of |Delta r|.

The `c` value and initial conditions are read/derived to match the apsis side
bit-for-bit: `c` is parsed from the apsis CSV header; the ICs use the same f64
expressions as the Rust example.

Run (in a Linux env — REBOUNDx does not build on Windows/MSVC):
    pip install rebound reboundx
    python reboundx_side.py                       # gr on
    python reboundx_side.py --no-1pn \
        --apsis-csv out/apsis_kepler.csv \
        --output out/rebound_kepler.csv           # gr off (control)

Protocol notebook: paper/notebooks/2026-05-29-reboundx-parity-gr.md
"""

# rebound/reboundx are untyped C extensions (no py.typed/stubs); silence the
# unavoidable strict-mode noise for this file only.
# pyright: reportMissingImports=false, reportMissingTypeStubs=false, reportUnknownMemberType=false, reportAttributeAccessIssue=false, reportUnknownVariableType=false, reportUnknownArgumentType=false

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import rebound
import reboundx

# ── Protocol constants (mirror the Rust example exactly) ────────────────── #

A_MERCURY: float = 0.387098
E_MERCURY: float = 0.20563
M_SUN: float = 1.0
M_MERCURY: float = 1.660114e-7
N_ORBITS: int = 500


def main() -> int:
    args = parse_args()
    apsis_csv = Path(args.apsis_csv).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sample_times = read_apsis_sample_times(apsis_csv)
    c_value = read_apsis_c(apsis_csv)
    if len(sample_times) != N_ORBITS + 1:
        print(f"ERROR: apsis CSV has {len(sample_times)} samples; expected {N_ORBITS + 1}",
              file=sys.stderr)
        return 1

    # ── Initial conditions (COM frame, mirror the Rust f64 expressions) ───── #
    r_peri = A_MERCURY * (1.0 - E_MERCURY)
    v_peri = ((1.0 + E_MERCURY) / (A_MERCURY * (1.0 - E_MERCURY))) ** 0.5
    m_total = M_SUN + M_MERCURY
    sun_x = -(M_MERCURY / m_total) * r_peri
    sun_vy = -(M_MERCURY / m_total) * v_peri
    mercury_x = (M_SUN / m_total) * r_peri
    mercury_vy = (M_SUN / m_total) * v_peri

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.dt = sample_times[1] * 1.0e-3 if len(sample_times) > 1 else 1.0e-3
    sim.exact_finish_time = 1
    sim.add(m=M_SUN, x=sun_x, y=0.0, z=0.0, vx=0.0, vy=sun_vy, vz=0.0)
    sim.add(m=M_MERCURY, x=mercury_x, y=0.0, z=0.0, vx=0.0, vy=mercury_vy, vz=0.0)

    if not args.no_1pn:
        rebx = reboundx.Extras(sim)
        gr = rebx.load_force("gr")
        gr.params["c"] = c_value
        rebx.add_force(gr)

    rows = []
    for orbit, t_target in enumerate(sample_times):
        if t_target > 0.0:
            sim.integrate(t_target)
        b0, b1 = sim.particles[0], sim.particles[1]
        try:
            e_total = sim.energy()
        except AttributeError:
            e_total = sim.calculate_energy()
        rows.append((orbit, sim.t, b0.x, b0.y, b0.vx, b0.vy, b1.x, b1.y, b1.vx, b1.vy, e_total))

    mode = "Newtonian (gr off, control)" if args.no_1pn else "REBOUNDx gr (single dominant mass)"
    with output_path.open("w", newline="", encoding="utf-8") as f:
        f.write("# REBOUNDx parity — Sun-Mercury 1PN — REBOUND side\n")
        f.write("# protocol: paper/notebooks/2026-05-29-reboundx-parity-gr.md\n")
        f.write(f"# integrator: IAS15 (REBOUND); force: {mode}\n")
        f.write(f"# rebound_version: {rebound.__version__}, reboundx_version: {reboundx.__version__}\n")
        f.write(f"# c={c_value:.18e}\n")
        f.write("orbit,t,x0,y0,vx0,vy0,x1,y1,vx1,vy1,e_total\n")
        for row in rows:
            f.write(str(row[0]) + "".join(f",{v:.18e}" for v in row[1:]) + "\n")

    print(f"wrote {len(rows)} samples to {output_path} (no_1pn={args.no_1pn})", flush=True)
    return 0


def read_apsis_sample_times(path: Path) -> list[float]:
    if not path.exists():
        raise FileNotFoundError(
            f"apsis CSV not found at {path}. Run the apsis side first: "
            f"`cargo run --release --example reboundx_parity_gr -p apsis-1pn`."
        )
    times: list[float] = []
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            times.append(float(row["t"]))
    return times


def read_apsis_c(path: Path) -> float:
    """Parse the `# c=...` header line so c matches the apsis side bit-for-bit."""
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("# c="):
                return float(line.split("=", 1)[1].strip())
            if not line.startswith("#"):
                break
    raise ValueError(f"no '# c=' header line in {path}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="REBOUND + REBOUNDx gr side of the Sun-Mercury 1PN parity.")
    p.add_argument("--apsis-csv", default="out/apsis.csv")
    p.add_argument("--output", default="out/rebound.csv")
    p.add_argument("--no-1pn", action="store_true", help="omit the gr force (Newtonian control)")
    return p.parse_args()


if __name__ == "__main__":
    sys.exit(main())
