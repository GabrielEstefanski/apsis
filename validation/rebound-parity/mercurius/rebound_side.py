"""REBOUND parity — Mercurius, REBOUND MERCURIUS side.

Mirror of `crates/apsis/examples/rebound_parity_mercurius.rs`. Runs the
canonical Sun + 4 outer planets + 1 Jupiter-crossing test particle
scenario under REBOUND's MERCURIUS integrator for 10^4 years, sampling
state and conservation diagnostics at the same target times the apsis
side actually landed at.

For parity, REBOUND lands at apsis's *actual* sample times via
`exact_finish_time = 1` (apsis is fixed-step, so its sample times are
exactly k * 1 yr, but the orchestration is structured the same way as
the existing Kepler / Pythagorean / figure-8 scenarios so adaptive
integrators on the apsis side are handled identically without the
harness needing a special branch).

Run:
    python rebound_side.py
    python rebound_side.py --apsis-csv ./out/apsis.csv --output ./out/rebound.csv

Protocol notebook:
    paper/notebooks/2026-05-13-rebound-parity-mercurius.md

Constants below mirror those in the Rust example
(`crates/apsis/examples/rebound_parity_mercurius.rs`). A change here is
a protocol change — update the notebook in lockstep.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import rebound

# ── Protocol constants ──────────────────────────────────────────────────── #

M_SUN: float = 1.0

M_JUPITER: float = 9.55e-4
M_SATURN: float = 2.86e-4
M_URANUS: float = 4.37e-5
M_NEPTUNE: float = 5.15e-5
M_TEST: float = 1.0e-9

A_JUPITER: float = 5.20
A_SATURN: float = 9.58
A_URANUS: float = 19.18
A_NEPTUNE: float = 30.07

A_TEST: float = 4.20
E_TEST: float = 0.40
I_TEST: float = 0.05

N_YEARS: int = 10_000
DT_YEARS: float = 0.01
ALPHA_HILL: float = 3.0

# G in solar AU-year units. Same SI constants and same f64 association
# as apsis's `UnitSystem::solar().g()`:
#   G_code = G_SI · mass_kg · time_s · time_s / (length_m · length_m · length_m)
# Match the multiply-then-divide order so the result agrees to f64
# precision with the apsis side.
G_SI: float = 6.67430e-11
AU_M: float = 1.495978707e11
YR_S: float = 3.15576e7
GM_SUN_SI: float = 1.3271244e20  # IAU 2015 nominal; apsis units.rs GM_SUN_SI
MSUN_KG: float = GM_SUN_SI / G_SI
G_SOLAR: float = G_SI * MSUN_KG * YR_S * YR_S / (AU_M * AU_M * AU_M)


def main() -> int:
    args = parse_args()

    apsis_csv = Path(args.apsis_csv).resolve()
    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Discover apsis's actual sample times ────────────────────────────── #
    sample_times = read_apsis_sample_times(apsis_csv)
    if len(sample_times) != N_YEARS + 1:
        print(
            f"ERROR: apsis CSV has {len(sample_times)} samples; "
            f"expected {N_YEARS + 1}",
            file=sys.stderr,
        )
        return 1
    if sample_times[0] != 0.0:
        print(
            f"ERROR: first apsis sample time should be 0; got {sample_times[0]}",
            file=sys.stderr,
        )
        return 1

    # ── Build REBOUND simulation ────────────────────────────────────────── #
    sim = rebound.Simulation()
    sim.G = G_SOLAR
    sim.integrator = "mercurius"
    sim.dt = DT_YEARS
    sim.ri_mercurius.r_crit_hill = ALPHA_HILL
    # Force REBOUND to land exactly at the requested t. Mercurius is a
    # fixed-step outer integrator (its inner IAS15 sub-integration
    # happens within the outer dt window), but exact_finish_time keeps
    # the harness orchestration identical to the IAS15-on-IAS15 parity
    # scenarios elsewhere.
    sim.exact_finish_time = 1

    # ── Initial conditions (mirror Rust side bit-for-bit) ───────────────── #
    #
    # Same heliocentric construction, same circular-velocity formula,
    # same vis-viva at periapsis. The only frame transformation is the
    # COM shift via `sim.move_to_com()` after building, which is what the
    # Rust side's `com_shift` helper also does (deterministically
    # reproducible across implementations).
    sim.add(
        m=M_SUN,
        x=0.0, y=0.0, z=0.0,
        vx=0.0, vy=0.0, vz=0.0,
    )
    add_circular_planet(sim, M_JUPITER, A_JUPITER, 0.0)
    add_circular_planet(sim, M_SATURN, A_SATURN, math.pi / 2.0)
    add_circular_planet(sim, M_URANUS, A_URANUS, math.pi)
    add_circular_planet(sim, M_NEPTUNE, A_NEPTUNE, 1.5 * math.pi)
    add_eccentric_inclined_planet(sim, M_TEST, A_TEST, E_TEST, I_TEST)

    sim.move_to_com()

    n_bodies = sim.N

    # ── Integrate to each target time, recording state ──────────────────── #
    rows: list[tuple] = []
    for year, t_target in enumerate(sample_times):
        if t_target > 0.0:
            sim.integrate(t_target)

        body_state: list[float] = []
        for p in sim.particles:
            body_state.extend(
                [p.x, p.y, p.z, p.vx, p.vy, p.vz]
            )
        try:
            e_total = sim.energy()  # REBOUND >= 4
        except AttributeError:
            e_total = sim.calculate_energy()  # REBOUND 3.x fallback

        # Total z-component of angular momentum: Σ m (x vy − y vx).
        # REBOUND has `sim.angular_momentum()` (returns a vec3); use the
        # arithmetic form for explicit determinism across REBOUND
        # versions.
        lz_total = sum(
            p.m * (p.x * p.vy - p.y * p.vx) for p in sim.particles
        )

        rows.append((year, sim.t, body_state, e_total, lz_total))

    # ── Write CSV with the apsis-side schema ────────────────────────────── #
    with output_path.open("w", newline="") as f:
        f.write("# REBOUND parity — Mercurius — REBOUND side\n")
        f.write("# protocol: paper/notebooks/2026-05-13-rebound-parity-mercurius.md\n")
        f.write(f"# integrator: MERCURIUS (REBOUND), r_crit_hill={ALPHA_HILL}\n")
        f.write(f"# rebound_version: {rebound.__version__}\n")
        f.write(f"# units: solar AU-year (G = {G_SOLAR:.18e})\n")
        f.write(
            f"# n_years={N_YEARS}, dt_years={DT_YEARS}, n_bodies={n_bodies}\n"
        )
        f.write(_header(n_bodies) + "\n")
        for year, t, body_state, e_total, lz_total in rows:
            row_str = f"{year},{t:.18e}"
            for v in body_state:
                row_str += f",{v:.18e}"
            row_str += f",{e_total:.18e},{lz_total:.18e}"
            f.write(row_str + "\n")

    print(f"wrote {len(rows)} samples to {output_path}", flush=True)
    return 0


# ── Initial-condition helpers (mirror Rust side) ────────────────────────── #


def add_circular_planet(
    sim: "rebound.Simulation",
    mass: float,
    a: float,
    nu: float,
) -> None:
    """Place a circular planet at heliocentric distance `a`, true anomaly
    `nu`, in the (x, y) plane around the Sun at the origin. Mirrors
    `circular_planet` in the Rust example."""
    v_c = math.sqrt(G_SOLAR * M_SUN / a)
    x = a * math.cos(nu)
    y = a * math.sin(nu)
    vx = -v_c * math.sin(nu)
    vy = v_c * math.cos(nu)
    sim.add(m=mass, x=x, y=y, z=0.0, vx=vx, vy=vy, vz=0.0)


def add_eccentric_inclined_planet(
    sim: "rebound.Simulation",
    mass: float,
    a: float,
    e: float,
    i: float,
) -> None:
    """Eccentric inclined planet starting at periapsis. Position in the
    orbit plane is `(r_peri, 0)`; tangent velocity is `v_peri`. Orbit
    plane rotated by inclination `i` around the x-axis: position stays
    `(r_peri, 0, 0)`; velocity becomes `(0, v_peri·cos(i),
    v_peri·sin(i))`. Mirrors `eccentric_inclined_planet` in the Rust
    example."""
    r_peri = a * (1.0 - e)
    v_peri = math.sqrt(G_SOLAR * M_SUN / a * (1.0 + e) / (1.0 - e))
    sim.add(
        m=mass,
        x=r_peri, y=0.0, z=0.0,
        vx=0.0, vy=v_peri * math.cos(i), vz=v_peri * math.sin(i),
    )


# ── CSV I/O ─────────────────────────────────────────────────────────────── #


def _header(n_bodies: int) -> str:
    parts = ["year", "t"]
    for i in range(n_bodies):
        parts += [f"x{i}", f"y{i}", f"z{i}", f"vx{i}", f"vy{i}", f"vz{i}"]
    parts += ["e_total", "lz_total"]
    return ",".join(parts)


def read_apsis_sample_times(path: Path) -> list[float]:
    """Read the `t` column of the apsis CSV, skipping comment lines."""
    if not path.exists():
        raise FileNotFoundError(
            f"apsis CSV not found at {path}. Run the apsis side first via "
            f"`cargo run --release --example rebound_parity_mercurius -p apsis`."
        )
    times: list[float] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            times.append(float(row["t"]))
    return times


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="REBOUND MERCURIUS side of the Mercurius parity comparison.",
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
