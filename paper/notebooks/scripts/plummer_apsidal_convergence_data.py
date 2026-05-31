"""Generate the §3.2 convergence-figure data: geometric softened-Plummer apsidal
precession deviation from the exact full-potential oracle, vs softening ε.

Two series, both per RADIAL period — the oracle's own convention, so no
rate/period division and no osculating elements enter:

  - closed form (leading O(ε²) secular term)  — computed here
  - apsis IAS15 geometric measurement          — READ from the Rust sweep CSV

against the exact apsidal-angle quadrature oracle. Writes |1 − X/oracle|.

NO pinned literals: apsis comes from `cargo run --example softened_plummer_sweep`
(geometric apsidal angle per radial period, periapsis-passage detection); the
oracle and closed form are evaluated from the quadrature module.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

# Sibling import from this script's own directory (robust to CWD).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plummer_apsidal_quadrature import (
    closed_form_per_orbit_rad,
    precession_per_orbit_rad,
)

DATA = Path(__file__).resolve().parents[2] / "figures" / "data"
APSIS_CSV = DATA / "apsis_softened_sweep.csv"
OUT = DATA / "plummer_apsidal_convergence.csv"


def read_apsis() -> dict[float, float]:
    """eps -> apsis geometric apsidal precession per radial period (rad)."""
    out: dict[float, float] = {}
    with APSIS_CSV.open(encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#") or row[0] == "eps":
                continue
            out[float(row[0])] = float(row[1])
    return out


def main() -> None:
    apsis = read_apsis()

    rows = []
    for eps in sorted(apsis):
        oracle = precession_per_orbit_rad(eps)
        closed = closed_form_per_orbit_rad(eps)
        ap = apsis[eps]
        dev_closed = abs(closed / oracle - 1.0)
        dev_apsis = abs(ap / oracle - 1.0)
        rows.append((eps, oracle, closed, ap, dev_closed, dev_apsis))
        print(f"  eps={eps:.4e}  dev_closed={dev_closed:.3e}  dev_apsis={dev_apsis:.3e}")

    with OUT.open("w", newline="", encoding="utf-8") as fh:
        fh.write("# Softened-Plummer apsidal precession: deviation from exact quadrature oracle\n")
        fh.write("# all per radial period (oracle convention); no pinned literals\n")
        fh.write(f"# apsis: {APSIS_CSV.name} (IAS15, geometric apsidal angle per radial period)\n")
        fh.write("# closed form + exact oracle from plummer_apsidal_quadrature.py\n")
        w = csv.writer(fh)
        w.writerow(["eps", "oracle", "closed", "apsis", "dev_closed", "dev_apsis"])
        w.writerows(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
