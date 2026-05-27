"""Generate the Mercury 1PN long-horizon perihelion precession figure.

Reads `paper/figures/data/mercury_1pn_long_horizon_ias15.csv` (a frozen
snapshot from `validation/mercury-1pn-long-horizon/out/ias15.csv`),
computes the unwrapped perihelion-orientation drift `Δω(t)`, and plots
it against the closed-form Schwarzschild test-particle prediction
`Δω = 6π / (c² · a · (1 − e²)) · N_orbits`.

Output: ``paper/figures/mercury_1pn_long_horizon.pdf``.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

# ── Protocol constants (mirror crates/apsis-1pn + validation harness) ──── #

A: float = 0.387098
E: float = 0.20563

C_SI: float = 299_792_458.0
AU_SI: float = 149_597_870_700.0
YEAR_S: float = 365.25 * 86_400.0
TWO_PI: float = 2.0 * math.pi
C_SOLAR_UNITS: float = C_SI * (YEAR_S / TWO_PI) / AU_SI

DELTA_OMEGA_GR_PER_ORBIT: float = 6.0 * math.pi / (
    C_SOLAR_UNITS * C_SOLAR_UNITS * A * (1.0 - E * E)
)

RAD_TO_ARCSEC: float = 180.0 / math.pi * 3600.0


# ── IO ─────────────────────────────────────────────────────────────────── #


def load_omega(path: Path) -> tuple[list[int], list[float]]:
    """Return ``(orbit_numbers, omega_osc_rad)`` from a harness CSV."""
    orbits: list[int] = []
    omegas: list[float] = []
    with path.open() as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for row in reader:
            orbits.append(int(row["orbit"]))
            omegas.append(float(row["omega_osc"]))
    return orbits, omegas


def unwrap(angles: list[float]) -> list[float]:
    """Remove ``2π`` jumps from an angle trajectory."""
    if not angles:
        return []
    out = [angles[0]]
    for a in angles[1:]:
        prev = out[-1]
        d = a - prev
        while d > math.pi:
            d -= 2.0 * math.pi
        while d < -math.pi:
            d += 2.0 * math.pi
        out.append(prev + d)
    return out


# ── Figure ─────────────────────────────────────────────────────────────── #


def render(csv_path: Path, output_path: Path) -> None:
    orbits, omega_osc = load_omega(csv_path)
    omega_unwrapped = unwrap(omega_osc)
    delta_omega_rad = [w - omega_unwrapped[0] for w in omega_unwrapped]
    predicted_rad = [DELTA_OMEGA_GR_PER_ORBIT * k for k in orbits]

    measured_arcsec = [r * RAD_TO_ARCSEC for r in delta_omega_rad]
    predicted_arcsec = [r * RAD_TO_ARCSEC for r in predicted_rad]
    residual_mas = [(m - p) * 1000.0 for m, p in zip(measured_arcsec, predicted_arcsec)]

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(6.5, 5.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )

    ax_top.plot(orbits, predicted_arcsec, linestyle="--", linewidth=1.0,
                color="0.4", alpha=0.8,
                label="Schwarzschild GR prediction (analytical)", zorder=1)
    ax_top.plot(orbits, measured_arcsec, linewidth=1.5,
                color="C0", label="apsis IAS15 + 1PN (measured)", zorder=2)
    ax_top.set_ylabel(r"$\Delta\omega$ [arcsec]")
    ax_top.legend(loc="upper left", frameon=False)
    ax_top.grid(True, alpha=0.3)
    ax_top.text(
        0.97, 0.05,
        "curves overlap at this scale —\nsee residual in lower panel",
        transform=ax_top.transAxes, fontsize=8, alpha=0.7,
        ha="right", va="bottom", style="italic",
    )

    n_orbits = orbits[-1]
    measured_end = measured_arcsec[-1]
    predicted_end = predicted_arcsec[-1]
    rel_err = abs(measured_end - predicted_end) / abs(predicted_end)

    ax_bot.axhline(0.0, color="0.3", linestyle="--", linewidth=1.0,
                   label="zero residual (perfect agreement)")
    ax_bot.plot(orbits, residual_mas, linewidth=1.0, color="C1",
                label="residual: measured $-$ GR")
    ax_bot.set_xlabel("orbit number")
    ax_bot.set_ylabel("residual [mas]")
    ax_bot.legend(loc="upper left", frameon=False, fontsize=8)
    ax_bot.grid(True, alpha=0.3)
    # Leave headroom above the residual curve so the drift annotation
    # sits in empty space.
    residual_max = max(abs(min(residual_mas)), abs(max(residual_mas)))
    ax_bot.set_ylim(top=residual_max * 1.45)
    ax_bot.text(
        0.98, 0.95,
        f"{rel_err * 1e6:.0f} ppm cumulative drift",
        transform=ax_bot.transAxes, ha="right", va="top",
        fontsize=8, alpha=0.85,
        bbox=dict(facecolor="white", edgecolor="none", pad=2.0),
    )

    fig.suptitle(
        f"Mercury perihelion precession — {n_orbits} orbits, "
        f"agreement within {rel_err * 1e6:.0f} ppm",
        fontsize=11,
    )

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    render(
        csv_path=root / "data" / "mercury_1pn_long_horizon_ias15.csv",
        output_path=root / "mercury_1pn_long_horizon.pdf",
    )
