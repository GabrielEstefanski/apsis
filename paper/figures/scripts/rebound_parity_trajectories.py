"""Configuration-space parity across three canonical N-body scenarios:
Kepler e=0.5 (regular), Figure-8 choreography (regular, symmetric),
and Pythagorean three-body (chaotic). Three-panel row, each panel
plots apsis IAS15 trajectory overlaid by REBOUND IAS15.

Frozen CSV snapshots in `paper/figures/data/rebound_parity_*`. Refresh
by re-running the harnesses under `validation/rebound-parity/` and
copying their `out/{apsis,rebound,comparison}.json` outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "rebound_parity_trajectories.pdf"


def sci(v: float) -> str:
    """Format as LaTeX m×10ⁿ with a two-decimal mantissa."""
    exp = int(np.floor(np.log10(abs(v))))
    mant = v / 10 ** exp
    if round(mant, 2) >= 10.0:  # 9.999 rounds to 10.00 at 2dp -> renormalise
        exp += 1
        mant /= 10.0
    return rf"{mant:.2f}\times10^{{{exp}}}"


def cross_energy(apsis: pd.DataFrame, rebound: pd.DataFrame) -> float:
    """Max over the run of |E_apsis − E_rebound| / |E_0|."""
    m = min(len(apsis), len(rebound))
    e0 = abs(float(apsis["e_total"].iloc[0]))
    ea = np.asarray(apsis["e_total"])[:m]
    er = np.asarray(rebound["e_total"])[:m]
    return float(np.max(np.abs(ea - er) / e0))


def load(scenario: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = pd.read_csv(
        DATA / f"rebound_parity_{scenario}_apsis.csv", comment="#", encoding="latin-1"
    )
    r = pd.read_csv(
        DATA / f"rebound_parity_{scenario}_rebound.csv", comment="#", encoding="latin-1"
    )
    return a, r


def kepler_ellipse_xy(
    a: float, e: float, n: int = 360
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    r = a * (1.0 - e * e) / (1.0 + e * np.cos(theta))
    return r * np.cos(theta), r * np.sin(theta)


def annotate(ax: Axes, text: str) -> None:
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.85),
    )


def panel_kepler(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    # apsis/rebound are the once-per-orbit gate samples (the |Δr| headline over
    # 100 orbits); the dense trace traces the ellipse for the visual.
    ex, ey = kepler_ellipse_xy(1.0, 0.5)
    ax.plot(ex, ey, color="grey", linewidth=0.6, alpha=0.5, label="analytical")
    trace_a, trace_r = load("kepler_trace")
    ax.plot(trace_a["x1"], trace_a["y1"], color="#1f77b4", linewidth=3.0, alpha=0.35, label="apsis")
    # REBOUND over one orbit only: it retraces the same ellipse each orbit, and
    # overlapping passes would fill the dashes into a solid line.
    one = len(trace_r) // 2 + 1
    ax.plot(
        trace_r["x1"][:one],
        trace_r["y1"][:one],
        color="black",
        linewidth=1.1,
        linestyle="--",
        label="REBOUND",
    )
    dr = np.hypot(
        np.asarray(apsis["x1"]) - np.asarray(rebound["x1"]),
        np.asarray(apsis["y1"]) - np.asarray(rebound["y1"]),
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Kepler $e=0.5$ — 100 orbits", fontsize=10)
    ax.legend(loc="lower left", fontsize=7, framealpha=0.85)
    annotate(ax, rf"cross-impl $|\Delta\mathbf{{r}}|_{{\max}} = {sci(dr.max())}$")


def panel_figure8(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    # apsis is the coloured solid base; REBOUND overlays as a dashed line (solid
    # model vs dashed reference). The three choreography bodies trace one curve,
    # so REBOUND is drawn once (body 0); three phase-offset dashed curves would
    # interleave into a solid.
    shades = ["#1f77b4", "#4a90d9", "#7fb3e3"]
    for b in range(3):
        ax.plot(apsis[f"x{b}"], apsis[f"y{b}"], color=shades[b], linewidth=3.0, alpha=0.35)
    # REBOUND over one period only: body 0 traces the full figure-8 each period,
    # and overlapping passes would fill the dashes into a solid line.
    one = len(rebound) // 10 + 1
    ax.plot(rebound["x0"][:one], rebound["y0"][:one], color="black", linewidth=1.1, linestyle="--")
    apsis_proxy = Line2D([], [], color=shades[0], linewidth=3.0, alpha=0.35, label="apsis")
    reb_proxy = Line2D([], [], color="black", linestyle="--", linewidth=1.1, label="REBOUND")
    ax.legend(handles=[apsis_proxy, reb_proxy], loc="lower right", fontsize=7, framealpha=0.85)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Figure-8 choreography — 10 periods", fontsize=10)
    annotate(ax, rf"cross-impl $|\Delta E|/|E_0| = {sci(cross_energy(apsis, rebound))}$")


def panel_pythagorean(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    colors = ["C0", "C1", "C2"]
    for b in range(3):
        ax.plot(
            apsis[f"x{b}"],
            apsis[f"y{b}"],
            color=colors[b],
            linewidth=1.0,
            label=f"body {b}",
        )
    for b in range(3):
        ax.plot(
            rebound[f"x{b}"],
            rebound[f"y{b}"],
            color="black",
            linewidth=1.0,
            linestyle="--",
            alpha=0.7,
            label="REBOUND" if b == 0 else None,
        )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Pythagorean (Burrau 1913) — $T=70$", fontsize=10)
    ax.legend(loc="lower right", fontsize=7, framealpha=0.85)
    annotate(
        ax,
        rf"cross-impl $|\Delta E|/|E_0| = {sci(cross_energy(apsis, rebound))}$"
        + "\n(chaotic; both at REBOUND-class floor)",
    )


def main() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 6.0))
    panel_kepler(axes[0], *load("kepler"))
    panel_figure8(axes[1], *load("figure8"))
    panel_pythagorean(axes[2], *load("pythagorean"))
    fig.suptitle("REBOUND parity: configuration-space trajectories", fontsize=11)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    # Pin metadata date so re-runs are byte-identical (matplotlib stamps
    # CreationDate=now by default, which otherwise dirties git on every regen).
    fig.savefig(OUT, format="pdf", bbox_inches="tight", metadata={"CreationDate": None})
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
