"""Configuration-space parity across three canonical N-body scenarios.

Kepler e=0.5 and the figure-8 choreography: apsis and REBOUND IAS15 agree to
machine precision, so each trajectory is drawn once (the apsis curve) and REBOUND
is shown as decimated reference markers riding that curve over a single pass —
the samples land on the line because the cross-implementation separation stays
far below plot resolution (the f64 floor for the figure-8; residual phase drift
for Kepler). The Pythagorean three-body is chaotic: the trajectories genuinely
diverge, so both are drawn (apsis solid, REBOUND dashed) and the curves visibly
drift apart.

Frozen CSV snapshots in `paper/figures/data/rebound_parity_*`. Refresh by
re-running the harnesses under `validation/rebound-parity/` and copying their
`out/{apsis,rebound}.csv` outputs (the Kepler panel also uses the dense
`rebound_parity_kepler_trace_*` produced by the example's `--trace-output`).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "rebound_parity_trajectories.pdf"

APSIS = "#1f4e79"  # model under validation: solid coloured line
REBOUND_C = "#2b2b2b"  # reference oracle: charcoal, constant role across panels
BODY = ["#1f4e79", "#c97a1e", "#3a7d44"]  # colourblind-safe per-body apsis (Pythagorean)


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


def annotate(ax: Axes, text: str) -> None:
    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.9),
    )


def panel_kepler(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    # apsis/rebound are the once-per-orbit gate samples (the |Δr| headline over
    # 100 orbits); the dense trace draws one clean orbit and REBOUND rides it.
    trace_a, trace_r = load("kepler_trace")
    one = trace_a["orbit"] == 0
    # The orbit is periodic; the dense trace samples [0, P) and omits the closing
    # arc, so repeat the first sample to join the ellipse.
    xa, ya = np.asarray(trace_a["x1"][one]), np.asarray(trace_a["y1"][one])
    ax.plot(
        np.append(xa, xa[0]), np.append(ya, ya[0]), color=APSIS, linewidth=1.7,
        solid_capstyle="round", zorder=2, label="apsis",
    )
    idx = np.arange(0, int(one.sum()), 8)
    ax.plot(
        np.asarray(trace_r["x1"][one])[idx], np.asarray(trace_r["y1"][one])[idx],
        linestyle="none", marker="o", markerfacecolor="none", markeredgecolor=REBOUND_C,
        markeredgewidth=1.1, markersize=6, zorder=3, label="REBOUND (samples)",
    )
    ax.plot(0, 0, marker="*", color="goldenrod", markersize=11, linestyle="none",
            markeredgecolor="#7a5c00", markeredgewidth=0.5, label="primary (focus)")
    dr = np.hypot(
        np.asarray(apsis["x1"]) - np.asarray(rebound["x1"]),
        np.asarray(apsis["y1"]) - np.asarray(rebound["y1"]),
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.16)  # breathing room so the ellipse sits smaller in the frame
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Kepler $e=0.5$ — 100 orbits", fontsize=10)
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9)
    annotate(ax, rf"cross-impl $|\Delta\mathbf{{r}}|_{{\max}} = {sci(dr.max())}$")


def panel_figure8(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    # The three bodies chase each other on one shared curve, so it is drawn once
    # (apsis body 0 over a single period) with the t=0 body positions marked;
    # REBOUND rides it as decimated samples.
    one = len(apsis) // 10 + 1
    ax.plot(
        apsis["x0"][:one], apsis["y0"][:one], color=APSIS, linewidth=1.7,
        solid_capstyle="round", zorder=2, label="apsis (1 curve, 3 bodies)",
    )
    for b in range(3):
        ax.plot(
            apsis[f"x{b}"].iloc[0], apsis[f"y{b}"].iloc[0], marker="o", color=APSIS,
            markersize=5, linestyle="none", zorder=4,
        )
    idx = np.arange(0, one, 10)
    ax.plot(
        np.asarray(rebound["x0"])[idx], np.asarray(rebound["y0"])[idx], linestyle="none",
        marker="o", markerfacecolor="none", markeredgecolor=REBOUND_C, markeredgewidth=1.1,
        markersize=6, zorder=3, label="REBOUND (samples)",
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Figure-8 choreography — 10 periods", fontsize=10)
    ax.legend(loc="lower center", fontsize=7, framealpha=0.9)
    annotate(ax, rf"cross-impl $|\Delta E|/|E_0| = {sci(cross_energy(apsis, rebound))}$")


def panel_pythagorean(ax: Axes, apsis: pd.DataFrame, rebound: pd.DataFrame) -> None:
    # Chaotic: the trajectories diverge, so both are drawn (apsis solid per body,
    # REBOUND dashed) and the curves visibly drift apart.
    for b in range(3):
        ax.plot(apsis[f"x{b}"], apsis[f"y{b}"], color=BODY[b], linewidth=0.9, zorder=2,
                label=f"body {b}")
    for b in range(3):
        ax.plot(
            rebound[f"x{b}"], rebound[f"y{b}"], color=REBOUND_C, linewidth=0.9,
            linestyle=(0, (4, 3)), alpha=0.85, zorder=3,
            label="REBOUND" if b == 0 else None,
        )
    ax.set_aspect("equal", adjustable="datalim")
    # Frame the interaction region (tangle + binary recoil) and let the ejected
    # body's ray exit the top edge; the bottom tracks the binary's recoil depth
    # so the crop follows the physics, not the arbitrary integration cutoff.
    ymin = min(float(np.min(apsis[f"y{b}"])) for b in range(3))
    ax.set_ylim(ymin - 0.5, -ymin)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Pythagorean (Burrau 1913) — $T=70$", fontsize=10)
    ax.legend(loc="lower right", fontsize=6, framealpha=0.9)
    annotate(
        ax,
        rf"cross-impl $|\Delta E|/|E_0| = {sci(cross_energy(apsis, rebound))}$"
        + "\n(chaotic; REBOUND-class floor)",
    )


def main() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.6))
    panel_kepler(axes[0], *load("kepler"))
    panel_figure8(axes[1], *load("figure8"))
    panel_pythagorean(axes[2], *load("pythagorean"))
    fig.suptitle("Configuration-space trajectories: apsis vs REBOUND IAS15", fontsize=11)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    # Pin metadata date so re-runs are byte-identical (matplotlib stamps
    # CreationDate=now by default, which otherwise dirties git on every regen).
    fig.savefig(OUT, format="pdf", bbox_inches="tight", metadata={"CreationDate": None})
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
