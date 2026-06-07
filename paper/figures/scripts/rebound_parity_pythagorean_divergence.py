"""Cross-implementation trajectory divergence for the chaotic Pythagorean
three-body (Burrau 1913). apsis and REBOUND IAS15 start from identical conditions;
max_b |Δr_b|(t) climbs from the round-off floor through an exponential
(positive-Lyapunov) ramp to the ejection, while the cross-implementation energy
difference stays at the REBOUND-class chaotic floor (~4e-10) — the divergence is
chaos, not a solver disagreement. Companion to the regular-orbit panels, where the
two codes agree to machine precision.

Frozen CSV snapshot in `paper/figures/data/rebound_parity_pythagorean_*`. Refresh
by re-running the harness under `validation/rebound-parity/pythagorean/` and
copying its `out/{apsis,rebound}.csv` outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "rebound_parity_pythagorean_divergence.pdf"

APSIS = "#1f4e79"
LEVELS = (1e-8, 1e-6, 1e-3, 1e-1)  # decade anchors marking the Lyapunov ramp


def sci(v: float) -> str:
    """Format as LaTeX m×10ⁿ with a two-decimal mantissa."""
    exp = int(np.floor(np.log10(abs(v))))
    mant = v / 10 ** exp
    if round(mant, 2) >= 10.0:  # 9.999 rounds to 10.00 at 2dp -> renormalise
        exp += 1
        mant /= 10.0
    return rf"{mant:.2f}\times10^{{{exp}}}"


def load(scenario: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    a = pd.read_csv(
        DATA / f"rebound_parity_{scenario}_apsis.csv", comment="#", encoding="latin-1"
    )
    r = pd.read_csv(
        DATA / f"rebound_parity_{scenario}_rebound.csv", comment="#", encoding="latin-1"
    )
    return a, r


def main() -> None:
    apsis, rebound = load("pythagorean")
    m = min(len(apsis), len(rebound))
    t = np.asarray(apsis["t"])[:m]
    sep = np.maximum.reduce([
        np.hypot(
            np.asarray(apsis[f"x{b}"])[:m] - np.asarray(rebound[f"x{b}"])[:m],
            np.asarray(apsis[f"y{b}"])[:m] - np.asarray(rebound[f"y{b}"])[:m],
        )
        for b in range(3)
    ])
    e0 = abs(float(apsis["e_total"].iloc[0]))
    cross_e = float(
        np.max(np.abs(np.asarray(apsis["e_total"])[:m] - np.asarray(rebound["e_total"])[:m]) / e0)
    )

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.semilogy(t, np.clip(sep, 1e-16, None), color=APSIS, linewidth=1.4)
    ax.set_ylim(1e-16, 1.0)

    # Decade anchors: first time the separation crosses each level, labelled at
    # the foot of the axis so they clear the curve and the peak callout.
    for lv in LEVELS:
        i = int(np.argmax(sep >= lv))
        if sep[i] >= lv:
            ax.axvline(t[i], color="grey", linewidth=0.6, linestyle=":", alpha=0.6)
            ax.text(t[i], 1.6e-16, rf"$t={t[i]:.1f}$", fontsize=7, color="grey",
                    ha="center", va="bottom")

    ipk = int(np.argmax(sep))
    ax.plot(t[ipk], sep[ipk], marker="o", color=APSIS, markersize=5, zorder=4)
    ax.annotate(
        rf"peak ${sci(sep[ipk])}$ @ $t={t[ipk]:.1f}$ (ejection)",
        xy=(t[ipk], sep[ipk]), xytext=(38, 3e-12), fontsize=8, ha="left",
        arrowprops=dict(arrowstyle="->", color="grey", lw=0.7),
    )

    ax.set_xlabel("t")
    ax.set_ylabel(r"cross-impl separation  $\max_b|\Delta\mathbf{r}_b|$")
    ax.set_title(
        "Pythagorean (Burrau 1913): apsis vs REBOUND IAS15 diverge (chaotic)", fontsize=10
    )
    ax.text(
        0.03, 0.97,
        "shared ICs; trajectories diverge exponentially\n"
        rf"energy held at $|\Delta E|/|E_0| = {sci(cross_e)}$",
        transform=ax.transAxes, fontsize=8, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.85),
    )
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    # Pin metadata date so re-runs are byte-identical (matplotlib stamps
    # CreationDate=now by default, which otherwise dirties git on every regen).
    fig.savefig(OUT, format="pdf", bbox_inches="tight", metadata={"CreationDate": None})
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
