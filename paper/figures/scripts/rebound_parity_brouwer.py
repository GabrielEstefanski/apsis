"""Long-horizon energy parity: retrograde Kepler e=0.5 over 10⁴
orbits, |E(t)−E₀|/|E₀| in log-log, both apsis IAS15 and REBOUND
IAS15 on the same axes plus a √N reference line (Brouwer 1937).

Frozen CSV snapshot in `paper/figures/data/rebound_parity_retrograde_*`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "rebound_parity_brouwer.pdf"


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
    apsis, rebound = load("retrograde")
    fig, ax = plt.subplots(figsize=(7.0, 4.0))

    e0_a = float(apsis["e_total"].iloc[0])
    e0_r = float(rebound["e_total"].iloc[0])
    rel_a = np.abs(np.asarray(apsis["e_total"]) - e0_a) / np.abs(e0_a)
    rel_r = np.abs(np.asarray(rebound["e_total"]) - e0_r) / np.abs(e0_r)
    orbit_a = np.asarray(apsis["orbit"])
    orbit_r = np.asarray(rebound["orbit"])
    mask_a = orbit_a >= 1
    mask_r = orbit_r >= 1

    ax.loglog(orbit_a[mask_a], rel_a[mask_a], "C0-", linewidth=0.9, label="apsis")
    ax.loglog(orbit_r[mask_r], rel_r[mask_r], "C1:", linewidth=1.3, alpha=0.7, label="REBOUND")
    n_vals = np.array([1, 10000])
    floor = float(np.median(rel_a[mask_a]))
    sqrt_n = floor * np.sqrt(n_vals / np.median(orbit_a[mask_a]))
    ax.loglog(n_vals, sqrt_n, "k--", linewidth=0.8, alpha=0.5, label=r"$\sqrt{N}$ (Brouwer 1937)")

    m = min(len(apsis), len(rebound))
    ea = np.asarray(apsis["e_total"])[:m]
    er = np.asarray(rebound["e_total"])[:m]
    cross_max = float(np.max(np.abs(ea - er) / abs(e0_r)))

    ax.set_xlabel("orbit")
    ax.set_ylabel(r"$|E(t) - E_0|\,/\,|E_0|$")
    ax.set_title(r"Retrograde Kepler $e=0.5$ — $10^4$ orbits", fontsize=10)
    ax.text(
        0.03,
        0.97,
        rf"max cross-impl $|\Delta E|/|E_0| = {sci(cross_max)}$",
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.85),
    )
    ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    # Pin metadata date so re-runs are byte-identical (matplotlib stamps
    # CreationDate=now by default, which otherwise dirties git on every regen).
    fig.savefig(OUT, format="pdf", bbox_inches="tight", metadata={"CreationDate": None})
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
