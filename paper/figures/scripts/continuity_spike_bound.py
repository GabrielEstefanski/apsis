"""Continuity counter-test energy-spike figure (paper §3.3).

Top panel: each R_c-crossing's measured ``|ΔE/E|`` against the closed-form
jump-bound. Bottom panel: the separation ``r(t)`` crossing R_c — aligned in
time with the spikes above, so the spike-per-crossing bijection is visible.
Reads ``data/continuity_spike_bound.csv`` and ``data/continuity_trajectory.csv``.
Output: ``continuity_spike_bound.pdf``.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load_crossings(path: Path) -> list[dict[str, float]]:
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(lines)]


def load_trajectory(path: Path) -> tuple[list[float], list[float]]:
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    reader = csv.DictReader(lines)
    t, r = [], []
    for row in reader:
        t.append(float(row["t"]))
        r.append(float(row["r"]))
    return t, r


def render(data_dir: Path, output_path: Path) -> None:
    crossings = load_crossings(data_dir / "continuity_spike_bound.csv")
    t_traj, r_traj = load_trajectory(data_dir / "continuity_trajectory.csv")

    t_cross = [c["t_cross"] for c in crossings]
    spike = [c["spike_rel"] for c in crossings]
    bound = max(c["bound_rel"] for c in crossings)

    fig, (ax_e, ax_r) = plt.subplots(
        2, 1, figsize=(6.5, 4.8), sharex=True,
        gridspec_kw={"height_ratios": [2, 1.1], "hspace": 0.1},
    )

    # ── Top: spike magnitudes vs the analytic bound ──────────────────────
    ax_e.axhline(
        bound, color="C3", linestyle="--", linewidth=1.2, zorder=2,
        label=r"jump-bound $\Delta F\,v_{\rm cross}\,\delta t/|E_0|$",
    )
    ax_e.scatter(
        t_cross, spike, s=42, color="C0", edgecolor="white", linewidth=0.6, zorder=3,
        label=r"spike at each $R_c$ crossing",
    )
    ax_e.set_yscale("log")
    ax_e.set_ylim(2e-6, 1e-3)
    ax_e.set_ylabel(r"$|\Delta E / E|$")
    ax_e.grid(True, which="major", alpha=0.3)
    ax_e.legend(loc="lower center", frameon=False, fontsize=8, ncol=2)

    # ── Bottom: separation crossing R_c, aligned with the spikes above ───
    ax_r.plot(t_traj, r_traj, color="0.5", linewidth=0.9, zorder=2, label=r"separation $r(t)$")
    ax_r.axhline(1.0, color="0.5", linestyle=":", linewidth=1.1, zorder=1, label=r"$R_c$")
    for tc in t_cross:
        ax_r.axvline(tc, color="0.8", linewidth=0.6, zorder=0)
    ax_r.set_ylabel(r"separation $r/a$")
    ax_r.set_xlabel("simulation time")
    ax_r.set_xlim(left=0.0)
    ax_r.set_ylim(0.3, 2.7)
    ax_r.grid(True, which="major", alpha=0.3)
    ax_r.legend(loc="upper center", frameon=False, fontsize=8, ncol=2)

    fig.suptitle(
        "Continuity counter-test: energy-error spikes stay within the analytic bound",
        fontsize=10,
    )

    # Pin metadata date so re-runs are byte-identical.
    fig.savefig(output_path, bbox_inches="tight", metadata={"CreationDate": None})
    plt.close(fig)
    print(f"wrote {output_path}")


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    render(data_dir=root / "data", output_path=root / "continuity_spike_bound.pdf")
