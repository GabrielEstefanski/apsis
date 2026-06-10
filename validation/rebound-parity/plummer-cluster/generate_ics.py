"""Seeded Plummer-model IC generator for the cluster parity protocol.

Sampling: Aarseth, Hénon & Wielen (1974) — inverse-CDF radius from a uniform
mass fraction, isotropic directions, velocity modulus by rejection from
g(q) = q²(1−q²)^{7/2} against the local escape speed. Sampled at G = M = a = 1,
then centred (one pass: subtract mass-weighted mean position and velocity) and
rescaled per realisation to N-body units: PE = −1/2, KE = 1/4 (E = −1/4,
exact virial ratio at t = 0 for the unsoftened potential).

The committed CSV is the canonical artefact; this script is its provenance.
ε is computed here (0.98 · N^(−0.26) Plummer scale lengths; Athanassoula et
al. 2000) and embedded in the header so every consumer parses one source.

Run:
    python generate_ics.py --n 256
    python generate_ics.py --n 1000
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

SEED = 20260609
A_PLUMMER = 3.0 * math.pi / 16.0  # scale radius in N-body units (W = -3*pi*GM^2/32a = -1/2)


def eps_protocol(n: int) -> float:
    return 0.98 * float(n) ** -0.26 * A_PLUMMER


def iso_dir(rng: np.random.Generator) -> np.ndarray:
    z = 2.0 * rng.random() - 1.0
    phi = 2.0 * math.pi * rng.random()
    s = math.sqrt(1.0 - z * z)
    return np.array([s * math.cos(phi), s * math.sin(phi), z])


def sample_plummer(n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    pos = np.empty((n, 3))
    vel = np.empty((n, 3))
    for i in range(n):
        m = rng.random()
        while m <= 0.0:
            m = rng.random()
        r = (m ** (-2.0 / 3.0) - 1.0) ** -0.5
        pos[i] = iso_dir(rng) * r
        v_esc = math.sqrt(2.0) * (1.0 + r * r) ** -0.25
        while True:
            x4 = rng.random()
            x5 = rng.random()
            if 0.1 * x5 < x4 * x4 * (1.0 - x4 * x4) ** 3.5:
                break
        vel[i] = iso_dir(rng) * (x4 * v_esc)
    return pos, vel


def potential_energy(masses: np.ndarray, pos: np.ndarray) -> float:
    n = len(masses)
    pe = 0.0
    for i in range(n):
        d = pos[i + 1 :] - pos[i]
        r = np.sqrt(np.einsum("ij,ij->i", d, d))
        pe -= float(np.sum(masses[i] * masses[i + 1 :] / r))
    return pe


def kinetic_energy(masses: np.ndarray, vel: np.ndarray) -> float:
    return 0.5 * float(np.sum(masses * np.einsum("ij,ij->i", vel, vel)))


def main() -> int:
    p = argparse.ArgumentParser(description="Plummer cluster IC generator (protocol artefact).")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--output", default=None)
    args = p.parse_args()

    n = args.n
    out = Path(args.output) if args.output else Path(f"ics_n{n}.csv")
    rng = np.random.default_rng(args.seed)

    masses = np.full(n, 1.0 / n)
    pos, vel = sample_plummer(n, rng)

    pos -= np.average(pos, axis=0, weights=masses)
    vel -= np.average(vel, axis=0, weights=masses)

    pe = potential_energy(masses, pos)
    pos *= pe / -0.5
    ke = kinetic_energy(masses, vel)
    vel *= math.sqrt(0.25 / ke)

    pe = potential_energy(masses, pos)
    ke = kinetic_energy(masses, vel)
    com = np.average(pos, axis=0, weights=masses)
    ptot = np.sum(masses[:, None] * vel, axis=0)
    r_sorted = np.sort(np.sqrt(np.einsum("ij,ij->i", pos, pos)))
    r_half = float(r_sorted[n // 2 - 1])

    assert abs(pe + 0.5) < 1e-12, f"PE after rescale: {pe}"
    assert abs(ke - 0.25) < 1e-12, f"KE after rescale: {ke}"
    assert float(np.linalg.norm(com)) < 1e-13, f"COM residual: {com}"
    assert float(np.linalg.norm(ptot)) < 1e-13, f"P residual: {ptot}"
    assert 0.5 < r_half < 1.1, f"half-mass radius {r_half} far from Plummer 0.77"

    eps = eps_protocol(n)
    with out.open("w", newline="\n") as f:
        f.write(
            "# Plummer cluster ICs — protocol:"
            " paper/notebooks/2026-06-09-rebound-parity-plummer-cluster.md\n"
        )
        f.write(f"# n={n}\n")
        f.write(f"# seed={args.seed}\n")
        f.write(f"# eps={eps!r}\n")
        f.write(f"# pe={pe!r} ke={ke!r}\n")
        f.write(
            f"# com_residual={float(np.linalg.norm(com))!r}"
            f" p_residual={float(np.linalg.norm(ptot))!r}\n"
        )
        f.write(f"# r_half={r_half!r}\n")
        f.write("body,m,x,y,z,vx,vy,vz\n")
        for i in range(n):
            f.write(
                f"{i},{float(masses[i])!r},{float(pos[i,0])!r},{float(pos[i,1])!r},{float(pos[i,2])!r},"
                f"{float(vel[i,0])!r},{float(vel[i,1])!r},{float(vel[i,2])!r}\n"
            )
    print(f"wrote {n} bodies to {out} (eps={eps:.6e}, r_half={r_half:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
