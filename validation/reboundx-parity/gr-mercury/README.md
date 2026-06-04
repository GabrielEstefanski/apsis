# REBOUNDx parity — Sun–Mercury 1PN (apsis-1pn vs REBOUNDx `gr`)

First parity check against a **REBOUNDx effect** (not just the base REBOUND
integrator): apsis-1pn's first post-Newtonian operator versus REBOUNDx's
`gr` force, on the Sun–Mercury orbit.

The two implement **different 1PN formulations** — apsis-1pn is the
test-particle Schwarzschild form applied pairwise (inertial frame); REBOUNDx
`gr` is single-dominant-mass in Jacobi coordinates with an iterative solve.
For Sun–Mercury (mass ratio ~1.7×10⁻⁷) both reduce to the test-particle
limit, so this measures the formulation/gauge difference, **not** bit-parity.

## Result (frozen reference in `out/`)

| | apsis | REBOUNDx | analytic |
|---|---|---|---|
| apsidal precession (″/century) | +42.9783 | +42.9784 | +42.9824 |

- **apsis vs REBOUNDx precession: 7×10⁻⁷** (gated ≤ 2×10⁻⁵).
- Both vs analytic Schwarzschild: ~9×10⁻⁵ (reported, not gated here —
  accuracy is owned by `crates/apsis-1pn/tests/mercury_precession_gate.rs`).
- Osculating invariants (a, e, h) and cross-implementation energy agree at the
  ULP floor (~7×10⁻¹⁵), gated ≤ 1×10⁻¹³.
- 1PN-off control: everything at the ULP floor, precession ≈ 0 — validates the
  harness independent of the 1PN physics.

## Run (Linux only)

REBOUNDx does **not** build on Windows/MSVC (C99 VLAs in `gr_full.c`); use WSL,
a container, or Linux CI. Install reboundx fresh against its co-installed
rebound (the cached-wheel RPATH to `librebound` breaks on reuse):

```bash
python -m venv .venv && . .venv/bin/activate
pip install --no-cache-dir -r requirements.txt
python run.py        # apsis (cargo) -> reboundx_side.py -> compare.py, gr + control
```

`run.py` needs both `cargo` and the reboundx venv in one Linux environment.
The apsis side (`reboundx_parity_gr` example) runs anywhere cargo runs; only
the REBOUNDx side is Linux-bound.

## Files

- `crates/apsis-1pn/examples/reboundx_parity_gr.rs` — apsis side (CSV dump).
- `reboundx_side.py` — REBOUND + REBOUNDx `gr` side.
- `compare.py` — gated comparator (orbital invariants + secular precession).
- `out/` — frozen reference CSVs + `comparison_{gr,control}.json`.
- Protocol: [`paper/notebooks/2026-05-29-reboundx-parity-gr.md`](../../../paper/notebooks/2026-05-29-reboundx-parity-gr.md).

## Scope

Confirms apsis-1pn matches an independent code in its **valid (test-particle)
regime**. It does not stress comparable-mass post-Newtonian dynamics — the
apsis-1pn formulation is not valid there and warns at registration.
