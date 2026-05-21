# Python examples

Runnable scripts demonstrating each federated operator and their
composition.

## Setup

```bash
# from repository root
python -m venv .venv
. .venv/bin/activate    # or .venv\Scripts\activate on Windows
pip install maturin
maturin develop --release
```

## Run

```bash
python python-examples/mercury_perihelion.py
python python-examples/radiation_dust.py
python python-examples/central_precession.py
python python-examples/federation_composition.py
python python-examples/implicit_midpoint.py
```

## What each example shows

| script                        | operator                     | demonstrates                                                                       |
|-------------------------------|------------------------------|------------------------------------------------------------------------------------|
| `mercury_perihelion.py`       | `apsis.gr`                   | Mercury perihelion advance matches the GR closed form over 100 orbits              |
| `radiation_dust.py`           | `apsis.radiation`            | Dust grain at beta = 0.1 sees effective gravity reduced by (1 - beta), Burns 1979  |
| `central_precession.py`       | `apsis.central`              | Regime-based (`from_raw`) + observable-inversion (`from_apsidal_rate`, round-trip) |
| `federation_composition.py`   | `apsis.gr` + `apsis.central` | Two operators registered on one System compose additively                          |
| `implicit_midpoint.py`        | `apsis.gr` + integrator      | Mercury 1PN under the symplectic A-stable ImplicitMidpoint integrator              |

Each script asserts on its measurable claim and exits non-zero if
the claim breaks, so `python python-examples/<name>.py` doubles as a
contract test.

## See also

Rust counterparts under `crates/apsis/examples/` and
`crates/apsis-1pn/examples/`. The cross-implementation parity
portfolio against REBOUND lives at `crates/apsis/examples/rebound_parity_*.rs`.
