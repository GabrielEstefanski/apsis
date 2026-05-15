# apsis (Python)

*Python bindings for the [apsis](../../README.md) N-body simulation library.*

A researcher-first Python API for setting up and integrating
gravitational systems. The underlying numerical work — adaptive
high-order integration, conservation diagnostics, post-Newtonian
corrections — is performed in Rust by the `apsis` core crate; this
package is a thin façade that exposes that surface to Python in an
idiomatic, kwargs-friendly form.

> *Status: pre-release scaffolding (`v0.1.0` alpha). The package
> compiles and imports; the user-facing API surface is being landed
> incrementally over the next phase. The headline Mercury-perihelion
> reproduction in five lines is the deliberate target of the first
> tagged release.*

---

## Quickstart (target API, landing incrementally)

Reproduce Mercury's perihelion precession at the General-Relativistic
prediction:

```python
import apsis

sys = apsis.mercury_with_gr(orbits=500)
sys.run()
print(f"{sys.precession_rate():.3f} arcsec/century")  # ≈ 42.98
```

Or build a custom system explicitly:

```python
import apsis

sun = apsis.Body.star(mass=1.0)
mercury = (apsis.Body.rocky(mass=3e-6)
              .at(0.307, 0.0)
              .with_velocity(0.0, 1.98))

sys = apsis.System(
    bodies=[sun, mercury],
    integrator="ias15",
    dt=1e-3,
)
sys.integrate_for(100.0)

print(f"dE/E = {sys.energy_delta:.3e}")
print(f"|L| drift = {sys.lz_delta:.3e}")
```

Sample a trajectory ready for `matplotlib`:

```python
traj = sys.sample(duration=10.0, n_samples=1000)
# traj.t : np.ndarray, shape (n_samples,)
# traj.x : np.ndarray, shape (n_samples, n_bodies)
# traj.y : np.ndarray, shape (n_samples, n_bodies)
import matplotlib.pyplot as plt
plt.plot(traj.x[:, 1], traj.y[:, 1])  # body 1's path
plt.show()
```

---

## Installation

When the first wheels land on PyPI:

```bash
pip install apsis
```

Wheels are built for `manylinux2014_x86_64`, `manylinux2014_aarch64`,
`x86_64-apple-darwin`, `aarch64-apple-darwin` (M-series Macs), and
`x86_64-pc-windows-msvc`, against the stable Python ABI from CPython
3.9 onward. A single wheel per platform serves all supported Python
versions.

### Building from source (development)

A Rust toolchain (`rustup install stable`, ≥ 1.85) and a Python
development environment are required. From the workspace root:

```bash
pip install maturin
maturin develop --manifest-path crates/apsis-py/Cargo.toml
python -c "import apsis; print(apsis.__version__)"
```

`maturin develop` builds the Rust extension and installs it into the
current virtual environment in editable mode. The resulting
`apsis._native` module is loaded by `python/apsis/__init__.py` at
import time.

---

## Architecture

The package is a façade. The Rust core crate `apsis` is the source of
truth for every numerical decision; the binding layer here is
restricted to:

1. Type translation across the FFI boundary (Rust enums ↔ Python
   strings/enum classes; `Vec<Body>` ↔ Python list; `f64` ↔ Python
   `float`).
2. Argument validation at the boundary (rejecting unknown integrator
   names with `ValueError`; rejecting negative masses with
   `ValueError`; rejecting empty body lists where the underlying API
   would panic).
3. Re-presentation of return values in shapes researchers reach for
   (NumPy arrays for trajectories; `dataclass`-style accessors for
   metrics).

If a feature requires logic, it goes in the `apsis` crate. The
binding is not the place to compose new behaviour — the source-level
documentation of `crates/apsis-py/src/lib.rs` makes this an
explicit reviewer rule.

---

## Validation

The Rust core is validated by:

- the **Mercury perihelion precession** test (4.4 ppm of the
  General-Relativistic prediction over 500 orbits), and
- the **cross-implementation parity portfolio** against REBOUND
  on canonical scenarios (Kepler `e = 0.5` over 100 orbits at 1–3 ULP;
  the Chenciner–Montgomery figure-8 over 10 periods at 1 ULP across
  twelve gated metrics).

The Python bindings inherit this validation by construction: every
wrapped API call delegates to the same Rust functions exercised by
those test suites. The Python-side test suite (`tests/`) covers the
binding surface itself — that kwargs translate to the right Rust
calls, that invalid arguments raise clean Python exceptions, and
that NumPy arrays returned by sampling APIs match the corresponding
Rust trajectory data bit-for-bit.

See the workspace [`README.md`](../../README.md) for the project's
research positioning and the lab-notebook trail under
[`docs/experiments/`](../../docs/experiments/) for the methodological
record behind each validation result.
