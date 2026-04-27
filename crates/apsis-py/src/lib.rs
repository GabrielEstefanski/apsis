//! Python bindings for the [`apsis`] N-body simulation library.
//!
//! # Façade-only invariant
//!
//! This crate is a **thin** binding layer over [`apsis`]. It exists to
//! translate values across the Rust↔Python boundary and to validate
//! arguments at that boundary; it does **not** implement physics,
//! integration, conservation diagnostics, or any algorithmic logic of
//! its own. If a feature requires logic, it belongs in the [`apsis`]
//! crate, not here.
//!
//! Concretely, every PR to this crate must satisfy:
//!
//! 1. The Rust types exposed via [`pyo3`] are wrappers around
//!    [`apsis`] types (newtypes, references, or owned copies). They
//!    do not redeclare the underlying state.
//! 2. Each `#[pymethods]` body either delegates to a single [`apsis`]
//!    public-API call or performs argument validation followed by
//!    such a delegation. No multi-step physics composition lives
//!    here.
//! 3. Integrator behaviour, conservation invariants, and the
//!    public-API contract are validated by the parent crate's test
//!    suite (`crates/apsis/tests/`, `crates/apsis-1pn/tests/`, and
//!    the cross-implementation parity portfolio under
//!    `validation/rebound-parity/`). This crate's tests
//!    (`tests/test_*.py`) cover only the binding surface — that
//!    Python kwargs translate to the right Rust call and that
//!    Python-side type errors are raised at the boundary, not
//!    inside the integrator.
//!
//! Reviewers: any PR that pushes new logic into this crate without
//! a matching change in [`apsis`] is a façade violation. Reject it
//! and ask the contributor to land the logic in the core first.
//!
//! # Module layout
//!
//! The Rust side is organised by concern, one wrapper per file:
//!
//! - [`body`]: `Body` Python class — point-mass kinematics, softening,
//!   and material classification, with the nine material factories
//!   (`star`, `rocky`, `gas_giant`, ...) and the immutable fluent
//!   builder (`at`, `with_velocity`, `with_density`, `unsoftened`).
//! - [`integrator`]: `IntegratorKind` enum exposed to Python under
//!   upper-case acronym names (`IAS15`, `YOSHIDA4`, ...) plus the
//!   string-slug normalisation [`integrator::resolve`] used by every
//!   wrapper that takes an `integrator=` kwarg.
//! - [`system`]: `System` Python class — orchestration with kwargs
//!   constructor, run-loop verbs (`step`, `integrate_for`,
//!   `integrate_until`), and read-only diagnostic properties
//!   (`t`, `bodies`, `energy`, `energy_delta`, ...).
//! - [`trajectory`]: `Trajectory` Python class — dense NumPy-backed
//!   record returned by `System.sample`, with shape-`(n_samples,)`
//!   `t` / `energy` axes and shape-`(n_samples, n_bodies)` `x` / `y`
//!   / `vx` / `vy` axes ready for `matplotlib`.
//! - [`convert`]: shared boundary helpers (error formatting, 2-vector
//!   parsing, slug normalisation). Owned by no single wrapper; called
//!   from all of them.
//!
//! Each module owns one [`#[pyclass]`](pyo3::pyclass) (or one
//! cohesive group of related classes) and exposes a `pub(crate)
//! register` function that is called from the [`#[pymodule]`](pyo3::pymodule)
//! entry point [`_native`] below. Adding a new class is a single-file
//! addition plus one line of registration here; nothing else changes.

use pyo3::prelude::*;

mod body;
mod convert;
mod integrator;
mod system;
mod trajectory;
mod units;

/// `apsis._native`: the Rust-built extension module.
///
/// User-facing imports flow through `python/apsis/__init__.py`, which
/// re-exports selectively from this module. The `_native` namespace is
/// considered private; researchers should `import apsis` and access
/// classes through the package façade so the Rust/Python split
/// remains an implementation detail.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    body::register(m)?;
    integrator::register(m)?;
    system::register(m)?;
    trajectory::register(m)?;
    units::register(m)?;
    Ok(())
}
