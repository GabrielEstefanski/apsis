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
//! - [`body`]: `Body` Python class — point-mass kinematics + softening.
//! - [`integrator`]: `IntegratorKind` and `TemplateKind` enums, plus
//!   the string-to-enum normalisation used by every kwargs-style
//!   constructor.
//! - [`system`]: `System` Python class — orchestration, builders,
//!   run methods, accessors.
//! - [`trajectory`]: `Trajectory` Python class — dense sampling
//!   results returned as NumPy arrays.
//! - [`perturbations`]: pre-bundled `PerturbationForce` instances
//!   (currently `OnePN`) so researchers can reproduce headline
//!   results without consuming the trait-level API.
//!
//! Each module owns one [`#[pyclass]`](pyo3::pyclass) (or one
//! cohesive group of related classes) and re-exports it through
//! [`register`] into the top-level [`#[pymodule]`](pyo3::pymodule)
//! [`_native`]. Adding a new class is a single-file addition plus
//! one line of [`register`] wiring; nothing else in this crate
//! changes.

use pyo3::prelude::*;

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
    Ok(())
}
