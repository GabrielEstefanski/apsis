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
//! Each module owns one `#[pyclass]` and exposes a `pub(crate) register`
//! called from the [`_native`] entry point below.
//!
//! - [`body`] — `Body` class with the nine material factories and fluent builder
//! - [`integrator`] — `IntegratorKind` enum + slug parser
//! - [`perturbation`] — `add_hamiltonian_perturbation` boundary helper
//!   that pulls a `Box<dyn HamiltonianOperator>` out of the user-facing
//!   pure-Python `apsis.Perturbation` wrapper (see `python/apsis/__init__.py`)
//! - [`stats`] — `Stats` and `AdaptiveStats` frozen diagnostic snapshots
//! - [`system`] — `System` orchestrator (constructor, run loop, diagnostics)
//! - [`trajectory`] — `Trajectory` NumPy-backed return value of `System.sample`
//! - [`units`] — `UnitSystem` class + `apsis.units` submodule of singletons
//! - [`convert`] — boundary helpers (error formatting, 2-vector parsing, slugify)

use pyo3::prelude::*;

mod body;
mod convert;
mod integrator;
mod perturbation;
mod stats;
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
    stats::register(m)?;
    system::register(m)?;
    trajectory::register(m)?;
    units::register(m)?;
    Ok(())
}
