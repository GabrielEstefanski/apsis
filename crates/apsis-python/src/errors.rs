//! Typed Python exceptions for apsis core domain errors.
//!
//! Exposes [`UnitSystemMismatchError`] — raised when
//! `System.add_*_perturbation` receives an operator built for a
//! different `UnitSystem` than the `System`. The Rust core returns
//! [`apsis::physics::integrator::UnitSystemMismatch`]; this module
//! converts it into a typed exception callers can `except`.

use apsis::physics::integrator::UnitSystemMismatch;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(
    apsis,
    UnitSystemMismatchError,
    PyException,
    "Operator's UnitSystem disagrees with the System's UnitSystem at \
     registration time. Carries `operator`, `operator_units`, and \
     `system_units` attributes for programmatic recovery; the message \
     is the human-readable Display impl from the Rust core."
);

/// Convert the Rust [`UnitSystemMismatch`] into a [`UnitSystemMismatchError`]
/// with structured fields preserved as Python attributes (alongside the
/// human-readable message).
pub(crate) fn unit_system_mismatch_to_pyerr(e: UnitSystemMismatch) -> PyErr {
    Python::with_gil(|py| {
        let err = UnitSystemMismatchError::new_err(format!("{e}"));
        // Best-effort attribute attachment: walk the live exception
        // instance and stash the structured fields. If anything fails
        // (e.g. an exception class without writable attrs) the message
        // alone still carries the diagnosis, so the user-facing failure
        // mode is "exception with full message" rather than "no error".
        let instance = err.value(py).clone();
        let _ = instance.setattr("operator", e.operator);
        let _ = instance.setattr("operator_units", format!("{}", e.operator_units));
        let _ = instance.setattr("system_units", format!("{}", e.system_units));
        err
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("UnitSystemMismatchError", m.py().get_type::<UnitSystemMismatchError>())?;
    Ok(())
}
