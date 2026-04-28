//! Perturbation transport — extracts a `Box<dyn PerturbationForce>` out
//! of an `apsis.Perturbation` Python object's underlying `PyCapsule`.
//!
//! The user-facing `Perturbation` class is defined in
//! `python/apsis/__init__.py` (pure Python). This module is just the
//! `add_perturbation` boundary: pull the capsule attribute, hand it to
//! [`apsis_py_core::take_box_from_capsule`], forward the boxed trait
//! object to the core. See `apsis-py-core` for the rationale on why
//! the wrapper is pure-Python rather than a pyclass.

use apsis::physics::integrator::PerturbationForce;
use apsis_py_core::take_box_from_capsule;
use pyo3::prelude::*;
use pyo3::types::PyCapsule;

use crate::convert::value_error;

/// Pull the boxed perturbation out of a Python object the user passed
/// to `System.add_perturbation`. Accepts anything carrying a
/// `_capsule` attribute that holds a valid apsis perturbation capsule
/// — the canonical case is the `apsis.Perturbation` class defined in
/// `python/apsis/__init__.py`.
pub(crate) fn take_perturbation_from_python(
    p: &Bound<'_, PyAny>,
) -> PyResult<Box<dyn PerturbationForce>> {
    let capsule_attr = p.getattr("_capsule").map_err(|_| {
        value_error(
            "perturbation",
            "expected an apsis.Perturbation instance (object with a `_capsule` attribute)",
        )
    })?;
    let capsule: &Bound<'_, PyCapsule> = capsule_attr.downcast::<PyCapsule>().map_err(|_| {
        value_error(
            "perturbation",
            "the `_capsule` attribute is not a PyCapsule; this Perturbation was not built \
             by an apsis-compatible perturbation crate",
        )
    })?;
    take_box_from_capsule(capsule)
}
