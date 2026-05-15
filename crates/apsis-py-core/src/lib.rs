//! Capsule transport for `Box<dyn HamiltonianOperator>` across the
//! Python FFI boundary, plus duck-type extractors for apsis types.

use std::ffi::CStr;
use std::sync::Mutex;

use apsis::physics::integrator::HamiltonianOperator;
use apsis::units::UnitSystem;
use pyo3::exceptions::PyValueError;
use pyo3::ffi::c_str;
use pyo3::prelude::*;
use pyo3::types::{PyCapsule, PyCapsuleMethods};

const CAPSULE_NAME: &CStr = c_str!("apsis_perturbation_box_v3");

type CapsulePayload = Mutex<Option<Box<dyn HamiltonianOperator>>>;

/// Wrap a freshly-built operator into a `PyCapsule`.
/// Ownership transfers to Python.
pub fn box_into_capsule(
    py: Python<'_>,
    inner: Box<dyn HamiltonianOperator>,
) -> PyResult<Bound<'_, PyCapsule>> {
    let payload: CapsulePayload = Mutex::new(Some(inner));
    PyCapsule::new(py, payload, Some(CAPSULE_NAME.to_owned()))
}

/// Extract the boxed operator from a capsule. Single-consume.
pub fn take_box_from_capsule(
    capsule: &Bound<'_, PyCapsule>,
) -> PyResult<Box<dyn HamiltonianOperator>> {
    if capsule.name()?.map(|n| n != CAPSULE_NAME).unwrap_or(true) {
        return Err(PyValueError::new_err(
            "perturbation: capsule type tag does not match apsis_perturbation_box_v3",
        ));
    }

    // SAFETY: capsule constructed with `CapsulePayload` (verified via
    // type-tag check above) and the payload is `Send + Sync + 'static`.
    let payload: &CapsulePayload = unsafe { capsule.reference::<CapsulePayload>() };
    let mut guard = payload
        .lock()
        .map_err(|_| PyValueError::new_err("perturbation: capsule mutex poisoned"))?;

    guard.take().ok_or_else(|| {
        PyValueError::new_err(
            "perturbation: this Perturbation has already been attached to a System; \
             construct a fresh instance for each add_perturbation call",
        )
    })
}

/// Extract an apsis [`UnitSystem`] from a Python object exposing
/// `length_scale_si`, `time_scale_si`, `mass_scale_si` attributes.
pub fn extract_unit_system(units: &Bound<'_, PyAny>) -> PyResult<UnitSystem> {
    let l = units.getattr("length_scale_si")?.extract::<f64>()?;
    let t = units.getattr("time_scale_si")?.extract::<f64>()?;
    let m = units.getattr("mass_scale_si")?.extract::<f64>()?;
    UnitSystem::custom(l, t, m).map_err(|e| {
        PyValueError::new_err(format!("units: failed to construct UnitSystem: {e}"))
    })
}
