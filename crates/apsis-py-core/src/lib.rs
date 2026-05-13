//! Shared FFI helpers for the apsis Python binding crates.
//!
//! Cross-extension type sharing in PyO3 is fragile: when two `cdylib`s
//! both register the same `#[pyclass]`, each gets its own Python class
//! object and `isinstance` checks fail at the boundary. The workaround
//! used here is the canonical CPython one — transport an opaque
//! `Box<dyn HamiltonianOperator>` via `PyCapsule` (a typed C pointer
//! Python knows how to lifecycle-manage) and present it to the user
//! as a pure-Python `Perturbation` class defined once in
//! `apsis/__init__.py`. The capsule lives inside the Python wrapper;
//! the Rust types here are just the safe enter / safe exit hooks.
//!
//! Crate layout: `rlib`-only on purpose, no `#[pymodule]`, so multiple
//! `cdylib`s can link the helpers without colliding on
//! `PyInit__native`.
//!
//! Non-conservative operators are not yet exposed across the FFI; when
//! they are, a parallel `_nc_v1` capsule will live alongside this one.

use std::ffi::CStr;
use std::sync::Mutex;

use apsis::physics::integrator::HamiltonianOperator;
use pyo3::exceptions::PyValueError;
use pyo3::ffi::c_str;
use pyo3::prelude::*;
use pyo3::types::{PyCapsule, PyCapsuleMethods};

/// Capsule type tag. Bumping the suffix is the breaking-change marker
/// when the underlying transport contract changes — every consumer
/// has to recompile against the new name.
///
/// - `_v2`: `accumulate` migrated from `&mut [(f64, f64)]` to
///   `&mut [apsis::math::Vec3]` for the 3D port.
/// - `_v3`: payload narrowed from the unified `PerturbationForce` trait
///   to the conservation-segregated `HamiltonianOperator` trait
///   (force + closed-form Hamiltonian); non-conservative operators
///   reserved for a future parallel capsule.
const CAPSULE_NAME: &CStr = c_str!("apsis_perturbation_box_v3");

/// Capsule payload. The `Mutex<Option<...>>` shape gives us:
/// - `Send + Sync` so PyO3's `PyCapsule::new` accepts it
/// - single-consume semantics via `Option::take`
/// - safe access from the Rust side (capsule contents are `'static`)
type CapsulePayload = Mutex<Option<Box<dyn HamiltonianOperator>>>;

/// Wrap a freshly-built `Box<dyn HamiltonianOperator>` into a `PyCapsule`
/// the user-facing `apsis.Perturbation` Python class consumes.
///
/// Ownership transfers to Python: if the capsule is dropped without
/// being passed through [`take_box_from_capsule`], the inner box is
/// freed by the capsule's destructor.
pub fn box_into_capsule(
    py: Python<'_>,
    inner: Box<dyn HamiltonianOperator>,
) -> PyResult<Bound<'_, PyCapsule>> {
    let payload: CapsulePayload = Mutex::new(Some(inner));
    PyCapsule::new(py, payload, Some(CAPSULE_NAME.to_owned()))
}

/// Pull the boxed operator back out of a capsule produced by
/// [`box_into_capsule`]. Single-consume: subsequent calls on the same
/// capsule return `Err` so a Python user who reuses the same
/// `Perturbation` value across `System.add_perturbation` calls sees a
/// clear error rather than a use-after-free.
pub fn take_box_from_capsule(
    capsule: &Bound<'_, PyCapsule>,
) -> PyResult<Box<dyn HamiltonianOperator>> {
    if capsule.name()?.map(|n| n != CAPSULE_NAME).unwrap_or(true) {
        return Err(PyValueError::new_err(
            "perturbation: capsule type tag does not match apsis_perturbation_box_v3; \
             this object did not come from an apsis-compatible perturbation crate \
             (rebuild against apsis ≥ 0.3 for the operator-split contract)",
        ));
    }

    // SAFETY: the capsule was constructed with `CapsulePayload` as the
    // payload type (verified via the type-tag check above) and the
    // payload is `Send + Sync + 'static`, so dereferencing it as `&CapsulePayload`
    // for the duration of the `take` call is sound.
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
