//! Operator submodules registered under `apsis._native.<name>`.

use apsis::physics::integrator::HamiltonianOperator;
use apsis_py_core::box_into_capsule;
use pyo3::prelude::*;

#[cfg(feature = "gr")]
pub mod gr;

/// Wrap a freshly-built boxed operator into an `apsis.Perturbation`
/// instance ready for `System.add_hamiltonian_perturbation(...)`.
pub(crate) fn wrap_in_apsis_perturbation(
    py: Python<'_>,
    inner: Box<dyn HamiltonianOperator>,
    label: &str,
) -> PyResult<PyObject> {
    let capsule = box_into_capsule(py, inner)?;
    let apsis = py.import("apsis")?;
    let perturbation_cls = apsis.getattr("Perturbation")?;
    Ok(perturbation_cls.call1((capsule, label))?.into())
}
