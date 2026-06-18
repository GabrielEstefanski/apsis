//! `apsis._native.central` — central-potential perturbations
//! (Tamayo 2020, observable-inversion), wrapping [`apsis_central::CentralForce`].

use apsis_central::CentralForce;
use apsis_py_core::extract_unit_system;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::wrap_in_apsis_perturbation;
use crate::body::PyBody;

/// Central-potential perturbation parameterised by `(a_central, γ)`.
/// Constructed via the named factories below; each returns an
/// `apsis.Perturbation` ready for `System.add_hamiltonian_perturbation`.
#[pyclass(module = "apsis.central", name = "CentralForce")]
pub struct PyCentralForce;

#[pymethods]
impl PyCentralForce {
    /// Construct from explicit `a_central` and `γ`. `source` is the
    /// body whose central potential is modified.
    #[staticmethod]
    #[pyo3(signature = (*, source, a_central, gamma, units))]
    fn from_raw(
        py: Python<'_>,
        source: usize,
        a_central: f64,
        gamma: f64,
        units: &Bound<'_, PyAny>,
    ) -> PyResult<PyObject> {
        let rust_units = extract_unit_system(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(CentralForce::from_raw(source, a_central, gamma, rust_units)),
            "CentralForce(from_raw)",
        )
    }

    /// Construct by inverting a desired apsidal rate. `bodies` is the
    /// body sequence the `System` will be built from (used to compute
    /// the target's mean motion); `target` must be on a bound orbit
    /// around `source`. Accepts any iterable of `Body`.
    #[staticmethod]
    #[pyo3(signature = (*, source, target, omega_dot, gamma, bodies, units))]
    fn from_apsidal_rate(
        py: Python<'_>,
        source: usize,
        target: usize,
        omega_dot: f64,
        gamma: f64,
        bodies: &Bound<'_, PyAny>,
        units: &Bound<'_, PyAny>,
    ) -> PyResult<PyObject> {
        let rust_units = extract_unit_system(units)?;
        let rust_bodies: Vec<apsis::domain::body::Body> = bodies
            .try_iter()?
            .map(|item| Ok::<_, PyErr>(item?.extract::<PyRef<'_, PyBody>>()?.inner.clone()))
            .collect::<PyResult<_>>()?;
        let force = CentralForce::from_apsidal_rate(
            source,
            target,
            omega_dot,
            gamma,
            &rust_bodies,
            rust_units,
        )
        .map_err(|e| PyValueError::new_err(format!("CentralForce.from_apsidal_rate: {e}")))?;
        wrap_in_apsis_perturbation(py, Box::new(force), "CentralForce(from_apsidal_rate)")
    }
}

pub(crate) fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent.py();
    let m = PyModule::new(py, "central")?;
    m.add_class::<PyCentralForce>()?;
    parent.add_submodule(&m)?;
    py.import("sys")?.getattr("modules")?.set_item("apsis._native.central", &m)?;
    Ok(())
}
