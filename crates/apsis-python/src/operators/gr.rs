//! `apsis._native.gr` — first post-Newtonian Schwarzschild correction
//! (Anderson et al. 1975), wrapping [`apsis_1pn::PostNewtonian1PN`].

use apsis_1pn::PostNewtonian1PN;
use apsis_py_core::extract_unit_system;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::wrap_in_apsis_perturbation;

/// First post-Newtonian Schwarzschild correction. Constructed via the
/// named factories below; each returns an `apsis.Perturbation` ready
/// for `System.add_hamiltonian_perturbation(...)`.
#[pyclass(module = "apsis.gr", name = "PostNewtonian1PN")]
pub struct PyPostNewtonian1PN;

#[pymethods]
impl PyPostNewtonian1PN {
    /// 1PN with `c` derived from the supplied `apsis.UnitSystem`.
    #[staticmethod]
    #[pyo3(signature = (*, units))]
    fn for_units(py: Python<'_>, units: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        let rust_units = extract_unit_system(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::for_units(rust_units)),
            "PostNewtonian1PN(for_units)",
        )
    }

    /// 1PN with an explicit `c` value, pinned to `units` for the
    /// `System` registration check.
    #[staticmethod]
    #[pyo3(signature = (*, c, units))]
    fn from_raw_c(py: Python<'_>, c: f64, units: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if !c.is_finite() || c <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "c: expected a strictly positive finite float, got {c}"
            )));
        }
        let rust_units = extract_unit_system(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::from_raw_c(c, rust_units)),
            &format!("PostNewtonian1PN(c={c})"),
        )
    }
}

pub(crate) fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let gr = PyModule::new(parent.py(), "gr")?;
    gr.add_class::<PyPostNewtonian1PN>()?;
    gr.add("C_SOLAR_UNITS", apsis_1pn::C_SOLAR_UNITS)?;
    parent.add_submodule(&gr)?;
    Ok(())
}
