//! `apsis._native.radiation` — radiation pressure (Burns et al. 1979),
//! wrapping [`apsis_radiation::RadiationPressure`].

use apsis_py_core::extract_unit_system;
use apsis_radiation::RadiationPressure;
use pyo3::prelude::*;

use super::wrap_in_apsis_perturbation;

/// Radiation pressure as a fractional reduction of central gravity per
/// receiver, indexed by `β` (Burns 1979). Returns an
/// `apsis.Perturbation` ready for `System.add_hamiltonian_perturbation`.
#[pyclass(module = "apsis.radiation", name = "RadiationPressure")]
pub struct PyRadiationPressure;

#[pymethods]
impl PyRadiationPressure {
    /// Construct from an explicit `β` per body. `source` is the
    /// radiating body's index; `betas[i]` is the β applied to body `i`
    /// (set `betas[source] = 0.0`). `units` pins the operator to a
    /// `UnitSystem` for the `System` registration check.
    #[staticmethod]
    #[pyo3(signature = (*, source, betas, units))]
    fn from_raw_betas(
        py: Python<'_>,
        source: usize,
        betas: Vec<f64>,
        units: &Bound<'_, PyAny>,
    ) -> PyResult<PyObject> {
        let rust_units = extract_unit_system(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(RadiationPressure::from_raw_betas(source, betas, rust_units)),
            "RadiationPressure(from_raw_betas)",
        )
    }
}

pub(crate) fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent.py();
    let m = PyModule::new(py, "radiation")?;
    m.add_class::<PyRadiationPressure>()?;
    parent.add_submodule(&m)?;
    py.import("sys")?.getattr("modules")?.set_item("apsis._native.radiation", &m)?;
    Ok(())
}
