//! Python binding for [`apsis_1pn`].
//!
//! **This crate proves that the apsis perturbation extension model is
//! preserved across both Rust and Python boundaries** — without
//! duplicating physics, breaking ownership semantics, or requiring
//! kernel modification. A `Box<dyn HamiltonianOperator>` constructed
//! in [`apsis_1pn`] crosses into Python via a typed
//! [`PyCapsule`](pyo3::types::PyCapsule) (transport defined in
//! [`apsis_py_core`]), travels in the pure-Python `apsis.Perturbation`
//! wrapper, and is unwrapped at `System.add_hamiltonian_perturbation`
//! back into Rust. The 1PN formula itself is implemented exactly once,
//! in [`apsis_1pn`]; this crate is plumbing only.
//!
//! Treat this crate as the **template** when writing Python bindings
//! for new perturbation crates: every factory below is a one-liner
//! built on [`apsis_py_core::box_into_capsule`].
//!
//! See [`README`](https://github.com/gabrielbragaestefanski/apsis/tree/master/crates/apsis-1pn-py)
//! for the full extension-contract specification, the critical kernel
//! precondition, and the rationale.
//!
//! # ⚠ Critical precondition
//!
//! Attaching 1PN to a softened-gravity system **invalidates the
//! physical model**. For Mercury-like orbits, the numerical apsidal
//! precession from Plummer softening alone is ~2 × 10³ larger than
//! the relativistic signal *and inverts its sign* — energy and
//! angular momentum stay conserved at machine precision while the
//! trajectory is physically wrong. **This is not a numerical error —
//! it is a model violation.** Pass `exact_gravity=True` or call
//! `Body.<material>(...).unsoftened()`; a violation emits a structured
//! warning at registration.
//!
//! # Use
//!
//! ```python
//! import apsis
//! import apsis_1pn
//!
//! sun = apsis.Body.star(mass=1.0).unsoftened()
//! mercury = (apsis.Body.rocky(mass=1.66e-7)
//!            .at((0.387, 0.0))
//!            .with_velocity((0.0, 1.61))
//!            .unsoftened())
//!
//! sys = apsis.System(
//!     bodies=[sun, mercury],
//!     units=apsis.units.SOLAR,
//!     integrator="ias15",
//!     dt=1e-3,
//!     exact_gravity=True,
//! )
//! sys.add_hamiltonian_perturbation(apsis_1pn.PostNewtonian1PN.solar_units())
//! ```

use apsis_1pn::PostNewtonian1PN;
use apsis_py_core::box_into_capsule;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Construct an `apsis.Perturbation` instance from a freshly-built
/// boxed operator. Imports `apsis.Perturbation` once per call; the
/// cost is negligible compared to the enclosing simulation work.
fn wrap_in_apsis_perturbation(
    py: Python<'_>,
    inner: Box<dyn apsis::physics::integrator::HamiltonianOperator>,
    label: &str,
) -> PyResult<PyObject> {
    let capsule = box_into_capsule(py, inner)?;
    let apsis = py.import("apsis")?;
    let perturbation_cls = apsis.getattr("Perturbation")?;
    let result = perturbation_cls.call1((capsule, label))?;
    Ok(result.into())
}

/// First post-Newtonian gravitational correction (Schwarzschild,
/// test-particle form, applied pairwise).
///
/// Use the named factories to construct an instance in the appropriate
/// unit system; the result is a fully-formed `apsis.Perturbation` ready
/// to attach via `System.add_hamiltonian_perturbation(...)`.
///
/// # Kernel preconditions
///
/// 1PN is derived around the bit-exact Newtonian potential. Attaching
/// it to a softened-gravity system substitutes a different unperturbed
/// Hamiltonian whose apsidal precession alone is ~2 × 10³ larger than
/// the 1PN signal for a Mercury-like orbit, silently inverting the
/// sign of the measured precession. Either pass `exact_gravity=True`
/// to `apsis.System(...)` or call `Body.<material>(...).unsoftened()`
/// on every body.
#[pyclass(module = "apsis_1pn", name = "PostNewtonian1PN")]
pub struct PyPostNewtonian1PN;

/// Duck-type extraction of an `apsis.UnitSystem` Python object into the
/// Rust [`apsis::units::UnitSystem`]. apsis-1pn-py does not depend on
/// apsis-py directly (each is its own cdylib), so we reach the L/T/M
/// scales through method calls and reconstruct via `UnitSystem::custom`.
fn unit_system_from_python(units: &Bound<'_, PyAny>) -> PyResult<apsis::units::UnitSystem> {
    let l = units.call_method0("length_scale_si")?.extract::<f64>()?;
    let t = units.call_method0("time_scale_si")?.extract::<f64>()?;
    let m = units.call_method0("mass_scale_si")?.extract::<f64>()?;
    apsis::units::UnitSystem::custom(l, t, m).map_err(|e| {
        PyValueError::new_err(format!(
            "units: failed to construct UnitSystem from the supplied object: {e}"
        ))
    })
}

#[pymethods]
impl PyPostNewtonian1PN {
    // ── Named-regime constructors (Pattern A) ─────────────────────────────────

    /// 1PN calibrated for the simulator's canonical solar-system units
    /// (`G = 1`, length = AU, mass = M_sun, time chosen so that
    /// `c = `[`apsis_1pn::C_SOLAR_UNITS`]).
    #[staticmethod]
    fn solar_units(py: Python<'_>) -> PyResult<PyObject> {
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::solar_units()),
            "PostNewtonian1PN(solar_units)",
        )
    }

    /// 1PN with `c` derived from a supplied `apsis.UnitSystem`. The
    /// recommended path for non-solar unit choices — the user picks
    /// the unit system, `c` is computed exactly so the relativistic
    /// correction stays consistent with the rest of the integration.
    #[staticmethod]
    #[pyo3(signature = (*, units))]
    fn for_units(py: Python<'_>, units: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        let rust_units = unit_system_from_python(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::for_units(rust_units)),
            "PostNewtonian1PN(for_units)",
        )
    }

    // ── Raw escape (with optional validation) ─────────────────────────────────

    /// 1PN with an explicit speed of light in the caller's unit system.
    /// No validation — use [`from_raw_c_validated`](Self::from_raw_c_validated)
    /// for a cross-checked construction, or [`for_units`](Self::for_units)
    /// to skip raw input entirely.
    #[staticmethod]
    #[pyo3(signature = (*, c))]
    fn from_raw_c(py: Python<'_>, c: f64) -> PyResult<PyObject> {
        if !c.is_finite() || c <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "c: expected a strictly positive finite float, got {c}"
            )));
        }
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::from_raw_c(c)),
            &format!("PostNewtonian1PN(c={c})"),
        )
    }

    /// 1PN with an explicit `c`, cross-checked against the `c` derived
    /// from a supplied `apsis.UnitSystem`. Raises `ValueError` when the
    /// relative error exceeds `1e-9` — the protection against silent
    /// unit-mismatch errors when `c` comes from an external source.
    #[staticmethod]
    #[pyo3(signature = (*, c, units))]
    fn from_raw_c_validated(
        py: Python<'_>,
        c: f64,
        units: &Bound<'_, PyAny>,
    ) -> PyResult<PyObject> {
        if !c.is_finite() || c <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "c: expected a strictly positive finite float, got {c}"
            )));
        }
        let rust_units = unit_system_from_python(units)?;
        let pn = PostNewtonian1PN::from_raw_c_validated(c, rust_units)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        wrap_in_apsis_perturbation(py, Box::new(pn), &format!("PostNewtonian1PN(c={c} validated)"))
    }
}

/// `apsis_1pn._native`: the Rust-built extension module.
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<PyPostNewtonian1PN>()?;
    m.add("C_SOLAR_UNITS", apsis_1pn::C_SOLAR_UNITS)?;
    Ok(())
}
