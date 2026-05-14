//! Python binding for [`apsis_1pn`].
//!
//! A `Box<dyn HamiltonianOperator>` constructed in [`apsis_1pn`] crosses
//! into Python via a typed [`PyCapsule`](pyo3::types::PyCapsule)
//! (transport defined in [`apsis_py_core`]), travels in the pure-Python
//! `apsis.Perturbation` wrapper, and is unwrapped at
//! `System.add_hamiltonian_perturbation` back into Rust. The 1PN formula
//! itself lives in [`apsis_1pn`]; this crate is plumbing only.
//!
//! Each factory below is a one-liner
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
//! `Body.<material>(...)`; a violation emits a structured
//! warning at registration.
//!
//! # Use
//!
//! ```python
//! import apsis
//! import apsis_1pn
//!
//! sun = apsis.Body.star(mass=1.0)
//! mercury = (apsis.Body.rocky(mass=1.66e-7)
//!            .at((0.387, 0.0))
//!            .with_velocity((0.0, 1.61))
//!            )
//!
//! sys = apsis.System(
//!     bodies=[sun, mercury],
//!     units=apsis.units.SOLAR_CANONICAL,
//!     integrator="ias15",
//!     dt=1e-3,
//!     exact_gravity=True,
//! )
//! # Same UnitSystem on both sides — registration check passes.
//! sys.add_hamiltonian_perturbation(
//!     apsis_1pn.PostNewtonian1PN.for_units(units=apsis.units.SOLAR_CANONICAL),
//! )
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
/// to `apsis.System(...)` or call `Body.<material>(...)`
/// on every body.
#[pyclass(module = "apsis_1pn", name = "PostNewtonian1PN")]
pub struct PyPostNewtonian1PN;

/// Duck-type extraction of an `apsis.UnitSystem` Python object into the
/// Rust [`apsis::units::UnitSystem`]. apsis-1pn-py does not depend on
/// apsis-py directly (each is its own cdylib), so we reach the L/T/M
/// scales through attribute access and reconstruct via `UnitSystem::custom`.
///
/// `length_scale_si` / `time_scale_si` / `mass_scale_si` are exposed by
/// `apsis.UnitSystem` as `#[getter]` properties (not methods), so this
/// uses `getattr` rather than `call_method0` — the latter would resolve
/// to `getattr(...)()` and fail with `'float' object is not callable`.
fn unit_system_from_python(units: &Bound<'_, PyAny>) -> PyResult<apsis::units::UnitSystem> {
    let l = units.getattr("length_scale_si")?.extract::<f64>()?;
    let t = units.getattr("time_scale_si")?.extract::<f64>()?;
    let m = units.getattr("mass_scale_si")?.extract::<f64>()?;
    apsis::units::UnitSystem::custom(l, t, m).map_err(|e| {
        PyValueError::new_err(format!(
            "units: failed to construct UnitSystem from the supplied object: {e}"
        ))
    })
}

#[pymethods]
impl PyPostNewtonian1PN {
    // ── Named-regime constructor (Pattern A) ──────────────────────────────────

    /// 1PN with `c` derived from the supplied `apsis.UnitSystem`. The
    /// recommended constructor — pass the same `UnitSystem` used to
    /// build the `System`, and the relativistic correction stays
    /// consistent with the rest of the integration. The `System`
    /// registration check panics on mismatch, so unit-system confusion
    /// cannot survive `add_hamiltonian_perturbation`.
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

    // ── Raw escape ────────────────────────────────────────────────────────────

    /// 1PN with an explicit `c` value, pinned to the supplied
    /// `apsis.UnitSystem`. No cross-check between `c` and `units` —
    /// `c` is taken as given. The `System` registration check still
    /// validates that `units` matches the `System`'s own unit system
    /// and panics on mismatch.
    ///
    /// Use when `c` is computed by neighbouring code (so cross-checking
    /// against `units` is redundant), or for hypothetical experiments
    /// where `c` is intentionally non-physical. Prefer
    /// [`for_units`](Self::for_units) for normal physics.
    #[staticmethod]
    #[pyo3(signature = (*, c, units))]
    fn from_raw_c(py: Python<'_>, c: f64, units: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if !c.is_finite() || c <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "c: expected a strictly positive finite float, got {c}"
            )));
        }
        let rust_units = unit_system_from_python(units)?;
        wrap_in_apsis_perturbation(
            py,
            Box::new(PostNewtonian1PN::from_raw_c(c, rust_units)),
            &format!("PostNewtonian1PN(c={c})"),
        )
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
