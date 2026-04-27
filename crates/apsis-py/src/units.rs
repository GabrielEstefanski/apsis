//! PyO3 wrapper for [`apsis::units::UnitSystem`].
//!
//! The Python surface mirrors the Rust core bit-for-bit: the named
//! factories surface as module-level singletons (`apsis.units.SOLAR`,
//! `apsis.units.SI`, ...), and a [`UnitSystem`] class exposes the
//! same accessors and explicit conversions. Construction of arbitrary
//! unit systems goes through the `custom()` classmethod, which forwards
//! to [`apsis::units::UnitSystem::custom`] and translates a
//! [`UnitError`] into a `ValueError` at the FFI boundary.
//!
//! # Façade-only
//!
//! Nothing here implements unit logic. The conversions, the derived
//! `g`, the SI constants, and the validation rules all live in
//! [`apsis::units`]; this file is a translation layer.
//!
//! # Why module-level singletons (not a `UnitsKind` enum)
//!
//! `apsis.units.SOLAR` is a `UnitSystem` instance, not an enum
//! variant. A researcher who imports `from apsis.units import SOLAR`
//! gets a fully-formed object with `.g`, `.length_scale_si`, and
//! every conversion method ready to call — no enum-to-instance
//! resolution step in the way. The downside (allocating one
//! `UnitSystem` per named factory at import time) is trivial since
//! the type is `Copy` and 3 × `f64` + 3 × `&'static str` per slot.

use apsis::units::{UnitError as CoreUnitError, UnitSystem as CoreUnitSystem};
use pyo3::prelude::*;

use crate::convert::value_error;

/// Closed system of units for length, time, and mass.
///
/// Construct via one of the named factories (:meth:`solar`, :meth:`si`,
/// :meth:`canonical`, :meth:`henon`, :meth:`cgs`) or :meth:`custom`.
/// Once chosen, a ``UnitSystem`` is immutable — there is no setter
/// for any of its fields, and every `.replace_*` style call would
/// return a new instance rather than mutate.
///
/// All physics inputs (positions, velocities, masses, dt) passed to a
/// :class:`System` constructed against this unit system are interpreted
/// in its canonical units. The wrapper performs no dimensional checking
/// — passing a value in the wrong unit is a silent physical error,
/// matching REBOUND's convention.
#[pyclass(module = "apsis", name = "UnitSystem", frozen)]
#[derive(Clone, Copy)]
pub(crate) struct PyUnitSystem {
    pub(crate) inner: CoreUnitSystem,
}

impl PyUnitSystem {
    /// Internal constructor used by the named-factory module attributes
    /// and by [`crate::system::PySystem`] when accepting a Python
    /// `UnitSystem` value.
    pub(crate) fn from_core(inner: CoreUnitSystem) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyUnitSystem {
    /// Hénon-style canonical N-body units: ``G = 1`` by construction.
    ///
    /// Length and time are nominally ``1`` in SI; the mass scale absorbs
    /// the ``1/G_SI`` factor required to make the derived ``G`` exactly
    /// one. This is the implicit default of REBOUND when no units are
    /// specified, and the convention of stellar-dynamics literature.
    #[staticmethod]
    fn canonical() -> Self {
        Self::from_core(CoreUnitSystem::canonical())
    }

    /// Alias for :meth:`canonical` using the literature name. The
    /// returned object is identical; pick the name that reads best
    /// in the surrounding code.
    #[staticmethod]
    fn henon() -> Self {
        Self::from_core(CoreUnitSystem::henon())
    }

    /// SI units: metre, second, kilogram. ``G ≈ 6.674e-11``.
    #[staticmethod]
    fn si() -> Self {
        Self::from_core(CoreUnitSystem::si())
    }

    /// Solar-system canonical units: AU, year, solar mass.
    /// Derived ``G ≈ 39.478`` (the IAU approximation of ``4π²``).
    #[staticmethod]
    fn solar() -> Self {
        Self::from_core(CoreUnitSystem::solar())
    }

    /// CGS units: centimetre, second, gram.
    /// ``G ≈ 6.674e-8`` (cm³ g⁻¹ s⁻²).
    #[staticmethod]
    fn cgs() -> Self {
        Self::from_core(CoreUnitSystem::cgs())
    }

    /// Build a custom unit system from explicit SI scales.
    ///
    /// All three arguments must be strictly positive and finite.
    /// Zero, negative, infinite, or NaN values raise ``ValueError``
    /// at the boundary because they would otherwise produce a
    /// non-finite ``g`` that explodes inside the integrator.
    #[staticmethod]
    #[pyo3(signature = (*, length_m, time_s, mass_kg))]
    fn custom(length_m: f64, time_s: f64, mass_kg: f64) -> PyResult<Self> {
        match CoreUnitSystem::custom(length_m, time_s, mass_kg) {
            Ok(u) => Ok(Self::from_core(u)),
            Err(CoreUnitError::InvalidLength(v)) => Err(value_error(
                "length_m",
                format!("expected a strictly positive finite float, got {v}"),
            )),
            Err(CoreUnitError::InvalidTime(v)) => Err(value_error(
                "time_s",
                format!("expected a strictly positive finite float, got {v}"),
            )),
            Err(CoreUnitError::InvalidMass(v)) => Err(value_error(
                "mass_kg",
                format!("expected a strictly positive finite float, got {v}"),
            )),
        }
    }

    // ── Scale accessors ──────────────────────────────────────────────────

    /// SI metres per code-unit length. ``solar().length_scale_si``
    /// returns one astronomical unit in metres.
    #[getter]
    fn length_scale_si(&self) -> f64 {
        self.inner.length_scale_si()
    }

    /// SI seconds per code-unit time. ``solar().time_scale_si``
    /// returns one Julian year in seconds.
    #[getter]
    fn time_scale_si(&self) -> f64 {
        self.inner.time_scale_si()
    }

    /// SI kilograms per code-unit mass. ``solar().mass_scale_si``
    /// returns one solar mass in kilograms.
    #[getter]
    fn mass_scale_si(&self) -> f64 {
        self.inner.mass_scale_si()
    }

    /// Newtonian gravitational constant in this system's canonical
    /// units. Always derived from ``G_SI · M · T² / L³`` — never read
    /// from a hardcoded literature value.
    #[getter]
    fn g(&self) -> f64 {
        self.inner.g()
    }

    // ── Display labels ───────────────────────────────────────────────────

    /// Display symbol for the length axis (``"AU"``, ``"m"``, ``"cm"``, ...).
    #[getter]
    fn length_label(&self) -> &'static str {
        self.inner.length_label()
    }

    /// Display symbol for the time axis (``"yr"``, ``"s"``, ...).
    #[getter]
    fn time_label(&self) -> &'static str {
        self.inner.time_label()
    }

    /// Display symbol for the mass axis (``"Msun"``, ``"kg"``, ``"g"``, ...).
    #[getter]
    fn mass_label(&self) -> &'static str {
        self.inner.mass_label()
    }

    // ── Explicit conversions ─────────────────────────────────────────────

    /// Convert a length expressed in this system's canonical units to SI metres.
    fn length_to_si(&self, x: f64) -> f64 {
        self.inner.length_to_si(x)
    }

    /// Convert a length expressed in SI metres to this system's canonical units.
    fn length_from_si(&self, x: f64) -> f64 {
        self.inner.length_from_si(x)
    }

    /// Convert a duration expressed in this system's canonical units to SI seconds.
    fn time_to_si(&self, x: f64) -> f64 {
        self.inner.time_to_si(x)
    }

    /// Convert a duration expressed in SI seconds to this system's canonical units.
    fn time_from_si(&self, x: f64) -> f64 {
        self.inner.time_from_si(x)
    }

    /// Convert a mass expressed in this system's canonical units to SI kilograms.
    fn mass_to_si(&self, x: f64) -> f64 {
        self.inner.mass_to_si(x)
    }

    /// Convert a mass expressed in SI kilograms to this system's canonical units.
    fn mass_from_si(&self, x: f64) -> f64 {
        self.inner.mass_from_si(x)
    }

    // ── Equality & display ───────────────────────────────────────────────
    //
    // `PartialEq` on the core type compares only SI scales (labels are
    // metadata) — `__eq__` mirrors that contract: two `solar()`
    // instances are equal, and a `custom()` with identical scales is
    // also equal to `solar()` even though their labels differ.

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        // Hash the three SI scales as bit patterns. Same equivalence
        // class as `__eq__` (labels ignored). NaN scales are rejected
        // at construction so `is_nan` is impossible here.
        let mut h: u64 = 0;
        for v in [
            self.inner.length_scale_si(),
            self.inner.time_scale_si(),
            self.inner.mass_scale_si(),
        ] {
            h = h.wrapping_mul(0x100000001b3).wrapping_add(v.to_bits());
        }
        h
    }

    fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}

/// Register the `UnitSystem` class and the named-factory module-level
/// singletons (`apsis._native.units.SOLAR`, ...) so a Python user
/// writes `apsis.units.SOLAR` rather than `apsis.UnitSystem.solar()`.
///
/// Each singleton is a `Py<PyUnitSystem>` constructed at module-init
/// time. They're cheap (`Copy` underlying type) and stable across the
/// lifetime of the interpreter.
pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let units_module = PyModule::new(py, "units")?;
    units_module.add_class::<PyUnitSystem>()?;

    units_module.add("CANONICAL", PyUnitSystem::canonical())?;
    units_module.add("HENON", PyUnitSystem::henon())?;
    units_module.add("SI", PyUnitSystem::si())?;
    units_module.add("SOLAR", PyUnitSystem::solar())?;
    units_module.add("CGS", PyUnitSystem::cgs())?;

    // SI constants — single source of truth, mirrored from the Rust core.
    units_module.add("G_SI", apsis::units::G_SI)?;
    units_module.add("AU_M", apsis::units::AU_M)?;
    units_module.add("YR_S", apsis::units::YR_S)?;
    units_module.add("MSUN_KG", apsis::units::MSUN_KG)?;
    units_module.add("CM_M", apsis::units::CM_M)?;
    units_module.add("G_KG", apsis::units::G_KG)?;

    m.add_submodule(&units_module)?;
    // Also expose the class at the top level so `apsis.UnitSystem`
    // works for direct-construction call sites.
    m.add_class::<PyUnitSystem>()?;
    Ok(())
}
