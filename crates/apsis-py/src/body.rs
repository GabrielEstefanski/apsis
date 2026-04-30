//! Python-side wrapper of [`apsis::domain::body::Body`].
//!
//! The wrapper exposes a researcher-first API: nine material factories
//! (`Body.star`, `Body.rocky`, `Body.gas_giant`, ...) that mirror the
//! corresponding constructors in [`apsis::domain::body::Body`], a
//! kwargs-only signature on each factory so position and velocity
//! never depend on argument order, and an immutable fluent-builder
//! tail (`at`, `with_velocity`, `with_density`, `unsoftened`) for the
//! cases where chaining reads more naturally than a single call site.
//!
//! # Façade-only invariant
//!
//! Every `#[pymethods]` body in this file delegates to a single call on
//! [`apsis::domain::body::Body`] — the wrapper translates types at the
//! boundary and never composes physics. The set of valid materials,
//! the default softening rule, the density model, and the body-state
//! invariants are all owned by the core crate; this module is the
//! Python-shaped door into them.
//!
//! # Why builders return new bodies
//!
//! [`apsis::domain::body::Body`] is `Copy`, so each builder method here
//! produces a fresh `PyBody` rather than mutating the receiver. This
//! gives Python users the same value-semantics they get from
//! NumPy scalar operations or `dataclasses.replace`: chaining is safe,
//! aliasing is not surprising, and a body passed into a `System`
//! constructor is not retroactively mutated by later calls.
//!
//! ```text
//!   sun = Body.star(mass=1.0)
//!   far_sun = sun.at((10.0, 0.0))
//!   assert sun.position == (0.0, 0.0)         # unchanged
//!   assert far_sun.position == (10.0, 0.0)    # new instance
//! ```

use apsis::domain::body::Body as CoreBody;
use apsis::domain::materials::Material;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::convert::{value_error, xyz_triple};

/// Point-mass body with kinematics, mass, softening, and material class.
///
/// Bodies are constructed through one of the nine material factories
/// (`Body.star`, `Body.rocky`, ...) — each chooses sensible defaults for
/// density, softening, and any visual/physical properties tied to the
/// material class. Position and velocity always default to zero; pass
/// them as kwargs (`position=(x, y)`, `velocity=(vx, vy)`) at the
/// factory or via the fluent builder methods (`at`, `with_velocity`).
///
/// All builder methods return a fresh `Body` — bodies are immutable
/// by convention on the Python side, matching the value-semantics of
/// the underlying `apsis::domain::body::Body` (which is `Copy`). A
/// body handed to a `System` constructor is therefore not affected by
/// any subsequent builder call on its original handle.
///
/// Examples:
///
/// ```python
/// import apsis
///
/// # One-liner with kwargs
/// mercury = apsis.Body.rocky(
///     mass=3e-6,
///     position=(0.307, 0.0),
///     velocity=(0.0, 1.98),
/// )
///
/// # Equivalent fluent form
/// mercury = (apsis.Body.rocky(mass=3e-6)
///            .at((0.307, 0.0))
///            .with_velocity((0.0, 1.98)))
///
/// # Switch off softening (exact 1/r gravity for this body)
/// sun = apsis.Body.star(mass=1.0).unsoftened()
/// ```
#[pyclass(module = "apsis", name = "Body", frozen)]
#[derive(Clone, Copy)]
pub(crate) struct PyBody {
    pub(crate) inner: CoreBody,
}

impl PyBody {
    /// Build a `PyBody` of `material` with kwargs-driven kinematics.
    /// Single point of construction for every factory below; per-factory
    /// methods are one-liner wrappers that name a specific material.
    fn build(
        material: Material,
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        if !mass.is_finite() {
            return Err(value_error("mass", format!("expected a finite float, got {mass}")));
        }
        if mass <= 0.0 {
            return Err(value_error(
                "mass",
                format!("expected a strictly positive value, got {mass}"),
            ));
        }

        let mut inner = CoreBody::of(mass, material);

        if let Some(obj) = position {
            let (x, y, z) = xyz_triple("position", obj)?;
            inner.x = x;
            inner.y = y;
            inner.z = z;
        }
        if let Some(obj) = velocity {
            let (vx, vy, vz) = xyz_triple("velocity", obj)?;
            inner.vx = vx;
            inner.vy = vy;
            inner.vz = vz;
        }
        if let Some(eps) = softening {
            if !eps.is_finite() || eps < 0.0 {
                return Err(value_error(
                    "softening",
                    format!("expected a finite non-negative float, got {eps}"),
                ));
            }
            inner.softening = eps;
        }

        Ok(Self { inner })
    }
}

#[pymethods]
impl PyBody {
    // ── Material factories ───────────────────────────────────────────────

    /// Main-sequence luminous body. Default density and luminous material.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn star(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Star, mass, position, velocity, softening)
    }

    /// Brown dwarf — sub-stellar, deuterium-burning regime.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn brown_dwarf(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::BrownDwarf, mass, position, velocity, softening)
    }

    /// White dwarf — compact stellar remnant.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn white_dwarf(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::WhiteDwarf, mass, position, velocity, softening)
    }

    /// Gas giant — Jupiter-class hydrogen/helium envelope.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn gas_giant(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Gas, mass, position, velocity, softening)
    }

    /// Ice giant — Neptune-class water/methane envelope.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn ice_giant(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::IceGiant, mass, position, velocity, softening)
    }

    /// Rocky body — terrestrial planet or large rocky satellite.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn rocky(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Rocky, mass, position, velocity, softening)
    }

    /// Icy body — water-dominated composition (outer satellites, KBOs).
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn icy(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Icy, mass, position, velocity, softening)
    }

    /// Asteroid — rocky minor body.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn asteroid(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Asteroid, mass, position, velocity, softening)
    }

    /// Comet — volatile-rich minor body.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None, softening=None))]
    fn comet(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
        softening: Option<f64>,
    ) -> PyResult<Self> {
        Self::build(Material::Comet, mass, position, velocity, softening)
    }

    // ── Fluent builder ───────────────────────────────────────────────────
    // Each method takes a single 2-tuple (no positional-swap risk) and
    // returns a fresh `Body`; the Python view is value-typed.

    /// Place the body at `position = (x, y)` or `(x, y, z)`. Returns a
    /// new `Body`. A 2-element sequence is treated as planar input with
    /// `z = 0`.
    fn at(&self, position: &Bound<'_, PyAny>) -> PyResult<Self> {
        let (x, y, z) = xyz_triple("position", position)?;
        let mut inner = self.inner;
        inner.x = x;
        inner.y = y;
        inner.z = z;
        Ok(Self { inner })
    }

    /// Set the body's velocity to `(vx, vy)` or `(vx, vy, vz)`. Returns
    /// a new `Body`. A 2-element sequence is treated as planar input
    /// with `vz = 0`.
    fn with_velocity(&self, velocity: &Bound<'_, PyAny>) -> PyResult<Self> {
        let (vx, vy, vz) = xyz_triple("velocity", velocity)?;
        let mut inner = self.inner;
        inner.vx = vx;
        inner.vy = vy;
        inner.vz = vz;
        Ok(Self { inner })
    }

    /// Override the material-default density. Physical radius is
    /// recomputed from the new value at the call site (delegated to
    /// `apsis::domain::body::Body::with_density`). Returns a new `Body`.
    fn with_density(&self, density: f64) -> PyResult<Self> {
        if !density.is_finite() || density <= 0.0 {
            return Err(value_error(
                "density",
                format!("expected a strictly positive finite float, got {density}"),
            ));
        }
        Ok(Self { inner: self.inner.with_density(density) })
    }

    /// Drop this body's Plummer softening to zero, restoring the exact
    /// $1/r$ potential for every interaction it participates in. Use
    /// when measuring a deviation-from-Kepler signal whose magnitude
    /// is below the apsidal precession introduced by the default
    /// material-scaled softening (post-Newtonian, $J_2$ oblateness,
    /// tidal dissipation). Returns a new `Body`.
    fn unsoftened(&self) -> Self {
        Self { inner: self.inner.unsoftened() }
    }

    // ── Read-only properties ─────────────────────────────────────────────

    /// Body mass in simulation units.
    #[getter]
    fn mass(&self) -> f64 {
        self.inner.mass
    }

    /// Position as a 3-tuple $(x, y, z)$. For bodies confined to the
    /// `xy`-plane the third component is zero by default.
    #[getter]
    fn position(&self) -> (f64, f64, f64) {
        (self.inner.x, self.inner.y, self.inner.z)
    }

    /// Velocity as a 3-tuple $(v_x, v_y, v_z)$.
    #[getter]
    fn velocity(&self) -> (f64, f64, f64) {
        (self.inner.vx, self.inner.vy, self.inner.vz)
    }

    /// $x$-component of position. Convenience for plotting.
    #[getter]
    fn x(&self) -> f64 {
        self.inner.x
    }

    /// $y$-component of position. Convenience for plotting.
    #[getter]
    fn y(&self) -> f64 {
        self.inner.y
    }

    /// $z$-component of position. Convenience for plotting.
    #[getter]
    fn z(&self) -> f64 {
        self.inner.z
    }

    /// $x$-component of velocity. Convenience for plotting.
    #[getter]
    fn vx(&self) -> f64 {
        self.inner.vx
    }

    /// $y$-component of velocity. Convenience for plotting.
    #[getter]
    fn vy(&self) -> f64 {
        self.inner.vy
    }

    /// $z$-component of velocity. Convenience for plotting.
    #[getter]
    fn vz(&self) -> f64 {
        self.inner.vz
    }

    /// Plummer softening length $\epsilon$. Pairwise softening is
    /// combined in quadrature: $\epsilon_{ij}^2 = (\epsilon_i^2 +
    /// \epsilon_j^2)/2$.
    #[getter]
    fn softening(&self) -> f64 {
        self.inner.softening
    }

    /// Bulk density of the body. Drives the physical radius through
    /// $r = (3m/4\pi\rho)^{1/3}$.
    #[getter]
    fn density(&self) -> f64 {
        self.inner.density
    }

    /// Physical radius derived from mass and density. Independent of
    /// any softening calibration applied by the system.
    #[getter]
    fn radius(&self) -> f64 {
        self.inner.physical_radius
    }

    /// Material class as a canonical slug (e.g. `"star"`, `"rocky"`,
    /// `"gas_giant"`). The slug round-trips with the factory names:
    /// `Body.star(...).material == "star"`.
    #[getter]
    fn material(&self) -> &'static str {
        material_slug(self.inner.material)
    }

    /// Bolometric luminosity in internal energy / time units. Stays at
    /// zero for non-luminous materials and for luminous bodies that
    /// have not been processed by the radiation pipeline.
    #[getter]
    fn luminosity(&self) -> f64 {
        self.inner.luminosity
    }

    fn __repr__(&self) -> String {
        format!(
            "Body(material={:?}, mass={}, position=({}, {}, {}), velocity=({}, {}, {}), softening={})",
            material_slug(self.inner.material),
            self.inner.mass,
            self.inner.x,
            self.inner.y,
            self.inner.z,
            self.inner.vx,
            self.inner.vy,
            self.inner.vz,
            self.inner.softening,
        )
    }
}

/// Map a core `Material` variant to its canonical Python slug.
/// Compile-checked exhaustive — a new variant in
/// `apsis::domain::materials` will surface as a missing-arm error here
/// at the next `cargo build` and force a deliberate decision on the
/// Python-facing name. That deliberate-decision property is the entire
/// reason this lives in the binding rather than in the core: the Python
/// slug is part of the binding's API contract, not of the physics.
fn material_slug(m: Material) -> &'static str {
    match m {
        Material::Star => "star",
        Material::BrownDwarf => "brown_dwarf",
        Material::WhiteDwarf => "white_dwarf",
        Material::Gas => "gas_giant",
        Material::IceGiant => "ice_giant",
        Material::Rocky => "rocky",
        Material::Icy => "icy",
        Material::Asteroid => "asteroid",
        Material::Comet => "comet",
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBody>()?;
    Ok(())
}
