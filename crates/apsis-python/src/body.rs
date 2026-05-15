//! Python-side wrapper of [`apsis::domain::body::Body`].
//!
//! The wrapper exposes a researcher-first API: nine preset factories
//! (`Body.star`, `Body.rocky`, `Body.gas_giant`, ...) that mirror the
//! corresponding constructors in [`apsis::domain::body::Body`], a
//! kwargs-only signature on each factory so position and velocity
//! never depend on argument order, and an immutable fluent-builder
//! tail (`at`, `with_velocity`, `with_density`) for the cases where
//! chaining reads more naturally than a single call site.
//!
//! # Façade-only invariant
//!
//! Every `#[pymethods]` body in this file delegates to a single call on
//! [`apsis::domain::body::Body`] — the wrapper translates types at the
//! boundary and never composes physics. The set of valid presets, the
//! density model, and the body-state
//! invariants are all owned by the core crate; this module is the
//! Python-shaped door into them.
//!
//! # Material slug as binding-layer tag
//!
//! [`apsis::domain::body::Body`] no longer carries a runtime material
//! taxonomy field — physical defaults are applied once at construction
//! by the preset and never referenced again. The Python wrapper still
//! exposes a `body.material` slug for ergonomic introspection (`"star"`,
//! `"rocky"`, ...), tracked locally by [`PyBody::slug`] and propagated
//! through the fluent builder methods. The slug is a binding-layer
//! convenience, not a core-crate concept.
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
use apsis::domain::body_preset::{self, BodyPreset};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::convert::{value_error, xyz_triple};

/// Point-mass body with kinematics, mass, and a binding-layer material
/// slug.
///
/// Bodies are constructed through one of the nine preset factories
/// (`Body.star`, `Body.rocky`, ...) — each chooses sensible defaults
/// for density and any visual/physical properties tied to the preset.
/// Position and velocity always default to zero; pass them as kwargs
/// (`position=(x, y)`, `velocity=(vx, vy)`) at the factory or via the
/// fluent builder methods (`at`, `with_velocity`).
///
/// All builder methods return a fresh `Body` — bodies are immutable
/// by convention on the Python side, matching the value-semantics of
/// the underlying `apsis::domain::body::Body` (which is `Copy`).
///
/// Examples:
///
/// ```python
/// import apsis
///
/// mercury = apsis.Body.rocky(
///     mass=3e-6,
///     position=(0.307, 0.0),
///     velocity=(0.0, 1.98),
/// )
/// ```
#[pyclass(module = "apsis", name = "Body", frozen)]
#[derive(Clone, Copy)]
pub(crate) struct PyBody {
    pub(crate) inner: CoreBody,
    /// Construction-time tag exposed via [`material`](Self::material).
    /// Tracked alongside the core body so the Python `body.material`
    /// property keeps returning a stable slug after chained builder
    /// calls — pure binding-layer convenience.
    pub(crate) slug: &'static str,
}

impl PyBody {
    /// Build a `PyBody` from a `BodyPreset` with kwargs-driven
    /// kinematics. Single point of construction for every factory
    /// below; per-factory methods are one-liner wrappers that name a
    /// specific preset.
    fn build(
        preset: &'static BodyPreset,
        slug: &'static str,
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
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

        let mut inner = CoreBody::from_preset(preset, mass);

        if let Some(obj) = position {
            let (x, y, z) = xyz_triple("position", obj)?;
            inner.pos_x = x;
            inner.pos_y = y;
            inner.pos_z = z;
        }
        if let Some(obj) = velocity {
            let (vx, vy, vz) = xyz_triple("velocity", obj)?;
            inner.vel_x = vx;
            inner.vel_y = vy;
            inner.vel_z = vz;
        }

        Ok(Self { inner, slug })
    }
}

#[pymethods]
impl PyBody {
    // ── Preset factories ─────────────────────────────────────────────────

    /// Main-sequence luminous body. Default density and luminous preset.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn star(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::STAR, "star", mass, position, velocity)
    }

    /// Brown dwarf — sub-stellar, deuterium-burning regime.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn brown_dwarf(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::BROWN_DWARF, "brown_dwarf", mass, position, velocity)
    }

    /// White dwarf — compact stellar remnant.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn white_dwarf(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::WHITE_DWARF, "white_dwarf", mass, position, velocity)
    }

    /// Gas giant — Jupiter-class hydrogen/helium envelope.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn gas_giant(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::GAS, "gas_giant", mass, position, velocity)
    }

    /// Ice giant — Neptune-class water/methane envelope.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn ice_giant(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::ICE_GIANT, "ice_giant", mass, position, velocity)
    }

    /// Rocky body — terrestrial planet or large rocky satellite.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn rocky(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::ROCKY, "rocky", mass, position, velocity)
    }

    /// Icy body — water-dominated composition (outer satellites, KBOs).
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn icy(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::ICY, "icy", mass, position, velocity)
    }

    /// Asteroid — rocky minor body.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn asteroid(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::ASTEROID, "asteroid", mass, position, velocity)
    }

    /// Comet — volatile-rich minor body.
    #[staticmethod]
    #[pyo3(signature = (*, mass, position=None, velocity=None))]
    fn comet(
        mass: f64,
        position: Option<&Bound<'_, PyAny>>,
        velocity: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Self::build(&body_preset::COMET, "comet", mass, position, velocity)
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
        inner.pos_x = x;
        inner.pos_y = y;
        inner.pos_z = z;
        Ok(Self { inner, slug: self.slug })
    }

    /// Set the body's velocity to `(vx, vy)` or `(vx, vy, vz)`. Returns
    /// a new `Body`. A 2-element sequence is treated as planar input
    /// with `vz = 0`.
    fn with_velocity(&self, velocity: &Bound<'_, PyAny>) -> PyResult<Self> {
        let (vx, vy, vz) = xyz_triple("velocity", velocity)?;
        let mut inner = self.inner;
        inner.vel_x = vx;
        inner.vel_y = vy;
        inner.vel_z = vz;
        Ok(Self { inner, slug: self.slug })
    }

    /// Override the preset-default density. Physical radius is
    /// recomputed from the new value at the call site (delegated to
    /// `apsis::domain::body::Body::with_density`). Returns a new `Body`.
    fn with_density(&self, density: f64) -> PyResult<Self> {
        if !density.is_finite() || density <= 0.0 {
            return Err(value_error(
                "density",
                format!("expected a strictly positive finite float, got {density}"),
            ));
        }
        Ok(Self { inner: self.inner.with_density(density), slug: self.slug })
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
        (self.inner.pos_x, self.inner.pos_y, self.inner.pos_z)
    }

    /// Velocity as a 3-tuple $(v_x, v_y, v_z)$.
    #[getter]
    fn velocity(&self) -> (f64, f64, f64) {
        (self.inner.vel_x, self.inner.vel_y, self.inner.vel_z)
    }

    /// $x$-component of position. Convenience for plotting.
    #[getter]
    fn x(&self) -> f64 {
        self.inner.pos_x
    }

    /// $y$-component of position. Convenience for plotting.
    #[getter]
    fn y(&self) -> f64 {
        self.inner.pos_y
    }

    /// $z$-component of position. Convenience for plotting.
    #[getter]
    fn z(&self) -> f64 {
        self.inner.pos_z
    }

    /// $x$-component of velocity. Convenience for plotting.
    #[getter]
    fn vx(&self) -> f64 {
        self.inner.vel_x
    }

    /// $y$-component of velocity. Convenience for plotting.
    #[getter]
    fn vy(&self) -> f64 {
        self.inner.vel_y
    }

    /// $z$-component of velocity. Convenience for plotting.
    #[getter]
    fn vz(&self) -> f64 {
        self.inner.vel_z
    }

    /// Bulk density of the body. Drives the physical radius through
    /// $r = (3m/4\pi\rho)^{1/3}$.
    #[getter]
    fn density(&self) -> f64 {
        self.inner.density
    }

    /// Physical radius derived from mass and density.
    #[getter]
    fn radius(&self) -> f64 {
        self.inner.physical_radius
    }

    /// Construction-time preset slug (e.g. `"star"`, `"rocky"`,
    /// `"gas_giant"`). Round-trips with the factory names:
    /// `Body.star(...).material == "star"`. The slug is a binding-
    /// layer convenience; the underlying core body holds no material
    /// taxonomy field.
    #[getter]
    fn material(&self) -> &'static str {
        self.slug
    }

    /// Bolometric luminosity in solar luminosities. Set at construction
    /// time by luminous presets ([`star`], [`brown_dwarf`],
    /// [`white_dwarf`]); zero for non-luminous classes.
    #[getter]
    fn luminosity(&self) -> f64 {
        self.inner.luminosity
    }

    /// Radiation-pressure receiver coefficient `Q_pr`. Positive on
    /// radiation receivers (asteroids, comets, icy grains); zero on
    /// emitters and large planets.
    #[getter]
    fn q_pr(&self) -> f64 {
        self.inner.q_pr
    }

    fn __repr__(&self) -> String {
        format!(
            "Body(material={:?}, mass={}, position=({}, {}, {}), velocity=({}, {}, {}))",
            self.slug,
            self.inner.mass,
            self.inner.pos_x,
            self.inner.pos_y,
            self.inner.pos_z,
            self.inner.vel_x,
            self.inner.vel_y,
            self.inner.vel_z,
        )
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyBody>()?;
    Ok(())
}
