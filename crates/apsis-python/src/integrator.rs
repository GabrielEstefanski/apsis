//! Python-side wrapper of [`apsis::physics::integrator::IntegratorKind`].
//!
//! Researchers naming an integrator at the API boundary use either an
//! enum constant (`apsis.IntegratorKind.IAS15`) or the canonical
//! string slug (`"ias15"`). Both forms resolve to the same underlying
//! `apsis` enum variant via [`resolve`]; no other path constructs an
//! `IntegratorKind` from user input.
//!
//! The Python enum names use upper-case acronyms (`IAS15`, `YOSHIDA4`,
//! `VELOCITY_VERLET`, `WISDOM_HOLMAN`); the string slugs (`"ias15"`,
//! `"yoshida4"`, ...) mirror the core `IntegratorKind`'s `FromStr`, so
//! the same string works in `run.toml`, on the CLI, and at the Python
//! kwarg.

use std::str::FromStr;

use apsis::physics::integrator::IntegratorKind as CoreIntegratorKind;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::convert::{slugify, value_error};

/// Numerical integration scheme applied to the simulation's body state.
///
/// Each variant corresponds 1:1 to a Rust integrator implementation in
/// [`apsis::physics::integrator`]. Variants are chosen at construction
/// time and never reconfigured silently; switching integrator
/// mid-simulation is a deliberate `System` mutation (not yet exposed to
/// Python).
///
/// See the project's [`docs/integrator.md`] for the per-integrator
/// contract (execution profile, force-model determinism requirement,
/// selection rubric).
///
/// [`docs/integrator.md`]: https://github.com/GabrielEstefanski/gravity-sim/blob/develop/docs/integrator.md
#[pyclass(eq, eq_int, frozen, hash, module = "apsis", name = "IntegratorKind")]
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub(crate) enum IntegratorKind {
    /// 15th-order adaptive Gauss-Radau integrator (Rein & Spiegel 2015).
    /// Machine-precision energy conservation; off-line precision
    /// runs. Pairs exclusively with direct $O(N^2)$ gravity.
    #[pyo3(name = "IAS15")]
    Ias15,

    /// 4th-order symplectic composition (Yoshida 1990). Default for
    /// interactive playback at any $N$; bounded per-step wall time.
    #[pyo3(name = "YOSHIDA4")]
    Yoshida4,

    /// 2nd-order symplectic leapfrog. Lowest cost; energy oscillates
    /// at order $\Delta t^2$ around the initial value.
    #[pyo3(name = "VELOCITY_VERLET")]
    VelocityVerlet,

    /// Mixed-variable symplectic map for hierarchies with one
    /// dominant primary (Wisdom & Holman 1991). Analytic Kepler
    /// drift + numerical perturbation kick.
    #[pyo3(name = "WISDOM_HOLMAN")]
    WisdomHolman,

    /// Hybrid symplectic close-encounter integrator (Rein et al. 2019).
    /// Wisdom-Holman outer step with a K-weighted planet-planet kick;
    /// IAS15 sub-integrates the (1−K)-weighted residual on the
    /// encountering subset over the same outer interval. Requires a
    /// hierarchical mass distribution.
    #[pyo3(name = "MERCURIUS")]
    Mercurius,

    /// Wisdom-Holman split with compensated summation on per-step
    /// position and velocity accumulators (Rein & Tamayo 2015).
    /// Reduces round-off envelope from $O(N \cdot \varepsilon)$ to
    /// $O(\sqrt{N} \cdot \varepsilon)$, unlocking long-horizon
    /// planetary integration. Same hierarchical-mass requirement
    /// as Wisdom-Holman.
    #[pyo3(name = "WHFAST")]
    WHFast,

    /// Single-stage Gauss-Legendre implicit symplectic method
    /// (Hairer-Lubich-Wanner 2006, Chapter II.1.4). A-stable, time-
    /// symmetric, no central-mass dominance assumption — accepts BH
    /// binaries, equal-mass triples, particle clouds. Not L-stable.
    #[pyo3(name = "IMPLICIT_MIDPOINT")]
    ImplicitMidpoint,
}

#[pymethods]
impl IntegratorKind {
    fn __repr__(&self) -> String {
        format!("IntegratorKind.{}", self.name())
    }

    fn __str__(&self) -> String {
        self.slug().into()
    }

    /// Canonical Python enum name (`"IAS15"`, `"YOSHIDA4"`, ...).
    #[getter]
    fn name(&self) -> &'static str {
        match self {
            Self::Ias15 => "IAS15",
            Self::Yoshida4 => "YOSHIDA4",
            Self::VelocityVerlet => "VELOCITY_VERLET",
            Self::WisdomHolman => "WISDOM_HOLMAN",
            Self::Mercurius => "MERCURIUS",
            Self::WHFast => "WHFAST",
            Self::ImplicitMidpoint => "IMPLICIT_MIDPOINT",
        }
    }

    /// Lower-case canonical slug used by the core `FromStr` impl and
    /// by config files (`"ias15"`, `"yoshida4"`, ...).
    #[getter]
    fn slug(&self) -> &'static str {
        self.into_core().slug()
    }
}

impl IntegratorKind {
    pub(crate) fn from_core(core: CoreIntegratorKind) -> Self {
        match core {
            CoreIntegratorKind::Ias15 => Self::Ias15,
            CoreIntegratorKind::Yoshida4 => Self::Yoshida4,
            CoreIntegratorKind::VelocityVerlet => Self::VelocityVerlet,
            CoreIntegratorKind::WisdomHolman => Self::WisdomHolman,
            CoreIntegratorKind::Mercurius => Self::Mercurius,
            CoreIntegratorKind::WHFast => Self::WHFast,
            CoreIntegratorKind::ImplicitMidpoint => Self::ImplicitMidpoint,
        }
    }

    pub(crate) fn into_core(self) -> CoreIntegratorKind {
        match self {
            Self::Ias15 => CoreIntegratorKind::Ias15,
            Self::Yoshida4 => CoreIntegratorKind::Yoshida4,
            Self::VelocityVerlet => CoreIntegratorKind::VelocityVerlet,
            Self::WisdomHolman => CoreIntegratorKind::WisdomHolman,
            Self::Mercurius => CoreIntegratorKind::Mercurius,
            Self::WHFast => CoreIntegratorKind::WHFast,
            Self::ImplicitMidpoint => CoreIntegratorKind::ImplicitMidpoint,
        }
    }
}

/// Resolve a Python value (`IntegratorKind` instance or string slug) to
/// the core `apsis` enum variant. Raises `ValueError` at the boundary
/// for unknown strings or wrong types — researchers see a typo at call
/// time, not after a partial integration completes.
///
/// This is the single entry point through which user-facing strings
/// become `apsis::physics::integrator::IntegratorKind`. Other modules
/// in the binding (`crate::system`) call this rather than re-doing
/// the parse, so the set of accepted spellings is uniform.
pub(crate) fn resolve(obj: &Bound<'_, PyAny>) -> PyResult<CoreIntegratorKind> {
    if let Ok(kind) = obj.extract::<IntegratorKind>() {
        return Ok(kind.into_core());
    }
    if let Ok(s) = obj.extract::<String>() {
        let slug = slugify(&s);
        return CoreIntegratorKind::from_str(&slug).map_err(|err| {
            value_error(
                "integrator",
                format!(
                    "{err} (got {s:?}; accepted: 'ias15', 'yoshida4', \
                     'velocity_verlet', 'wisdom_holman', 'whfast', \
                     'mercurius', 'implicit_midpoint', or any \
                     IntegratorKind variant)"
                ),
            )
        });
    }
    Err(value_error(
        "integrator",
        format!(
            "expected an IntegratorKind variant or a string slug, got {}",
            obj.get_type().name().map(|s| s.to_string()).unwrap_or_else(|_| "<?>".into()),
        ),
    ))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<IntegratorKind>()?;
    Ok(())
}
