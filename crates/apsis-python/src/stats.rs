//! Frozen diagnostic-snapshot classes returned by `System.stats` /
//! `System.adaptive_stats`.
//!
//! Two classes mirror the Rust core's split:
//! - [`PyStats`] always-present scalars that every integrator carries
//!   (time, steps, current dt, energy / angular-momentum drift, etc.)
//! - [`PyAdaptiveStats`] controller counters that only adaptive
//!   integrators (currently IAS15) populate; `None` from the binding
//!   side for fixed-step integrators rather than a struct of zeros.

use apsis::core::system::System as CoreSystem;
use apsis::physics::integrator::IntegratorKind as CoreIntegratorKind;
use apsis::physics::integrator::traits::AdaptiveStats as CoreAdaptiveStats;
use pyo3::prelude::*;

use crate::integrator::IntegratorKind as PyIntegratorKind;

fn integrator_slug(kind: CoreIntegratorKind) -> &'static str {
    match kind {
        CoreIntegratorKind::Ias15 => "ias15",
        CoreIntegratorKind::Yoshida4 => "yoshida4",
        CoreIntegratorKind::VelocityVerlet => "velocity_verlet",
        CoreIntegratorKind::WisdomHolman => "wisdom_holman",
        CoreIntegratorKind::Mercurius => "mercurius",
        CoreIntegratorKind::WHFast => "whfast",
        CoreIntegratorKind::ImplicitMidpoint => "implicit_midpoint",
    }
}

/// Snapshot of the simulation's cumulative scalar diagnostics.
///
/// Returned by `System.stats`; immutable. A researcher prints
/// `print(sys.stats)` at the end of a run and gets every headline
/// number without composing several property reads.
#[pyclass(module = "apsis", name = "Stats", frozen)]
#[derive(Clone)]
pub(crate) struct PyStats {
    #[pyo3(get)]
    t: f64,
    #[pyo3(get)]
    steps: u64,
    #[pyo3(get)]
    dt: f64,
    #[pyo3(get)]
    energy: f64,
    #[pyo3(get)]
    energy_drift: Option<f64>,
    #[pyo3(get)]
    abs_energy_drift: f64,
    #[pyo3(get)]
    kinetic_energy: f64,
    #[pyo3(get)]
    potential_energy: f64,
    #[pyo3(get)]
    lz: f64,
    #[pyo3(get)]
    lz_drift: Option<f64>,
    #[pyo3(get)]
    abs_lz_drift: f64,
    #[pyo3(get)]
    integrator: PyIntegratorKind,
    #[pyo3(get)]
    force_evaluations: u64,

    integrator_kind: CoreIntegratorKind,
}

impl PyStats {
    pub(crate) fn from_system(sys: &CoreSystem) -> Self {
        let m = sys.metrics();
        Self {
            t: m.t,
            steps: m.steps,
            dt: m.dt,
            energy: m.total_energy,
            energy_drift: m.rel_energy_error,
            abs_energy_drift: m.abs_energy_error,
            kinetic_energy: m.kinetic,
            potential_energy: m.potential,
            lz: m.angular_momentum_z,
            lz_drift: m.rel_angular_momentum_error,
            abs_lz_drift: m.abs_angular_momentum_error,
            integrator: PyIntegratorKind::from_core(m.integrator_kind),
            force_evaluations: m.steps * (m.integrator_kind.force_evals_per_step() as u64),
            integrator_kind: m.integrator_kind,
        }
    }
}

#[pymethods]
impl PyStats {
    fn __repr__(&self) -> String {
        let de = match self.energy_drift {
            Some(rel) => format!("{:.3e}", rel),
            None => format!("{:.3e} (abs, |E0|~0)", self.abs_energy_drift),
        };
        let dl = match self.lz_drift {
            Some(rel) => format!("{:.3e}", rel),
            None => format!("{:.3e} (abs, |Lz0|~0)", self.abs_lz_drift),
        };
        format!(
            "Stats(t={:.6}, steps={}, dt={:.3e}, dE={}, dLz={}, integrator={:?})",
            self.t,
            self.steps,
            self.dt,
            de,
            dl,
            integrator_slug(self.integrator_kind),
        )
    }
}

/// Snapshot of an adaptive integrator's controller counters.
///
/// Returned by `System.adaptive_stats` for IAS15 (and any future
/// adaptive scheme); `None` from the binding for fixed-step
/// integrators since they don't run a controller.
///
/// Each counter is monotonically non-decreasing across the run.
/// Per-step rates (e.g. `picard_stagnations / substeps`) are the
/// usual diagnostic — sustained values ≫ 0 flag controller stress.
#[pyclass(module = "apsis", name = "AdaptiveStats", frozen)]
#[derive(Clone)]
pub(crate) struct PyAdaptiveStats {
    #[pyo3(get)]
    substeps: u64,
    #[pyo3(get)]
    rejections: u64,
    #[pyo3(get)]
    rejections_picard: u64,
    #[pyo3(get)]
    rejections_truncation: u64,
    #[pyo3(get)]
    picard_iters: u64,
    #[pyo3(get)]
    picard_stagnations: u64,
    #[pyo3(get)]
    shrink_grow_cycles: u64,
    #[pyo3(get)]
    degraded: u64,
}

impl PyAdaptiveStats {
    pub(crate) fn from_core(s: CoreAdaptiveStats) -> Self {
        Self {
            substeps: s.substeps,
            rejections: s.rejections,
            rejections_picard: s.rejections_picard,
            rejections_truncation: s.rejections_truncation,
            picard_iters: s.picard_iters,
            picard_stagnations: s.picard_stagnations,
            shrink_grow_cycles: s.shrink_grow_cycles,
            degraded: s.degraded,
        }
    }
}

#[pymethods]
impl PyAdaptiveStats {
    fn __repr__(&self) -> String {
        format!(
            "AdaptiveStats(substeps={}, rejections={} (picard={}, trunc={}), \
             picard_iters={}, picard_stagnations={}, shrink_grow_cycles={}, degraded={})",
            self.substeps,
            self.rejections,
            self.rejections_picard,
            self.rejections_truncation,
            self.picard_iters,
            self.picard_stagnations,
            self.shrink_grow_cycles,
            self.degraded,
        )
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyStats>()?;
    m.add_class::<PyAdaptiveStats>()?;
    Ok(())
}
