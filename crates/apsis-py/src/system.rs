//! Python-side wrapper of [`apsis::core::system::System`].
//!
//! The wrapper exposes a researcher-first API for the orchestrator:
//! a kwargs-only constructor that names every dial a user is likely
//! to set (`bodies`, `integrator`, `dt`, ...), simple verb methods
//! for the run loop (`step`, `integrate_for`, `integrate_until`),
//! and read-only properties for the cheap-to-query state
//! (`t`, `bodies`, `energy`, `energy_delta`, ...).
//!
//! # Façade-only invariant
//!
//! Every `#[pymethods]` body delegates to one or two calls on
//! [`apsis::core::system::System`]; nothing here implements a step,
//! a force evaluation, or a conservation diagnostic of its own. The
//! integrator-specific behaviour, the force-model determinism rule
//! enforced by `set_integrator`, and the orchestrator's hook
//! discipline are all owned by the core crate. This module exists
//! to translate Python kwargs into those calls and to lift the
//! results back through the FFI boundary in shapes that read like
//! research code rather than like Rust types.
//!
//! # What is not (yet) wrapped
//!
//! Hooks, snapshot serialisation, headless run-config loading, and
//! custom Python-defined perturbations are intentionally outside
//! Phase 1 — each requires either Python-side callbacks (and
//! therefore a careful PyO3 GIL story) or a serialisation contract
//! that hasn't yet stabilised on the Rust side. They land in
//! follow-up commits, not in this scaffolding.

use apsis::core::system::System as CoreSystem;
use apsis::physics::integrator::IntegratorKind as CoreIntegratorKind;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::body::PyBody;
use crate::convert::value_error;
use crate::integrator::{IntegratorKind as PyIntegratorKind, resolve as resolve_integrator};
use crate::trajectory::PyTrajectory;
use crate::units::PyUnitSystem;

/// Orchestrator for the simulation: bodies, chosen integrator, and
/// the run loop. Bodies are passed at construction time; the
/// integrator and step size are kwargs that the researcher decides
/// up front and that don't typically change mid-run.
///
/// `System` holds a mutable simulation state: each call to `step`,
/// `integrate_for`, or `integrate_until` advances bodies, accumulates
/// adaptive-controller statistics, and updates the energy /
/// angular-momentum baselines that drive the conservation accessors.
/// The cheap read-only properties (`t`, `energy`, `energy_delta`,
/// `lz`, `lz_delta`, ...) are O(1) reads of cached state — calling
/// them at any cadence has no effect on the simulation.
///
/// Example — Kepler two-body orbit at high eccentricity, integrated
/// for one century in canonical units:
///
/// ```python
/// import apsis
///
/// sun = apsis.Body.star(mass=1.0).unsoftened()
/// mercury = (apsis.Body.rocky(mass=3e-6)
///            .at((0.307, 0.0))
///            .with_velocity((0.0, 1.98))
///            .unsoftened())
///
/// sys = apsis.System(
///     bodies=[sun, mercury],
///     integrator="ias15",
///     dt=1e-3,
/// )
///
/// sys.integrate_for(100.0)
///
/// print(f"|dE/E_0| = {abs(sys.energy_delta):.3e}")  # ~1e-15
/// ```
#[pyclass(module = "apsis", name = "System", unsendable)]
pub(crate) struct PySystem {
    inner: CoreSystem,
}

#[pymethods]
impl PySystem {
    /// Construct a system from an explicit body list, a unit system,
    /// and integrator settings. All arguments are kwargs-only — there
    /// is no positional form, so a researcher reading the call site
    /// sees every dial named.
    ///
    /// Arguments:
    ///
    /// - `bodies`: list of [`Body`] instances. May be empty for the
    ///   degenerate case (`integrate_for` will then advance time
    ///   without changing anything), but typically holds two or more.
    /// - `units`: which [`UnitSystem`] interprets the body state and
    ///   the time step. Mandatory and immutable — pick from the named
    ///   factories (`apsis.units.SOLAR`, `apsis.units.SI`,
    ///   `apsis.units.CANONICAL`, `apsis.units.HENON`, `apsis.units.CGS`)
    ///   or build one with `UnitSystem.custom(...)`. There is no
    ///   default; the unit system is part of the simulation's
    ///   physical contract and must be stated explicitly.
    /// - `integrator`: which integrator to drive the run loop with.
    ///   Accepts an [`IntegratorKind`] variant or a canonical string
    ///   slug (`"ias15"`, `"yoshida4"`, `"velocity_verlet"`,
    ///   `"wisdom_holman"`).
    /// - `dt`: initial time step in the chosen unit system's
    ///   canonical time. For self-adaptive integrators (IAS15) this
    ///   is a hint that the controller mutates as it runs; for
    ///   fixed-step schemes it is the exact step.
    /// - `epsilon`: target relative truncation error per substep,
    ///   used by IAS15. `None` (default) keeps the integrator's
    ///   built-in default ($10^{-9}$, per Rein & Spiegel 2015).
    ///   Ignored by fixed-step integrators.
    /// - `exact_gravity`: drop the Plummer softening on every body in
    ///   one call. Equivalent to building the bodies with
    ///   `Body.<material>(...).unsoftened()` but applies system-wide
    ///   without per-body chaining. Default `False`.
    #[new]
    #[pyo3(signature = (*, bodies, units, integrator, dt, epsilon=None, exact_gravity=false))]
    fn new(
        bodies: Vec<PyBody>,
        units: PyUnitSystem,
        integrator: &Bound<'_, PyAny>,
        dt: f64,
        epsilon: Option<f64>,
        exact_gravity: bool,
    ) -> PyResult<Self> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(value_error(
                "dt",
                format!("expected a strictly positive finite float, got {dt}"),
            ));
        }
        if let Some(eps) = epsilon {
            if !eps.is_finite() || eps <= 0.0 {
                return Err(value_error(
                    "epsilon",
                    format!("expected a strictly positive finite float, got {eps}"),
                ));
            }
        }

        let kind: CoreIntegratorKind = resolve_integrator(integrator)?;
        let body_vec = bodies.into_iter().map(|b| b.inner).collect::<Vec<_>>();

        let mut sys = CoreSystem::new(body_vec, units.inner)
            .with_integrator(kind)
            .with_dt(dt);
        if exact_gravity {
            sys = sys.with_exact_gravity();
        }
        if let Some(eps) = epsilon {
            sys.set_ias15_epsilon(eps);
        }

        Ok(Self { inner: sys })
    }

    // ── Run loop ─────────────────────────────────────────────────────────
    //
    // Each method advances the simulation; nothing here reads back
    // state. The cheap-property getters below are the read path.

    /// Advance the simulation by exactly one integrator step. For
    /// fixed-step schemes (Velocity Verlet, Yoshida 4, Wisdom-Holman)
    /// this consumes the configured `dt`; for adaptive IAS15 it
    /// consumes one controller-chosen sub-step (typically smaller
    /// than `dt` near close encounters, larger in smooth regions).
    fn step(&mut self) {
        self.inner.step();
    }

    /// Advance the simulation by `duration` time units relative to
    /// the current `t`. For adaptive integrators, the actual final
    /// `t` may slightly exceed `current_t + duration` by at most one
    /// sub-step (the loop exits as soon as the threshold is crossed).
    fn integrate_for(&mut self, duration: f64) -> PyResult<u64> {
        if !duration.is_finite() || duration < 0.0 {
            return Err(value_error(
                "duration",
                format!("expected a finite non-negative float, got {duration}"),
            ));
        }
        Ok(self.inner.integrate_for(duration))
    }

    /// Advance the simulation until `t >= t_end`. No-op when the
    /// current `t` already meets that condition. Returns the number
    /// of integrator steps executed during the call.
    fn integrate_until(&mut self, t_end: f64) -> PyResult<u64> {
        if !t_end.is_finite() {
            return Err(value_error(
                "t_end",
                format!("expected a finite float, got {t_end}"),
            ));
        }
        Ok(self.inner.integrate_until(t_end))
    }

    /// Integrate forward by `duration` time units while recording the
    /// state at `n_samples` evenly spaced times, returning a
    /// [`Trajectory`](crate::trajectory::PyTrajectory) of NumPy arrays
    /// ready for `matplotlib`, `pandas`, or any other Python-side
    /// analysis tool.
    ///
    /// The first sample is taken before any integration runs, so
    /// `traj.t[0]` is the system's current `t` at the start of the
    /// call and `traj.x[0, :]` matches the bodies' positions at that
    /// instant. The last sample is taken after integration to
    /// `start_t + duration`, so `traj.t[-1] >= start_t + duration`
    /// (the integrator may overshoot the final target by at most one
    /// adaptive sub-step under IAS15). Intermediate samples are spaced
    /// uniformly in target time; the actual `traj.t` values are
    /// monotonically non-decreasing but inherit the same overshoot
    /// behaviour as `integrate_until`.
    ///
    /// Sampling **advances the system state**: after the call,
    /// `sys.t == traj.t[-1]` and the bodies hold the configuration
    /// recorded at the last sample. To preserve the pre-sample state,
    /// hold a snapshot of the bodies before calling.
    ///
    /// Arguments:
    ///
    /// - `duration`: total time advanced during the sampling window.
    ///   Must be strictly positive.
    /// - `n_samples`: number of samples to record, including the
    ///   initial state. Must be ≥ 1.
    ///
    /// Memory cost is `n_samples × (4 × n_bodies + 2)` `f64` values —
    /// 32 MB for `n_samples = 1000` and `n_bodies = 1000`, well below
    /// any realistic physical-research scenario. For sampling regimes
    /// large enough to exceed RAM, prefer a streaming consumer (which
    /// is not yet exposed; file an issue if this becomes a real
    /// constraint).
    fn sample(
        &mut self,
        py: Python<'_>,
        duration: f64,
        n_samples: usize,
    ) -> PyResult<PyTrajectory> {
        if !duration.is_finite() || duration <= 0.0 {
            return Err(value_error(
                "duration",
                format!("expected a strictly positive finite float, got {duration}"),
            ));
        }
        if n_samples == 0 {
            return Err(value_error("n_samples", "expected at least 1, got 0"));
        }

        let n_bodies = self.inner.bodies().len();
        let t_start = self.inner.t();
        let stride = if n_samples == 1 {
            // Degenerate single-sample run: just record current state.
            // Stride is unused but must be defined.
            0.0
        } else {
            duration / (n_samples - 1) as f64
        };

        let mut t_buf = Vec::with_capacity(n_samples);
        let mut x_buf = Vec::with_capacity(n_samples * n_bodies);
        let mut y_buf = Vec::with_capacity(n_samples * n_bodies);
        let mut vx_buf = Vec::with_capacity(n_samples * n_bodies);
        let mut vy_buf = Vec::with_capacity(n_samples * n_bodies);
        let mut energy_buf = Vec::with_capacity(n_samples);

        // Prime the energy / angular-momentum cache so the pre-integration
        // sample carries a real K + U rather than the construction-time
        // zero. After the first integrator step, both `step()` and the
        // adaptive controller maintain the cache themselves.
        self.inner.refresh_energy_diagnostics();

        record_state(
            &self.inner,
            &mut t_buf,
            &mut x_buf,
            &mut y_buf,
            &mut vx_buf,
            &mut vy_buf,
            &mut energy_buf,
        );

        for i in 1..n_samples {
            let t_target = t_start + stride * i as f64;
            self.inner.integrate_until(t_target);
            record_state(
                &self.inner,
                &mut t_buf,
                &mut x_buf,
                &mut y_buf,
                &mut vx_buf,
                &mut vy_buf,
                &mut energy_buf,
            );
        }

        PyTrajectory::build(py, t_buf, x_buf, y_buf, vx_buf, vy_buf, energy_buf, n_samples, n_bodies)
    }

    // ── State & diagnostics (O(1) properties) ────────────────────────────

    /// Current simulation time in simulation units.
    #[getter]
    fn t(&self) -> f64 {
        self.inner.t()
    }

    /// Number of completed integrator steps since construction.
    #[getter]
    fn steps(&self) -> u64 {
        self.inner.steps()
    }

    /// Snapshot of the body list at the current state. Returns a fresh
    /// Python list — mutating it does not affect the simulation, and
    /// the bodies themselves are immutable on the Python side. Stable
    /// across calls only if the simulation is not stepped between
    /// them.
    #[getter]
    fn bodies(&self) -> Vec<PyBody> {
        self.inner.bodies().iter().copied().map(|b| PyBody { inner: b }).collect()
    }

    /// Total mechanical energy at the most recently completed step,
    /// $E = K + U$. Cached during the step itself (no body sum at
    /// read time).
    #[getter]
    fn energy(&self) -> f64 {
        self.inner.energy()
    }

    /// Relative energy drift, $\delta E / E_0 = (E - E_0) / |E_0|$.
    /// `E_0` is captured at the first step of the simulation; this is
    /// the cheapest single-number diagnostic of integrator quality.
    #[getter]
    fn energy_delta(&self) -> f64 {
        self.inner.energy_delta()
    }

    /// Kinetic energy at the most recent step,
    /// $K = \sum_i \tfrac{1}{2} m_i |\mathbf{v}_i|^2$.
    #[getter]
    fn kinetic_energy(&self) -> f64 {
        self.inner.kinetic_energy()
    }

    /// Potential energy at the most recent step,
    /// $U = -\sum_{i<j} G m_i m_j / |\mathbf{r}_i - \mathbf{r}_j|$.
    /// Includes any active perturbations whose contribution maps to a
    /// potential.
    #[getter]
    fn potential_energy(&self) -> f64 {
        self.inner.potential_energy()
    }

    /// Total $z$-component of angular momentum,
    /// $L_z = \sum_i m_i (x_i v_{y,i} - y_i v_{x,i})$. The 2D
    /// simulator's only conserved component of $\mathbf{L}$.
    #[getter]
    fn lz(&self) -> f64 {
        self.inner.lz()
    }

    /// Relative angular-momentum drift,
    /// $\delta L_z / L_{z,0} = (L_z - L_{z,0}) / |L_{z,0}|$.
    /// Falls back to absolute drift when `L_{z,0}` is below an
    /// internal numerical threshold (figure-8-style choreographies
    /// where total angular momentum is zero by construction).
    #[getter]
    fn lz_delta(&self) -> f64 {
        self.inner.lz_delta()
    }

    /// Current controller time step. For self-adaptive integrators
    /// (IAS15) this tracks the controller's `dt_next` proposal — it
    /// changes from step to step as the geometry demands. For
    /// fixed-step schemes this stays at the value passed to
    /// `dt` in the constructor.
    #[getter]
    fn dt(&self) -> f64 {
        self.inner.metrics().dt
    }

    /// Which integrator is driving this system, as an
    /// [`IntegratorKind`] variant. Stable across `step` /
    /// `integrate_for` calls; mutated only by an explicit
    /// `set_integrator` (not yet exposed to Python).
    #[getter]
    fn integrator(&self) -> PyIntegratorKind {
        PyIntegratorKind::from_core(self.inner.metrics().integrator_kind)
    }

    /// The unit system the simulation was constructed against.
    ///
    /// Frozen for the lifetime of the system — there is no setter,
    /// and no integration step can mutate it. Read it to log the
    /// run's physical contract, compare to a saved baseline, or
    /// drive a downstream conversion (e.g. plot positions in metres
    /// when the simulation runs in `solar()` units).
    #[getter]
    fn units(&self) -> PyUnitSystem {
        PyUnitSystem::from_core(*self.inner.units())
    }

    // ── Mutators ─────────────────────────────────────────────────────────

    /// Translate every body so the system's centre of mass is at
    /// the origin. Idempotent when the COM is already there to
    /// within $10^{-14}$. Routed through the active integrator's
    /// compensation buffers, so IAS15's `csx` invariant is preserved
    /// across the shift.
    fn recenter_com(&mut self) {
        self.inner.recenter_com();
    }

    fn __repr__(&self) -> String {
        let m = self.inner.metrics();
        format!(
            "System(integrator={:?}, n_bodies={}, t={}, dt={}, dE/E0={:.3e})",
            kind_slug(m.integrator_kind),
            self.inner.bodies().len(),
            m.t,
            m.dt,
            m.rel_energy_error,
        )
    }
}

fn kind_slug(kind: CoreIntegratorKind) -> &'static str {
    match kind {
        CoreIntegratorKind::Ias15 => "ias15",
        CoreIntegratorKind::Yoshida4 => "yoshida4",
        CoreIntegratorKind::VelocityVerlet => "velocity_verlet",
        CoreIntegratorKind::WisdomHolman => "wisdom_holman",
    }
}

/// Append one row of state to the trajectory buffers.
///
/// The 1-D buffers (`t_buf`, `energy_buf`) get one entry; the 2-D
/// buffers (`x_buf`, `y_buf`, `vx_buf`, `vy_buf`) get one entry per
/// body in row-major (sample-major) order so the resulting NumPy
/// `Array2::from_shape_vec((n_samples, n_bodies), ...)` reshapes
/// correctly without further transposition.
///
/// Pulled out into a dedicated helper rather than inlined into the
/// sample loop because the same recording pattern repeats once for
/// the initial-state sample and once per integration target — DRY
/// at the wrapper level even though each individual call is a few
/// `push` statements.
fn record_state(
    sys: &CoreSystem,
    t_buf: &mut Vec<f64>,
    x_buf: &mut Vec<f64>,
    y_buf: &mut Vec<f64>,
    vx_buf: &mut Vec<f64>,
    vy_buf: &mut Vec<f64>,
    energy_buf: &mut Vec<f64>,
) {
    t_buf.push(sys.t());
    energy_buf.push(sys.energy());
    for b in sys.bodies() {
        x_buf.push(b.x);
        y_buf.push(b.y);
        vx_buf.push(b.vx);
        vy_buf.push(b.vy);
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySystem>()?;
    Ok(())
}
