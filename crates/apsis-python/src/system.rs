//! Python-side wrapper of [`apsis::core::system::System`].
//!
//! Kwargs-only constructor, verb-method run loop (`step`,
//! `integrate_for`, `integrate_until`, `sample`), and O(1) read-only
//! diagnostics. Every `#[pymethods]` body delegates to the core crate;
//! no physics or conservation logic lives here.

use apsis::core::system::System as CoreSystem;
use apsis::physics::integrator::IntegratorKind as CoreIntegratorKind;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::body::PyBody;
use crate::convert::value_error;
use crate::integrator::{IntegratorKind as PyIntegratorKind, resolve as resolve_integrator};
use crate::perturbation::take_perturbation_from_python;
use crate::stats::{PyAdaptiveStats, PyStats};
use crate::trajectory::{PyTrajectory, TrajectoryBuffers};
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
/// sun = apsis.Body.star(mass=1.0)
/// mercury = (apsis.Body.rocky(mass=3e-6)
///            .at((0.307, 0.0))
///            .with_velocity((0.0, 1.98))
///            )
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
    /// - `mercurius_alpha`: Hill-radius multiplier for the Mercurius
    ///   close-encounter changeover. `None` (default) keeps the
    ///   integrator's built-in $\alpha = 3$, matching REBOUND. Ignored
    ///   by non-Mercurius integrators.
    #[new]
    #[pyo3(signature = (*, bodies, units, integrator, dt, epsilon=None, mercurius_alpha=None))]
    fn new(
        bodies: Vec<PyBody>,
        units: PyUnitSystem,
        integrator: &Bound<'_, PyAny>,
        dt: f64,
        epsilon: Option<f64>,
        mercurius_alpha: Option<f64>,
    ) -> PyResult<Self> {
        if !dt.is_finite() || dt <= 0.0 {
            return Err(value_error(
                "dt",
                format!("expected a strictly positive finite float, got {dt}"),
            ));
        }
        if let Some(eps) = epsilon
            && (!eps.is_finite() || eps <= 0.0)
        {
            return Err(value_error(
                "epsilon",
                format!("expected a strictly positive finite float, got {eps}"),
            ));
        }
        if let Some(alpha) = mercurius_alpha
            && (!alpha.is_finite() || alpha < 0.0)
        {
            return Err(value_error(
                "mercurius_alpha",
                format!("expected a non-negative finite float, got {alpha}"),
            ));
        }

        let kind: CoreIntegratorKind = resolve_integrator(integrator)?;
        let body_vec = bodies.into_iter().map(|b| b.inner).collect::<Vec<_>>();

        let mut sys = CoreSystem::new(body_vec, units.inner).with_integrator(kind).with_dt(dt);
        if let Some(eps) = epsilon {
            sys.set_ias15_epsilon(eps);
        }
        if let Some(alpha) = mercurius_alpha {
            sys.set_mercurius_alpha(alpha);
        }

        Ok(Self { inner: sys })
    }

    // ── Run loop ─────────────────────────────────────────────────────────

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
            return Err(value_error("t_end", format!("expected a finite float, got {t_end}")));
        }
        Ok(self.inner.integrate_until(t_end))
    }

    /// Close any attached records and fire each hook's lifecycle-end
    /// callback. Idempotent; called automatically when the System is
    /// garbage-collected. Call explicitly for deterministic close.
    fn finish(&mut self) {
        self.inner.finish();
    }

    /// Record the system state at a set of target times, returning a
    /// [`Trajectory`](crate::trajectory::PyTrajectory) of NumPy arrays
    /// ready for `matplotlib`, `pandas`, or any other Python-side
    /// analysis tool.
    ///
    /// # Two invocation forms
    ///
    /// **Explicit times** (primary):
    ///
    /// ```python
    /// import numpy as np
    /// traj = sys.sample(times=np.linspace(0.0, 100.0, 1024))
    /// traj = sys.sample(times=np.logspace(-3, 2, 200))
    /// traj = sys.sample(times=[0.0, 1.0, 10.0, 100.0])
    /// ```
    ///
    /// The `times` argument accepts any 1-D sequence of floats — NumPy
    /// arrays, Python lists, tuples. The simulator integrates forward
    /// (using `integrate_until`) to each target time in order and
    /// records the state. There is **no interpolation**: each row of
    /// the returned trajectory is the integrator's actual output at
    /// (or just past) the requested time, with overshoot ≤ 1 adaptive
    /// sub-step under IAS15.
    ///
    /// **Evenly spaced** (convenience):
    ///
    /// ```python
    /// traj = sys.sample(duration=10.0, n_samples=128)
    /// ```
    ///
    /// Equivalent to passing
    /// `times=np.linspace(sys.t, sys.t + duration, n_samples)`.
    ///
    /// Pass exactly one form: either `times=` alone, or both
    /// `duration=` and `n_samples=` together. Mixing the two forms
    /// raises `ValueError` at the FFI boundary.
    ///
    /// # Validation
    ///
    /// `times` must be:
    /// - non-empty (at least one sample),
    /// - finite (no `NaN` or `±∞`),
    /// - monotonically non-decreasing (the simulator integrates
    ///   forward only — moving backwards would require either
    ///   reversibility, which not every integrator provides, or a
    ///   reset, which silently destroys state),
    /// - `times[0] >= sys.t` for the same reason.
    ///
    /// Each rule is reported with its offending index so a researcher
    /// debugging a notebook sees the exact element that broke the
    /// contract.
    ///
    /// # Side effects
    ///
    /// Sampling **advances the system state**: after the call,
    /// `sys.t == traj.t[-1]` and the bodies hold the configuration
    /// recorded at the last sample. The energy / angular-momentum
    /// cache is primed so `traj.energy[0]` carries a real `K + U`
    /// even when sampling starts from a freshly-constructed system.
    ///
    #[pyo3(signature = (*, times=None, duration=None, n_samples=None))]
    fn sample(
        &mut self,
        py: Python<'_>,
        times: Option<&Bound<'_, PyAny>>,
        duration: Option<f64>,
        n_samples: Option<usize>,
    ) -> PyResult<PyTrajectory> {
        let resolved_times: Vec<f64> = match (times, duration, n_samples) {
            (Some(arr), None, None) => extract_times(arr)?,
            (None, Some(dur), Some(n)) => build_evenly_spaced_times(self.inner.t(), dur, n)?,
            (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
                return Err(value_error(
                    "sample",
                    "pass either times= or both duration= and n_samples=, not both",
                ));
            },
            (None, Some(_), None) | (None, None, Some(_)) => {
                return Err(value_error(
                    "sample",
                    "duration= and n_samples= must be passed together",
                ));
            },
            (None, None, None) => {
                return Err(value_error(
                    "sample",
                    "pass either times= or both duration= and n_samples=",
                ));
            },
        };

        validate_times(&resolved_times, self.inner.t())?;

        let n_samples = resolved_times.len();
        let n_bodies = self.inner.bodies().len();

        let mut buf = TrajectoryBuffers {
            t: Vec::with_capacity(n_samples),
            x: Vec::with_capacity(n_samples * n_bodies),
            y: Vec::with_capacity(n_samples * n_bodies),
            z: Vec::with_capacity(n_samples * n_bodies),
            vx: Vec::with_capacity(n_samples * n_bodies),
            vy: Vec::with_capacity(n_samples * n_bodies),
            vz: Vec::with_capacity(n_samples * n_bodies),
            energy: Vec::with_capacity(n_samples),
            dt: Vec::with_capacity(n_samples),
            energy_drift: Vec::with_capacity(n_samples),
            lz_drift: Vec::with_capacity(n_samples),
        };

        // Pre-step cache is zero until the first integrator step; without
        // priming, `traj.energy[0]` would be the construction-time zero.
        self.inner.refresh_energy_diagnostics();

        for &t_target in &resolved_times {
            self.inner.integrate_until(t_target);
            record_state(&self.inner, &mut buf);
        }

        PyTrajectory::build(py, buf, n_samples, n_bodies)
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
    ///
    /// Round-tripped bodies report `material == "body"` because the
    /// core `Body` carries no preset reference; the construction-time
    /// slug is binding-layer state that does not survive the System
    /// round-trip. Use the slug on freshly-constructed bodies (before
    /// they enter a System) when material introspection matters.
    #[getter]
    fn bodies(&self) -> Vec<PyBody> {
        self.inner.bodies().iter().copied().map(|b| PyBody { inner: b, slug: "body" }).collect()
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
    /// Frozen — no setter, and integration cannot mutate it.
    #[getter]
    fn units(&self) -> PyUnitSystem {
        PyUnitSystem::from_core(*self.inner.units())
    }

    // ── Adaptive controller counters ─────────────────────────────────────
    //
    // Zero for fixed-step integrators (Velocity Verlet, Yoshida 4,
    // Wisdom-Holman) — those expose nothing to count, so a single number
    // per dimension keeps the API uniform across integrator choices.
    // For richer per-counter access on adaptive integrators, prefer
    // `sys.adaptive_stats` which returns `None` for fixed-step instead
    // of conflating "zero" with "not applicable".

    /// Accepted IAS15 sub-steps. `0` for fixed-step integrators.
    #[getter]
    fn substeps(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.substeps).unwrap_or(0)
    }

    /// Total IAS15 step rejections (controller shrunk `dt` and retried).
    /// `0` for fixed-step integrators.
    #[getter]
    fn step_rejections(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.rejections).unwrap_or(0)
    }

    /// IAS15 Picard predictor–corrector early-exits via the stagnation
    /// guard. Sustained `picard_stagnations / substeps ≫ 0` indicates
    /// the warmstart is biasing the predictor outside its convergence
    /// basin. `0` for fixed-step integrators.
    #[getter]
    fn picard_stagnations(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.picard_stagnations).unwrap_or(0)
    }

    /// Number of "shrink → grow" reversals in the IAS15 controller's
    /// `dt_next`. Reveals controller chatter; healthy smooth runs see
    /// `shrink_grow_cycles / substeps ≈ 0`. `0` for fixed-step.
    #[getter]
    fn shrink_grow_cycles(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.shrink_grow_cycles).unwrap_or(0)
    }

    /// Cumulative Picard iterations across all IAS15 attempts (accepted
    /// and rejected). `0` for fixed-step integrators.
    #[getter]
    fn picard_iters(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.picard_iters).unwrap_or(0)
    }

    /// Accepted IAS15 sub-steps that hit the `DT_MIN` floor or the
    /// step deadline without meeting tolerance. Should be `0` in
    /// healthy scenes; `0` for fixed-step integrators.
    #[getter]
    fn degraded_steps(&self) -> u64 {
        self.inner.adaptive_stats().map(|s| s.degraded).unwrap_or(0)
    }

    /// Estimated total force evaluations since construction. Computed
    /// as `steps × integrator.force_evals_per_step()`; for IAS15 the
    /// per-step factor (14) is an amortised average over the
    /// Gauss-Radau stages and the typical Picard iteration count.
    /// Use `picard_iters` for the exact IAS15 cost driver.
    #[getter]
    fn force_evaluations(&self) -> u64 {
        let m = self.inner.metrics();
        m.steps * (m.integrator_kind.force_evals_per_step() as u64)
    }

    /// Full snapshot of cumulative diagnostics as a [`Stats`] object.
    /// O(1); reads cached counters with no integration side effects.
    #[getter]
    fn stats(&self) -> PyStats {
        PyStats::from_system(&self.inner)
    }

    /// Adaptive-integrator counters as a frozen [`AdaptiveStats`] object,
    /// or `None` for fixed-step integrators where no controller state
    /// exists.
    #[getter]
    fn adaptive_stats(&self) -> Option<PyAdaptiveStats> {
        self.inner.adaptive_stats().map(PyAdaptiveStats::from_core)
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

    /// Attach a Hamiltonian-class perturbation constructed by a
    /// downstream binding crate (`apsis_1pn`, future J2 / tidal
    /// packages). The perturbation is consumed by the call — the same
    /// `Perturbation` instance cannot be attached twice; build a fresh
    /// one for each system.
    ///
    /// Raises ``apsis.UnitSystemMismatchError`` when the perturbation
    /// was constructed for a different ``UnitSystem`` than the
    /// ``System``'s own. The exception carries ``operator``,
    /// ``operator_units``, and ``system_units`` attributes so callers
    /// can decide policy (log, skip, swap, fall back).
    ///
    /// Kernel-requirement violations (e.g. attaching a 1PN correction
    /// on top of a softened kernel) emit structured warnings on
    /// ``stderr`` — they are non-fatal and do not raise. The default
    /// kernel is exact ``NewtonKernel`` (ε = 0); the warning only fires
    /// when the caller explicitly opts into a softened kernel.
    ///
    /// Non-conservative operators (drag, radiation reaction) travel in a
    /// separate capsule; the Python wrapper for them is not yet exposed.
    fn add_hamiltonian_perturbation(&mut self, perturbation: &Bound<'_, PyAny>) -> PyResult<()> {
        let boxed = take_perturbation_from_python(perturbation)?;
        // Result error variant is `Box<UnitSystemMismatch>` to keep
        // the Result enum small per clippy::result_large_err; deref
        // before passing to the PyErr converter.
        self.inner
            .add_hamiltonian_perturbation(boxed)
            .map_err(|e| crate::errors::unit_system_mismatch_to_pyerr(*e))
    }

    // ── Records ──────────────────────────────────────────────────────────

    /// Attach an Apsis Record writer to this system. Subsequent
    /// ``step()`` calls write to ``path``; the file is closed (with a
    /// trailer) when the System is dropped or the record hook is
    /// otherwise destroyed.
    ///
    /// Snapshot cadence — at most one of these keyword arguments:
    ///
    /// - ``every_steps``: emit a Snapshot when ``steps() % N == 0``
    /// - ``every_time``: emit a Snapshot when sim time crosses a
    ///   multiple of ``Δt``
    /// - ``dense``: emit a Snapshot every step (debug mode)
    ///
    /// Diagnostic cadence (optional, independent) — at most one of:
    ///
    /// - ``diagnostics_every_steps``: emit a Diagnostic frame (ΔE/E,
    ///   ΔLz/Lz) when ``steps() % N == 0``
    /// - ``diagnostics_every_time``: emit a Diagnostic frame when sim
    ///   time crosses a multiple of ``Δt``
    ///
    /// Default: bookends + events only (initial + final + collisions /
    /// escapes); no diagnostic frames.
    ///
    /// The header is gathered from the System's current state
    /// (operators, kernel, units, integrator, seed) and includes a
    /// BLAKE3 hash of the workspace ``Cargo.lock``.
    #[pyo3(signature = (
        path,
        *,
        seed = None,
        every_steps = None,
        every_time = None,
        dense = false,
        diagnostics_every_steps = None,
        diagnostics_every_time = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn attach_record(
        &mut self,
        path: String,
        seed: Option<u64>,
        every_steps: Option<u32>,
        every_time: Option<f64>,
        dense: bool,
        diagnostics_every_steps: Option<u32>,
        diagnostics_every_time: Option<f64>,
    ) -> PyResult<()> {
        use apsis::records::{
            DiagnosticCadence, RecordHook, RecordPolicy, provenance::header_from_system,
        };
        let policy = match (every_steps, every_time, dense) {
            (None, None, false) => RecordPolicy::BookendsAndEvents,
            (Some(n), None, false) => RecordPolicy::EveryNSteps(n),
            (None, Some(dt), false) => RecordPolicy::EveryTime(dt),
            (None, None, true) => RecordPolicy::Dense,
            _ => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "attach_record: set at most one of every_steps / every_time / dense",
                ));
            },
        };
        let diagnostics = match (diagnostics_every_steps, diagnostics_every_time) {
            (None, None) => DiagnosticCadence::Off,
            (Some(n), None) => DiagnosticCadence::EveryNSteps(n),
            (None, Some(dt)) => DiagnosticCadence::EveryTime(dt),
            (Some(_), Some(_)) => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "attach_record: set at most one of diagnostics_every_steps / diagnostics_every_time",
                ));
            },
        };
        let seed = seed.unwrap_or(self.inner.seed());
        let header = header_from_system(&self.inner, seed, None).map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("apsis record header: {e}"))
        })?;
        let hook = RecordHook::with_header(&path, header, policy)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("apsis record open: {e}")))?
            .with_diagnostics(diagnostics);
        self.inner.hooks_mut().register(0, Box::new(hook));
        Ok(())
    }

    // ── Provenance ───────────────────────────────────────────────────────

    /// Reference list for the registered operator stack. Returns one
    /// dictionary per citation, in registration order, with keys:
    ///
    /// - ``crate_name`` (str)
    /// - ``crate_version`` (str)
    /// - ``commit_hash`` (str | None) — full SHA when the implementing
    ///   crate was built from a git checkout, ``None`` otherwise
    /// - ``doi`` (str | None) — bare DOI suffix (e.g. ``10.1086/153180``)
    /// - ``bibtex`` (str) — full BibTeX entry / entries for the
    ///   underlying paper(s)
    ///
    /// Operators that don't publish a citation (test fakes, internal
    /// tooling, default ``None``) are silently skipped. The returned
    /// list is empty when no operators are registered, or when none of
    /// them publish a citation.
    ///
    /// Stable across runs given the same operator stack — diff two
    /// outputs to confirm the dependency graph stayed bit-equal.
    fn citations(&self, py: Python<'_>) -> PyResult<PyObject> {
        let cites = self.inner.citations();
        let list = pyo3::types::PyList::empty(py);
        for c in cites {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("crate_name", c.crate_name)?;
            dict.set_item("crate_version", c.crate_version)?;
            dict.set_item("commit_hash", c.commit_hash)?;
            dict.set_item("doi", c.doi)?;
            dict.set_item("bibtex", c.bibtex)?;
            list.append(dict)?;
        }
        Ok(list.into())
    }

    /// Render the registered operator stack's citations as a
    /// human-readable provenance block. Suitable for paper supplementary
    /// material or for embedding in a snapshot file.
    ///
    /// The layout is identical to the Rust ``System::provenance()``
    /// renderer — call sites can diff two outputs (across runs, across
    /// machines, across language bindings) to confirm the dependency
    /// graph stayed bit-equal.
    fn provenance(&self) -> String {
        self.inner.provenance()
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
        CoreIntegratorKind::Mercurius => "mercurius",
        CoreIntegratorKind::WHFast => "whfast",
        CoreIntegratorKind::ImplicitMidpoint => "implicit_midpoint",
    }
}

/// Extract a `Vec<f64>` from any Python sequence-of-floats (NumPy array,
/// list, tuple). Owned buffer so the sampling loop can take a slice.
fn extract_times(obj: &Bound<'_, PyAny>) -> PyResult<Vec<f64>> {
    obj.extract::<Vec<f64>>().map_err(|_| {
        value_error(
            "times",
            format!(
                "expected a 1-D sequence of floats (NumPy array, list, tuple), got {}",
                obj.get_type().name().map(|s| s.to_string()).unwrap_or_else(|_| "<?>".into()),
            ),
        )
    })
}

/// Build `np.linspace(t_start, t_start + duration, n_samples)` as the
/// convenience-form expansion of `sample(duration=, n_samples=)`.
fn build_evenly_spaced_times(t_start: f64, duration: f64, n_samples: usize) -> PyResult<Vec<f64>> {
    if !duration.is_finite() || duration <= 0.0 {
        return Err(value_error(
            "duration",
            format!("expected a strictly positive finite float, got {duration}"),
        ));
    }
    if n_samples == 0 {
        return Err(value_error("n_samples", "expected at least 1, got 0"));
    }

    let stride = if n_samples == 1 { 0.0 } else { duration / (n_samples - 1) as f64 };
    Ok((0..n_samples).map(|i| t_start + stride * i as f64).collect())
}

/// Validate a `times` array: non-empty, finite, monotonically
/// non-decreasing, `times[0] >= current_t`. Equality is allowed — fine
/// `np.linspace` can hit float-equality ties at consecutive indices.
fn validate_times(times: &[f64], current_t: f64) -> PyResult<()> {
    if times.is_empty() {
        return Err(value_error("times", "expected at least one sample time, got an empty array"));
    }
    for (i, &t) in times.iter().enumerate() {
        if !t.is_finite() {
            return Err(value_error("times", format!("times[{i}] = {t} is not finite")));
        }
    }
    for i in 1..times.len() {
        if times[i] < times[i - 1] {
            return Err(value_error(
                "times",
                format!(
                    "times must be monotonically non-decreasing; \
                     times[{}] = {} < times[{}] = {}",
                    i,
                    times[i],
                    i - 1,
                    times[i - 1]
                ),
            ));
        }
    }
    if times[0] < current_t {
        return Err(value_error(
            "times",
            format!(
                "first sample time {} is before the system's current t = {}; \
                 cannot integrate backwards",
                times[0], current_t
            ),
        ));
    }
    Ok(())
}

/// Push one row of state into the flat trajectory buffers in
/// sample-major (row-major) order so the consumer's `Array2::from_shape_vec`
/// reshape is straight `(n_samples, n_bodies)` without a transpose.
fn record_state(sys: &CoreSystem, buf: &mut TrajectoryBuffers) {
    let m = sys.metrics();
    buf.t.push(m.t);
    buf.energy.push(m.total_energy);
    buf.dt.push(m.dt);
    buf.energy_drift.push(m.rel_energy_error);
    buf.lz_drift.push(m.rel_angular_momentum_error);
    for b in sys.bodies() {
        buf.x.push(b.pos_x);
        buf.y.push(b.pos_y);
        buf.z.push(b.pos_z);
        buf.vx.push(b.vel_x);
        buf.vy.push(b.vel_y);
        buf.vz.push(b.vel_z);
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PySystem>()?;
    Ok(())
}
