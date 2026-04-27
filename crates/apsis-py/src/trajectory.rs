//! Dense trajectory recorder backed by NumPy arrays.
//!
//! [`PyTrajectory`] is what [`crate::system::PySystem::sample`] returns:
//! a fixed-size record of positions, velocities, total energy, and
//! sample times, materialised once into NumPy and handed out as
//! zero-copy `Bound<'py, PyArrayN>` views. Field-named struct, not a
//! 6-tuple, so IDE autocomplete works and there's no positional pitfall.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;

/// Dense recording of a simulation interval, returned by
/// [`System.sample`](crate::system::PySystem::sample).
///
/// All arrays are NumPy `ndarray`s materialised once at construction
/// time. The 1-D arrays (`t`, `energy`) have shape `(n_samples,)`;
/// the 2-D arrays (`x`, `y`, `vx`, `vy`) have shape `(n_samples,
/// n_bodies)` with the body index on the second axis. A researcher
/// plotting body $k$'s trajectory does
/// `plt.plot(traj.x[:, k], traj.y[:, k])`; plotting the energy
/// drift is `plt.plot(traj.t, traj.energy)`.
///
/// `Trajectory` is immutable once constructed — there is no mutator
/// method on the Python side, and the underlying arrays are stored
/// behind shared references so reading them is side-effect-free.
#[pyclass(module = "apsis", name = "Trajectory", frozen)]
pub(crate) struct PyTrajectory {
    t: Py<PyArray1<f64>>,
    x: Py<PyArray2<f64>>,
    y: Py<PyArray2<f64>>,
    vx: Py<PyArray2<f64>>,
    vy: Py<PyArray2<f64>>,
    energy: Py<PyArray1<f64>>,
    dt: Py<PyArray1<f64>>,
    energy_drift: Py<PyArray1<f64>>,
    lz_drift: Py<PyArray1<f64>>,
    n_samples: usize,
    n_bodies: usize,
}

pub(crate) struct TrajectoryBuffers {
    pub t: Vec<f64>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub energy: Vec<f64>,
    pub dt: Vec<f64>,
    pub energy_drift: Vec<f64>,
    pub lz_drift: Vec<f64>,
}

impl PyTrajectory {
    /// Materialise a `Trajectory` from row-major flat buffers. Crate-private
    /// — Python users only get a `Trajectory` via `System.sample`, so the
    /// shape invariants stay enforced at the single producer.
    pub(crate) fn build(
        py: Python<'_>,
        b: TrajectoryBuffers,
        n_samples: usize,
        n_bodies: usize,
    ) -> PyResult<Self> {
        let to_2d = |label: &str, data: Vec<f64>| -> PyResult<Py<PyArray2<f64>>> {
            Array2::from_shape_vec((n_samples, n_bodies), data)
                .map(|arr| arr.into_pyarray(py).unbind())
                .map_err(|e| {
                    crate::convert::value_error(
                        label,
                        format!("internal shape error in trajectory buffer: {e}"),
                    )
                })
        };

        Ok(Self {
            t: b.t.into_pyarray(py).unbind(),
            x: to_2d("x", b.x)?,
            y: to_2d("y", b.y)?,
            vx: to_2d("vx", b.vx)?,
            vy: to_2d("vy", b.vy)?,
            energy: b.energy.into_pyarray(py).unbind(),
            dt: b.dt.into_pyarray(py).unbind(),
            energy_drift: b.energy_drift.into_pyarray(py).unbind(),
            lz_drift: b.lz_drift.into_pyarray(py).unbind(),
            n_samples,
            n_bodies,
        })
    }
}

#[pymethods]
impl PyTrajectory {
    /// Number of recorded samples (length of the time axis).
    #[getter]
    fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Number of bodies tracked in this trajectory. Equals the body
    /// count of the originating `System` at the time `sample` was
    /// called.
    #[getter]
    fn n_bodies(&self) -> usize {
        self.n_bodies
    }

    /// Sample times in simulation units. Shape `(n_samples,)`. The
    /// first entry is the system's `t` at the start of the sampling
    /// call; the last entry is approximately `start_t + duration`
    /// (within one adaptive sub-step for IAS15).
    #[getter]
    fn t<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.t.bind(py).clone()
    }

    /// Body $x$-coordinates. Shape `(n_samples, n_bodies)`.
    #[getter]
    fn x<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.x.bind(py).clone()
    }

    /// Body $y$-coordinates. Shape `(n_samples, n_bodies)`.
    #[getter]
    fn y<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.y.bind(py).clone()
    }

    /// Body $x$-velocities. Shape `(n_samples, n_bodies)`.
    #[getter]
    fn vx<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.vx.bind(py).clone()
    }

    /// Body $y$-velocities. Shape `(n_samples, n_bodies)`.
    #[getter]
    fn vy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.vy.bind(py).clone()
    }

    /// Total mechanical energy at each sample. Shape `(n_samples,)`.
    /// Useful for plotting the conservation diagnostic
    /// `(traj.energy - traj.energy[0]) / abs(traj.energy[0])` vs
    /// `traj.t`.
    #[getter]
    fn energy<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.energy.bind(py).clone()
    }

    /// Controller time-step at each sample. Shape `(n_samples,)`.
    /// Constant for fixed-step integrators; for IAS15 it traces the
    /// adaptive `dt_next` history — `plt.semilogy(traj.t, traj.dt)`
    /// reveals the controller's response to close encounters.
    #[getter]
    fn dt<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.dt.bind(py).clone()
    }

    /// Relative energy drift `δE / E₀` at each sample. Shape `(n_samples,)`.
    /// First entry is zero by definition (baseline). Already normalised —
    /// `plt.semilogy(traj.t, np.abs(traj.energy_drift))` is the standard
    /// conservation diagnostic plot without further arithmetic.
    #[getter]
    fn energy_drift<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.energy_drift.bind(py).clone()
    }

    /// Relative angular-momentum drift `δLz / Lz₀` at each sample.
    /// Shape `(n_samples,)`. First entry zero by baseline; falls back
    /// to absolute drift when `Lz₀` is below the core's numerical
    /// threshold (figure-8-like configurations with zero net angular
    /// momentum).
    #[getter]
    fn lz_drift<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lz_drift.bind(py).clone()
    }

    fn __repr__(&self) -> String {
        format!("Trajectory(n_samples={}, n_bodies={})", self.n_samples, self.n_bodies)
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTrajectory>()?;
    Ok(())
}
