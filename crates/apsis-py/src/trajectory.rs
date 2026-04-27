//! Dense trajectory recorder backed by NumPy arrays.
//!
//! The `Trajectory` Python class is the result of a single
//! [`crate::system::PySystem::sample`] call: a fixed-size record of
//! the simulation state at evenly spaced sample times, with positions,
//! velocities, total energy, and the sample times themselves laid out
//! as NumPy arrays ready for `matplotlib`, `pandas`, or any other
//! Python-side analysis tool.
//!
//! # Why pre-built NumPy arrays
//!
//! The alternative — holding raw `Vec<f64>` on the Rust side and
//! converting on each property access — would force a researcher who
//! plots `traj.x[:, 0]` and `traj.x[:, 1]` to pay two array
//! constructions for what is conceptually one read. A `Trajectory` is
//! produced once, materialised into NumPy at construction time, and
//! handed out as zero-copy `Bound<'py, PyArrayN>` views thereafter.
//! Memory cost is the same; the API call cost drops from O(n_samples
//! × n_bodies) per access to O(1).
//!
//! # Why a struct rather than a 6-tuple
//!
//! A `(t, x, y, vx, vy, energy)` return type would force every caller
//! to remember the order, which is exactly the positional-argument
//! pitfall called out in the binding's design discussion. The
//! `Trajectory` struct names every field, makes IDE autocomplete
//! discover them, and surfaces `n_samples` and `n_bodies` as
//! first-class properties so a researcher does not need to
//! `traj.t.shape[0]` every other line.

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
    n_samples: usize,
    n_bodies: usize,
}

impl PyTrajectory {
    /// Materialise a `Trajectory` from the row-major flat buffers
    /// populated by the sampling loop. Each 2-D array is reshaped
    /// from a `Vec<f64>` of length `n_samples × n_bodies`; the 1-D
    /// arrays come straight from `Vec<f64>` of length `n_samples`.
    ///
    /// This is the only entry point the binding offers for building a
    /// `Trajectory` — Python users construct one exclusively via
    /// `System.sample(...)`. Keeping construction internal preserves
    /// the invariant that every field has consistent dimensions.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build(
        py: Python<'_>,
        t: Vec<f64>,
        x: Vec<f64>,
        y: Vec<f64>,
        vx: Vec<f64>,
        vy: Vec<f64>,
        energy: Vec<f64>,
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
            t: t.into_pyarray(py).unbind(),
            x: to_2d("x", x)?,
            y: to_2d("y", y)?,
            vx: to_2d("vx", vx)?,
            vy: to_2d("vy", vy)?,
            energy: energy.into_pyarray(py).unbind(),
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

    fn __repr__(&self) -> String {
        format!("Trajectory(n_samples={}, n_bodies={})", self.n_samples, self.n_bodies)
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTrajectory>()?;
    Ok(())
}
