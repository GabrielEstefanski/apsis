//! `apsis.records` — read-only Python view of `.apsis` records.
//!
//! Writing is intentionally Rust-only in v0.1 (the writer's `Header` includes
//! a BLAKE3 hash of `Cargo.lock` that the Rust binary computes at file
//! creation time). Python consumers open existing records, inspect provenance,
//! and iterate events / dense snapshots.

use apsis::records::{Record, frame};
use pyo3::prelude::*;

#[pyclass(module = "apsis", name = "Record", frozen)]
pub(crate) struct PyRecord {
    inner: Record,
}

#[pymethods]
impl PyRecord {
    #[new]
    fn new(path: String) -> PyResult<Self> {
        let inner = Record::open(&path).map_err(|e| {
            pyo3::exceptions::PyIOError::new_err(format!("apsis record open failed: {e}"))
        })?;
        Ok(Self { inner })
    }

    /// Raw TOML header as a string. Self-describing; parse with `tomllib`
    /// in Python 3.11+ for structured access.
    #[getter]
    fn header(&self) -> PyResult<String> {
        self.inner.header().to_toml().map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("header re-serialise: {e}"))
        })
    }

    /// List of event tuples in time order. Each tuple is one of:
    /// - `("collision", t, body_a, body_b, distance)`
    /// - `("escape",    t, body, radius)`
    fn events(&self, py: Python<'_>) -> PyResult<PyObject> {
        let list = pyo3::types::PyList::empty(py);
        let iter = self
            .inner
            .events()
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("events: {e}")))?;
        for ev in iter {
            let ev =
                ev.map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("event read: {e}")))?;
            match ev {
                frame::Event::Collision { t, body_a, body_b, distance } => {
                    list.append(("collision", t, body_a, body_b, distance))?;
                },
                frame::Event::Escape { t, body, radius } => {
                    list.append(("escape", t, body, radius))?;
                },
            }
        }
        Ok(list.into())
    }

    /// Number of dense snapshots in the record (initial bookend + per-policy
    /// snapshots + final bookend).
    fn snapshot_count(&self) -> PyResult<usize> {
        let iter = self
            .inner
            .dense()
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("dense: {e}")))?;
        Ok(iter.count())
    }

    /// Diagnostic frames in time order. Each tuple is
    /// `(t, d_energy_rel, d_lz_rel)`. Empty when the record was written
    /// without a diagnostic cadence.
    fn diagnostics(&self, py: Python<'_>) -> PyResult<PyObject> {
        let list = pyo3::types::PyList::empty(py);
        let iter = self
            .inner
            .diagnostics()
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(format!("diagnostics: {e}")))?;
        for d in iter {
            let d = d.map_err(|e| {
                pyo3::exceptions::PyIOError::new_err(format!("diagnostic read: {e}"))
            })?;
            list.append((d.t, d.d_energy_rel, d.d_lz_rel))?;
        }
        Ok(list.into())
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRecord>()?;
    Ok(())
}
