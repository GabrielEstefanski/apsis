//! Boundary helpers shared by the binding wrappers: error construction,
//! 3-vector kwarg parsing, and string-to-enum normalisation.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Build a `ValueError` with a `"<context>: <detail>"` message.
pub(crate) fn value_error(context: &str, detail: impl AsRef<str>) -> PyErr {
    PyValueError::new_err(format!("{context}: {}", detail.as_ref()))
}

/// Parse a `position`/`velocity` kwarg: a 2-element `(x, y)` (planar,
/// `z = 0`) or a 3-element `(x, y, z)` sequence.
pub(crate) fn xyz_triple(field: &str, obj: &Bound<'_, PyAny>) -> PyResult<(f64, f64, f64)> {
    if let Ok(triple) = obj.extract::<(f64, f64, f64)>() {
        return Ok(triple);
    }
    if let Ok((x, y)) = obj.extract::<(f64, f64)>() {
        return Ok((x, y, 0.0));
    }
    Err(value_error(
        field,
        format!(
            "expected a 2- or 3-element sequence of floats, got {}",
            obj.get_type().name().map(|s| s.to_string()).unwrap_or_else(|_| "<?>".into()),
        ),
    ))
}

/// Normalise a string kwarg to the lowercase slug `apsis`'s `FromStr`
/// expects: trims and lowercases ASCII; non-ASCII passes through (so
/// e.g. `trappist-1` matches on its hyphenated form).
pub(crate) fn slugify(input: &str) -> String {
    input.trim().to_ascii_lowercase()
}
