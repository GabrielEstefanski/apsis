//! Boundary helpers for the Rust↔Python translation layer.
//!
//! This module owns the small set of utilities every binding wrapper
//! reaches for: error construction with researcher-friendly messages,
//! 3-vector parsing for position/velocity kwargs, and the canonical
//! string-to-enum normalisation policy applied uniformly across
//! [`crate::body`], [`crate::system`], and elsewhere.
//!
//! Centralising these here is what keeps the per-wrapper modules thin:
//! a `#[pymethods]` body that already delegates to a single [`apsis`]
//! call only needs `convert::xyz_triple(...)?` and `convert::value_error(...)`
//! at the boundary, never its own ad-hoc parsing or formatting.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// Build a `ValueError` with a structured message in the form
/// `"<context>: <detail>"`. Used at the boundary whenever a Python
/// argument fails domain validation (unknown enum variant, wrong
/// dimensionality, negative mass, etc.).
///
/// Researchers reading the traceback see the class/function they
/// called as the context and the specific failure as the detail —
/// shorter than a Rust `Result<_, ApsisError>` chain and immediately
/// actionable in a notebook.
pub(crate) fn value_error(context: &str, detail: impl AsRef<str>) -> PyErr {
    PyValueError::new_err(format!("{context}: {}", detail.as_ref()))
}

/// Extract a 3-component vector from a Python object, accepting either
/// a 2-element sequence `(x, y)` — treated as planar input with
/// implicit `z = 0` — or a 3-element sequence `(x, y, z)`.
///
/// Used for `position` and `velocity` kwargs across the binding.
/// Keeping the parser in one place ensures both `Body` and `System`
/// reject the same set of malformed inputs the same way; there is no
/// per-wrapper drift in what "a position vector" means. Planar Python
/// code (`position=(0.5, 0.0)`) keeps working without change;
/// out-of-plane scenarios opt in by passing the third component.
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

/// Normalise a Python string to the lowercase canonical slug expected
/// by `apsis`'s `FromStr` impls. Trims whitespace and lowercases ASCII;
/// non-ASCII characters pass through untouched so a future preset
/// named e.g. `"trappist-1"` can be matched on its hyphenated form.
///
/// This is the single normalisation policy applied to every string-typed
/// enum kwarg in the binding (`integrator`, `template`, ...). Variants
/// of capitalisation that researchers naturally reach for at a REPL
/// — `"IAS15"`, `"Ias15"`, `"ias15"` — all collapse to the same slug
/// before lookup, so the binding is liberal in what it accepts without
/// the underlying enum table needing to know about Python conventions.
pub(crate) fn slugify(input: &str) -> String {
    input.trim().to_ascii_lowercase()
}
