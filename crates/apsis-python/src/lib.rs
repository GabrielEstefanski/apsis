//! `apsis._native` — PyO3 cdylib for the apsis Python distribution.

use pyo3::prelude::*;

mod body;
mod convert;
mod errors;
mod integrator;
mod operators;
mod perturbation;
mod records;
mod stats;
mod system;
mod trajectory;
mod units;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    body::register(m)?;
    errors::register(m)?;
    integrator::register(m)?;
    records::register(m)?;
    stats::register(m)?;
    system::register(m)?;
    trajectory::register(m)?;
    units::register(m)?;

    #[cfg(feature = "gr")]
    operators::gr::register(m)?;
    #[cfg(feature = "radiation")]
    operators::radiation::register(m)?;
    #[cfg(feature = "central")]
    operators::central::register(m)?;

    Ok(())
}
