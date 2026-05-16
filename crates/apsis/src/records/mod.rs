//! Apsis Record — binary reproducibility certificate for a simulation run.
//!
//! See `docs/adr/011-apsis-record.md` for the architectural decision and
//! the full format specification.

pub mod frame;
pub mod header;
pub mod hook;
pub mod policy;
pub mod provenance;
pub mod reader;

pub use header::Header;
pub use hook::{FORMAT_VER, MAGIC, RecordHook};
pub use policy::RecordPolicy;
pub use reader::{Record, RecordError};

#[cfg(test)]
mod tests;
