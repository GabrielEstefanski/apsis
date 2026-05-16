//! Apsis Record — binary reproducibility certificate for a simulation run.
//!
//! See `docs/adr/011-apsis-record.md` for the architectural decision and
//! the full format specification.

pub mod frame;
pub mod header;
pub mod hook;
pub mod policy;
pub mod provenance;

pub use header::Header;
pub use hook::{FORMAT_VER, MAGIC, RecordHook};
pub use policy::RecordPolicy;

// Submodules added in subsequent commits:
// pub mod reader;
