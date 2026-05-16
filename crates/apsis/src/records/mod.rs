//! Apsis Record — binary reproducibility certificate for a simulation run.
//!
//! See `docs/adr/011-apsis-record.md` for the architectural decision and
//! the full format specification.

pub mod frame;
pub mod policy;

pub use policy::RecordPolicy;

// Submodules added in subsequent commits:
// pub mod header;
// pub mod provenance;
// pub mod hook;
// pub mod reader;
