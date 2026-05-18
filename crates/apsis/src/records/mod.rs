//! Apsis Record — binary reproducibility certificate for a simulation run.
//!
//! See `docs/adr/011-apsis-record.md` for the architectural decision and
//! the full format specification.

pub mod diff;
pub mod format;
pub mod frame;
pub mod header;
pub mod hook;
pub mod policy;
pub mod provenance;
pub mod reader;
pub mod resume;

pub use diff::{FrameStreamDiff, HeaderChange, RecordDiff};
pub use format::{FORMAT_VER, MAGIC};
pub use header::Header;
pub use hook::RecordHook;
pub use policy::{DiagnosticCadence, RecordPolicy};
pub use reader::{Record, RecordError};
pub use resume::{RestoreError, restore_into};

#[cfg(test)]
mod tests;
