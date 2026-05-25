//! File-format constants for Apsis Records.
//!
//! These constants are written by [`crate::records::hook::RecordHook`]
//! and read by [`crate::records::reader::Record`]; living in a neutral
//! module keeps the reader from depending on the writer for what is
//! really a format-level contract.

/// File-format version embedded in the prefix after the magic bytes.
/// Bumping requires the `tests::schema_version` pin + an ADR update.
/// v0.2 (ADR-012): added `Diagnostic` frame kind, `ResumeState` frame
/// kind, and `kernel.exactness` / `kernel.continuity` /
/// `apsis.rustc_version` / `apsis.generated_by` header fields.
pub const FORMAT_VER: u16 = 2;

/// Four-byte file magic identifying an Apsis Record.
pub const MAGIC: &[u8; 4] = b"APSR";
