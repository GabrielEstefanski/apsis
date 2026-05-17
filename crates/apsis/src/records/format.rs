//! File-format constants for Apsis Records.
//!
//! These constants are written by [`crate::records::hook::RecordHook`]
//! and read by [`crate::records::reader::Record`]; living in a neutral
//! module keeps the reader from depending on the writer for what is
//! really a format-level contract.

/// File-format version embedded in the prefix after the magic bytes.
/// Bumping requires the `tests::schema_version` pin + an ADR update.
pub const FORMAT_VER: u16 = 1;

/// Four-byte file magic identifying an Apsis Record.
pub const MAGIC: &[u8; 4] = b"APSR";
