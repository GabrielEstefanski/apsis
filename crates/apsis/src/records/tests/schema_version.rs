//! Pin the format version. Bumping this constant requires (a) a new ADR
//! documenting the format change, and (b) updating this assert.

use crate::records::format::FORMAT_VER;

#[test]
fn format_version_pinned_to_2() {
    assert_eq!(
        FORMAT_VER, 2,
        "FORMAT_VER bumped — add a new ADR documenting the format change \
         and update this pin."
    );
}
