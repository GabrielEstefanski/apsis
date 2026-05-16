//! Capture the build-time git commit hash so apsis records can pin the
//! source state of the core crate in their header (`apsis.git_sha`).
//!
//! Sets `APSIS_GIT_COMMIT` to the full SHA when this crate is built from
//! a git checkout, with a `-dirty` suffix when the working tree carries
//! uncommitted changes. Sets the empty string when `git rev-parse` is
//! unavailable or the source is not a git working tree (tarball
//! distribution, vendored source, sandboxed CI). The runtime treats an
//! empty hash as "no commit known".
//!
//! # Why the dirty check (elevation over operator-crate convention)
//!
//! The operator crates (apsis-1pn, apsis-radiation, apsis-central) capture
//! `APSIS_<CRATE>_GIT_COMMIT` for `Operator::citation` without a dirty
//! check. That is acceptable for human-facing citation strings, but the
//! Apsis Record header carries `apsis.git_sha` as part of a paper-grade
//! reproducibility claim: a reviewer holding `{record, Cargo.lock}` must
//! be able to identify the source state bit-exactly. A SHA captured from
//! a dirty working tree silently misidentifies the source. The `-dirty`
//! suffix makes the discrepancy visible. Aligning the operator crates to
//! the same standard is a future cleanup tracked separately.
//!
//! Re-runs only when `HEAD` or the index changes, not on every source
//! edit, so the incremental build cost is one extra `git rev-parse` +
//! one `git diff-index` per relevant change.

use std::process::Command;

fn main() {
    let mut sha = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_default();

    if !sha.is_empty() {
        let dirty = Command::new("git")
            .args(["diff-index", "--quiet", "HEAD", "--"])
            .status()
            .map(|s| !s.success())
            .unwrap_or(false);
        if dirty {
            sha.push_str("-dirty");
        }
    }

    println!("cargo:rustc-env=APSIS_GIT_COMMIT={sha}");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=../../.git/index");
    println!("cargo:rerun-if-changed=build.rs");
}
