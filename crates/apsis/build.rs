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
//! # Rerun policy
//!
//! The dirty-state check has to inspect the workdir, not just `.git/`
//! state, because plain edits to source files don't touch `.git/index`
//! until they are staged. The build script forces a rerun on every
//! cargo invocation by declaring a nonexistent sentinel as a
//! `rerun-if-changed` input — cargo emits a one-line warning about the
//! missing file and reruns the script. Two `git` invocations per build
//! is cheaper than a silently stale `-dirty` suffix in a paper-grade
//! record header.

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

    let rustc = Command::new(std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into()))
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_default();

    println!("cargo:rustc-env=APSIS_GIT_COMMIT={sha}");
    println!("cargo:rustc-env=APSIS_RUSTC_VERSION={rustc}");
    println!("cargo:rerun-if-changed=NONEXISTENT_BUILD_FORCE_RERUN");
}
