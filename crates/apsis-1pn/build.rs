//! Capture the build-time git commit hash so [`apsis_1pn::PostNewtonian1PN`]'s
//! [`Operator::citation`](apsis::physics::integrator::Operator::citation)
//! can pin the implementation to a specific source state.
//!
//! - Sets `APSIS_1PN_GIT_COMMIT` to the full SHA when this crate is built
//!   from a git checkout (the common case).
//! - Sets `APSIS_1PN_GIT_COMMIT` to the empty string when `git rev-parse`
//!   is unavailable or the source is not a git working tree (tarball
//!   distribution, vendored source, sandboxed CI). The runtime treats an
//!   empty hash as "no commit known" and renders the citation without a
//!   commit line.
//!
//! Re-runs only when `HEAD` changes, not on every source edit, so the
//! incremental build cost is one extra `git rev-parse` per branch
//! switch / commit.

use std::process::Command;

fn main() {
    let hash = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_default();

    println!("cargo:rustc-env=APSIS_1PN_GIT_COMMIT={hash}");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=build.rs");
}
