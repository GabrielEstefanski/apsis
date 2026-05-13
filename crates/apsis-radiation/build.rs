//! Capture the build-time git commit hash so this crate's
//! [`Operator::citation`](apsis::physics::integrator::Operator::citation)
//! impls can pin the implementation to a specific source state.
//!
//! Sets `APSIS_RADIATION_GIT_COMMIT` to the full SHA when built from a
//! git checkout, and to the empty string otherwise (tarball, vendored
//! source, sandboxed CI). The runtime treats an empty hash as "no
//! commit known" and renders the citation without a commit line.
//!
//! Re-runs only when `HEAD` changes, not on every source edit.

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

    println!("cargo:rustc-env=APSIS_RADIATION_GIT_COMMIT={hash}");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=build.rs");
}
