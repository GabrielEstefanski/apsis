//! Capture the build-time git commit hash for [`Operator::citation`].
//! Sets `APSIS_CENTRAL_GIT_COMMIT` to the full SHA when built from a
//! git checkout, empty otherwise (tarball, vendored source). Re-runs
//! only when `HEAD` changes.

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

    println!("cargo:rustc-env=APSIS_CENTRAL_GIT_COMMIT={hash}");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=build.rs");
}
