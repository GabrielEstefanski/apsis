//! Operator citation provenance.
//!
//! Each registered operator declares the paper / derivation it
//! implements via [`Operator::citation`](crate::physics::integrator::Operator::citation).
//! [`crate::core::system::System::citations`] aggregates them across
//! the registered stack, and [`crate::core::system::System::provenance`]
//! renders a human-readable block suitable for paper supplementary
//! material or for embedding into a snapshot file.
//!
//! # Federation-thesis alignment
//!
//! The federation model treats every perturbation crate as a citable
//! scientific artifact: it has a paper, a DOI, a versioned Cargo
//! dependency, and (when the build is captured from a git checkout)
//! a commit hash. `Citation` carries those four pieces as static
//! data attached to the operator type, so the simulation's full
//! reference list can be read off the operator stack at runtime —
//! the dependency graph IS the references list.
//!
//! ## Reproducibility envelope
//!
//! `crate_version` + `commit_hash` together pin the implementation
//! to a specific source state. Two runs that report identical
//! provenance blocks ran the same Rust code, modulo platform-level
//! f64 variance the apsis::contract notes as out of scope.
//! `commit_hash` is `Option` because not every build comes from a
//! git checkout (CI from tarball, vendored source, etc.); the
//! supplying crate's `build.rs` decides how to populate it.

/// Reference card for an operator's underlying physics + the source
/// state of the implementing crate.
///
/// All fields are `&'static str` so the citation is zero-cost at
/// runtime and trivially embeddable in `Box<dyn Operator>`. Carry
/// rich BibTeX text in `bibtex` rather than splitting into structured
/// fields — paper.md / supplementary material consumers want the
/// raw entry, and operators with multiple references concatenate
/// them in one string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Citation {
    /// BibTeX entry (or entries, concatenated). Operators that
    /// derive from multiple references include all of them here.
    pub bibtex: &'static str,

    /// DOI of the primary reference, when available. Use the bare
    /// suffix (e.g. `"10.1093/mnras/stv1257"`), not a full URL.
    /// `None` for textbook references, internal documentation, or
    /// pre-publication code.
    pub doi: Option<&'static str>,

    /// Cargo package name implementing the operator (typically the
    /// crate's `CARGO_PKG_NAME` at build time).
    pub crate_name: &'static str,

    /// Crate version at build time (typically `CARGO_PKG_VERSION`).
    /// Pinned in the consumer's `Cargo.lock`; `crate_name @
    /// crate_version` is sufficient to reproduce the operator's
    /// behaviour bit-for-bit on a single platform.
    pub crate_version: &'static str,

    /// Git commit hash at build time, when the build is from a git
    /// checkout. `None` for tarball / vendored / non-git builds.
    /// Operators that want to expose this populate it in their crate's
    /// `build.rs` via `cargo:rustc-env=`.
    pub commit_hash: Option<&'static str>,
}

impl Citation {
    /// Construct a `Citation` with the calling crate's
    /// `CARGO_PKG_NAME` and `CARGO_PKG_VERSION` filled in. Convenience
    /// for the common case where the crate that defines the operator
    /// also publishes it. Use this in macro form so the env vars are
    /// resolved at the operator's crate compile site, not apsis core's.
    ///
    /// ```ignore
    /// # use apsis::physics::integrator::Citation;
    /// // In your operator crate:
    /// fn my_citation() -> Citation {
    ///     Citation {
    ///         bibtex: "@article{...}",
    ///         doi: Some("10.xxxx/yyyy"),
    ///         crate_name: env!("CARGO_PKG_NAME"),
    ///         crate_version: env!("CARGO_PKG_VERSION"),
    ///         commit_hash: option_env!("MY_CRATE_GIT_COMMIT"),
    ///     }
    /// }
    /// ```
    ///
    /// Kept as documentation rather than a constructor function so the
    /// `env!` macro expansion happens in the consumer crate's compile
    /// context, capturing that crate's name and version (not apsis
    /// core's).
    pub const PROVENANCE_RECIPE: &'static str = "see Citation rustdoc";
}

/// Render a list of citations as a human-readable multi-line text
/// block. Format is deterministic so consumers can diff two
/// provenance blocks across runs to confirm the dependency graph
/// stayed bit-equal.
///
/// Layout:
/// ```text
/// Provenance (3 operators):
///
///   apsis-1pn 0.1.0 (commit f2d8e91)
///     DOI: 10.1007/BF00769986
///     @article{anderson1975, ...}
///
///   apsis-j2 0.1.0
///     DOI: 10.xxxx/yyyy
///     @article{...}
/// ```
pub fn render_provenance(citations: &[Citation]) -> String {
    if citations.is_empty() {
        return "Provenance: no operators with citations registered.\n".to_string();
    }
    let mut out = format!(
        "Provenance ({} operator{}):\n\n",
        citations.len(),
        if citations.len() == 1 { "" } else { "s" }
    );
    for c in citations {
        let header = match c.commit_hash {
            Some(h) => {
                format!("  {} {} (commit {})\n", c.crate_name, c.crate_version, short_commit(h))
            },
            None => format!("  {} {}\n", c.crate_name, c.crate_version),
        };
        out.push_str(&header);
        if let Some(doi) = c.doi {
            out.push_str(&format!("    DOI: {doi}\n"));
        }
        for line in c.bibtex.lines() {
            out.push_str(&format!("    {line}\n"));
        }
        out.push('\n');
    }
    out
}

fn short_commit(hash: &str) -> &str {
    // First 7 chars is the standard short SHA. Empty string if the
    // build set commit_hash to "".
    hash.get(..7).unwrap_or(hash)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(crate_name: &'static str, commit: Option<&'static str>) -> Citation {
        Citation {
            bibtex: "@article{sample, year={2026}}",
            doi: Some("10.0000/sample"),
            crate_name,
            crate_version: "0.1.0",
            commit_hash: commit,
        }
    }

    #[test]
    fn empty_provenance_says_no_citations() {
        let s = render_provenance(&[]);
        assert!(s.contains("no operators"));
    }

    #[test]
    fn provenance_pluralizes() {
        let one = render_provenance(&[sample("apsis-1pn", None)]);
        assert!(one.contains("(1 operator):"));
        let two = render_provenance(&[sample("apsis-1pn", None), sample("apsis-j2", None)]);
        assert!(two.contains("(2 operators):"));
    }

    #[test]
    fn provenance_includes_short_commit_when_present() {
        let s = render_provenance(&[sample("apsis-1pn", Some("f2d8e91abcdef1234567890"))]);
        assert!(s.contains("commit f2d8e91"));
        // Long hash should not leak into the formatted output.
        assert!(!s.contains("f2d8e91abcdef1234567890"));
    }

    #[test]
    fn provenance_omits_commit_when_none() {
        let s = render_provenance(&[sample("apsis-1pn", None)]);
        assert!(!s.contains("commit "));
        assert!(s.contains("apsis-1pn 0.1.0"));
    }

    #[test]
    fn provenance_includes_doi_and_bibtex() {
        let s = render_provenance(&[sample("apsis-1pn", None)]);
        assert!(s.contains("DOI: 10.0000/sample"));
        assert!(s.contains("@article{sample"));
    }

    #[test]
    fn short_commit_handles_under_7_chars() {
        // Defensive: a build script might set a short hash already.
        assert_eq!(short_commit("abc"), "abc");
        assert_eq!(short_commit("1234567890abcdef"), "1234567");
        assert_eq!(short_commit(""), "");
    }
}
