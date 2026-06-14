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

    /// One-sentence description of the operator's physics, used as
    /// the `note` field of the `@software` entry emitted by
    /// [`crate::core::system::System::cite`].
    pub description: Option<&'static str>,

    /// Canonical source-repository URL, used as the `url` field of
    /// the `@software` entry emitted by
    /// [`crate::core::system::System::cite`]. Per-crate rather than a
    /// workspace constant so a future spinoff crate stays honest
    /// without touching core.
    pub url: Option<&'static str>,

    /// Author / maintainer name used as the `author` field of the
    /// `@software` entry emitted by
    /// [`crate::core::system::System::cite`]. BibTeX convention
    /// `"Surname, G. B."`; multiple authors join with ` and `.
    pub author: Option<&'static str>,
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
    ///         description: Some("short one-sentence summary of the physics"),
    ///         url: Some("https://github.com/your-org/your-crate"),
    ///         author: Some("Lastname, F. M."),
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

/// Short BLAKE3 prefix used in `cite()` output — first 4 hex chars,
/// "...", last 4 hex chars. Keeps the `note` line under a single
/// terminal width while remaining distinctive across forks.
fn short_lock_hash(hash: &str) -> String {
    if hash.len() <= 11 {
        return hash.to_string();
    }
    format!("{}...{}", &hash[..4], &hash[hash.len() - 4..])
}

/// Render the registered operator stack as a `@software` BibTeX
/// block suitable for direct inclusion in a paper's `.bib`. One
/// entry per `(Citation, KernelRequirements)` tuple, in the order
/// the caller supplies — the caller owns dedupe by `crate_name`
/// (see [`crate::core::system::System::cite`]).
///
/// Entry shape per Anderson et al. (1975) / Will (1993) / Burns et
/// al. (1979) / Tamayo et al. (2019) operator stack:
///
/// ```bibtex
/// @software{apsis-1pn_0.1.0,
///   title   = {apsis-1pn},
///   version = {0.1.0},
///   commit  = {f2d8e91},
///   url     = {https://github.com/GabrielEstefanski/apsis},
///   note    = {First-post-Newtonian Schwarzschild correction.
///              Cargo.lock blake3: 7f2a...e3c1;
///              kernel_requirements: exact_and_smooth},
/// }
/// ```
///
/// `commit` is omitted when `Citation::commit_hash` is `None`;
/// `url` is omitted when `Citation::url` is `None`; the description
/// half of `note` is omitted when `Citation::description` is `None`.
pub fn render_cite_block(
    entries: &[(Citation, crate::physics::gravity::kernel::KernelRequirements)],
    lock_blake3: &str,
) -> String {
    let mut out = String::new();
    let lock_short = short_lock_hash(lock_blake3);
    for (i, (c, req)) in entries.iter().enumerate() {
        if i > 0 {
            out.push('\n');
        }
        out.push_str(&format!(
            "@software{{{name}_{ver},\n",
            name = c.crate_name,
            ver = c.crate_version
        ));
        if let Some(author) = c.author {
            out.push_str(&format!("  author  = {{{}}},\n", bibtex_escape(author)));
        }
        out.push_str(&format!("  title   = {{{}}},\n", bibtex_escape(c.crate_name)));
        out.push_str(&format!("  version = {{{}}},\n", bibtex_escape(c.crate_version)));
        if let Some(hash) = c.commit_hash {
            out.push_str(&format!("  commit  = {{{}}},\n", short_commit(hash)));
        }
        if let Some(url) = c.url {
            out.push_str(&format!("  url     = {{{}}},\n", bibtex_escape(url)));
        }
        let req_slug = kernel_requirements_slug(req);
        let desc_prefix = c
            .description
            .map(|d| format!("{}.\n             ", bibtex_escape(d)))
            .unwrap_or_default();
        out.push_str(&format!(
            "  note    = {{{desc_prefix}Cargo.lock blake3: {lock_short};\n             \
             kernel_requirements: {req_slug}}},\n"
        ));
        out.push_str("}\n");
    }
    out
}

/// Escape `{` and `}` so a value containing braces does not close the
/// surrounding BibTeX field early. Returns the input unchanged when
/// no escape is needed — common case for crate names, versions, URLs.
fn bibtex_escape(s: &str) -> std::borrow::Cow<'_, str> {
    if !s.contains(['{', '}']) {
        return std::borrow::Cow::Borrowed(s);
    }
    let mut out = String::with_capacity(s.len() + 4);
    for ch in s.chars() {
        match ch {
            '{' => out.push_str("\\{"),
            '}' => out.push_str("\\}"),
            c => out.push(c),
        }
    }
    std::borrow::Cow::Owned(out)
}

fn kernel_requirements_slug(req: &crate::physics::gravity::kernel::KernelRequirements) -> String {
    use crate::physics::gravity::kernel::{Continuity, Exactness};
    let e = req.required_exactness.map(|e| match e {
        Exactness::Exact => "exact",
        Exactness::Softened => "softened",
        Exactness::Modified => "modified",
    });
    let c = req.min_continuity.map(|c| match c {
        Continuity::C0 => "c0",
        Continuity::C1 => "c1",
        Continuity::C2 => "c2",
        Continuity::Smooth => "smooth",
    });
    match (e, c) {
        (None, None) => "unconstrained".to_string(),
        (Some(e), None) => e.to_string(),
        (None, Some(c)) => c.to_string(),
        (Some(e), Some(c)) => format!("{e}_and_{c}"),
    }
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
            description: Some("sample operator for tests"),
            url: Some("https://example.invalid/sample"),
            author: Some("Test, A."),
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

    // ── render_cite_block ────────────────────────────────────────────────────

    use crate::physics::gravity::kernel::{Continuity, Exactness, KernelRequirements};

    fn req_exact_and_smooth() -> KernelRequirements {
        KernelRequirements {
            required_exactness: Some(Exactness::Exact),
            min_continuity: Some(Continuity::Smooth),
        }
    }

    fn req_none() -> KernelRequirements {
        KernelRequirements { required_exactness: None, min_continuity: None }
    }

    fn full_citation() -> Citation {
        Citation {
            bibtex: "@article{anderson1975, ...}",
            doi: Some("10.1086/153779"),
            crate_name: "apsis-1pn",
            crate_version: "0.1.0",
            commit_hash: Some("f2d8e91abcdef1234567890"),
            description: Some("First-post-Newtonian Schwarzschild correction"),
            url: Some("https://github.com/GabrielEstefanski/apsis"),
            author: Some("Estefanski, G. B."),
        }
    }

    #[test]
    fn cite_block_empty_when_no_entries() {
        let block = render_cite_block(&[], "deadbeef".repeat(8).as_str());
        assert!(block.is_empty());
    }

    /// Locks the note-field wrap so a future collapse to single line
    /// surfaces here, not as a PDF overflow on re-render.
    #[test]
    fn cite_block_wraps_note_field_with_thirteen_space_indent() {
        let block =
            render_cite_block(&[(full_citation(), req_exact_and_smooth())], &"a".repeat(64));
        assert!(block.contains("\n             Cargo.lock blake3: "));
        assert!(block.contains("\n             kernel_requirements: "));
    }

    #[test]
    fn cite_block_renders_full_entry_with_paper_md_shape() {
        let block = render_cite_block(
            &[(full_citation(), req_exact_and_smooth())],
            "7f2a000000000000000000000000000000000000000000000000000000003c1",
        );
        assert!(block.starts_with("@software{apsis-1pn_0.1.0,\n"));
        assert!(block.contains("  author  = {Estefanski, G. B.},\n"));
        assert!(block.contains("  title   = {apsis-1pn},\n"));
        assert!(block.contains("  version = {0.1.0},\n"));
        assert!(block.contains("  commit  = {f2d8e91},\n"));
        assert!(block.contains("  url     = {https://github.com/GabrielEstefanski/apsis},\n"));
        assert!(block.contains("First-post-Newtonian Schwarzschild correction."));
        assert!(block.contains("kernel_requirements: exact_and_smooth"));
        assert!(block.contains("Cargo.lock blake3: 7f2a...03c1"));
        assert!(block.trim_end().ends_with('}'));
    }

    #[test]
    fn cite_block_omits_commit_when_none() {
        let mut c = full_citation();
        c.commit_hash = None;
        let block = render_cite_block(&[(c, req_none())], &"a".repeat(64));
        assert!(!block.contains("commit  ="));
    }

    #[test]
    fn cite_block_omits_url_when_none() {
        let mut c = full_citation();
        c.url = None;
        let block = render_cite_block(&[(c, req_none())], &"a".repeat(64));
        assert!(!block.contains("url     ="));
    }

    #[test]
    fn cite_block_note_skips_description_when_none() {
        let mut c = full_citation();
        c.description = None;
        let block = render_cite_block(&[(c, req_none())], &"a".repeat(64));
        assert!(block.contains("Cargo.lock blake3:"));
        assert!(!block.contains("First-post-Newtonian"));
        // No leading ". " from the dropped description.
        assert!(!block.contains(". Cargo.lock"));
    }

    #[test]
    fn cite_block_renders_two_entries_separated_by_blank_line() {
        let c1 = full_citation();
        let mut c2 = full_citation();
        c2.crate_name = "apsis-radiation";
        let block =
            render_cite_block(&[(c1, req_exact_and_smooth()), (c2, req_none())], &"a".repeat(64));
        assert!(block.contains("@software{apsis-1pn_0.1.0,"));
        assert!(block.contains("@software{apsis-radiation_0.1.0,"));
        assert!(block.contains("kernel_requirements: exact_and_smooth"));
        assert!(block.contains("kernel_requirements: unconstrained"));
        // Two entries → two open braces of @software{ }.
        assert_eq!(block.matches("@software{").count(), 2);
    }

    #[test]
    fn kernel_requirements_slug_covers_every_combo() {
        use crate::physics::gravity::kernel::{Continuity, Exactness, KernelRequirements};
        assert_eq!(kernel_requirements_slug(&req_none()), "unconstrained");
        assert_eq!(kernel_requirements_slug(&req_exact_and_smooth()), "exact_and_smooth");
        assert_eq!(
            kernel_requirements_slug(&KernelRequirements {
                required_exactness: Some(Exactness::Exact),
                min_continuity: None
            }),
            "exact"
        );
        assert_eq!(
            kernel_requirements_slug(&KernelRequirements {
                required_exactness: None,
                min_continuity: Some(Continuity::C1)
            }),
            "c1"
        );
        assert_eq!(
            kernel_requirements_slug(&KernelRequirements {
                required_exactness: Some(Exactness::Modified),
                min_continuity: Some(Continuity::C0)
            }),
            "modified_and_c0"
        );
    }

    #[test]
    fn short_lock_hash_truncates_long_hex() {
        let h = short_lock_hash("7f2a4c8e1d9b3a5fa0e3c1b2d4f6a8c0e2d4b6f8a0c2e4b6d8f0a2c4e6d8f001");
        assert_eq!(h, "7f2a...f001");
    }

    #[test]
    fn short_lock_hash_keeps_short_input_verbatim() {
        assert_eq!(short_lock_hash("abc"), "abc");
        assert_eq!(short_lock_hash("1234567890a"), "1234567890a");
    }

    /// `}` and `{` in description/url do not close the surrounding
    /// BibTeX field early. Guards against a third-party operator
    /// shipping a Citation with brace-bearing prose.
    #[test]
    fn cite_block_escapes_braces_in_description_and_url() {
        let c = Citation {
            bibtex: "@misc{x}",
            doi: None,
            crate_name: "apsis-fake",
            crate_version: "0.1.0",
            commit_hash: None,
            description: Some("contains } and { chars"),
            url: Some("https://x.example/{owner}/repo"),
            author: None,
        };
        let block = render_cite_block(&[(c, req_none())], &"a".repeat(64));
        assert!(block.contains("contains \\} and \\{ chars"));
        assert!(block.contains("https://x.example/\\{owner\\}/repo"));
    }

    #[test]
    fn bibtex_escape_passes_clean_input_through_without_alloc() {
        let clean = "no special chars";
        let escaped = bibtex_escape(clean);
        assert!(matches!(escaped, std::borrow::Cow::Borrowed(_)));
        assert_eq!(escaped, clean);
    }

    #[test]
    fn bibtex_escape_handles_braces() {
        assert_eq!(bibtex_escape("a{b}c"), "a\\{b\\}c");
    }
}
