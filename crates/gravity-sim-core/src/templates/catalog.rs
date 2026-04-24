//! Legacy string-keyed catalog. Derived from [`TemplateKind`] so there is a
//! single source of truth for the preset list — add a variant to the enum
//! and this table picks it up automatically at runtime.
//!
//! New code should prefer [`TemplateKind`] directly; this table exists for
//! the interactive app's template picker and for string-keyed config paths.

use crate::templates::{Template, TemplateKind, category::TemplateCategory};
use std::sync::LazyLock;

/// Runtime view of a preset: name, category, and the [`TemplateKind`]
/// discriminant used to dispatch to the builder.
pub struct TemplateEntry {
    pub name: &'static str,
    pub category: TemplateCategory,
    pub kind: TemplateKind,
}

impl TemplateEntry {
    /// Build the template for this entry using the given seed.
    pub fn build(&self, seed: u64) -> Template {
        self.kind.build(seed)
    }
}

/// All built-in preset entries, in the canonical order of [`TemplateKind::all`].
///
/// Lazy-initialised on first access; the entries are derived from the enum,
/// so adding a `TemplateKind` variant extends this table automatically.
pub static TEMPLATES: LazyLock<Vec<TemplateEntry>> = LazyLock::new(|| {
    TemplateKind::all()
        .iter()
        .map(|&kind| TemplateEntry { name: kind.name(), category: kind.category(), kind })
        .collect()
});
