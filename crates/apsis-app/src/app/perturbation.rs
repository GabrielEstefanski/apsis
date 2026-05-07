//! Bundled perturbation registry — the federation seam consumed by
//! the interactive shell.
//!
//! Each downstream perturbation crate publishes a
//! [`PerturbationDescriptor`] value; this module exposes the static
//! list of descriptors compiled into the current build. Adding a new
//! plugin (e.g. `apsis-j2`, `apsis-tidal`) requires:
//!
//! 1. A `Cargo.toml` dependency on the plugin crate.
//! 2. One entry in [`bundled_descriptors`].
//!
//! No `match` arms, no enum variants, no string-typed dispatch — the
//! shell never names a concrete perturbation type.

use apsis::physics::integrator::PerturbationDescriptor;

/// One slot in the user-facing perturbation registry.
///
/// `descriptor` carries the plugin metadata (name, description, kernel
/// preconditions, builder); `enabled` reflects the user's checkbox state
/// in the Config tab.
pub struct PerturbationCatalogEntry {
    pub descriptor: Box<dyn PerturbationDescriptor>,
    pub enabled: bool,
}

/// Descriptors for every perturbation plugin compiled into the current
/// `apsis-app` build. The list is whitelisted at compile time — runtime
/// plugin loading remains a separate, optional feature (Issue #30
/// federation path 1).
pub fn bundled_descriptors() -> Vec<Box<dyn PerturbationDescriptor>> {
    vec![Box::new(apsis_1pn::Descriptor)]
}

/// Build the default catalog (every bundled descriptor, all disabled).
pub fn default_catalog() -> Vec<PerturbationCatalogEntry> {
    bundled_descriptors()
        .into_iter()
        .map(|descriptor| PerturbationCatalogEntry { descriptor, enabled: false })
        .collect()
}
