//! Free functions used by `System` that don't require `&self`.

use crate::domain::body::Body;

/// Default prefix when a body is added without an explicit name and
/// no preset hint is available (e.g. via [`System::add_body`]). Spawn
/// UIs and template loaders that know which preset produced the body
/// should pass that preset's `display_name` through
/// [`System::add_named_body`] instead.
pub(crate) const DEFAULT_NAME_PREFIX: &str = "Body";

/// Generate an auto-name `"<prefix> N"` for a new body given the
/// names already in use. Counts existing names that start with the
/// requested prefix and appends `N + 1`.
pub(crate) fn auto_name(prefix: &str, existing: &[String]) -> String {
    let count = existing.iter().filter(|n| n.starts_with(prefix)).count() + 1;
    format!("{prefix} {count}")
}

/// Resolve a final body name from an optional explicit value and a
/// fallback prefix used when the explicit value is missing or blank.
pub(crate) fn resolved_name(
    explicit: Option<String>,
    fallback_prefix: &str,
    existing: &[String],
) -> String {
    explicit
        .map(|name| name.trim().to_owned())
        .filter(|name| !name.is_empty())
        .unwrap_or_else(|| auto_name(fallback_prefix, existing))
}

/// Compute the minimum pairwise separation and maximum effective softening
/// length over all body pairs.
///
/// Skipped (returns sentinels) when N < 2 or N > `N_CLOSENESS_THRESHOLD`,
/// to keep overhead bounded for large asteroid-belt simulations.
pub(crate) fn compute_closeness(bodies: &[Body]) -> (f64, f64) {
    const N_CLOSENESS_THRESHOLD: usize = 512;

    if bodies.len() < 2 || bodies.len() > N_CLOSENESS_THRESHOLD {
        return (f64::MAX, 0.0);
    }

    let mut r_min = f64::MAX;
    let mut soft_max = 0.0_f64;

    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let r = (dx * dx + dy * dy).sqrt();
            if r < r_min {
                r_min = r;
            }
            let eps2_ij = (bodies[i].softening * bodies[i].softening
                + bodies[j].softening * bodies[j].softening)
                * 0.5;
            let eps_ij = eps2_ij.sqrt();
            if eps_ij > soft_max {
                soft_max = eps_ij;
            }
        }
    }

    (r_min, soft_max)
}
