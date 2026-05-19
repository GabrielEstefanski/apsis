//! Free functions used by `System` that don't require `&self`.

use crate::domain::body::Body;

/// Default prefix when a body is added without an explicit name and
/// no preset hint is available (e.g. via [`System::add_body`]). Callers
/// that know which preset produced the body should pass that preset's
/// `display_name` through [`System::add_named_body`] instead.
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

/// Minimum pairwise separation across all body pairs.
///
/// Skipped (returns `f64::MAX`) when N < 2 or N > `N_CLOSENESS_THRESHOLD`,
/// to keep overhead bounded for large asteroid-belt simulations.
///
/// Distances are 3D (`dx² + dy² + dz²`); a previous 2D-only implementation
/// silently understated `r_min` for any inclined or out-of-plane pair.
pub(crate) fn compute_closeness(bodies: &[Body]) -> f64 {
    const N_CLOSENESS_THRESHOLD: usize = 512;

    if bodies.len() < 2 || bodies.len() > N_CLOSENESS_THRESHOLD {
        return f64::MAX;
    }

    let mut r_min = f64::MAX;

    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dx = bodies[i].pos_x - bodies[j].pos_x;
            let dy = bodies[i].pos_y - bodies[j].pos_y;
            let dz = bodies[i].pos_z - bodies[j].pos_z;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_min {
                r_min = r;
            }
        }
    }

    r_min
}
