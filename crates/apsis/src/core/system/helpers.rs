//! Free functions used by `System` that don't require `&self`.

use crate::domain::body::Body;

const MASS_TO_SOLAR: f64 = 1.0;
const RADIUS_TO_SOLAR: f64 = 1.0 / 0.00465;
const L_SUN: f64 = 1.0;

pub(crate) fn mass_to_solar() -> f64 {
    MASS_TO_SOLAR
}
pub(crate) fn radius_to_solar() -> f64 {
    RADIUS_TO_SOLAR
}
pub(crate) fn l_sun() -> f64 {
    L_SUN
}

/// Generate an auto-name for a new body given existing names.
/// Counts existing names that start with the material prefix and appends N+1.
pub(crate) fn auto_name(
    material: crate::domain::materials::Material,
    existing: &[String],
) -> String {
    let prefix = material.display_name();
    let count = existing.iter().filter(|n| n.starts_with(prefix)).count() + 1;
    format!("{prefix} {count}")
}

pub(crate) fn resolved_name(
    explicit: Option<String>,
    material: crate::domain::materials::Material,
    existing: &[String],
) -> String {
    explicit
        .map(|name| name.trim().to_owned())
        .filter(|name| !name.is_empty())
        .unwrap_or_else(|| auto_name(material, existing))
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
            let dx = bodies[i].x - bodies[j].x;
            let dy = bodies[i].y - bodies[j].y;
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
