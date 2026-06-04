use crate::domain::body::Body;
use crate::templates::Template;

/// Instantiate a template at the origin while preserving explicit body
/// names. Bodies without an explicit name fall back to the preset's
/// `display_name` so the system-level auto-numbering produces
/// `"Rocky 1"`, `"Comet 1"`, etc.
pub fn instantiate(template: &Template) -> Vec<Body> {
    instantiate_at(template, 0.0, 0.0)
}

/// Instantiate a template with its mass-weighted centroid translated
/// to `(cx, cy, 0)` in the simulation frame.
///
/// Body positions in the template are relative to an arbitrary
/// origin; this function shifts the whole system so its CoM lands at
/// the requested 2D drop point. Z is preserved per body (templates
/// with non-zero z keep their out-of-plane structure intact).
pub fn instantiate_at(template: &Template, cx: f64, cy: f64) -> Vec<Body> {
    let total_mass: f64 = template.bodies.iter().map(|t| t.mass).sum();

    let (com_x, com_y) = if total_mass > 0.0 {
        template.bodies.iter().fold((0.0, 0.0), |(ax, ay), t| {
            let [px, py, _pz] = t.position.unwrap_or([0.0, 0.0, 0.0]);
            (ax + t.mass * px, ay + t.mass * py)
        })
    } else {
        (0.0, 0.0)
    };
    let com_x = com_x / total_mass;
    let com_y = com_y / total_mass;

    let dx = cx - com_x;
    let dy = cy - com_y;

    template
        .bodies
        .iter()
        .map(|t| {
            let [px, py, pz] = t.position.unwrap_or([0.0, 0.0, 0.0]);
            let [vx, vy, vz] = t.velocity;

            let mut b = Body::from_preset(t.preset, t.mass)
                .at_3d(px + dx, py + dy, pz)
                .with_velocity_3d(vx, vy, vz);
            if let Some(class) = t.class_override {
                b = b.with_class(class);
            }

            // Per-body density override: templates that quote real
            // bodies (Earth, Sun, Jupiter, …) supply published values
            // here, replacing the preset's heuristic EOS so
            // `physical_radius` matches NASA fact-sheet values.
            if let Some(rho) = t.density {
                b = b.with_density(rho);
            }

            // Per-body Bond-albedo override — same rationale: real
            // bodies override the preset placeholder so the
            // photometry pipeline computes the published apparent
            // magnitude, not a class-typical guess.
            if let Some(a) = t.albedo {
                b = b.with_albedo(a);
            }

            // Fall back to the preset's display name when the
            // template author didn't pick one explicitly. The
            // system-level numerator turns repeated prefixes into
            // `"Asteroid 1"`, `"Asteroid 2"`, etc.
            let name =
                t.name.map(str::to_owned).unwrap_or_else(|| t.preset.display_name.to_owned());

            b.with_name(name)
        })
        .collect()
}
