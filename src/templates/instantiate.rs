use crate::domain::body::{Body, NamedBody};
use crate::templates::Template;

/// Instantiate a template at the origin while preserving explicit body names.
pub fn instantiate(template: &Template) -> Vec<NamedBody> {
    instantiate_at(template, 0.0, 0.0)
}

/// Instantiate a template with its mass-weighted centroid (CoM) at `(cx, cy)`.
///
/// Body positions in the template are relative to an arbitrary origin; this
/// function translates the whole system so its CoM lands exactly at the
/// requested world position.
pub fn instantiate_at(template: &Template, cx: f64, cy: f64) -> Vec<NamedBody> {
    let total_mass: f64 = template.bodies.iter().map(|t| t.mass).sum();

    let (com_x, com_y) = if total_mass > 0.0 {
        template.bodies.iter().fold((0.0, 0.0), |(ax, ay), t| {
            let [px, py] = t.position.unwrap_or([0.0, 0.0]);
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
            let [px, py] = t.position.unwrap_or([0.0, 0.0]);

            let mut b =
                Body::new(px + dx, py + dy, t.velocity[0], t.velocity[1], t.mass, t.material);

            NamedBody { body: b, name: t.name.map(str::to_owned) }
        })
        .collect()
}
