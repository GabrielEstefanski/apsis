use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, builders::circular_orbit},
};

pub fn simple_system() -> Template {
    let mut bodies = Vec::new();

    bodies.push(TemplateBody {
        mass: 1.0,
        radius: 0.05,
        position: None,
        velocity: [0.0, 0.0],
        material: Material::Star,
    });

    let (pos, vel) = circular_orbit(1.0, 0.5, 0.0);

    bodies.push(TemplateBody {
        mass: 0.001,
        radius: 0.01,
        position: None,
        velocity: vel,
        material: Material::Rocky,
    });

    Template {
        name: "Simple System",
        bodies,
        scale: 1.0,
    }
}
