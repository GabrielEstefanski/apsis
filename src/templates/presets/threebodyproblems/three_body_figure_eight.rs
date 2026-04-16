use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_figure_eight() -> Template {
    Template {
        name: "Three Body - Figure Eight",
        description: "Three equal-mass bodies in the classic figure-eight periodic orbit.",
        bodies: vec![
            TemplateBody {
                name: Some("Body 1"),
                mass: 1.0,
                position: Some([-0.97000436, 0.24308753]),
                velocity: [0.4662036850, 0.4323657300],
                material: Material::Rocky,
                spin: 0.0,
            },
            TemplateBody {
                name: Some("Body 2"),
                mass: 1.0,
                position: Some([0.97000436, -0.24308753]),
                velocity: [0.4662036850, 0.4323657300],
                material: Material::Rocky,
                spin: 0.0,
            },
            TemplateBody {
                name: Some("Body 3"),
                mass: 1.0,
                position: Some([0.0, 0.0]),
                velocity: [-0.93240737, -0.86473146],
                material: Material::Rocky,
                spin: 0.0,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}
