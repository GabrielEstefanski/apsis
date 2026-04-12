use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody},
};

pub fn three_body_chaotic_ejection() -> Template {
    let eps = 0.001;

    Template {
        name: "Three Body - Chaotic Ejection",
        description: "Three equal-mass bodies in a chaotic configuration leading to ejection.",
        bodies: vec![
            TemplateBody {
                mass: 1.0,
                position: Some([-1.0, 0.0]),
                velocity: [0.3, 0.4],
                material: Material::Rocky,
                spin: 0.0,
            },
            TemplateBody {
                mass: 1.0,
                position: Some([1.0, 0.0]),
                velocity: [0.3, 0.4 + eps],
                material: Material::Rocky,
                spin: 0.0,
            },
            TemplateBody {
                mass: 1.0,
                position: Some([0.0, 0.1]),
                velocity: [-0.6, -0.8],
                material: Material::Rocky,
                spin: 0.0,
            },
        ],
        display_scale: 1.0,
        suggested_dt: Some(0.001),
    }
}
