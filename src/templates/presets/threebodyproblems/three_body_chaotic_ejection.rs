use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody},
};

pub fn three_body_chaotic_ejection() -> Template {
    let eps = 0.001;

    Template {
        name: "Three Body - Chaotic Ejection",
        bodies: vec![
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([-1.0, 0.0]),
                velocity: [0.3, 0.4],
                material: Material::Rocky,
            },
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([1.0, 0.0]),
                velocity: [0.3, 0.4 + eps],
                material: Material::Rocky,
            },
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([0.0, 0.1]),
                velocity: [-0.6, -0.8],
                material: Material::Rocky,
            },
        ],
        scale: 1.0,
    }
}
