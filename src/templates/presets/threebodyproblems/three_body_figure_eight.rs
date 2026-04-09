use crate::{
    domain::materials::Material,
    templates::{Template, TemplateBody},
};

pub fn three_body_figure_eight() -> Template {
    Template {
        name: "Three Body - Figure Eight",
        bodies: vec![
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([-0.97000436, 0.24308753]),
                velocity: [0.4662036850, 0.4323657300],
                material: Material::Rocky,
            },
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([0.97000436, -0.24308753]),
                velocity: [0.4662036850, 0.4323657300],
                material: Material::Rocky,
            },
            TemplateBody {
                mass: 1.0,
                radius: 0.006,
                position: Some([0.0, 0.0]),
                velocity: [-0.93240737, -0.86473146],
                material: Material::Rocky,
            },
        ],
        scale: 1.0,
    }
}
