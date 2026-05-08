use crate::{
    domain::body_preset,
    templates::{Template, TemplateBody, UnitSystem},
};

pub fn three_body_figure_eight(_seed: u64) -> Template {
    Template {
        name: "Three Body - Figure Eight",
        description: "Three equal-mass bodies in the classic figure-eight periodic orbit.",
        bodies: vec![
            TemplateBody {
                name: Some("Body 1"),
                mass: 1.0,
                position: Some([-0.97000436, 0.24308753, 0.0]),
                velocity: [0.4662036850, 0.4323657300, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("Body 2"),
                mass: 1.0,
                position: Some([0.97000436, -0.24308753, 0.0]),
                velocity: [0.4662036850, 0.4323657300, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
            TemplateBody {
                name: Some("Body 3"),
                mass: 1.0,
                position: Some([0.0, 0.0, 0.0]),
                velocity: [-0.93240737, -0.86473146, 0.0],
                class_override: None,
                preset: &body_preset::ROCKY,
                density: None,
                albedo: None,
            },
        ],
        display_scale: 1.0,
        orbital_up: None,
        default_view_distance: None,
        suggested_dt: Some(0.001),
        units: UnitSystem::dimensionless(),
    }
}
