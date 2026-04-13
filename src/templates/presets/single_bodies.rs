use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody};

pub fn star() -> Template {
    Template {
        name: "Star",
        description: "Single star at rest.",
        bodies: vec![TemplateBody {
            name: Some("Star"),
            mass: 1.0,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Star,
            spin: 0.0,
        }],
        display_scale: 1.0,
        suggested_dt: None,
    }
}

pub fn brown_dwarf() -> Template {
    Template {
        name: "Brown Dwarf",
        description: "Single brown dwarf at rest.",
        bodies: vec![TemplateBody {
            name: Some("Brown Dwarf"),
            mass: 0.04,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::BrownDwarf,
            spin: 0.0,
        }],
        display_scale: 1.0,
        suggested_dt: None,
    }
}

pub fn gas_giant() -> Template {
    Template {
        name: "Gas Giant",
        description: "Single gas giant (Jupiter-like) at rest.",
        bodies: vec![TemplateBody {
            name: Some("Gas Giant"),
            mass: 9.5e-4, // ≈ Jupiter mass in solar units
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Gas,
            spin: 0.0,
        }],
        display_scale: 1.0,
        suggested_dt: None,
    }
}

pub fn rocky_planet() -> Template {
    Template {
        name: "Rocky Planet",
        description: "Single rocky planet (Earth-like) at rest.",
        bodies: vec![TemplateBody {
            name: Some("Rocky Planet"),
            mass: 3.0e-6, // ≈ Earth mass in solar units
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Rocky,
            spin: 0.0,
        }],
        display_scale: 1.0,
        suggested_dt: None,
    }
}
