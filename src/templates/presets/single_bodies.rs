use crate::domain::materials::Material;
use crate::templates::{Template, TemplateBody};

pub fn star() -> Template {
    Template {
        name: "Star",
        bodies: vec![TemplateBody {
            mass: 1.0,
            radius: 0.05,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Star,
        }],
        scale: 1.0,
    }
}

pub fn brown_dwarf() -> Template {
    Template {
        name: "Brown Dwarf",
        bodies: vec![TemplateBody {
            mass: 0.04,
            radius: 0.025,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::BrownDwarf,
        }],
        scale: 1.0,
    }
}

pub fn gas_giant() -> Template {
    Template {
        name: "Gas Giant",
        bodies: vec![TemplateBody {
            // Jupiter ≈ 9.5e-4 M_sun
            mass: 9.5e-4,
            radius: 0.008,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Gas,
        }],
        scale: 1.0,
    }
}

pub fn rocky_planet() -> Template {
    Template {
        name: "Rocky Planet",
        bodies: vec![TemplateBody {
            // Earth ≈ 3.0e-6 M_sun
            mass: 3.0e-6,
            radius: 0.002,
            position: Some([0.0, 0.0]),
            velocity: [0.0, 0.0],
            material: Material::Rocky,
        }],
        scale: 1.0,
    }
}
