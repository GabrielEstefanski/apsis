use crate::domain::body_preset;
use crate::templates::{Template, TemplateBody, UnitSystem};

pub fn star(_seed: u64) -> Template {
    Template {
        name: "Star",
        description: "Single star at rest.",
        bodies: vec![TemplateBody {
            name: Some("Star"),
            mass: 1.0,
            position: Some([0.0, 0.0, 0.0]),
            velocity: [0.0, 0.0, 0.0],
            class_override: None,
            preset: &body_preset::STAR,
        }],
        display_scale: 1.0,
        suggested_dt: None,
        units: UnitSystem::solar_au(),
    }
}

pub fn brown_dwarf(_seed: u64) -> Template {
    Template {
        name: "Brown Dwarf",
        description: "Single brown dwarf at rest.",
        bodies: vec![TemplateBody {
            name: Some("Brown Dwarf"),
            mass: 0.04,
            position: Some([0.0, 0.0, 0.0]),
            velocity: [0.0, 0.0, 0.0],
            class_override: None,
            preset: &body_preset::BROWN_DWARF,
        }],
        display_scale: 1.0,
        suggested_dt: None,
        units: UnitSystem::solar_au(),
    }
}

pub fn gas_giant(_seed: u64) -> Template {
    Template {
        name: "Gas Giant",
        description: "Single gas giant (Jupiter-like) at rest.",
        bodies: vec![TemplateBody {
            name: Some("Gas Giant"),
            mass: 9.5e-4, // ≈ Jupiter mass in solar units
            position: Some([0.0, 0.0, 0.0]),
            velocity: [0.0, 0.0, 0.0],
            class_override: None,
            preset: &body_preset::GAS,
        }],
        display_scale: 1.0,
        suggested_dt: None,
        units: UnitSystem::solar_au(),
    }
}

pub fn rocky_planet(_seed: u64) -> Template {
    Template {
        name: "Rocky Planet",
        description: "Single rocky planet (Earth-like) at rest.",
        bodies: vec![TemplateBody {
            name: Some("Rocky Planet"),
            mass: 3.0e-6, // ≈ Earth mass in solar units
            position: Some([0.0, 0.0, 0.0]),
            velocity: [0.0, 0.0, 0.0],
            class_override: None,
            preset: &body_preset::ROCKY,
        }],
        display_scale: 1.0,
        suggested_dt: None,
        units: UnitSystem::solar_au(),
    }
}
