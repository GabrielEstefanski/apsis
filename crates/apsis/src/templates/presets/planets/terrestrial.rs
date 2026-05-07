pub fn terrestrial_planets() -> Template {
    let bodies = vec![
        // Mercúrio
        TemplateBody {
            mass: 0.055,
            radius: 0.38,
            position: None,
            velocity: [0.0, 0.0],
            class_override: None,
            material: Material::Rocky,
        },
        // Vênus
        TemplateBody {
            mass: 0.815,
            radius: 0.95,
            position: [-0.2, 0.0],
            velocity: [0.0, 0.0],
            class_override: None,
            material: Material::Rocky,
        },
        // Terra
        TemplateBody {
            mass: 1.0,
            radius: 1.0,
            position: [0.2, 0.0],
            velocity: [0.0, 0.0],
            class_override: None,
            material: Material::Rocky,
        },
        // Marte
        TemplateBody {
            mass: 0.107,
            radius: 0.53,
            position: [0.6, 0.0],
            velocity: [0.0, 0.0],
            class_override: None,
            material: Material::Rocky,
        },
    ];

    Template {
        name: "Terrestrial Planets",
        bodies,
        scale: 1.0,
    }
}
