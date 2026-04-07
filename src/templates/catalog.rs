use crate::templates::{Template, category::TemplateCategory, presets::*};

pub struct TemplateEntry {
    pub name: &'static str,
    pub category: TemplateCategory,
    pub build: fn() -> Template,
}

pub const TEMPLATES: &[TemplateEntry] = &[
    // ── Single bodies ──────────────────────────────────────────────────────── //
    TemplateEntry { name: "Star",         category: TemplateCategory::Bodies,  build: star         },
    TemplateEntry { name: "Brown Dwarf",  category: TemplateCategory::Bodies,  build: brown_dwarf  },
    TemplateEntry { name: "Gas Giant",    category: TemplateCategory::Bodies,  build: gas_giant    },
    TemplateEntry { name: "Rocky Planet", category: TemplateCategory::Bodies,  build: rocky_planet },
    // ── Multi-body systems ─────────────────────────────────────────────────── //
    TemplateEntry { name: "Simple",         category: TemplateCategory::Systems, build: simple_system      },
    TemplateEntry { name: "Star + Giant",   category: TemplateCategory::Systems, build: star_gas_giant     },
    TemplateEntry { name: "Binary Stars",   category: TemplateCategory::Systems, build: binary_star        },
    TemplateEntry { name: "Star + Comp.",   category: TemplateCategory::Systems, build: star_companion     },
    TemplateEntry { name: "Inner System",   category: TemplateCategory::Systems, build: inner_solar_system },
];
