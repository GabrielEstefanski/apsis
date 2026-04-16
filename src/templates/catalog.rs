use crate::templates::{
    Template,
    category::TemplateCategory,
    presets::{hierachical::simple_three_body, *},
};

pub struct TemplateEntry {
    pub name: &'static str,
    pub category: TemplateCategory,
    pub build: fn() -> Template,
}

pub const TEMPLATES: &[TemplateEntry] = &[
    // ── Single bodies ──────────────────────────────────────────────────────── //
    TemplateEntry { name: "Star", category: TemplateCategory::Bodies, build: star },
    TemplateEntry { name: "Brown Dwarf", category: TemplateCategory::Bodies, build: brown_dwarf },
    TemplateEntry { name: "Gas Giant", category: TemplateCategory::Bodies, build: gas_giant },
    TemplateEntry { name: "Rocky Planet", category: TemplateCategory::Bodies, build: rocky_planet },
    // ── Multi-body systems ─────────────────────────────────────────────────── //
    TemplateEntry { name: "Binary Stars", category: TemplateCategory::Systems, build: binary_star },
    TemplateEntry {
        name: "Star + Comp.",
        category: TemplateCategory::Systems,
        build: star_companion,
    },
    TemplateEntry {
        name: "Solar System",
        category: TemplateCategory::Systems,
        build: solar_system,
    },
    TemplateEntry { name: "TRAPPIST-1", category: TemplateCategory::Systems, build: trappist_1 },
    TemplateEntry { name: "Kepler-36", category: TemplateCategory::Systems, build: kepler_36 },
    TemplateEntry {
        name: "Alpha Centauri AB",
        category: TemplateCategory::Systems,
        build: alpha_centauri_ab,
    },
    TemplateEntry { name: "HD 80606 System", category: TemplateCategory::Systems, build: hd_80606 },
    TemplateEntry {
        name: "Sun–Earth L4/L5",
        category: TemplateCategory::Systems,
        build: sun_earth_lagrange,
    },
    TemplateEntry {
        name: "Sun–Earth L1/L2/L3",
        category: TemplateCategory::Systems,
        build: sun_earth_unstable_lagrange,
    },
    TemplateEntry {
        name: "Jupiter Trojans",
        category: TemplateCategory::Systems,
        build: jupiter_trojans,
    },
    TemplateEntry {
        name: "Hierarchical",
        category: TemplateCategory::Systems,
        build: simple_three_body,
    },
    // ── Three-body problems ─────────────────────────────────────────────────── //
    TemplateEntry {
        name: "3-Body Chaotic Ejection",
        category: TemplateCategory::ThreeBodyProblems,
        build: three_body_chaotic_ejection,
    },
    TemplateEntry {
        name: "3-Body Figure Eight",
        category: TemplateCategory::ThreeBodyProblems,
        build: three_body_figure_eight,
    },
    TemplateEntry {
        name: "3-Body Lagrange Triangle",
        category: TemplateCategory::ThreeBodyProblems,
        build: three_body_lagrange_triangle,
    },
];
