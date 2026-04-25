//! Discovery tool — print every built-in preset, its category, and body count.
//!
//! Run with:
//!
//! ```text
//! cargo run --example list_templates
//! ```
//!
//! Useful for "what's available?" when writing a new script or wiring a
//! config file. Doubles as a smoke test that every `TemplateKind` variant
//! actually builds.

use apsis::templates::TemplateKind;

fn main() {
    let mut rows: Vec<(TemplateKind, &'static str, &'static str, usize)> = TemplateKind::all()
        .iter()
        .map(|&k| {
            let tpl = k.build(0);
            let category = match k.category() {
                apsis::templates::TemplateCategory::Bodies => "body",
                apsis::templates::TemplateCategory::Systems => "system",
                apsis::templates::TemplateCategory::ThreeBodyProblems => "3-body",
            };
            (k, k.name(), category, tpl.bodies.len())
        })
        .collect();

    // Stable display: group by category, alphabetical within.
    rows.sort_by_key(|(_, name, cat, _)| (*cat, *name));

    println!("{:<28} {:<8} {:>6}", "preset", "kind", "N");
    println!("{}", "─".repeat(50));
    for (_, name, cat, n) in &rows {
        println!("{:<28} {:<8} {:>6}", name, cat, n);
    }
    println!();
    println!("{} presets total.", rows.len());
}
