// ── App-crate lint allowances ────────────────────────────────────────────────
//
// This crate is **not part of the library's validated public surface** (see
// the workspace README). It tracks egui/wgpu's API as a moving target and
// accumulates deprecation warnings whenever the upstream crates update a
// method name. Fixing them is scheduled as part of the routine egui version
// bumps, not per-commit.
//
// The allowances here are therefore deliberate technical debt scoped to the
// UI layer. The core library and the 1PN plugin have the full `-D warnings`
// lint gate applied; only this crate relaxes it.
#![allow(
    deprecated,
    dead_code,
    unused_variables,
    unused_imports,
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::excessive_precision,
    clippy::needless_late_init,
    clippy::manual_range_contains,
    clippy::manual_map,
    clippy::field_reassign_with_default,
    clippy::collapsible_if,
    clippy::let_and_return,
    clippy::comparison_chain,
    clippy::useless_vec,
    clippy::redundant_closure,
    clippy::map_flatten,
    clippy::unnecessary_map_or,
    clippy::search_is_some,
    clippy::missing_safety_doc,
    clippy::needless_collect,
    clippy::single_match,
    clippy::single_match_else,
    clippy::or_fun_call,
    clippy::nonminimal_bool,
    clippy::needless_bool,
    clippy::similar_names,
    clippy::unnecessary_to_owned,
    clippy::unnecessary_sort_by,
    clippy::unnecessary_unwrap,
    clippy::redundant_pattern_matching,
    clippy::doc_overindented_list_items,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::explicit_iter_loop,
    clippy::iter_nth_zero,
    clippy::iter_without_into_iter,
    clippy::should_implement_trait,
    clippy::new_without_default,
    clippy::while_let_on_iterator,
    clippy::needless_borrows_for_generic_args,
    clippy::extend_with_drain,
    clippy::if_same_then_else
)]

//! Gravity simulator interactive shell — egui/wgpu UI on top of
//! [`apsis`].
//!
//! The app crate owns everything visual: the main event loop, camera
//! controls, panels, theme, and GPU-side rendering. Physics, integrators,
//! persistence, and scenario templates live in `apsis` and are
//! consumed read-only through the public API.

pub mod app;
pub mod render;

use crate::app::ui::SimulationApp;
use apsis::core::system::System;
use apsis::io;
use apsis::units::UnitSystem;

/// Entry point shared by the `apsis` binary.
///
/// Dispatches to headless batch mode when `--config <path>` is present;
/// otherwise launches the interactive eframe/egui GUI.
pub fn run() {
    // ── Headless batch mode: apsis --config run.toml ─────────────────
    let args: Vec<String> = std::env::args().collect();
    if let Some(pos) = args.iter().position(|a| a == "--config") {
        let path = args.get(pos + 1).unwrap_or_else(|| {
            eprintln!("error: --config requires a path argument");
            std::process::exit(1);
        });
        let cfg =
            io::run_config::RunConfig::from_file(std::path::Path::new(path)).unwrap_or_else(|e| {
                eprintln!("error: failed to parse {path}: {e}");
                std::process::exit(1);
            });
        if let Err(e) = io::headless::run(&cfg) {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
        return;
    }

    // ── Interactive GUI mode ───────────────────────────────────────────────
    let system = System::new(vec![], UnitSystem::canonical())
        .with_theta(0.6)
        .with_dt(1e-4)
        .with_max_depth(32);

    let app = SimulationApp::new(system);

    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "Gravity Simulator",
        native_options,
        Box::new(|cc| {
            let mut fonts = egui::FontDefinitions::default();
            egui_phosphor::add_to_fonts(&mut fonts, egui_phosphor::Variant::Regular);
            egui_phosphor::add_to_fonts(&mut fonts, egui_phosphor::Variant::Fill);
            cc.egui_ctx.set_fonts(fonts);
            Ok(Box::new(app))
        }),
    )
    .unwrap();
}
