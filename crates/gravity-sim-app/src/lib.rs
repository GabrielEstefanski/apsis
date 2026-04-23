//! Gravity simulator interactive shell — egui/wgpu UI on top of
//! [`gravity_sim_core`].
//!
//! The app crate owns everything visual: the main event loop, camera
//! controls, panels, theme, and GPU-side rendering. Physics, integrators,
//! persistence, and scenario templates live in `gravity-sim-core` and are
//! consumed read-only through the public API.

pub mod app;
pub mod render;

use crate::app::ui::SimulationApp;
use gravity_sim_core::core::system::System;
use gravity_sim_core::io;

/// Entry point shared by the `gravity-sim` binary.
///
/// Dispatches to headless batch mode when `--config <path>` is present;
/// otherwise launches the interactive eframe/egui GUI.
pub fn run() {
    // ── Headless batch mode: gravity-sim --config run.toml ─────────────────
    let args: Vec<String> = std::env::args().collect();
    if let Some(pos) = args.iter().position(|a| a == "--config") {
        let path = args.get(pos + 1).unwrap_or_else(|| {
            eprintln!("error: --config requires a path argument");
            std::process::exit(1);
        });
        let cfg = io::run_config::RunConfig::from_file(std::path::Path::new(path))
            .unwrap_or_else(|e| {
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
    let system = System::new(vec![])
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
