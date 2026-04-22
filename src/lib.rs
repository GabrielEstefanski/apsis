//! Gravity Simulator — 2D N-body gravitational simulation.
//!
//! This crate exposes a `lib` + `bin` dual target:
//!
//! * **Library** (this file): the simulation core (physics, integrators,
//!   domain types) and the interactive app shell. Consumed by the
//!   `gravity-sim` binary and by `benches/` for Criterion-driven
//!   performance harnesses with versioned numerical baselines.
//! * **Binary** (`src/main.rs`): a thin entry point that delegates to
//!   [`run`] for both headless batch (`--config`) and interactive modes.
//!
//! The split was introduced so benches can link against the integrator
//! code directly (Criterion requires a library target), while keeping
//! the binary footprint unchanged.

pub mod app;
pub mod core;
pub mod domain;
pub mod io;
pub mod physics;
pub mod render;
pub mod templates;

use crate::app::ui::SimulationApp;
use crate::core::system::System;

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
    let theta = 0.6;
    let dt = 1e-4;
    let max_depth = 32;
    let trail_every = 1;

    let system = System::new(vec![], theta, dt, max_depth, trail_every);

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
