mod app;
mod core;
mod domain;
mod io;
mod physics;
mod render;
mod templates;

use app::ui::SimulationApp;
use core::system::System;

fn main() {
    // ── Headless batch mode: gravity-sim --config run.toml ─────────────────────
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

    // ── Interactive GUI mode ───────────────────────────────────────────────────
    let theta = 0.6;
    let dt = 1e-4;
    let max_depth = 32;
    let trail_every = 1;

    let system = System::new(
        vec![],
        theta,
        dt,
        max_depth,
        trail_every,
    );

    let app = SimulationApp::new(system);

    let native_options = eframe::NativeOptions::default();
    eframe::run_native("Gravity Simulator", native_options, Box::new(|_| Ok(Box::new(app))))
        .unwrap();
}
