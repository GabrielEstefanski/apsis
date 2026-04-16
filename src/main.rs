mod app;
mod core;
mod physics;
mod render;
mod templates;

use app::ui::SimulationApp;
use core::system::System;

fn main() {
    // Initial simulation parameters
    let theta = 0.6; // Barnes–Hut accuracy
    let dt = 1e-4; // Fixed time step
    let max_depth = 32; // Barnes–Hut tree depth
    // trail_every is kept for API compatibility but is no longer used in step();
    // the trail is now recorded once per rendered frame via push_trail().
    let trail_every = 1;

    let system = System::new(
        vec![], // start empty (or plug templates later)
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
