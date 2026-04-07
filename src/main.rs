mod app;
mod core;
mod domain;
mod physics;
mod templates;

use app::ui::SimulationApp;
use core::adaptive::{DtAdaptationConfig, DtController, ThetaController};
use core::system::System;

fn main() {
    let system = System::new(
        12000,
        4,
        32,
        ThetaController::new(5e-4, 0.25, 1.1).with_initial_theta(0.6),
        DtController::new(DtAdaptationConfig {
            enabled: true,
            min_dt: 1e-5,
            max_dt: 1e-4,
            target_rel_energy_error: 5e-5,
            accel_epsilon: 0.04,
            grow_limit: 1.02,
            shrink_limit: 0.5,
            dt_slew_fraction: 0.2,
        }),
        20,
        vec![],
    );

    let app = SimulationApp::new(system);

    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "Gravity Simulator",
        native_options,
        Box::new(|_| Box::new(app)),
    )
    .unwrap();
}
