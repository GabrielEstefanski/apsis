//! Smoke test proving `egui_kittest` is wired up and can drive a widget
//! through the kittest accessibility tree. Not a behavioural test of
//! `SimulationApp` — those land alongside the F1-F6 UI redesign so they
//! can be written against the new layout rather than the legacy one.

use std::cell::Cell;

use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;

#[test]
fn harness_can_drive_a_button_click() {
    let clicked = Cell::new(false);
    let mut harness = Harness::new_ui(|ui| {
        if ui.button("press").clicked() {
            clicked.set(true);
        }
    });

    harness.get_by_label("press").click();
    harness.run();

    assert!(clicked.get(), "kittest harness failed to register click on labelled button");
}
