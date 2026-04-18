use crate::app::theme::{ACCENT, BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use eframe::egui::{self, Color32, RichText, Stroke};

struct Binding {
    keys: &'static str,
    description: &'static str,
}

const BINDINGS: &[(&str, &[Binding])] = &[
    (
        "SIMULATION",
        &[
            Binding { keys: "Space", description: "Play / Pause" },
            Binding { keys: "Ctrl+Z", description: "Undo last add / delete / edit" },
        ],
    ),
    (
        "CAMERA",
        &[
            Binding { keys: "F", description: "Fit all bodies into view" },
            Binding { keys: "Scroll", description: "Zoom in / out" },
            Binding { keys: "Drag", description: "Pan camera" },
        ],
    ),
    ("HELP", &[Binding { keys: "H", description: "Toggle this shortcuts guide" }]),
];

impl SimulationApp {
    pub(in crate::app) fn draw_shortcuts_modal(&mut self, ctx: &egui::Context) {
        if !self.show_shortcuts_modal {
            return;
        }

        let mut open = true;

        egui::Window::new("Keyboard Shortcuts")
            .id(egui::Id::new("shortcuts_modal"))
            .collapsible(false)
            .resizable(false)
            .min_width(320.0)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .open(&mut open)
            .show(ctx, |ui| {
                ui.set_width(320.0);

                ui.label(RichText::new("KEYBOARD SHORTCUTS").size(11.0).color(TEXT_PRI).strong());

                ui.add_space(2.0);
                ui.separator();
                ui.add_space(4.0);

                for (section_name, bindings) in BINDINGS {
                    ui.label(RichText::new(*section_name).size(9.5).color(TEXT_DIM).strong());
                    ui.add_space(2.0);

                    egui::Grid::new(*section_name).num_columns(2).spacing([12.0, 3.0]).show(
                        ui,
                        |ui| {
                            for b in *bindings {
                                // Key badge
                                ui.add(
                                    egui::Button::new(
                                        RichText::new(b.keys).size(10.0).color(ACCENT).monospace(),
                                    )
                                    .fill(Color32::TRANSPARENT)
                                    .stroke(Stroke::new(0.5, BORDER))
                                    .min_size(egui::vec2(80.0, 18.0)),
                                );
                                // Description
                                ui.label(RichText::new(b.description).size(10.0).color(TEXT_SEC));
                                ui.end_row();
                            }
                        },
                    );

                    ui.add_space(8.0);
                }

                ui.separator();
                ui.add_space(4.0);
                ui.label(
                    RichText::new("More shortcuts will be added here.")
                        .size(9.5)
                        .color(TEXT_DIM)
                        .italics(),
                );
            });

        if !open {
            self.show_shortcuts_modal = false;
        }
    }
}
