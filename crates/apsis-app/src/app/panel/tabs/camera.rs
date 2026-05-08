//! Camera tool — view framing and scale semantics.
//!
//! Pure presentation: camera settings never modify body state. Covers zoom,
//! scale-mode semantics, and quick-view actions (fit, zero COM).

use crate::app::theme::{
    BORDER, TEXT_DIM, TEXT_PRI, TEXT_SEC, primary_btn, secondary_btn, section,
};
use crate::app::ui::{SemanticScaleMode, SimulationApp};
use eframe::egui::{self, RichText};

impl SimulationApp {
    pub(super) fn panel_tab_camera(&mut self, ui: &mut egui::Ui) {
        ui.add_space(2.0);
        ui.label(RichText::new("Camera").size(13.0).color(TEXT_PRI).strong());
        ui.label(RichText::new("Framing & scale. Purely visual.").size(10.0).color(TEXT_DIM));

        section(ui, "SCALE MODE");
        for (mode, blurb) in [
            (SemanticScaleMode::Physical, "True physical sizes (often invisible)"),
            (SemanticScaleMode::Comparative, "Relative sizes preserved, clamped"),
            (SemanticScaleMode::Illustrative, "Readable sizes, perceptual"),
        ] {
            let active = self.semantic_scale_mode == mode;
            let col = if active { TEXT_PRI } else { TEXT_SEC };
            let fill =
                if active { crate::app::theme::ACCENT_DIM } else { egui::Color32::TRANSPARENT };
            let resp = ui.add(
                egui::Button::new(
                    RichText::new(format!("{}  {}", if active { "●" } else { "○" }, mode.label()))
                        .size(11.0)
                        .color(col),
                )
                .fill(fill)
                .stroke(egui::Stroke::new(0.5, BORDER))
                .min_size(egui::vec2(ui.available_width(), 24.0))
                .corner_radius(4.0),
            );
            if resp.clicked() {
                self.semantic_scale_mode = mode;
            }
            resp.on_hover_text(blurb);
        }

        section(ui, "QUICK ACTIONS");
        if primary_btn(ui, "Fit view") {
            self.fit_to_view();
        }
        ui.add_space(4.0);
        if secondary_btn(ui, "Zero centre-of-mass velocity") {
            self.system.zero_com_velocity();
        }
    }
}
