//! Light confirmation modal shown when the user selects a
//! `ExecutionProfile::Precision` integrator from the UI.
//!
//! # Why it exists
//!
//! Selecting IAS15 (today the only Precision-profile integrator)
//! has visible consequences: the real-time playbar is replaced with
//! the Precision Run panel, the force model is auto-reconfigured to
//! direct O(N²), and interactive framerate is no longer guaranteed.
//! These are not hidden behaviours — they surface through panel
//! transitions and structured events — but the user deserves a
//! single heads-up before any of it happens.
//!
//! The modal is intentionally light: one sentence of explanation, a
//! session-scoped "Don't show again" checkbox, and two buttons
//! (Cancel / Continue). It does not attempt to teach or to defend
//! the choice — that's what the panel's Setup view does once the
//! user accepts.
//!
//! # Session-scoped `dont-show-again`
//!
//! The preference lives on [`SimulationApp::precision_confirmation_session_skip`]
//! and is reset on app restart. Persistent preferences would need a
//! settings-file write path; that's not in scope for the confirmation
//! itself. If the user wants to silence the modal forever, polishing
//! the preference-persistence layer is follow-up work.

use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, PANEL_BG, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use gravity_sim_core::physics::integrator::IntegratorKind;
use gravity_sim_core::physics::integrator::traits::ExecutionProfile;
use eframe::egui::{self, Color32, RichText, Stroke};

impl SimulationApp {
    /// Central helper for any code path that wants to change the
    /// integrator. Routes through the confirmation modal when the
    /// target is a Precision-profile integrator and the user has
    /// not opted out for this session.
    ///
    /// This is the only function call sites should use; talking to
    /// `physics_cfg.integrator` and `self.system.set_integrator`
    /// directly bypasses the modal and defeats the purpose.
    pub(in crate::app) fn request_integrator_change(&mut self, kind: IntegratorKind) {
        if self.physics_cfg.integrator == kind {
            return;
        }
        let needs_confirmation = kind.execution_profile() == ExecutionProfile::Precision
            && !self.precision_confirmation_session_skip;

        if needs_confirmation {
            self.precision_confirmation_pending = Some(kind);
        } else {
            self.apply_integrator_change(kind);
        }
    }

    fn apply_integrator_change(&mut self, kind: IntegratorKind) {
        self.physics_cfg.integrator = kind;
        self.system.set_integrator(kind);
    }

    /// Draw the precision-confirmation modal if one is pending. Does
    /// nothing otherwise. Call from the main update loop after the
    /// primary chrome has been painted — the modal overlays it.
    pub(in crate::app) fn draw_precision_confirmation_modal(&mut self, ctx: &egui::Context) {
        let Some(kind) = self.precision_confirmation_pending else { return };

        // Captured into locals so we can decide after the closure
        // returns whether to apply / cancel / remember-don't-show.
        let mut clicked_continue = false;
        let mut clicked_cancel = false;
        let mut dont_show = self.precision_confirmation_session_skip;

        egui::Window::new("Precision Mode")
            .id(egui::Id::new("precision_confirmation_modal"))
            .collapsible(false)
            .resizable(false)
            .min_width(420.0)
            .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .stroke(Stroke::new(0.5, BORDER))
                    .inner_margin(egui::Margin::symmetric(20, 16))
                    .outer_margin(egui::Margin::same(0)),
            )
            .show(ctx, |ui| {
                ui.set_width(420.0);

                ui.label(
                    RichText::new("PRECISION MODE")
                        .size(10.0)
                        .color(TEXT_DIM)
                        .strong(),
                );
                ui.add_space(4.0);
                ui.label(
                    RichText::new(format!("Switch to {}?", kind.label()))
                        .size(13.0)
                        .color(TEXT_PRI)
                        .strong(),
                );
                ui.add_space(10.0);

                ui.label(
                    RichText::new(
                        "This is a precision integrator. It runs accurate but slow — \
                         per-step wall time is unbounded in stiff regimes and the \
                         simulation will not render at interactive framerate while a run \
                         is in progress.",
                    )
                    .size(11.0)
                    .color(TEXT_SEC),
                );
                ui.add_space(4.0);
                ui.label(
                    RichText::new(
                        "The Precision Run panel will replace the playback controls. \
                         Body editing (add, remove, drag) is queued during runs and \
                         applied when you commit.",
                    )
                    .size(11.0)
                    .color(TEXT_SEC),
                );

                ui.add_space(14.0);
                ui.horizontal(|ui| {
                    ui.checkbox(
                        &mut dont_show,
                        RichText::new("Don't show again this session").size(10.5).color(TEXT_SEC),
                    );
                });

                ui.add_space(12.0);
                ui.horizontal(|ui| {
                    // Cancel (secondary)
                    let cancel = egui::Button::new(
                        RichText::new("Cancel").size(11.5).color(TEXT_SEC),
                    )
                    .fill(Color32::TRANSPARENT)
                    .stroke(Stroke::new(0.5, BORDER))
                    .min_size(egui::vec2(90.0, 26.0));
                    if ui.add(cancel).clicked() {
                        clicked_cancel = true;
                    }

                    ui.with_layout(
                        egui::Layout::right_to_left(egui::Align::Center),
                        |ui| {
                            let cont = egui::Button::new(
                                RichText::new("Continue").size(11.5).color(TEXT_PRI),
                            )
                            .fill(ACCENT_DIM)
                            .stroke(Stroke::new(1.0, ACCENT))
                            .min_size(egui::vec2(110.0, 26.0));
                            if ui.add(cont).clicked() {
                                clicked_continue = true;
                            }
                        },
                    );
                });
            });

        self.precision_confirmation_session_skip = dont_show;

        if clicked_continue {
            self.precision_confirmation_pending = None;
            self.apply_integrator_change(kind);
        } else if clicked_cancel {
            self.precision_confirmation_pending = None;
        }
    }
}
