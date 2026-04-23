//! Precision Run bottom panel — replaces the playbar when a
//! precision run is active.
//!
//! # Layout
//!
//! ```text
//! ┌───────────────────────────────────────────────────────────────┐
//! │ [■■■■■■■■■■■■░░░░░░░░░░░░░░░░░░]  62%   T = 3.72 / 6.00       │  primary
//! │ State · Running         substeps 128 · dt 1.8e-3 · 512 step/s │  secondary
//! │ [Pause]  [Abort]                                   [Commit]   │  controls
//! └───────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Scope — Block B1
//!
//! This file lands the *container* and primary readouts. The
//! control buttons (`Pause`, `Abort`, `Commit`) are rendered as
//! disabled placeholders here so the shape is honest; Block B3
//! wires them to [`PrecisionRunController`] intent methods. The
//! auto-correction notice (B6) and deferred-command pending-count
//! (B4) will append to this same panel.
//!
//! The panel is drawn only while the controller's state is not
//! `Idle`; the real-time playbar owns the bottom strip otherwise.

use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, DANGER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::precision_run::{RunOutcome, RunState, Telemetry};
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke};

/// Panel height. Taller than the playbar (36 px) because the panel
/// stacks three logical rows: progress, metrics, controls.
pub const PRECISION_PANEL_HEIGHT: f32 = 96.0;

impl SimulationApp {
    /// Draw the Precision Run bottom panel. The caller is responsible
    /// for choosing between this and the real-time playbar based on
    /// controller state — see [`ui::update`](crate::app::ui).
    pub(in crate::app) fn draw_precision_panel(&mut self, ctx: &egui::Context) {
        // Snapshot state + telemetry once under the lock so subsequent
        // UI math does not re-lock on every read.
        let (state, telemetry, t_sim_now) = {
            let ctrl_arc = self.system.precision_controller();
            let ctrl = ctrl_arc.lock().unwrap();
            (ctrl.state(), ctrl.telemetry().clone(), self.system.t())
        };

        egui::Panel::bottom("precision_panel")
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .stroke(Stroke::new(0.5, BORDER))
                    .inner_margin(egui::Margin::symmetric(14, 8)),
            )
            .default_size(PRECISION_PANEL_HEIGHT)
            .min_size(PRECISION_PANEL_HEIGHT)
            .max_size(PRECISION_PANEL_HEIGHT)
            .resizable(false)
            .show(ctx, |ui| {
                precision_panel_content(ui, state, &telemetry, t_sim_now);
            });
    }
}

fn precision_panel_content(
    ui: &mut egui::Ui,
    state: RunState,
    telemetry: &Telemetry,
    t_sim_now: f64,
) {
    ui.spacing_mut().item_spacing.y = 4.0;

    // ── Row 1: progress bar + primary labels ─────────────────────────────────
    draw_progress_row(ui, state, telemetry, t_sim_now);

    // ── Row 2: secondary metrics ─────────────────────────────────────────────
    draw_metrics_row(ui, state, telemetry);

    // ── Row 3: controls (disabled placeholders in B1) ────────────────────────
    draw_controls_row(ui, state);
}

fn draw_progress_row(
    ui: &mut egui::Ui,
    state: RunState,
    telemetry: &Telemetry,
    t_sim_now: f64,
) {
    let (t_target, t_start) = run_bounds(state).unwrap_or((0.0, 0.0));
    let fraction = progress_fraction(state, telemetry, t_sim_now, t_target, t_start);

    ui.horizontal(|ui| {
        // Progress bar takes most of the width; the right-side label is
        // fixed-width and holds the percentage + T = a / b readout.
        let avail = ui.available_width();
        let label_width = 220.0;
        let bar_width = (avail - label_width - 8.0).max(80.0);

        let bar_rect = ui.allocate_space(egui::vec2(bar_width, 12.0)).1;
        let painter = ui.painter_at(bar_rect);
        // Track.
        painter.rect_filled(bar_rect, 2.0, ACCENT_DIM);
        // Fill.
        let fill_w = (bar_rect.width() * fraction).clamp(0.0, bar_rect.width());
        let fill_rect =
            egui::Rect::from_min_size(bar_rect.min, egui::vec2(fill_w, bar_rect.height()));
        painter.rect_filled(fill_rect, 2.0, progress_fill_color(state));
        // Border.
        painter.rect_stroke(bar_rect, 2.0, Stroke::new(0.5, BORDER), egui::StrokeKind::Inside);

        ui.add_space(8.0);
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            ui.add_sized(
                [label_width, 14.0],
                egui::Label::new(
                    RichText::new(progress_label(state, telemetry, t_sim_now, t_target, t_start))
                        .monospace()
                        .size(11.0)
                        .color(TEXT_PRI),
                ),
            );
        });
    });
}

fn draw_metrics_row(ui: &mut egui::Ui, state: RunState, telemetry: &Telemetry) {
    ui.horizontal(|ui| {
        // State label (small tag on the left).
        let (label_text, label_color) = state_tag(state);
        ui.label(
            RichText::new("STATE")
                .size(9.0)
                .color(TEXT_DIM)
                .strong(),
        );
        ui.label(
            RichText::new(label_text)
                .size(11.0)
                .color(label_color)
                .strong(),
        );

        ui.add_space(12.0);
        metric_inline(ui, "substeps", &format!("{}", telemetry.substeps));
        metric_inline(ui, "dt", &fmt_dt(telemetry.current_dt));
        metric_inline(
            ui,
            "throughput",
            &format!("{:.0} step/s", telemetry.substeps_per_second_window),
        );
        metric_inline(
            ui,
            "sim rate",
            &format!("{:.2e} t/s", telemetry.sim_time_per_second_window),
        );

        // Acceptance/rejections: only show when there has been at
        // least one attempt — empty runs would render "0 / 0" noise.
        if telemetry.substeps > 0 || telemetry.rejections_total() > 0 {
            let accept = telemetry.acceptance_rate() * 100.0;
            let color = if accept >= 95.0 {
                SUCCESS
            } else if accept >= 80.0 {
                TEXT_SEC
            } else {
                DANGER
            };
            metric_inline_colored(
                ui,
                "accept",
                &format!("{:.1}%", accept),
                color,
            );
        }

        // Floor-hit indicator — only surfaces if there are any. Hidden
        // otherwise to keep the row quiet on healthy runs.
        if telemetry.degraded > 0 {
            metric_inline_colored(
                ui,
                "floor",
                &format!("×{}", telemetry.degraded),
                DANGER,
            );
        }
    });
}

fn draw_controls_row(ui: &mut egui::Ui, state: RunState) {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 6.0;

        // Left-aligned: Pause / Abort.
        let can_pause = matches!(state, RunState::Running { .. });
        let can_resume = matches!(state, RunState::Paused { .. });
        let can_abort = matches!(
            state,
            RunState::Running { .. } | RunState::Pausing { .. } | RunState::Paused { .. }
        );
        let is_completed = matches!(state, RunState::Completed { .. });
        let is_pausing = matches!(state, RunState::Pausing { .. });

        // Pause / Resume toggle — one button, label follows state.
        let (pause_label, pause_enabled) = if can_resume {
            ("Resume", true)
        } else if is_pausing {
            ("Pausing…", false)
        } else {
            ("Pause", can_pause)
        };
        placeholder_control_btn(ui, pause_label, pause_enabled, false);

        placeholder_control_btn(ui, "Abort", can_abort, true);

        // Right-aligned: Commit (only meaningful on completed runs).
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            let commit_label = if is_completed { "Close" } else { "Commit" };
            placeholder_control_btn(ui, commit_label, is_completed, false);
        });
    });
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn run_bounds(state: RunState) -> Option<(f64, f64)> {
    match state {
        RunState::Running { t_target, t_start, .. }
        | RunState::Pausing { t_target, t_start, .. }
        | RunState::Paused { t_target, t_start, .. }
        | RunState::Aborting { t_target, t_start } => Some((t_target, t_start)),
        _ => None,
    }
}

fn progress_fraction(
    state: RunState,
    telemetry: &Telemetry,
    t_sim_now: f64,
    t_target: f64,
    t_start: f64,
) -> f32 {
    match state {
        RunState::Completed { .. } => telemetry.last_progress_fraction,
        RunState::Idle => 0.0,
        _ => {
            let span = (t_target - t_start).max(f64::MIN_POSITIVE);
            (((t_sim_now - t_start) / span).clamp(0.0, 1.0)) as f32
        }
    }
}

fn progress_fill_color(state: RunState) -> Color32 {
    match state {
        RunState::Completed { outcome: RunOutcome::Reached } => SUCCESS,
        RunState::Completed { outcome: RunOutcome::UserAborted } => DANGER,
        RunState::Completed { outcome: RunOutcome::Errored } => DANGER,
        RunState::Aborting { .. } => DANGER,
        _ => ACCENT,
    }
}

fn progress_label(
    state: RunState,
    telemetry: &Telemetry,
    t_sim_now: f64,
    t_target: f64,
    t_start: f64,
) -> String {
    match state {
        RunState::Completed { .. } => format!(
            "{:.1}%   T = {:.3} (final)",
            telemetry.last_progress_fraction * 100.0,
            t_sim_now
        ),
        RunState::Idle => "—".to_string(),
        _ => {
            let span = (t_target - t_start).max(f64::MIN_POSITIVE);
            let pct = (((t_sim_now - t_start) / span).clamp(0.0, 1.0)) * 100.0;
            format!("{:.1}%   T = {:.3} / {:.3}", pct, t_sim_now, t_target)
        }
    }
}

fn state_tag(state: RunState) -> (&'static str, Color32) {
    match state {
        RunState::Idle => ("Idle", TEXT_SEC),
        RunState::Running { .. } => ("Running", ACCENT),
        RunState::Pausing { .. } => ("Pausing…", TEXT_SEC),
        RunState::Paused { .. } => ("Paused", ACCENT),
        RunState::Aborting { .. } => ("Aborting…", DANGER),
        RunState::Completed { outcome: RunOutcome::Reached } => ("Completed", SUCCESS),
        RunState::Completed { outcome: RunOutcome::UserAborted } => ("Aborted", DANGER),
        RunState::Completed { outcome: RunOutcome::Errored } => ("Errored", DANGER),
    }
}

fn metric_inline(ui: &mut egui::Ui, label: &str, value: &str) {
    metric_inline_colored(ui, label, value, TEXT_PRI);
}

fn metric_inline_colored(ui: &mut egui::Ui, label: &str, value: &str, value_color: Color32) {
    ui.label(
        RichText::new(label)
            .size(9.0)
            .color(TEXT_DIM)
            .strong(),
    );
    ui.label(
        RichText::new(value)
            .size(11.0)
            .monospace()
            .color(value_color),
    );
    ui.add_space(4.0);
}

/// Button rendered with the panel's control-strip styling. `enabled` gates
/// interaction; `danger` paints the stroke and hover with the danger palette.
/// B1 keeps click handling out of scope — the click is consumed and
/// discarded until B3 wires it to the controller.
fn placeholder_control_btn(
    ui: &mut egui::Ui,
    label: &str,
    enabled: bool,
    danger: bool,
) {
    let stroke_color = if danger { DANGER } else { BORDER };
    let text_color = if enabled {
        if danger { DANGER } else { TEXT_PRI }
    } else {
        TEXT_DIM
    };
    let btn = egui::Button::new(RichText::new(label).size(11.0).color(text_color))
        .fill(Color32::TRANSPARENT)
        .stroke(Stroke::new(0.5, stroke_color))
        .min_size(egui::vec2(82.0, 24.0));

    ui.add_enabled(enabled, btn);
}

fn fmt_dt(dt: f64) -> String {
    if dt == 0.0 {
        "—".into()
    } else if dt.abs() >= 1e-3 && dt.abs() < 1e3 {
        format!("{:.3}", dt)
    } else {
        format!("{:.2e}", dt)
    }
}
