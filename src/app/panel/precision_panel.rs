//! Precision Run bottom panel — replaces the playbar whenever a
//! [`ExecutionProfile::Precision`] integrator is selected.
//!
//! # Two layouts, one panel
//!
//! The panel renders two distinct layouts depending on the
//! controller's state:
//!
//! * **Setup** (state `Idle`): the user has chosen a Precision
//!   integrator but no run is active yet. The panel shows the
//!   target-duration input and a Start button; the rest of the real-
//!   time transport is suppressed by design because starting a run
//!   is the one action that matters here.
//!
//! * **Run** (state not Idle): progress bar, telemetry row, and
//!   run-lifecycle controls (Pause/Resume, Abort, Commit/Close).
//!
//! Keeping both layouts in one panel preserves the user's spatial
//! expectation — the Precision panel always lives at the bottom of
//! the window — and avoids the whiplash of a modal popping up and
//! vanishing around the actual run.
//!
//! # Scope
//!
//! * B1 delivered the Run layout with disabled placeholder buttons.
//! * B3 (this file's current form) wires the buttons to
//!   [`PrecisionRunController`] intent methods, adds the Setup
//!   layout, and introduces the Start flow. Buttons are no longer
//!   placeholders — they mutate controller state.
//! * B4/B5/B6 will layer on: pending-command chip, notification
//!   feed, auto-correction notice. Those insertion points are
//!   called out below with `// TODO(Bn)` markers.

use crate::app::icons;
use crate::app::theme::{ACCENT, ACCENT_DIM, BORDER, DANGER, PANEL_BG, SUCCESS, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::precision_run::{RunOutcome, RunState, Telemetry};
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke};

/// Panel height when rendering the active-run layout. The Setup
/// layout reuses the same height so the bottom strip does not
/// resize when the run starts.
pub const PRECISION_PANEL_HEIGHT: f32 = 96.0;

/// Intent produced by the controls row. Applied after the egui
/// borrow on `self` is released so we do not hold a double
/// mutable borrow across the panel body.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ControlIntent {
    None,
    Start,
    Pause,
    Resume,
    Abort,
    Acknowledge,
    ClearQueue,
}

impl SimulationApp {
    /// Draw the Precision Run bottom panel. The caller
    /// ([`ui::update`](crate::app::ui)) decides when to invoke
    /// this vs. the real-time playbar.
    pub(in crate::app) fn draw_precision_panel(&mut self, ctx: &egui::Context) {
        let (state, telemetry, t_sim_now) = {
            let ctrl_arc = self.system.precision_controller();
            let ctrl = ctrl_arc.lock().unwrap();
            (ctrl.state(), ctrl.telemetry().clone(), self.system.t())
        };
        let pending = self.system.pending_edits_count();
        let force_is_direct = self.system.metrics().force_is_direct;

        let mut intent = ControlIntent::None;
        let duration_ref = &mut self.precision_run_duration;

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
            .show(ctx, |ui| match state {
                RunState::Idle => {
                    setup_content(
                        ui,
                        duration_ref,
                        t_sim_now,
                        pending,
                        force_is_direct,
                        &mut intent,
                    );
                }
                _ => {
                    run_content(ui, state, &telemetry, t_sim_now, pending, &mut intent);
                }
            });

        // Apply intent.
        if intent != ControlIntent::None {
            self.apply_precision_intent(intent, t_sim_now);
        }
    }

    fn apply_precision_intent(&mut self, intent: ControlIntent, t_sim_now: f64) {
        match intent {
            ControlIntent::None => {}
            ControlIntent::Start => {
                let ctrl_arc = self.system.precision_controller();
                let mut ctrl = ctrl_arc.lock().unwrap();
                let t_target = t_sim_now + self.precision_run_duration.max(0.0);
                ctrl.start(t_target, t_sim_now);
            }
            ControlIntent::Pause => self.system.precision_controller().lock().unwrap().request_pause(),
            ControlIntent::Resume => self.system.precision_controller().lock().unwrap().resume(),
            ControlIntent::Abort => self.system.precision_controller().lock().unwrap().request_abort(),
            ControlIntent::Acknowledge => {
                // Acknowledge bridges the Completed state back to
                // Idle; the outcome decides what happens with the
                // edits the user queued during the run.
                let outcome = {
                    let ctrl_arc = self.system.precision_controller();
                    let ctrl = ctrl_arc.lock().unwrap();
                    match ctrl.state() {
                        RunState::Completed { outcome } => Some(outcome),
                        _ => None,
                    }
                };
                match outcome {
                    Some(RunOutcome::Reached) => self.system.commit_pending_edits(),
                    Some(RunOutcome::UserAborted) | Some(RunOutcome::Errored) => {
                        self.system.clear_pending_edits();
                    }
                    None => {}
                }
                self.system.precision_controller().lock().unwrap().acknowledge();
            }
            ControlIntent::ClearQueue => self.system.clear_pending_edits(),
        }
    }
}

// ── Setup layout (state == Idle) ──────────────────────────────────────────────

fn setup_content(
    ui: &mut egui::Ui,
    duration: &mut f64,
    t_sim_now: f64,
    pending: usize,
    force_is_direct: bool,
    intent: &mut ControlIntent,
) {
    ui.spacing_mut().item_spacing.y = 6.0;

    ui.horizontal(|ui| {
        ui.label(
            RichText::new("PRECISION RUN")
                .size(10.0)
                .color(TEXT_DIM)
                .strong(),
        );
        ui.label(
            RichText::new("·")
                .size(10.0)
                .color(TEXT_DIM),
        );
        ui.label(
            RichText::new("ready")
                .size(11.0)
                .color(ACCENT)
                .strong(),
        );
        if force_is_direct {
            ui.label(RichText::new("·").size(10.0).color(TEXT_DIM));
            ui.label(
                RichText::new("force: direct O(N²) required for deterministic physics")
                    .size(10.0)
                    .color(TEXT_SEC),
            );
        } else {
            ui.label(
                RichText::new("— select a target simulation time and press Start")
                    .size(10.0)
                    .color(TEXT_SEC),
            );
        }

        if pending > 0 {
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                let clear = egui::Button::new(
                    RichText::new("Clear").size(10.0).color(TEXT_SEC),
                )
                .fill(Color32::TRANSPARENT)
                .stroke(Stroke::new(0.5, BORDER))
                .min_size(egui::vec2(54.0, 20.0));
                if ui.add(clear).clicked() {
                    *intent = ControlIntent::ClearQueue;
                }
                ui.add_space(6.0);
                ui.label(
                    RichText::new(format!("{} queued edit{}", pending, if pending == 1 { "" } else { "s" }))
                        .size(10.0)
                        .color(ACCENT),
                );
            });
        }
    });

    ui.horizontal(|ui| {
        ui.label(RichText::new("DURATION").size(9.0).color(TEXT_DIM).strong());
        // Read the current duration *before* constructing the DragValue
        // (which borrows `duration` mutably for the lifetime of the
        // builder) so `.speed(...)` can use the value without aliasing.
        let duration_speed = (*duration * 0.05).max(1e-3);
        let drag = egui::DragValue::new(duration)
            .speed(duration_speed)
            .range(1e-6..=1e9)
            .custom_formatter(|v, _| format_duration_compact(v))
            .custom_parser(|s| parse_duration_compact(s));
        ui.add(drag).on_hover_text(
            "Duration the run will advance the simulation by. The target \
             simulation time is resolved at Start as t + duration.",
        );

        ui.add_space(16.0);
        ui.label(RichText::new("TARGET").size(9.0).color(TEXT_DIM).strong());
        let t_target = t_sim_now + *duration;
        ui.label(
            RichText::new(format!("T = {:.3}", t_target))
                .size(11.0)
                .monospace()
                .color(TEXT_PRI),
        );
        ui.label(
            RichText::new(format!("(from {:.3})", t_sim_now))
                .size(10.0)
                .monospace()
                .color(TEXT_SEC),
        );

        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            let start_btn = egui::Button::new(
                RichText::new(format!("{}  Start", icons::PRECISION_START))
                    .size(12.0)
                    .color(TEXT_PRI),
            )
            .fill(ACCENT_DIM)
            .stroke(Stroke::new(1.0, ACCENT))
            .min_size(egui::vec2(110.0, 28.0));
            if ui.add(start_btn).clicked() {
                *intent = ControlIntent::Start;
            }
        });
    });

    // TODO(B4): queued-commands chip ("N pending changes · Clear") goes
    // here — Setup is one of the views where the chip is relevant
    // (pending changes enqueued during a prior run that is now
    // Acknowledged).
    // TODO(B6): auto-correction notice ("Force model switched to
    // direct O(N²)") also lands here if it fired since app start.
}

// ── Run layout (state != Idle) ────────────────────────────────────────────────

fn run_content(
    ui: &mut egui::Ui,
    state: RunState,
    telemetry: &Telemetry,
    t_sim_now: f64,
    pending: usize,
    intent: &mut ControlIntent,
) {
    ui.spacing_mut().item_spacing.y = 4.0;

    draw_progress_row(ui, state, telemetry, t_sim_now);
    draw_metrics_row(ui, state, telemetry, pending);
    draw_controls_row(ui, state, intent);
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
        let avail = ui.available_width();
        let label_width = 220.0;
        let bar_width = (avail - label_width - 8.0).max(80.0);

        let bar_rect = ui.allocate_space(egui::vec2(bar_width, 12.0)).1;
        let painter = ui.painter_at(bar_rect);
        painter.rect_filled(bar_rect, 2.0, ACCENT_DIM);
        let fill_w = (bar_rect.width() * fraction).clamp(0.0, bar_rect.width());
        let fill_rect =
            egui::Rect::from_min_size(bar_rect.min, egui::vec2(fill_w, bar_rect.height()));
        painter.rect_filled(fill_rect, 2.0, progress_fill_color(state));
        painter.rect_stroke(
            bar_rect,
            2.0,
            Stroke::new(0.5, BORDER),
            egui::StrokeKind::Inside,
        );

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

fn draw_metrics_row(ui: &mut egui::Ui, state: RunState, telemetry: &Telemetry, pending: usize) {
    ui.horizontal(|ui| {
        let (label_text, label_color) = state_tag(state);
        ui.label(RichText::new("STATE").size(9.0).color(TEXT_DIM).strong());
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

        if telemetry.substeps > 0 || telemetry.rejections_total() > 0 {
            let accept = telemetry.acceptance_rate() * 100.0;
            let color = if accept >= 95.0 {
                SUCCESS
            } else if accept >= 80.0 {
                TEXT_SEC
            } else {
                DANGER
            };
            metric_inline_colored(ui, "accept", &format!("{:.1}%", accept), color);
        }

        if telemetry.degraded > 0 {
            metric_inline_colored(
                ui,
                "floor",
                &format!("×{}", telemetry.degraded),
                DANGER,
            );
        }

        if pending > 0 {
            metric_inline_colored(ui, "queued", &format!("{}", pending), ACCENT);
        }
    });
}

fn draw_controls_row(ui: &mut egui::Ui, state: RunState, intent: &mut ControlIntent) {
    ui.horizontal(|ui| {
        ui.spacing_mut().item_spacing.x = 6.0;

        let is_running = matches!(state, RunState::Running { .. });
        let is_paused = matches!(state, RunState::Paused { .. });
        let is_pausing = matches!(state, RunState::Pausing { .. });
        let is_aborting = matches!(state, RunState::Aborting { .. });
        let is_completed = matches!(state, RunState::Completed { .. });
        let can_abort = is_running || is_pausing || is_paused;

        // Pause / Resume / Pausing… toggle.
        let (icon, label, enabled, produces) = if is_paused {
            (icons::PRECISION_RESUME, "Resume", true, ControlIntent::Resume)
        } else if is_pausing {
            (icons::PRECISION_PAUSE, "Pausing…", false, ControlIntent::None)
        } else if is_running {
            (icons::PRECISION_PAUSE, "Pause", true, ControlIntent::Pause)
        } else {
            (icons::PRECISION_PAUSE, "Pause", false, ControlIntent::None)
        };
        if control_btn(ui, icon, label, enabled, false).clicked() {
            *intent = produces;
        }

        if control_btn(
            ui,
            icons::PRECISION_ABORT,
            "Abort",
            can_abort && !is_aborting,
            true,
        )
        .clicked()
        {
            *intent = ControlIntent::Abort;
        }

        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            match state {
                RunState::Completed { outcome: RunOutcome::Reached } => {
                    if control_btn(ui, icons::PRECISION_COMMIT, "Commit", true, false).clicked() {
                        *intent = ControlIntent::Acknowledge;
                    }
                }
                RunState::Completed { .. } => {
                    if control_btn(ui, icons::PRECISION_CLOSE, "Close", true, false).clicked() {
                        *intent = ControlIntent::Acknowledge;
                    }
                }
                _ => {
                    control_btn(ui, icons::PRECISION_COMMIT, "Commit", false, false);
                }
            }
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

fn metric_inline_colored(
    ui: &mut egui::Ui,
    label: &str,
    value: &str,
    value_color: Color32,
) {
    ui.label(RichText::new(label).size(9.0).color(TEXT_DIM).strong());
    ui.label(RichText::new(value).size(11.0).monospace().color(value_color));
    ui.add_space(4.0);
}

fn control_btn(
    ui: &mut egui::Ui,
    icon: &str,
    label: &str,
    enabled: bool,
    danger: bool,
) -> egui::Response {
    let stroke_color = if danger { DANGER } else { BORDER };
    let text_color = if enabled {
        if danger { DANGER } else { TEXT_PRI }
    } else {
        TEXT_DIM
    };
    let btn = egui::Button::new(
        RichText::new(format!("{}  {}", icon, label))
            .size(11.0)
            .color(text_color),
    )
    .fill(Color32::TRANSPARENT)
    .stroke(Stroke::new(0.5, stroke_color))
    .min_size(egui::vec2(92.0, 24.0));

    ui.add_enabled(enabled, btn)
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

fn format_duration_compact(v: f64) -> String {
    if v.abs() >= 1e4 || (v.abs() < 1e-2 && v != 0.0) {
        format!("{:.3e}", v)
    } else {
        format!("{:.3}", v)
    }
}

fn parse_duration_compact(s: &str) -> Option<f64> {
    s.trim().parse::<f64>().ok()
}
