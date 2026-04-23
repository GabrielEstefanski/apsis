//! Notification Center modal — renders the buffer populated by
//! [`crate::app::notifications`].
//!
//! Accessed via the bell icon in the top bar. The modal is
//! session-scoped: closing it does not clear the buffer, and
//! reopening shows the same list in the same order.

use crate::app::icons;
use crate::app::notifications::{NotificationEntry, NotificationFilter};
use crate::app::theme::{ACCENT, BORDER, DANGER, PANEL_BG, TEXT_DIM, TEXT_PRI, TEXT_SEC};
use crate::app::ui::SimulationApp;
use crate::core::log::{Level, Source};
use eframe::egui::{self, Align, Color32, Layout, RichText, Stroke};
use std::time::SystemTime;

impl SimulationApp {
    pub(in crate::app) fn draw_notifications_panel(&mut self, ctx: &egui::Context) {
        if !self.show_notifications_panel {
            return;
        }

        let snapshot: Vec<NotificationEntry> = {
            let store = self.notifications.lock().unwrap();
            store.entries().cloned().collect()
        };
        let (total_info, total_warn, total_error) = {
            let store = self.notifications.lock().unwrap();
            (
                store.count_at_level(Level::Info),
                store.count_at_level(Level::Warn),
                store.count_at_level(Level::Error),
            )
        };

        let filter = &mut self.notifications_filter;

        let mut clear_requested = false;
        let mut close_requested = false;

        egui::Window::new("Notifications")
            .id(egui::Id::new("notifications_panel"))
            .collapsible(false)
            .resizable(true)
            .min_width(480.0)
            .min_height(320.0)
            .default_width(540.0)
            .default_height(420.0)
            .anchor(egui::Align2::RIGHT_TOP, egui::vec2(-12.0, 48.0))
            .frame(
                egui::Frame::NONE
                    .fill(PANEL_BG)
                    .stroke(Stroke::new(0.5, BORDER))
                    .inner_margin(egui::Margin::symmetric(16, 14)),
            )
            .show(ctx, |ui| {
                ui.set_min_width(480.0);

                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new("NOTIFICATIONS")
                            .size(10.0)
                            .color(TEXT_DIM)
                            .strong(),
                    );
                    ui.label(
                        RichText::new(format!("· {} total", snapshot.len()))
                            .size(10.0)
                            .color(TEXT_SEC),
                    );
                    ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("Close").size(10.5).color(TEXT_SEC),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, BORDER))
                                .min_size(egui::vec2(62.0, 22.0)),
                            )
                            .clicked()
                        {
                            close_requested = true;
                        }
                        ui.add_space(6.0);
                        if ui
                            .add(
                                egui::Button::new(
                                    RichText::new("Clear all").size(10.5).color(TEXT_SEC),
                                )
                                .fill(Color32::TRANSPARENT)
                                .stroke(Stroke::new(0.5, BORDER))
                                .min_size(egui::vec2(78.0, 22.0)),
                            )
                            .clicked()
                        {
                            clear_requested = true;
                        }
                    });
                });

                ui.add_space(8.0);

                ui.horizontal(|ui| {
                    filter_chip(ui, filter, NotificationFilter::All, "All", snapshot.len());
                    filter_chip(
                        ui,
                        filter,
                        NotificationFilter::Info,
                        "Info",
                        total_info,
                    );
                    filter_chip(
                        ui,
                        filter,
                        NotificationFilter::Warn,
                        "Warnings",
                        total_warn,
                    );
                    filter_chip(
                        ui,
                        filter,
                        NotificationFilter::Error,
                        "Errors",
                        total_error,
                    );
                });

                ui.add_space(6.0);
                ui.separator();
                ui.add_space(4.0);

                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        let filtered: Vec<&NotificationEntry> = snapshot
                            .iter()
                            .rev()
                            .filter(|e| filter.matches(e.event.level))
                            .collect();

                        if filtered.is_empty() {
                            ui.add_space(40.0);
                            ui.vertical_centered(|ui| {
                                ui.label(
                                    RichText::new("No notifications to show")
                                        .size(11.0)
                                        .color(TEXT_DIM)
                                        .italics(),
                                );
                            });
                            return;
                        }

                        for entry in filtered {
                            draw_entry(ui, entry);
                        }
                    });
            });

        if clear_requested {
            self.notifications.lock().unwrap().clear();
        }
        if close_requested {
            self.show_notifications_panel = false;
        }
    }
}

fn filter_chip(
    ui: &mut egui::Ui,
    state: &mut NotificationFilter,
    this: NotificationFilter,
    label: &str,
    count: usize,
) {
    let selected = *state == this;
    let text_color = if selected { TEXT_PRI } else { TEXT_SEC };
    let stroke = if selected { ACCENT } else { BORDER };
    let btn = egui::Button::new(
        RichText::new(format!("{} · {}", label, count))
            .size(10.0)
            .color(text_color),
    )
    .fill(Color32::TRANSPARENT)
    .stroke(Stroke::new(if selected { 1.0 } else { 0.5 }, stroke))
    .min_size(egui::vec2(96.0, 22.0));
    if ui.add(btn).clicked() {
        *state = this;
    }
}

fn draw_entry(ui: &mut egui::Ui, entry: &NotificationEntry) {
    let level_icon = match entry.event.level {
        Level::Info => icons::LEVEL_INFO,
        Level::Warn => icons::LEVEL_WARN,
        Level::Error => icons::LEVEL_ERROR,
    };
    let level_color = match entry.event.level {
        Level::Info => TEXT_SEC,
        Level::Warn => ACCENT,
        Level::Error => DANGER,
    };

    ui.horizontal(|ui| {
        ui.add_sized(
            [20.0, 18.0],
            egui::Label::new(RichText::new(level_icon).size(13.0).color(level_color)),
        );

        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new(entry.event.message)
                        .size(11.5)
                        .color(TEXT_PRI)
                        .strong(),
                );
                if entry.count > 1 {
                    ui.label(
                        RichText::new(format!("×{}", entry.count))
                            .size(10.0)
                            .monospace()
                            .color(ACCENT),
                    );
                }
                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                    ui.label(
                        RichText::new(format_timestamp(entry.last_at_system()))
                            .size(9.5)
                            .monospace()
                            .color(TEXT_DIM),
                    );
                });
            });

            ui.horizontal(|ui| {
                ui.label(
                    RichText::new(source_label(entry.event.source))
                        .size(9.0)
                        .color(TEXT_DIM)
                        .strong(),
                );
                if !entry.event.fields.is_empty() {
                    let preview = preview_fields(&entry.event.fields);
                    ui.label(
                        RichText::new(preview)
                            .size(10.0)
                            .monospace()
                            .color(TEXT_SEC),
                    );
                }
            });
        });
    });

    ui.add_space(3.0);
    ui.separator();
    ui.add_space(3.0);
}

fn source_label(source: Source) -> &'static str {
    source.label()
}

fn preview_fields(fields: &[(&'static str, String)]) -> String {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(96);
    const MAX: usize = 4;
    for (i, (k, v)) in fields.iter().take(MAX).enumerate() {
        if i > 0 {
            s.push_str("  ·  ");
        }
        let _ = write!(&mut s, "{}={}", k, v);
    }
    if fields.len() > MAX {
        let _ = write!(&mut s, "  · +{} more", fields.len() - MAX);
    }
    s
}

fn format_timestamp(t: SystemTime) -> String {
    match t.duration_since(SystemTime::UNIX_EPOCH) {
        Ok(d) => {
            let secs = d.as_secs();
            let h = (secs / 3600) % 24;
            let m = (secs / 60) % 60;
            let s = secs % 60;
            format!("{:02}:{:02}:{:02}", h, m, s)
        }
        Err(_) => String::from("--:--:--"),
    }
}

trait EntryExt {
    fn last_at_system(&self) -> SystemTime;
}

impl EntryExt for NotificationEntry {
    fn last_at_system(&self) -> SystemTime {
        // `NotificationEntry::last_at` is an `Instant`, which has no
        // wall-clock mapping; the event's `timestamp` (a `SystemTime`
        // set at publication) is the right source for human-readable
        // display. Pull the original event's timestamp for the first
        // appearance and rely on `last_at`'s monotonic ordering for
        // sort order. A small compromise: coalesced runs show the
        // *first* event's clock time rather than the last.
        self.event.timestamp
    }
}
