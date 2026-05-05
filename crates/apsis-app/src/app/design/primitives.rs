//! Design primitives — UI building blocks consuming [`super::tokens`].
//! Render layout, typography, spacing, and visual state; do not depend
//! on simulation domain types.

use std::collections::HashMap;
use std::time::{Duration, Instant};

use eframe::egui::{
    self, Align, Color32, FontFamily, FontId, Layout, Response, RichText, Sense, Stroke, TextStyle,
    Ui, Widget,
};

use super::tokens::{border, color, motion, space, typography};

// ── Section ──────────────────────────────────────────────────────────────────

/// Section block — uppercase tracked-out heading followed by content.
///
/// Heading uses Plex Sans Medium 11px, foreground secondary, with extra
/// letter spacing so it reads as structural label rather than text.
pub struct Section<'a> {
    heading: &'a str,
}

impl<'a> Section<'a> {
    pub fn new(heading: &'a str) -> Self {
        Self { heading }
    }

    /// Render the heading and run `content` inside the section body.
    /// The heading is indented by `space::S4` (12 px) to align with the
    /// header title above; field rows below indent another 4 px to read
    /// as the heading's subordinates.
    pub fn show<R>(self, ui: &mut Ui, content: impl FnOnce(&mut Ui) -> R) -> R {
        ui.add_space(space::S5);
        ui.horizontal(|ui| {
            ui.add_space(space::S4);
            let medium = FontFamily::Name(typography::font::SANS_MEDIUM.into());
            ui.label(
                RichText::new(self.heading.to_uppercase())
                    .font(FontId::new(typography::text::XS, medium))
                    .color(color::fg::SECONDARY),
            );
        });
        ui.add_space(space::S2);
        content(ui)
    }
}

// ── Hairline ─────────────────────────────────────────────────────────────────

/// 1 px horizontal divider in `accent::HAIRLINE` (premultiplied α 0.15).
/// The full-width version spans the available rect with `inset` margin
/// on each side.
pub fn hairline(ui: &mut Ui, inset: f32) {
    let height = border::width::HAIRLINE;
    let (rect, _) =
        ui.allocate_exact_size(egui::vec2(ui.available_width(), height), Sense::hover());
    let painter = ui.painter();
    let y = rect.center().y;
    painter.line_segment(
        [egui::pos2(rect.left() + inset, y), egui::pos2(rect.right() - inset, y)],
        Stroke::new(border::width::HAIRLINE, color::accent::HAIRLINE),
    );
}

// ── Subgroup ─────────────────────────────────────────────────────────────────

/// Subgroup block — a labelled group of rows (`Position`, `Velocity`,
/// etc.) with a thin vertical connector painted on the left tying the
/// children to the heading. The connector reinforces hierarchy without
/// adding decoration; it reads as structure, not ornament.
///
/// Sans 14 px label in foreground secondary, indented `space::S5` so it
/// sits flush with row labels. Connector is `accent::HAIRLINE` (1 px,
/// premultiplied α 0.15) drawn to the left of the rows.
pub struct Subgroup<'a> {
    label: &'a str,
}

impl<'a> Subgroup<'a> {
    pub fn new(label: &'a str) -> Self {
        Self { label }
    }

    /// Render the heading and run `content` for the children. After the
    /// children are placed, a vertical connector is painted from just
    /// below the heading down to the bottom of the children.
    pub fn show<R>(self, ui: &mut Ui, content: impl FnOnce(&mut Ui) -> R) -> R {
        ui.add_space(space::S2);
        ui.horizontal(|ui| {
            ui.add_space(space::S5);
            ui.label(
                RichText::new(self.label)
                    .font(FontId::new(typography::text::BASE, FontFamily::Proportional))
                    .color(color::fg::SECONDARY),
            );
        });
        let label_bottom = ui.cursor().min.y;
        let panel_left = ui.min_rect().left();
        let result = content(ui);
        let children_bottom = ui.cursor().min.y;

        let connector_x = panel_left + space::S5 - space::S1;
        let line_top = label_bottom + space::S1;
        let line_bottom = children_bottom - space::S1;
        if line_bottom > line_top {
            ui.painter().line_segment(
                [egui::pos2(connector_x, line_top), egui::pos2(connector_x, line_bottom)],
                Stroke::new(border::width::HAIRLINE, color::accent::HAIRLINE),
            );
        }
        result
    }
}

// ── FieldRow ─────────────────────────────────────────────────────────────────

/// A single labelled value row.
///
/// Layout: label (sans, foreground secondary, indented) on the left;
/// value (mono, tabular-nums, foreground primary) right-aligned; unit
/// (mono, foreground tertiary) trailing.
///
/// Optional [`flash`] state turns the value background a dim accent for
/// `motion::FAST` whenever the formatted string changes — see
/// [`FlashTracker`] for the rate-limited trigger.
pub struct FieldRow<'a> {
    label: &'a str,
    value: &'a str,
    unit: &'a str,
    indent: f32,
    flash: Option<FlashState>,
}

#[derive(Clone, Copy)]
struct FlashState {
    /// Wall-clock instant the most recent flash started; the row paints a
    /// dim accent rectangle behind the value while `now − started <
    /// motion::FAST`.
    started: Instant,
}

impl<'a> FieldRow<'a> {
    pub fn new(label: &'a str, value: &'a str, unit: &'a str) -> Self {
        Self { label, value, unit, indent: space::S5, flash: None }
    }

    /// Indent the label by an extra increment (used for subgroup items
    /// like `x`, `y`, `z` under a `Position` subgroup label).
    pub fn indented(mut self, levels: usize) -> Self {
        self.indent += space::S5 * levels as f32;
        self
    }

    /// Attach the most recent flash instant for this row, if any. The row
    /// renders the flash background when this value is recent.
    pub fn flash(mut self, started: Option<Instant>) -> Self {
        self.flash = started.map(|started| FlashState { started });
        self
    }

    fn flash_alpha(&self, now: Instant) -> Option<f32> {
        let state = self.flash?;
        let elapsed = now.saturating_duration_since(state.started);
        if elapsed >= motion::FAST {
            return None;
        }
        let progress = elapsed.as_secs_f32() / motion::FAST.as_secs_f32();
        Some(1.0 - progress)
    }
}

impl<'a> Widget for FieldRow<'a> {
    fn ui(self, ui: &mut Ui) -> Response {
        let row_height = typography::lh::DENSE;
        let total_w = ui.available_width();
        let (rect, response) =
            ui.allocate_exact_size(egui::vec2(total_w, row_height), Sense::hover());

        // Flash background — dim accent rectangle, fades over motion::FAST.
        if let Some(alpha) = self.flash_alpha(Instant::now()) {
            let painter = ui.painter();
            let base = color::accent::FLASH;
            let scaled = Color32::from_rgba_premultiplied(
                ((base.r() as f32) * alpha) as u8,
                ((base.g() as f32) * alpha) as u8,
                ((base.b() as f32) * alpha) as u8,
                ((base.a() as f32) * alpha) as u8,
            );
            painter.rect_filled(rect, border::radius::NONE, scaled);
        }

        let mut child = ui.new_child(
            egui::UiBuilder::new().max_rect(rect).layout(Layout::left_to_right(Align::Center)),
        );
        child.add_space(self.indent);
        child.label(
            RichText::new(self.label)
                .font(FontId::new(typography::text::BASE, FontFamily::Proportional))
                .color(color::fg::SECONDARY),
        );

        // Right-aligned value + unit using a nested right_to_left layout.
        child.with_layout(Layout::right_to_left(Align::Center), |ui| {
            if !self.unit.is_empty() {
                ui.label(
                    RichText::new(self.unit)
                        .font(FontId::new(typography::text::BASE, FontFamily::Monospace))
                        .color(color::fg::TERTIARY),
                );
                ui.add_space(space::S2);
            }
            ui.label(
                RichText::new(self.value)
                    .font(FontId::new(typography::text::BASE, FontFamily::Monospace))
                    .color(color::fg::PRIMARY),
            );
        });

        response
    }
}

// ── IconButton ───────────────────────────────────────────────────────────────

/// Action row with an optional Phosphor icon, a label, and an optional
/// keyboard shortcut hint trailing on the right.
///
/// Hover paints `color::bg::HOVER` behind the row. The widget returns
/// the underlying [`Response`] so callers can read `clicked()`.
pub struct IconButton<'a> {
    icon: Option<&'a str>,
    label: &'a str,
    shortcut: Option<&'a str>,
    danger: bool,
}

impl<'a> IconButton<'a> {
    pub fn new(label: &'a str) -> Self {
        Self { icon: None, label, shortcut: None, danger: false }
    }

    pub fn with_icon(mut self, icon: &'a str) -> Self {
        self.icon = Some(icon);
        self
    }

    pub fn with_shortcut(mut self, shortcut: &'a str) -> Self {
        self.shortcut = Some(shortcut);
        self
    }

    /// Mark the row as destructive — hover background turns dim red
    /// (`signal::ERROR` at low alpha).
    pub fn danger(mut self) -> Self {
        self.danger = true;
        self
    }
}

impl<'a> Widget for IconButton<'a> {
    fn ui(self, ui: &mut Ui) -> Response {
        let row_height = space::S6;
        let total_w = ui.available_width();
        let (rect, response) =
            ui.allocate_exact_size(egui::vec2(total_w, row_height), Sense::click());

        if response.hovered() {
            let painter = ui.painter();
            let bg = if self.danger {
                Color32::from_rgba_premultiplied(46, 22, 25, 26)
            } else {
                color::bg::HOVER
            };
            painter.rect_filled(rect, border::radius::SHARP, bg);
        }

        let mut child = ui.new_child(
            egui::UiBuilder::new().max_rect(rect).layout(Layout::left_to_right(Align::Center)),
        );
        child.add_space(space::S5);

        if let Some(icon) = self.icon {
            child.label(
                RichText::new(icon)
                    .font(FontId::new(typography::text::BASE, FontFamily::Proportional))
                    .color(color::fg::SECONDARY),
            );
            child.add_space(space::S2);
        }

        child.label(
            RichText::new(self.label)
                .font(FontId::new(typography::text::BASE, FontFamily::Proportional))
                .color(color::fg::PRIMARY),
        );

        if let Some(shortcut) = self.shortcut {
            child.with_layout(Layout::right_to_left(Align::Center), |ui| {
                ui.add_space(space::S5);
                ui.label(
                    RichText::new(shortcut)
                        .font(FontId::new(typography::text::SM, FontFamily::Monospace))
                        .color(color::fg::TERTIARY),
                );
            });
        }

        response
    }
}

// ── Bloomberg-flash tracker ──────────────────────────────────────────────────

/// Per-field state for the Bloomberg-flash visual signal.
///
/// The tracker remembers the most recent **formatted** value for each
/// field key. When a new formatted value differs from the cached one,
/// a flash starts — *unless* a flash for the same key fired within the
/// `min_interval`, in which case the trigger is suppressed and the cache
/// is silently updated. The interval clamp prevents continuously-evolving
/// fields (e.g. simulation time at 1000× rate) from producing a constant
/// flash that visually degenerates into noise.
#[derive(Debug)]
pub struct FlashTracker {
    last_value: HashMap<String, String>,
    last_flash: HashMap<String, Instant>,
    min_interval: Duration,
}

impl FlashTracker {
    /// Default minimum spacing between consecutive flashes for the same
    /// field — `500ms`. Calmer than 300ms while still surfacing every
    /// real value change; chosen after observing the Inspector demo at
    /// 1e6× sim rate, where 300ms produced perceptual jitter across
    /// adjacent rows. Tunable via [`Self::with_min_interval`].
    pub const DEFAULT_MIN_INTERVAL: Duration = Duration::from_millis(500);

    pub fn new() -> Self {
        Self {
            last_value: HashMap::new(),
            last_flash: HashMap::new(),
            min_interval: Self::DEFAULT_MIN_INTERVAL,
        }
    }

    pub fn with_min_interval(mut self, min_interval: Duration) -> Self {
        self.min_interval = min_interval;
        self
    }

    /// Record an observation. Returns the flash start instant to feed
    /// into [`FieldRow::flash`] — `Some(instant)` if a fresh flash should
    /// run now or is still in progress; `None` if the row should not
    /// flash this frame.
    ///
    /// The instant returned reflects the **most recent** flash trigger
    /// for this key, even when the current call is suppressed by the
    /// rate limit. That keeps the visual fade smooth instead of cutting
    /// when a new sub-interval observation arrives.
    pub fn observe(&mut self, key: &str, formatted: &str) -> Option<Instant> {
        let now = Instant::now();
        let changed = match self.last_value.get(key) {
            Some(prev) => prev != formatted,
            None => false,
        };
        self.last_value.insert(key.to_owned(), formatted.to_owned());

        if changed {
            let allow = match self.last_flash.get(key) {
                Some(&t) => now.saturating_duration_since(t) >= self.min_interval,
                None => true,
            };
            if allow {
                self.last_flash.insert(key.to_owned(), now);
            }
        }

        self.last_flash.get(key).copied()
    }
}

impl Default for FlashTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread::sleep;

    #[test]
    fn flash_tracker_first_observation_does_not_fire() {
        let mut t = FlashTracker::new();
        let first = t.observe("M", "+1.234e+10");
        assert!(first.is_none(), "first observation must not flash");
    }

    #[test]
    fn flash_tracker_unchanged_value_does_not_fire() {
        let mut t = FlashTracker::new();
        t.observe("M", "+1.234e+10");
        let again = t.observe("M", "+1.234e+10");
        assert!(again.is_none(), "unchanged value must not flash");
    }

    #[test]
    fn flash_tracker_changed_value_fires() {
        let mut t = FlashTracker::new();
        t.observe("M", "+1.234e+10");
        let now = t.observe("M", "+1.235e+10");
        assert!(now.is_some(), "changed value must flash");
    }

    #[test]
    fn flash_tracker_rate_limits_within_min_interval() {
        let mut t = FlashTracker::new().with_min_interval(Duration::from_millis(100));
        t.observe("M", "+1.000e+00");
        let first = t.observe("M", "+2.000e+00").unwrap();
        // Immediately change again — same flash instant should be returned
        // (the second change is suppressed by rate limit).
        let suppressed = t.observe("M", "+3.000e+00").unwrap();
        assert_eq!(first, suppressed, "rate limit must suppress within-window changes",);
    }

    #[test]
    fn flash_tracker_allows_new_flash_after_min_interval() {
        let mut t = FlashTracker::new().with_min_interval(Duration::from_millis(20));
        t.observe("M", "+1.000e+00");
        let first = t.observe("M", "+2.000e+00").unwrap();
        sleep(Duration::from_millis(40));
        let second = t.observe("M", "+3.000e+00").unwrap();
        assert!(second > first, "after the min interval, a new change must produce a fresh flash",);
    }

    #[test]
    fn flash_tracker_keys_are_independent() {
        let mut t = FlashTracker::new();
        t.observe("M", "+1.000e+00");
        t.observe("ν", "+2.000e+00");
        let m_change = t.observe("M", "+1.500e+00");
        let nu_unchanged = t.observe("ν", "+2.000e+00");
        assert!(m_change.is_some(), "M changed");
        assert!(nu_unchanged.is_none(), "ν unchanged");
    }
}
