//! Theme bridge — derives [`egui::Visuals`] and [`egui::Style`] from design
//! tokens.
//!
//! [`install`] swaps the global visuals and style on an [`egui::Context`].
//! Call once at app entry after egui has been initialised. Loading Plex Sans
//! and Plex Mono into [`egui::FontDefinitions`] is a separate concern handled
//! by the app's font setup.

use eframe::egui::{self, CornerRadius, FontFamily, FontId, Stroke, TextStyle, Visuals};

use super::tokens::{border, color, space, typography};

/// Install the design system globally on the given context.
pub fn install(ctx: &egui::Context) {
    ctx.set_visuals(visuals());
    ctx.set_style(style());
}

fn visuals() -> Visuals {
    let mut v = Visuals::dark();

    v.window_fill = color::bg::SURFACE;
    v.panel_fill = color::bg::SURFACE;
    v.faint_bg_color = color::bg::HOVER;
    v.extreme_bg_color = color::bg::CANVAS;

    v.widgets.noninteractive.bg_fill = color::bg::SURFACE;
    v.widgets.noninteractive.weak_bg_fill = color::bg::SURFACE;
    v.widgets.noninteractive.fg_stroke = Stroke::new(border::width::HAIRLINE, color::fg::TERTIARY);
    v.widgets.noninteractive.bg_stroke =
        Stroke::new(border::width::HAIRLINE, color::accent::HAIRLINE);
    v.widgets.noninteractive.corner_radius = corner_radius(border::radius::NONE);

    v.widgets.inactive.bg_fill = color::bg::RAISED;
    v.widgets.inactive.weak_bg_fill = color::bg::RAISED;
    v.widgets.inactive.fg_stroke = Stroke::new(border::width::HAIRLINE, color::fg::SECONDARY);
    v.widgets.inactive.bg_stroke = Stroke::new(border::width::HAIRLINE, color::bg::POPOVER);
    v.widgets.inactive.corner_radius = corner_radius(border::radius::SHARP);

    v.widgets.hovered.bg_fill = color::bg::HOVER;
    v.widgets.hovered.weak_bg_fill = color::bg::HOVER;
    v.widgets.hovered.fg_stroke = Stroke::new(border::width::HAIRLINE, color::fg::PRIMARY);
    v.widgets.hovered.bg_stroke = Stroke::new(border::width::HAIRLINE, color::accent::SOLID);
    v.widgets.hovered.corner_radius = corner_radius(border::radius::SHARP);

    v.widgets.active.bg_fill = color::accent::HAIRLINE;
    v.widgets.active.weak_bg_fill = color::accent::HAIRLINE;
    v.widgets.active.fg_stroke = Stroke::new(border::width::REGULAR, color::accent::SOLID);
    v.widgets.active.bg_stroke = Stroke::new(border::width::REGULAR, color::accent::SOLID);
    v.widgets.active.corner_radius = corner_radius(border::radius::SHARP);

    v.widgets.open.bg_fill = color::bg::POPOVER;
    v.widgets.open.weak_bg_fill = color::bg::POPOVER;
    v.widgets.open.fg_stroke = Stroke::new(border::width::REGULAR, color::fg::PRIMARY);
    v.widgets.open.bg_stroke = Stroke::new(border::width::HAIRLINE, color::accent::HAIRLINE);
    v.widgets.open.corner_radius = corner_radius(border::radius::SHARP);

    v.selection.bg_fill = color::accent::HAIRLINE;
    v.selection.stroke = Stroke::new(border::width::REGULAR, color::accent::SOLID);

    v.hyperlink_color = color::accent::SOLID;

    v.window_shadow = egui::epaint::Shadow::NONE;
    v.popup_shadow = egui::epaint::Shadow::NONE;

    v.window_corner_radius = corner_radius(border::radius::SHARP);
    v.menu_corner_radius = corner_radius(border::radius::SHARP);

    v
}

fn style() -> egui::Style {
    let mut s = egui::Style::default();

    s.spacing.item_spacing = egui::vec2(space::S2, space::S2);
    s.spacing.window_margin = margin(space::S4);
    s.spacing.button_padding = egui::vec2(space::S3, space::S2);
    s.spacing.menu_margin = margin(space::S2);
    s.spacing.indent = space::S4;
    s.spacing.interact_size.y = space::S6;
    s.spacing.scroll = scroll_style();

    s.text_styles = text_styles();

    s
}

fn text_styles() -> std::collections::BTreeMap<TextStyle, FontId> {
    use FontFamily::{Monospace, Proportional};
    let mut map = std::collections::BTreeMap::new();
    map.insert(TextStyle::Small, FontId::new(typography::text::XS, Proportional));
    map.insert(TextStyle::Body, FontId::new(typography::text::BASE, Proportional));
    map.insert(TextStyle::Button, FontId::new(typography::text::BASE, Proportional));
    map.insert(TextStyle::Heading, FontId::new(typography::text::LG, Proportional));
    map.insert(TextStyle::Monospace, FontId::new(typography::text::BASE, Monospace));
    map
}

fn scroll_style() -> egui::style::ScrollStyle {
    egui::style::ScrollStyle {
        bar_width: 6.0,
        handle_min_length: 32.0,
        bar_inner_margin: 4.0,
        bar_outer_margin: 0.0,
        floating_width: 10.0,
        floating_allocated_width: 12.0,
        ..Default::default()
    }
}

fn corner_radius(r: f32) -> CornerRadius {
    CornerRadius::same(r as u8)
}

fn margin(m: f32) -> egui::Margin {
    egui::Margin::same(m as i8)
}
