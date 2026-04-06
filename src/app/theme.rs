use eframe::egui;
use eframe::egui::{Color32, Pos2, RichText, Stroke};

pub const BG: Color32 = Color32::from_rgb(7, 7, 9);
pub const PANEL_BG: Color32 = Color32::from_rgb(11, 11, 14);
pub const BORDER: Color32 = Color32::from_rgb(28, 28, 36);
pub const ACCENT: Color32 = Color32::from_rgb(200, 200, 210);
pub const ACCENT_DIM: Color32 = Color32::from_rgb(36, 36, 48);
pub const TEXT_PRI: Color32 = Color32::from_rgb(210, 210, 215);
pub const TEXT_SEC: Color32 = Color32::from_rgb(85, 85, 100);
pub const TEXT_DIM: Color32 = Color32::from_rgb(42, 42, 55);
pub const DANGER: Color32 = Color32::from_rgb(190, 70, 70);
pub const SUCCESS: Color32 = Color32::from_rgb(75, 170, 110);

pub fn body_radius(mass: f64) -> f32 {
    (mass.powf(0.33) as f32 * 2.8).clamp(2.5, 20.0)
}

pub fn fmt_world(v: f32) -> String {
    if v >= 1_000.0 {
        format!("{:.0}k", v / 1_000.0)
    } else if v >= 1.0 {
        format!("{:.0}", v)
    } else if v >= 0.01 {
        format!("{:.2}", v)
    } else {
        format!("{:.1e}", v)
    }
}

pub fn sci(v: f64) -> String {
    format!("{:+.3e}", v)
}

pub fn fix4(v: f64) -> String {
    format!("{:.4}", v)
}

pub fn nice_grid_world(scale: f32) -> f32 {
    let raw = 70.0_f32 / scale;
    let exp = raw.log10().floor();
    let frac = raw / 10_f32.powf(exp);
    let nice = if frac < 1.5 {
        1.0
    } else if frac < 3.5 {
        2.0
    } else if frac < 7.5 {
        5.0
    } else {
        10.0
    };
    nice * 10_f32.powf(exp)
}

pub fn apply_visuals(ctx: &egui::Context) {
    let mut vis = egui::Visuals::dark();
    vis.window_fill = PANEL_BG;
    vis.panel_fill = PANEL_BG;
    vis.faint_bg_color = Color32::from_rgb(14, 14, 18);
    vis.extreme_bg_color = Color32::from_rgb(7, 7, 10);
    vis.widgets.noninteractive.bg_fill = Color32::from_rgb(16, 16, 20);
    vis.widgets.noninteractive.fg_stroke = Stroke::new(0.5, BORDER);
    vis.widgets.inactive.bg_fill = Color32::from_rgb(18, 18, 23);
    vis.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT_SEC);
    vis.widgets.hovered.bg_fill = Color32::from_rgb(30, 30, 40);
    vis.widgets.hovered.fg_stroke = Stroke::new(1.0, ACCENT);
    vis.widgets.active.bg_fill = ACCENT_DIM;
    vis.widgets.active.fg_stroke = Stroke::new(1.5, ACCENT);
    vis.selection.bg_fill = ACCENT_DIM;
    vis.selection.stroke = Stroke::new(1.0, ACCENT);
    ctx.set_visuals(vis);
}

pub fn section(ui: &mut egui::Ui, label: &str) {
    ui.add_space(12.0);
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(9.5).color(TEXT_DIM).strong());
        ui.add_space(4.0);
        let rect = ui.available_rect_before_wrap();
        ui.painter().line_segment(
            [
                Pos2::new(rect.left(), rect.center().y),
                Pos2::new(rect.right(), rect.center().y),
            ],
            Stroke::new(0.5, TEXT_DIM),
        );
    });
    ui.add_space(5.0);
}

pub fn metric(ui: &mut egui::Ui, label: &str, value: &str, color: Color32) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(label).size(11.0).color(TEXT_SEC));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(RichText::new(value).monospace().size(11.0).color(color));
        });
    });
}

pub fn field(ui: &mut egui::Ui, label: &str, val: &mut String) {
    ui.horizontal(|ui| {
        ui.add_sized(
            [32.0, 18.0],
            egui::Label::new(RichText::new(label).size(10.0).color(TEXT_SEC)),
        );
        ui.add(
            egui::TextEdit::singleline(val)
                .desired_width(f32::INFINITY)
                .font(egui::TextStyle::Monospace),
        );
    });
}

pub fn primary_btn(ui: &mut egui::Ui, label: &str) -> bool {
    ui.add(
        egui::Button::new(RichText::new(label).size(12.0).color(TEXT_PRI))
            .fill(ACCENT_DIM)
            .stroke(Stroke::new(1.0, BORDER))
            .min_size(egui::vec2(ui.available_width(), 26.0)),
    )
    .clicked()
}

pub fn secondary_btn(ui: &mut egui::Ui, label: &str) -> bool {
    ui.add(
        egui::Button::new(RichText::new(label).size(11.0).color(TEXT_SEC))
            .fill(Color32::TRANSPARENT)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(ui.available_width(), 22.0)),
    )
    .clicked()
}

pub fn template_btn(ui: &mut egui::Ui, label: &str) -> bool {
    ui.add(
        egui::Button::new(RichText::new(label).size(10.5).color(TEXT_SEC))
            .fill(Color32::TRANSPARENT)
            .stroke(Stroke::new(0.5, BORDER))
            .min_size(egui::vec2(ui.available_width(), 20.0)),
    )
    .clicked()
}
