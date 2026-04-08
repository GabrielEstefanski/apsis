use egui::{Color32, Pos2, Stroke};

use crate::render::RenderBackend;

pub struct EguiBackend<'a> {
    pub(crate) painter: &'a egui::Painter,
}

impl RenderBackend for EguiBackend<'_> {
    fn begin(&mut self) {}

    fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]) {
        self.painter.circle_filled(
            Pos2::new(pos[0], pos[1]),
            radius,
            Color32::from_rgb(color[0], color[1], color[2]),
        );
    }

    fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]) {
        self.painter.circle_stroke(
            Pos2::new(pos[0], pos[1]),
            radius,
            Stroke::new(
                width,
                Color32::from_rgba_premultiplied(color[0], color[1], color[2], color[3]),
            ),
        );
    }

    fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]) {
        self.painter.line_segment(
            [Pos2::new(from[0], from[1]), Pos2::new(to[0], to[1])],
            Stroke::new(
                width,
                Color32::from_rgba_premultiplied(color[0], color[1], color[2], color[3]),
            ),
        );
    }

    fn end(&mut self) {}
}
