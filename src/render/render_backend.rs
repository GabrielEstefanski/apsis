pub trait RenderBackend {
    fn begin(&mut self);

    fn draw_circle(&mut self, pos: [f32; 2], radius: f32, color: [u8; 3]);

    fn draw_circle_stroke(&mut self, pos: [f32; 2], radius: f32, width: f32, color: [u8; 4]);

    fn draw_line_segment(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: [u8; 4]);

    fn end(&mut self);
}
