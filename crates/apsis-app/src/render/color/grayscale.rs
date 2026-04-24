//! Linear grayscale. Useful for publication figures where colour is reserved
//! for a different channel (e.g. body category).

use super::colormap::Colormap;

pub struct Grayscale;

impl Colormap for Grayscale {
    fn id(&self) -> &'static str {
        "grayscale"
    }
    fn name(&self) -> &'static str {
        "Grayscale"
    }
    fn sample(&self, t: f32) -> [u8; 3] {
        let v = (t.clamp(0.0, 1.0) * 255.0).round() as u8;
        [v, v, v]
    }
}
