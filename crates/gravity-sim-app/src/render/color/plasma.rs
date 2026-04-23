//! Plasma — matplotlib perceptually-uniform purple→pink→yellow colormap.

use super::colormap::{Colormap, sample_stops};

pub struct Plasma;

const STOPS: &[(f32, [u8; 3])] = &[
    (0.000, [13, 8, 135]),
    (0.125, [75, 3, 161]),
    (0.250, [125, 3, 168]),
    (0.375, [168, 34, 150]),
    (0.500, [203, 70, 121]),
    (0.625, [229, 107, 93]),
    (0.750, [248, 148, 65]),
    (0.875, [253, 195, 40]),
    (1.000, [240, 249, 33]),
];

impl Colormap for Plasma {
    fn id(&self) -> &'static str {
        "plasma"
    }
    fn name(&self) -> &'static str {
        "Plasma"
    }
    fn sample(&self, t: f32) -> [u8; 3] {
        sample_stops(STOPS, t)
    }
}
