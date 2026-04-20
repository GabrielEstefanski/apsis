//! Inferno — matplotlib perceptually-uniform black→red→yellow colormap.

use super::colormap::{Colormap, sample_stops};

pub struct Inferno;

const STOPS: &[(f32, [u8; 3])] = &[
    (0.000, [0, 0, 4]),
    (0.125, [31, 12, 72]),
    (0.250, [85, 15, 109]),
    (0.375, [136, 34, 106]),
    (0.500, [186, 54, 85]),
    (0.625, [227, 89, 51]),
    (0.750, [249, 140, 10]),
    (0.875, [249, 201, 50]),
    (1.000, [252, 255, 164]),
];

impl Colormap for Inferno {
    fn id(&self) -> &'static str {
        "inferno"
    }
    fn name(&self) -> &'static str {
        "Inferno"
    }
    fn sample(&self, t: f32) -> [u8; 3] {
        sample_stops(STOPS, t)
    }
}
