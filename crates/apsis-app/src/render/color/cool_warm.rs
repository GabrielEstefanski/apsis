//! Cool→warm — diverging blue/red colormap, Universe-Sandbox-style hot/cold.
//!
//! Kenneth Moreland's variant of the classic COLORBREWER RdBu, shown to be
//! the most perceptually uniform diverging colormap by Moreland (2009).

use super::colormap::{Colormap, sample_stops};

pub struct CoolWarm;

const STOPS: &[(f32, [u8; 3])] = &[
    (0.000, [59, 76, 192]),
    (0.250, [124, 159, 249]),
    (0.500, [221, 221, 221]),
    (0.750, [244, 154, 123]),
    (1.000, [180, 4, 38]),
];

impl Colormap for CoolWarm {
    fn id(&self) -> &'static str {
        "cool_warm"
    }
    fn name(&self) -> &'static str {
        "Cool-Warm"
    }
    fn sample(&self, t: f32) -> [u8; 3] {
        sample_stops(STOPS, t)
    }
}
