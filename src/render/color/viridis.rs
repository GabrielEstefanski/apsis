//! Viridis — matplotlib's default perceptually-uniform colormap.
//! Reference: Stéfan van der Walt & Nathaniel Smith, 2015.

use super::colormap::{Colormap, sample_stops};

pub struct Viridis;

// 9-stop subsample of the 256-entry matplotlib table. Dense enough that the
// piecewise-linear error stays below a single 8-bit step.
const STOPS: &[(f32, [u8; 3])] = &[
    (0.000, [68, 1, 84]),
    (0.125, [72, 34, 115]),
    (0.250, [64, 67, 135]),
    (0.375, [52, 94, 141]),
    (0.500, [41, 121, 142]),
    (0.625, [32, 144, 140]),
    (0.750, [34, 167, 132]),
    (0.875, [68, 190, 112]),
    (1.000, [253, 231, 36]),
];

impl Colormap for Viridis {
    fn id(&self) -> &'static str {
        "viridis"
    }
    fn name(&self) -> &'static str {
        "Viridis"
    }
    fn sample(&self, t: f32) -> [u8; 3] {
        sample_stops(STOPS, t)
    }
}
