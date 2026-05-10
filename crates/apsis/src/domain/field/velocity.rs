//! Velocity-magnitude field.

use super::body_field::{BodyField, FieldContext};

pub(crate) struct VelocityMagnitudeField;

impl BodyField for VelocityMagnitudeField {
    fn id(&self) -> &'static str {
        "velocity"
    }
    fn name(&self) -> &'static str {
        "Velocity"
    }
    fn unit_label(&self) -> &'static str {
        "|v|"
    }
    fn sample(&self, i: usize, ctx: &FieldContext) -> f64 {
        // `(x² + y² + z²)` in fixed component order. For planar input
        // `vz = 0` the trailing term is an IEEE-754-exact zero
        // addition, so the colour-bar normalisation observed in 2D
        // scenes is preserved bit-for-bit by the 3D sampler.
        let b = &ctx.bodies[i];
        (b.vel_x * b.vel_x + b.vel_y * b.vel_y + b.vel_z * b.vel_z).sqrt()
    }
}
