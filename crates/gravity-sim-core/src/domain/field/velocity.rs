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
        let b = &ctx.bodies[i];
        (b.vx * b.vx + b.vy * b.vy).sqrt()
    }
}
