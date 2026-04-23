//! Kinetic energy per body: ½ m |v|².

use super::body_field::{BodyField, FieldContext};

pub(crate) struct KineticEnergyField;

impl BodyField for KineticEnergyField {
    fn id(&self) -> &'static str {
        "kinetic_energy"
    }
    fn name(&self) -> &'static str {
        "Kinetic energy"
    }
    fn unit_label(&self) -> &'static str {
        "K"
    }
    fn sample(&self, i: usize, ctx: &FieldContext) -> f64 {
        let b = &ctx.bodies[i];
        0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy)
    }
    fn prefers_log(&self) -> bool {
        true
    }
}
