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
        // Same component-order convention as `physics::energy::kinetic_energy`
        // so per-body samples sum (over `bodies`) to the system total.
        // For planar input `vz = 0` the trailing `+ vz²` is an
        // IEEE-754-exact zero addition.
        let b = &ctx.bodies[i];
        0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz)
    }
    fn prefers_log(&self) -> bool {
        true
    }
}
