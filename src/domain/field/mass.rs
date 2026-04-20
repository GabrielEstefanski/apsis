//! Mass field — spans many decades across realistic systems, so it prefers a
//! logarithmic normalizer by default.

use super::body_field::{BodyField, FieldContext};

pub struct MassField;

impl BodyField for MassField {
    fn id(&self) -> &'static str {
        "mass"
    }
    fn name(&self) -> &'static str {
        "Mass"
    }
    fn unit_label(&self) -> &'static str {
        "m"
    }
    fn sample(&self, i: usize, ctx: &FieldContext) -> f64 {
        ctx.bodies[i].mass
    }
    fn prefers_log(&self) -> bool {
        true
    }
}
