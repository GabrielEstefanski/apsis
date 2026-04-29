//! Acceleration-magnitude field. Reads the live acceleration slice the
//! physics thread publishes on the [`RenderState`](crate::core::physics_thread::RenderState).

use super::body_field::{BodyField, FieldContext};

pub(crate) struct AccelerationMagnitudeField;

impl BodyField for AccelerationMagnitudeField {
    fn id(&self) -> &'static str {
        "acceleration"
    }
    fn name(&self) -> &'static str {
        "Acceleration"
    }
    fn unit_label(&self) -> &'static str {
        "|a|"
    }
    fn sample(&self, i: usize, ctx: &FieldContext) -> f64 {
        // Accelerations may not be published yet (first frame after load).
        match ctx.accelerations.get(i) {
            Some(a) => a.length(),
            None => 0.0,
        }
    }
    fn prefers_log(&self) -> bool {
        true
    }
}
