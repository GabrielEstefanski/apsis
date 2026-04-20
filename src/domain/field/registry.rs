//! Registry of [`BodyField`] instances.
//!
//! The registry owns the fields (as `Box<dyn BodyField>`) and hands out
//! borrowed references by string ID. The app keeps a single instance and
//! looks up the active field every frame — cheap because resolution is a
//! linear scan over a handful of entries.

use super::body_field::BodyField;

pub struct FieldRegistry {
    entries: Vec<Box<dyn BodyField>>,
}

impl FieldRegistry {
    pub fn new() -> Self {
        Self { entries: Vec::new() }
    }

    /// Construct the registry populated with every built-in field.
    pub fn standard() -> Self {
        use super::acceleration::AccelerationMagnitudeField;
        use super::kinetic_energy::KineticEnergyField;
        use super::mass::MassField;
        use super::velocity::VelocityMagnitudeField;

        let mut r = Self::new();
        r.register(Box::new(VelocityMagnitudeField));
        r.register(Box::new(MassField));
        r.register(Box::new(AccelerationMagnitudeField));
        r.register(Box::new(KineticEnergyField));
        r
    }

    pub fn register(&mut self, field: Box<dyn BodyField>) {
        self.entries.push(field);
    }

    pub fn get(&self, id: &str) -> Option<&dyn BodyField> {
        self.entries
            .iter()
            .find(|f| f.id() == id)
            .map(|f| f.as_ref())
    }

    pub fn iter(&self) -> impl Iterator<Item = &dyn BodyField> {
        self.entries.iter().map(|f| f.as_ref())
    }
}

impl Default for FieldRegistry {
    fn default() -> Self {
        Self::standard()
    }
}
