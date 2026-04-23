//! [`BodyField`] — scalar field sampled at each body.
//!
//! A *field* answers the question "what scalar value does this body carry?"
//! for a given visualization channel. Velocity magnitude, mass, acceleration
//! magnitude and kinetic energy are the built-ins; new fields plug in by
//! implementing this trait and registering with
//! [`FieldRegistry`](super::registry::FieldRegistry).
//!
//! Design rationale (SPLASH / yt lineage): visualization backends in the
//! scientific N-body community compose an independent *field* with an
//! independent *normalizer* and *colormap*. Keeping the three orthogonal
//! allows log-scaled mass with viridis, linear velocity with cool-warm,
//! etc. — without a combinatorial enum.

use crate::domain::body::Body;

/// Runtime state a [`BodyField`] may consult while sampling.
///
/// All slices are guaranteed to have the same length as `bodies` so fields
/// can use the body index directly. Fields that don't need accelerations
/// can ignore that slice.
pub struct FieldContext<'a> {
    pub bodies: &'a [Body],
    /// Gravitational accelerations in internal units from the last physics
    /// step. Empty when the physics thread hasn't yet published any.
    pub accelerations: &'a [(f64, f64)],
    /// Simulation time (internal units).
    pub t: f64,
    /// Active gravitational-constant multiplier.
    pub g_factor: f64,
}

/// A scalar field sampled per body.
///
/// Implementations must be `Send + Sync` so the registry and the active
/// [`ColorView`](crate::render::color::ColorView) can be stored behind a
/// trait object and shared across threads without cloning.
pub trait BodyField: Send + Sync {
    /// Stable identifier used for registry lookup and snapshot round-trip.
    /// Must be unique, ASCII, no spaces.
    fn id(&self) -> &'static str;

    /// Human-readable label shown in UI dropdowns.
    fn name(&self) -> &'static str;

    /// Unit suffix shown next to colour-bar tick labels (e.g. `"m/s"`).
    /// Empty string if the field is dimensionless.
    fn unit_label(&self) -> &'static str;

    /// Samples the field at body `i`. Must not panic — return `0.0` if the
    /// value is undefined (e.g. acceleration before the first step).
    fn sample(&self, i: usize, ctx: &FieldContext) -> f64;

    /// Whether a log normalizer is the natural default for this field.
    /// Mass, luminosity and acceleration typically return `true`; velocity
    /// returns `false`. The UI uses this hint when auto-picking a normalizer.
    fn prefers_log(&self) -> bool {
        false
    }
}
