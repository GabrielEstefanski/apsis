//! Scalar fields sampled per body for visualization.
//!
//! Each [`BodyField`] turns a body into a scalar (velocity magnitude, mass,
//! acceleration magnitude, kinetic energy, …). Fields are composed with a
//! [`Normalizer`](crate::render::color::Normalizer) and a
//! [`Colormap`](crate::render::color::Colormap) to form a
//! [`ColorView`](crate::render::color::ColorView).
//!
//! Adding a new field is a three-line change:
//! 1. implement [`BodyField`] in a new file,
//! 2. register it inside [`FieldRegistry::standard`].

pub mod acceleration;
pub mod body_field;
pub mod kinetic_energy;
pub mod mass;
pub mod registry;
pub mod velocity;

pub use body_field::{BodyField, FieldContext};
pub use registry::FieldRegistry;
