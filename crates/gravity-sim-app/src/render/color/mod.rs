//! Data-driven colour pipeline.
//!
//! ```text
//!   BodyField  ──► Normalizer ──► Colormap ──► [u8; 3] per body
//! ```
//!
//! The three stages are independent trait objects stored in
//! dedicated registries. A `ColorView` is *not* itself a trait — it's a
//! small selection struct ([`ColorViewSelection`]) paired with the
//! evaluator function [`compute`]. This keeps the hot path monomorphic
//! over the resolved references and avoids double indirection.
//!
//! The UI stores `Option<ColorViewSelection>`: `None` means "use material
//! colours" (the pre-existing default). Fields that prefer log-scaling
//! (mass, acceleration) surface that preference through
//! [`BodyField::prefers_log`](gravity_sim_core::domain::field::BodyField::prefers_log)
//! so the UI can auto-pick a sensible normalizer when the field changes.

pub mod colormap;
pub mod cool_warm;
pub mod grayscale;
pub mod inferno;
pub mod normalizer;
pub mod plasma;
pub mod registry;
pub mod view;
pub mod viridis;

pub use colormap::Colormap;
pub use normalizer::Normalizer;
pub use registry::{ColormapRegistry, NormalizerRegistry};
pub use view::{ColorViewOutput, ColorViewSelection, compute};
