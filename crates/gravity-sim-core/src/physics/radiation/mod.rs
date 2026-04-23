//! Radiation pressure and Poynting–Robertson drag.
//!
//! # Module layout
//!
//! | File              | Contents                                         |
//! |-------------------|--------------------------------------------------|
//! | [`source`]        | [`RadiationSource`] — the emitting body          |
//! | [`params`]        | [`RadiationParams`] — per-body receiver params   |
//! | [`force`]         | Pure force kernels — no Body, no System          |
//! | [`perturbation`]  | [`RadiationField`] — plugs into the integrator   |

pub mod force;
pub mod params;
pub mod perturbation;
pub mod source;

pub use params::RadiationParams;
pub use perturbation::RadiationField;
pub use source::RadiationSource;
