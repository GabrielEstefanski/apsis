//! Simulation core — integrator orchestration, calibration, and collision physics.
//!
//! # Module layout
//!
//! | Module | Responsibility |
//! |---|---|
//! | [`system`]      | [`System`] struct: step loop, body CRUD, metrics |
//! | [`metrics`]     | [`Metrics`] data-transfer object (no logic) |
//! | [`calibration`] | COM management, softening/radius calibration |
//! | [`collision`]   | Inelastic merge, energy-based collision detection |
//! | [`adaptive`]    | [`DtController`], [`ThetaController`] |
//! | [`diagnostics`] | Per-step acceleration and jerk statistics |

pub mod calibration;
pub mod diagnostics;
pub mod metrics;
pub mod system;
pub mod trail_buffer;
