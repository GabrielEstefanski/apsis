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
//! | [`recorder`]    | [`Recorder`] struct: recording and replay of simulations |
//! | [`snapshot`]     | [`Snapshot`] struct: serialisable system state for recording and replay |

pub mod body;
pub mod calibration;
pub mod diagnostics;
pub mod materials;
pub mod metrics;
pub mod physics_thread;
pub mod recorder;
pub mod snapshot;
pub mod system;
pub mod trail_buffer;
