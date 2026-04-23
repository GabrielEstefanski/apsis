//! Simulation runtime — step loop orchestration, diagnostics, and adaptive control.
//!
//! This module contains *runtime* state and orchestration only. Domain entities
//! live in [`crate::domain`]; persistence and export formats live in
//! [`crate::io`]; rendering infrastructure lives in [`crate::render`].
//!
//! | Module             | Responsibility                                             |
//! |--------------------|------------------------------------------------------------|
//! | [`system`]         | [`System`]: central state, `step()` dispatch, body CRUD    |
//! | [`physics_thread`] | Background thread, `PhysicsHandle`, `RenderState` publish  |
//! | [`metrics`]        | [`Metrics`] data-transfer object                           |
//! | [`calibration`]    | COM recentering, softening/radius helpers                  |
//! | [`adaptive`]       | [`DtController`], [`ThetaController`]                      |
//! | [`diagnostics`]    | Per-step acceleration and jerk statistics                  |
//! | [`log`]            | Structured event bus + `warn_diag!`/`info_diag!` macros    |
//! | [`precision_run`]  | `PrecisionRunController` state machine for IAS15-class runs|

pub mod adaptive;
pub mod calibration;
pub mod diagnostics;
pub mod hooks;
pub mod log;
pub mod metrics;
pub mod physics_thread;
pub mod precision_run;
pub mod system;
pub mod trail;
