//! Simulation runtime — step loop orchestration, diagnostics, and adaptive control.
//!
//! This module contains *runtime* state and orchestration only. Domain entities
//! live in [`crate::domain`]; persistence and export formats live in
//! [`crate::io`]; rendering infrastructure lives in [`crate::render`].
//!
//! | Module             | Responsibility                                             |
//! |--------------------|------------------------------------------------------------|
//! | [`system`]         | [`System`]: central state, `step()` dispatch, body CRUD    |
//! | [`metrics`]        | [`Metrics`] data-transfer object                           |
//! | [`calibration`]    | COM recentering helpers                                    |
//! | [`adaptive`]       | [`DtController`], [`ThetaController`]                      |
//! | [`diagnostics`]    | Per-step acceleration and jerk statistics                  |
//! | [`log`]            | Structured event bus + `warn_diag!`/`info_diag!` macros    |
//! | [`hooks`]          | `SimHook` extension surface for instrumentation/recording  |

pub mod adaptive;
pub mod calibration;
pub mod diagnostics;
pub mod hooks;
pub mod log;
pub mod metrics;
pub mod system;
