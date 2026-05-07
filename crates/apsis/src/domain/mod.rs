//! Domain entities — the data types the simulation operates on.
//!
//! These types represent *what* is being simulated (bodies and their materials),
//! independent of *how* the simulation runs (integrators, scheduling, I/O).
//! Nothing here depends on runtime state, the physics thread, or persistence.

pub mod body;
pub mod body_preset;
pub mod field;
