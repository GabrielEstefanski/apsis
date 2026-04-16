//! Persistence and export formats.
//!
//! | Module       | Responsibility                                    |
//! |--------------|---------------------------------------------------|
//! | [`snapshot`] | Binary `.grav` save/load with schema versioning   |
//! | [`recorder`] | Scientific CSV export of trajectories and metrics |

pub mod recorder;
pub mod snapshot;
