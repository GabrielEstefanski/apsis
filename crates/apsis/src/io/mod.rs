//! Persistence, export, and batch-run formats.
//!
//! | Module       | Responsibility                                    |
//! |--------------|---------------------------------------------------|
//! | [`recorder`]   | Scientific CSV export of trajectories and metrics |
//! | [`run_config`] | Declarative `run.toml` configuration for headless runs |
//! | [`headless`]   | Headless batch runner (no GPU/window required) |

pub mod headless;
pub mod recorder;
pub mod run_config;
