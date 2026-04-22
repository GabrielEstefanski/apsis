//! Shared infrastructure for the IAS15 benchmark harness.
//!
//! Layered from pure-data to side-effectful:
//!
//! * [`scenarios`] — scenario specifications (initial conditions + window).
//!   No runtime dependencies; can be consumed from anywhere.
//! * [`metrics`] — data types describing what a scenario run produces.
//! * [`runner`] — the one place that actually instantiates `System` and
//!   runs a scenario to completion, turning a [`scenarios::ScenarioSpec`]
//!   into a [`metrics::ScenarioMetrics`].
//! * [`baseline`] — persistent numerical-regression gate. Parses
//!   `benches/baselines/ias15.toml`, compares scenario metrics against
//!   it, and can regenerate the file from a batch of recording runs.
//!
//! The entry point (`benches/ias15.rs`) composes these: it runs
//! validation (runner → baseline::check) before handing control to
//! Criterion for timing, or rewrites the baseline file when the
//! `IAS15_BENCH_UPDATE_BASELINE` env var is set.

pub mod baseline;
pub mod metrics;
pub mod runner;
pub mod scenarios;
