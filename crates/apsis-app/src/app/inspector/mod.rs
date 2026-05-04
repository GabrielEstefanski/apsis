//! Inspector — read-only consumer of [`data::InspectorData`], rendered
//! via design primitives.
//!
//! Module layout:
//! * [`unit_strategy`] — distance/time threshold logic
//! * [`format`] — numeric formatters keyed by [`format::QuantityType`]
//! * [`data`] — `InspectorData` container
//! * [`view`] — UI consumer rendering against primitives
//!
//! NaN values pass through to the formatters and surface as `—`; the
//! consumer applies no fallbacks of its own.

pub mod data;
pub mod format;
pub mod unit_strategy;
pub mod view;

pub use data::{
    ActionData, ActionKind, CameraRelativeData, EnergyData, Header, Identity, InspectorData,
    KinematicState, OrbitData, PerturbationData, PerturbationReadout, RelationKind, RelationsData,
};
pub use view::{InspectorState, show};
