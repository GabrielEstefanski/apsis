//! Inspector data container.
//!
//! [`InspectorData`] is the structured payload the Inspector view consumes.
//! Builders / wire-up code (PR-B) populate this from the live simulation
//! state. The Inspector itself is a pure consumer — see
//! `feedback_scientific_app_idiom.md` (R2): no recompute, no fallback,
//! NaN renders as `—`.

use eframe::egui::Color32;

/// Top-level Inspector payload. Each `Option` section is `Some` when the
/// corresponding domain data is available; the view renders only the
/// sections that are `Some`. Auto-shown sections (`perturbations`,
/// `camera_relative`) follow the same rule — `Some` means "ready and
/// physically active", `None` means "not applicable now".
#[derive(Debug, Clone)]
pub struct InspectorData {
    pub header: Header,
    pub identity: Identity,
    pub state: KinematicState,
    pub orbit: Option<OrbitData>,
    /// Hierarchical relationship surfaced when the body's hierarchical
    /// primary differs from the strongest attractor used in `orbit`
    /// (Moon → Earth case). Auto-hidden otherwise.
    pub relations: Option<RelationsData>,
    pub energy: Option<EnergyData>,
    pub perturbations: Vec<PerturbationData>,
    pub camera_relative: Option<CameraRelativeData>,
    pub actions: Vec<ActionData>,
}

#[derive(Debug, Clone)]
pub struct Header {
    pub name: String,
    pub breadcrumb: String,
    /// Body colour swatch (`shape::SWATCH` 8×8 px filled rect).
    pub swatch: Color32,
}

#[derive(Debug, Clone)]
pub struct Identity {
    pub mass_kg: f64,
    pub radius_m: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct KinematicState {
    pub position_m: [f64; 3],
    pub velocity_m_s: [f64; 3],
}

/// Orbit element block. The struct mirrors
/// [`apsis::physics::orbital::OrbitalElements`] but flattened into f64
/// scalars and degree angles ready for direct formatting. `q` and `Q`
/// are populated as plain values (NaN when undefined).
#[derive(Debug, Clone)]
pub struct OrbitData {
    pub primary_name: String,
    pub semi_major_axis_m: f64,
    pub eccentricity: f64,
    pub inclination_rad: f64,
    pub period_s: f64,
    pub lon_ascending_node_rad: f64,
    pub argument_of_pericenter_rad: f64,
    pub true_anomaly_rad: f64,
    pub mean_anomaly_rad: f64,
    pub eccentric_anomaly_rad: f64,
    pub pericenter_m: f64,
    pub apocenter_m: f64,
}

#[derive(Debug, Clone)]
pub struct EnergyData {
    pub kinetic_j: f64,
    pub potential_j: f64,
    pub specific_j: f64,
}

/// One row in the PERTURBATIONS section. Each registered perturbation
/// reports its name and a key-value pair describing its current
/// contribution to the integrator (advance rate for 1PN, etc.).
#[derive(Debug, Clone)]
pub struct PerturbationData {
    pub name: String,
    pub active: bool,
    pub readouts: Vec<PerturbationReadout>,
}

#[derive(Debug, Clone)]
pub struct PerturbationReadout {
    pub label: String,
    pub value_str: String,
    pub unit: String,
}

/// Relation type that bound `secondary_name` to `primary_name`. Recorded
/// here so the consumer can show the user *how* the relationship was
/// established (Hill sphere / energy fallback) rather than asserting
/// it as fact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RelationKind {
    /// Hill-sphere containment — the standard hierarchical-capture case.
    BoundHillSphere,
    /// Energy fallback — Hill-sphere check failed but the secondary is
    /// gravitationally bound to the primary (`ε < 0`).
    BoundEnergy,
}

impl RelationKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::BoundHillSphere => "Bound (Hill sphere)",
            Self::BoundEnergy => "Bound (energy)",
        }
    }
}

#[derive(Debug, Clone)]
pub struct RelationsData {
    /// Hierarchical primary's name (e.g. "Earth" for the Moon).
    pub primary_name: String,
    /// The selected body's role in this relationship — the secondary in
    /// the binary `(primary, secondary)`. Same as the inspector header
    /// title; surfaced explicitly so the relationship reads symmetrically.
    pub secondary_name: String,
    pub kind: RelationKind,
    /// Natural reference frame label, e.g. `"Barycentric (Earth–Moon)"`.
    /// Computational quantities elsewhere in the inspector are not
    /// rotated into this frame; the label documents what frame would be
    /// natural for further analysis.
    pub frame_label: String,
}

/// Inspector payload for a multi-body selection.
///
/// Sections that are semantically undefined for groups (Orbit, Energy,
/// Relations, Camera-relative) are intentionally absent. Only
/// mass-aggregate quantities and COM kinematics are reported.
#[derive(Debug, Clone)]
pub struct AggregateData {
    pub count: usize,
    /// Body names in selection order, for the breadcrumb listing under
    /// the header. The view truncates with an ellipsis when the list
    /// would exceed the available width.
    pub body_names: Vec<String>,
    pub total_mass_kg: f64,
    /// Centre-of-mass position in metres.
    pub com_m: [f64; 3],
    /// Centre-of-mass velocity in m/s.
    pub v_com_m_s: [f64; 3],
    /// `max(|r_i − COM|)` — maximum distance from COM to any selected body.
    pub bounding_radius_m: f64,
    pub actions: Vec<ActionData>,
}

#[derive(Debug, Clone)]
pub struct CameraRelativeData {
    pub distance_m: f64,
    /// Signed radial velocity in m/s — negative when approaching.
    pub radial_velocity_m_s: f64,
    pub apparent_size_arcsec: f64,
    pub off_axis_rad: f64,
}

#[derive(Debug, Clone)]
pub struct ActionData {
    pub label: String,
    pub icon: Option<String>,
    pub shortcut: Option<String>,
    pub kind: ActionKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActionKind {
    Neutral,
    Destructive,
}
