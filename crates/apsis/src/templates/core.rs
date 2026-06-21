//! Template body definitions for N-body simulation scenarios.
//!
//! A [`Template`] is a complete initial-condition specification: a named
//! collection of [`TemplateBody`] descriptors that the simulation engine
//! converts into live [`Body`] instances via [`Template::instantiate`].
//! `TemplateBody` is a *descriptor*, not a `Body`: it holds only what the
//! scenario author specifies; derived properties (density, physical
//! radius) come from the material model at instantiation.
//!
//! Velocities are always in the inertial simulation frame, never relative
//! to a parent body — hierarchical builders (planet + moon) add the
//! parent velocity explicitly; see the solar system template.

use crate::domain::body_preset::BodyPreset;

// ── TemplateBody ──────────────────────────────────────────────────────────────

/// Descriptor for one body in a simulation template. Derived
/// quantities are not stored:
///
/// | Derived quantity    | Computed by                            |
/// |---------------------|----------------------------------------|
/// | Bulk density        | `preset.density.density_at(mass)`      |
/// | Physical radius     | `radius_from_density_mass(ρ, m)`       |
/// | Display colour      | `preset.default_color`                 |
/// | Radiation `q_pr`    | `preset.default_q_pr`                  |
/// | Luminosity          | `preset.luminosity.compute(...)`       |
#[derive(Debug, Clone, Copy)]
pub struct TemplateBody {
    /// Display name; `None` falls back to the preset's `display_name`.
    pub name: Option<&'static str>,

    /// Mass [simulation mass units, e.g. M_☉].
    pub mass: f64,

    /// Construction preset (colour, q_pr, luminosity, density model).
    /// Held by `&'static` so the built-in catalogue is zero-cost;
    /// user-defined presets can be `Box::leak`'d into the same shape.
    pub preset: &'static BodyPreset,

    /// Initial position `[x, y, z]` [simulation length units]. `None`
    /// defers placement to the instantiation logic (COM correction,
    /// grid placement). 2D scenarios set `z = 0`.
    pub position: Option<[f64; 3]>,

    /// Initial velocity `[vx, vy, vz]` in the inertial simulation
    /// frame [length / time]. 2D scenarios set `vz = 0`.
    pub velocity: [f64; 3],

    /// Override the preset's [`default_class`](BodyPreset::default_class)
    /// when a body's role differs from its material — e.g. Earth's Moon
    /// is `ROCKY` but renders as [`BodyClass::Moon`].
    pub class_override: Option<crate::domain::body_preset::BodyClass>,

    /// Explicit density [simulation units], overriding the preset's
    /// density model. Real-body templates pass published values so
    /// `physical_radius` tracks the fact sheet instead of the preset's
    /// EOS bounds; heuristic bodies leave `None`.
    pub density: Option<f64>,

    /// Bond-albedo override. Real-body templates pass the published
    /// value; `None` inherits the preset's class-typical placeholder.
    pub albedo: Option<f64>,
}

impl TemplateBody {
    /// Body at rest; position deferred to the scenario builder.
    pub fn at_rest(mass: f64, preset: &'static BodyPreset) -> Self {
        Self {
            name: None,
            mass,
            preset,
            position: None,
            velocity: [0.0, 0.0, 0.0],
            class_override: None,
            density: None,
            albedo: None,
        }
    }

    /// Body with explicit position and velocity.
    pub fn with_state(
        mass: f64,
        preset: &'static BodyPreset,
        position: [f64; 3],
        velocity: [f64; 3],
    ) -> Self {
        Self {
            name: None,
            mass,
            preset,
            position: Some(position),
            velocity,
            class_override: None,
            density: None,
            albedo: None,
        }
    }
}

// ── UnitSystem ────────────────────────────────────────────────────────────────

/// Physical unit system used by a simulation template: the meaning of
/// the three base units (mass, length, time) plus optional SI conversion
/// factors. The factors let downstream tools reconstruct physical values
/// from raw CSV numbers; they are `None` for dimensionless presets
/// (figure-eight, Pythagorean, …).
#[derive(Debug, Clone, Copy)]
pub struct UnitSystem {
    /// Short label for UI and CSV metadata, e.g. `"AU / M_☉ / T_AU"`.
    pub label: &'static str,

    /// Human-readable mass unit name, e.g. `"M_☉"` or `"–"` for dimensionless.
    pub mass_unit: &'static str,

    /// Human-readable length unit name, e.g. `"AU"` or `"–"`.
    pub length_unit: &'static str,

    /// Human-readable time unit name, e.g. `"T_AU"`, `"yr"`, or `"–"`.
    pub time_unit: &'static str,

    /// Kilograms per simulation mass unit.  `None` for dimensionless presets.
    pub mass_to_kg: Option<f64>,

    /// Metres per simulation length unit.  `None` for dimensionless presets.
    pub length_to_m: Option<f64>,

    /// Seconds per simulation time unit.  `None` for dimensionless presets.
    pub time_to_s: Option<f64>,
}

impl UnitSystem {
    /// Solar-system units for physically-calibrated presets: M_☉, AU,
    /// and `T_AU = √(AU³ / (G_SI · M_☉)) ≈ 58.1 days` (G = 1 in these
    /// base units; Earth's period is 2π T_AU).
    pub const fn solar_au() -> Self {
        Self {
            label: "AU / M_☉ / T_AU",
            mass_unit: "M_☉",
            length_unit: "AU",
            time_unit: "T_AU",
            mass_to_kg: Some(1.989e30),
            length_to_m: Some(1.496e11),
            time_to_s: Some(5.022e6),
        }
    }

    /// G = 1, no physical anchor — for mathematically-defined scenarios
    /// (figure-eight, Pythagorean, Lagrange triangle). Comparing with
    /// observations requires choosing an explicit physical mapping.
    pub const fn dimensionless() -> Self {
        Self {
            label: "dimensionless (G = 1)",
            mass_unit: "–",
            length_unit: "–",
            time_unit: "–",
            mass_to_kg: None,
            length_to_m: None,
            time_to_s: None,
        }
    }
}

// ── Template ──────────────────────────────────────────────────────────────────

/// Complete initial-condition specification for one simulation
/// scenario; the engine converts it into live [`Body`] objects via
/// [`Template::instantiate`]. After instantiation the engine zeroes the
/// net linear momentum.
#[derive(Debug, Clone)]
pub struct Template {
    /// Short human-readable name shown in the scenario picker UI.
    pub name: &'static str,

    /// One-line description shown below the name in the UI.
    pub description: &'static str,

    /// Body descriptors.  Order is preserved through instantiation.
    pub bodies: Vec<TemplateBody>,

    /// Initial display scale [pixels per simulation length unit] —
    /// renderer hint only, no effect on physics.
    pub display_scale: f64,

    /// Orbital plane normal in world coordinates, normalised. The
    /// canvas uses this to orient the default camera so the scenario
    /// loads with the orbital plane visible. `None` falls back to
    /// `[0.0, 0.0, 1.0]` — the convention `state_from_elements` writes
    /// heliocentric ecliptic templates into.
    pub orbital_up: Option<[f64; 3]>,

    /// Suggested camera-to-pivot distance for the initial view, in
    /// world units. Set to frame the bodies the author considers
    /// "primary" (e.g. inner planets for solar_system) instead of the
    /// bounding sphere of all bodies, which collapses the interesting
    /// part to a dot when scales span orders of magnitude. `None`
    /// falls back to bounding-sphere fit.
    pub default_view_distance: Option<f64>,

    /// Suggested time-step [simulation time units]; `None` lets the
    /// integrator choose. Set for scenarios needing a specific cadence
    /// (close binaries).
    pub suggested_dt: Option<f64>,

    /// Suggested integrator; `None` = no preference. Consumers enforce
    /// `explicit user choice > template suggestion > app default`.
    pub suggested_integrator: Option<crate::physics::integrator::IntegratorKind>,

    /// Unit system — informational only (CSV export, literature
    /// comparison): [`UnitSystem::solar_au`] for physically-calibrated
    /// scenarios, [`UnitSystem::dimensionless`] for mathematical ones.
    pub units: UnitSystem,
}

impl Template {
    /// Total mass of all bodies in the template [simulation mass units].
    pub fn total_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }
}
