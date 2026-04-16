//! Template body definitions for N-body simulation scenarios.
//!
//! A [`Template`] is a complete initial-condition specification: a named
//! collection of [`TemplateBody`] descriptors that the simulation engine
//! converts into live [`Body`] instances via [`Template::instantiate`].
//!
//! ## Design
//!
//! `TemplateBody` is intentionally a *descriptor*, not a `Body`. It holds only
//! the quantities that the scenario author specifies explicitly; all derived
//! physical properties (density, physical radius, softening, moment of inertia)
//! are computed by `Body::new` at instantiation time using the material model.
//!
//! This separation means templates are:
//! - **Compact** — no redundant derived fields.
//! - **Material-consistent** — density and radius always agree with the EOS.
//! - **Forward-compatible** — adding fields to `Body` does not break templates.
//!
//! ## Position semantics
//!
//! `position` is `None` when the scenario builder wants the engine to place the
//! body automatically (e.g. centre-of-mass correction, random placement).
//! Most explicit scenarios set `Some(pos)`.
//!
//! ## Velocity semantics
//!
//! Velocities are always in the **inertial simulation frame**, not relative to
//! any parent body.  Scenario builders that construct hierarchical systems
//! (planet + moon) must add the parent velocity explicitly — see the solar
//! system template for the canonical pattern.

use crate::core::materials::Material;

// ── TemplateBody ──────────────────────────────────────────────────────────────

/// Descriptor for one body in a simulation template.
///
/// All fields that can be derived from [`Material`] and mass are omitted here
/// and computed at instantiation time by [`Body::new`].  Only the quantities
/// that cannot be derived — or that the scenario author wants to override —
/// are present.
///
/// ## What is *not* stored here
///
/// | Derived quantity    | Computed by          |
/// |---------------------|----------------------|
/// | Bulk density        | `density(material, mass)` |
/// | Physical radius     | `radius_from_mass_density(...)` |
/// | Softening length    | `default_softening(mass)` |
/// | Moment of inertia   | `default_moment_inertia(...)` |
/// | Display colour      | `material.props().base_color` |
#[derive(Debug, Clone, Copy)]
pub struct TemplateBody {
    /// Optional authored display name preserved at instantiation time.
    pub name: Option<&'static str>,

    /// Mass [simulation mass units, e.g. M_☉].
    pub mass: f64,

    /// Material class — determines density, colour, and collision behaviour.
    pub material: Material,

    /// Initial position [simulation length units].
    ///
    /// `None` defers placement to the instantiation logic (e.g. the engine
    /// applies a centre-of-mass correction or places bodies on a grid).
    pub position: Option<[f64; 2]>,

    /// Initial velocity in the inertial simulation frame [length / time].
    pub velocity: [f64; 2],

    /// Initial spin angular velocity ω_z [rad / time_unit].
    ///
    /// Positive = counter-clockwise.  Most bodies start at rest (`0.0`);
    /// set this for pre-spun bodies or post-collision fragments.
    pub spin: f64,
}

impl TemplateBody {
    /// Construct a body at rest with no spin.
    ///
    /// Convenience constructor for the common case where position and velocity
    /// will be filled in by the scenario builder.
    pub fn at_rest(mass: f64, material: Material) -> Self {
        Self {
            name: None,
            mass,
            material,
            position: None,
            velocity: [0.0, 0.0],
            spin: 0.0,
        }
    }

    /// Construct a body with explicit position and velocity, no spin.
    pub fn with_state(
        mass: f64,
        material: Material,
        position: [f64; 2],
        velocity: [f64; 2],
    ) -> Self {
        Self {
            name: None,
            mass,
            material,
            position: Some(position),
            velocity,
            spin: 0.0,
        }
    }
}

// ── UnitSystem ────────────────────────────────────────────────────────────────

/// Physical unit system used by a simulation template.
///
/// Declares the meaning of the three independent simulation base units —
/// mass, length, and time — together with optional SI conversion factors.
///
/// # Design
///
/// The simulator always runs with `G = 1 × g_factor` in internal units.
/// Different presets use the same numerical equations but assign different
/// physical interpretations to their numbers:
///
/// | Preset class       | Mass unit | Length unit | Time unit           |
/// |--------------------|-----------|-------------|---------------------|
/// | Solar system, etc. | M_☉       | AU          | T_AU ≈ 58.1 days    |
/// | Figure-eight, etc. | –         | –           | –  (dimensionless)  |
///
/// The SI conversion factors allow downstream tools (Python, Julia, REBOUND) to
/// reconstruct physical values from the raw CSV numbers without re-deriving
/// the unit mapping.  They are `None` for purely dimensionless presets.
///
/// # Derived units
///
/// Given `mass_unit`, `length_unit`, `time_unit`:
///
/// | Quantity              | Derived unit                               |
/// |-----------------------|--------------------------------------------|
/// | Velocity              | `length / time`                            |
/// | Energy                | `mass · length² / time²`                  |
/// | Angular momentum      | `mass · length² / time`                   |
/// | Specific ang. mom.    | `length² / time`                           |
/// | Specific energy       | `length² / time²`                         |
/// | Period                | `time`                                     |
/// | Semi-major axis       | `length`                                   |
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
    /// Solar-system units: mass in M_☉, length in AU, time in T_AU.
    ///
    /// The time unit T_AU is defined by setting G = 1 in these base units:
    ///
    /// ```text
    /// T_AU = √(AU³ / (G_SI · M_☉))
    ///      = √(3.348 × 10³³ / 1.327 × 10²⁰)  s
    ///      ≈ 5.022 × 10⁶ s  ≈ 58.1 days
    /// ```
    ///
    /// Equivalently, Earth's orbital period is 2π T_AU ≈ 365.25 days ≈ 1 year.
    ///
    /// This is the natural unit system for all physically-calibrated presets
    /// in this simulator (Solar System, TRAPPIST-1, Alpha Centauri, etc.).
    pub const fn solar_au() -> Self {
        Self {
            label: "AU / M_☉ / T_AU",
            mass_unit: "M_☉",
            length_unit: "AU",
            // T_AU = sqrt(AU³ / (G_SI · M_☉)) ≈ 5.022e6 s ≈ 58.1 days
            time_unit: "T_AU",
            mass_to_kg: Some(1.989e30),
            length_to_m: Some(1.496e11),
            time_to_s: Some(5.022e6),
        }
    }

    /// Dimensionless units: G = 1, no physical anchor.
    ///
    /// Used for mathematically-defined scenarios (figure-eight, Pythagorean
    /// three-body, Lagrange triangle) where the masses and distances are
    /// chosen to satisfy a specific mathematical condition rather than to
    /// match any physical system.  Results cannot be directly compared with
    /// observations without choosing an explicit physical mapping first.
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

    /// Returns `true` if this unit system has a physical SI mapping.
    #[inline]
    pub fn is_physical(&self) -> bool {
        self.mass_to_kg.is_some()
    }
}

// ── Template ──────────────────────────────────────────────────────────────────

/// Complete initial-condition specification for one simulation scenario.
///
/// A template is the *source of truth* for a scenario.  The simulation engine
/// calls [`Template::instantiate`] (or equivalent) to convert descriptors into
/// live [`Body`] objects with fully computed physical properties.
///
/// ## Centre-of-mass correction
///
/// After instantiation the engine should zero the total linear momentum so the
/// system does not drift off-screen.  This is done by subtracting the
/// mass-weighted mean velocity from every body:
///
/// ```text
/// v_cm = Σ(mᵢ · vᵢ) / Σmᵢ
/// vᵢ  ← vᵢ − v_cm
/// ```
///
/// ## Scale
///
/// `display_scale` is a hint to the renderer — it does not affect the physics.
/// It sets the initial pixels-per-AU (or equivalent) so the scenario fills the
/// viewport sensibly.
#[derive(Debug, Clone)]
pub struct Template {
    /// Short human-readable name shown in the scenario picker UI.
    pub name: &'static str,

    /// One-line description shown below the name in the UI.
    pub description: &'static str,

    /// Body descriptors.  Order is preserved through instantiation.
    pub bodies: Vec<TemplateBody>,

    /// Suggested initial display scale [pixels per simulation length unit].
    ///
    /// The renderer uses this as the default zoom level when the scenario is
    /// first loaded.  Users can zoom freely after that.
    pub display_scale: f64,

    /// Suggested simulation time-step [simulation time units].
    ///
    /// `None` lets the integrator choose an adaptive step.  Set this for
    /// scenarios that require a specific cadence (e.g. close binary stars that
    /// need a small fixed step for accuracy).
    pub suggested_dt: Option<f64>,

    /// Physical unit system used by this template.
    ///
    /// Declares the meaning of the simulation's base units (mass, length, time)
    /// and provides optional SI conversion factors.  This field is purely
    /// informational — it does not affect the physics — but it is essential for
    /// interpreting exported CSV data and for comparing results with the
    /// literature.
    ///
    /// Use [`UnitSystem::solar_au`] for physically-calibrated scenarios and
    /// [`UnitSystem::dimensionless`] for mathematically-defined ones.
    pub units: UnitSystem,
}

impl Template {
    /// Total mass of all bodies in the template [simulation mass units].
    pub fn total_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }

    /// Number of bodies.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Mass-weighted centre-of-velocity.
    ///
    /// The instantiation step should subtract this from every body's velocity
    /// to ensure the system has zero net linear momentum in the simulation
    /// frame.  Bodies with `position = None` are included in the momentum sum
    /// using `velocity` as-is.
    pub fn centre_of_momentum_velocity(&self) -> [f64; 2] {
        let total = self.total_mass();
        if total <= 0.0 {
            return [0.0, 0.0];
        }
        let vx = self
            .bodies
            .iter()
            .map(|b| b.mass * b.velocity[0])
            .sum::<f64>()
            / total;
        let vy = self
            .bodies
            .iter()
            .map(|b| b.mass * b.velocity[1])
            .sum::<f64>()
            / total;
        [vx, vy]
    }

    /// Centre of mass position, computed only over bodies with known positions.
    ///
    /// Returns `None` if no body has a known position.
    pub fn centre_of_mass(&self) -> Option<[f64; 2]> {
        let known: Vec<_> = self
            .bodies
            .iter()
            .filter(|b| b.position.is_some())
            .collect();
        if known.is_empty() {
            return None;
        }
        let total: f64 = known.iter().map(|b| b.mass).sum();
        if total <= 0.0 {
            return None;
        }
        let cx = known
            .iter()
            .map(|b| b.mass * b.position.unwrap()[0])
            .sum::<f64>()
            / total;
        let cy = known
            .iter()
            .map(|b| b.mass * b.position.unwrap()[1])
            .sum::<f64>()
            / total;
        Some([cx, cy])
    }
}
