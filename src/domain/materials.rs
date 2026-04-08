use std::f64::consts::PI;

/// Discrete material classes used to parameterise bulk properties such as
/// density, colour, restitution and fragmentation threshold.  This is a
/// **phenomenological model**: each variant represents a broad astrophysical
/// category rather than a specific chemical composition.
///
/// The material determines:
/// - How density scales with mass (gravitational compression model).
/// - The body's display colour (derived from surface/atmospheric type).
/// - The coefficient of restitution used in arcade-mode bounces.
/// - A multiplicative scaling factor on the Leinhardt-Stewart disruption
///   threshold Q*_D, making hard bodies (stars) harder to shatter and
///   volatile bodies (comets) easier.
///
/// Radius is always derived, never set directly:
///   r = (3m / 4πρ)^(1/3)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Material {
    // ── Planetary bodies ─────────────────────────────────────────────────── //
    /// Silicate / metallic body (terrestrial planets, metallic asteroids).
    Rocky,

    /// Volatile-rich body (icy moons, ocean worlds).
    Icy,

    /// Gas-dominated giant (Jupiter / Saturn analogues).
    Gas,

    /// Ice-dominated giant with a small rocky core (Neptune / Uranus analogues).
    IceGiant,

    // ── Small solar-system bodies ────────────────────────────────────────── //
    /// Rubble-pile or monolithic rock (C/S/M-type asteroids).
    Asteroid,

    /// Dirty snowball: high porosity, very volatile (comet nuclei).
    Comet,

    /// Diffuse ejecta cloud that aggregates unresolved impact dust.
    DustCloud,

    // ── Stellar objects ───────────────────────────────────────────────────── //
    /// Main-sequence star (F/G/K spectral types — Sun-like).
    Star,

    /// Sub-stellar object, below hydrogen-burning limit (~13–80 M_Jup).
    BrownDwarf,

    /// Degenerate stellar remnant; Earth-sized but solar-mass density.
    WhiteDwarf,
}

/// Material parameters controlling how density evolves with mass and how the
/// body behaves in collisions and fragmentation events.
///
/// ## Density model
/// ρ(m) = base_density · (1 + compressibility · log10(m))
/// clamped to [density_min, density_max].
///
/// This captures first-order gravitational compression without requiring a
/// full equation of state.
///
/// ## Collision / fragmentation parameters
/// - `restitution`       — coefficient of restitution (CoR) for arcade-mode
///   bounces; 0 = perfectly inelastic, 1 = elastic.
/// - `disruption_scale`  — multiplier on the Leinhardt-Stewart Q*_D threshold.
///   > 1 makes the body harder to shatter; < 1 makes it easier.
///
/// ## Visual
/// - `base_color` — representative [R, G, B] for the body type.
///   The renderer may tint this further (e.g. temperature glow for stars).
#[derive(Debug, Clone, Copy)]
pub struct MaterialProps {
    /// Baseline density for a body of mass ~1 in simulation units.
    pub base_density: f64,

    /// Strength of density increase with mass (dimensionless).
    pub compressibility: f64,

    /// Lower bound for density (prevents unrealistic expansion on merges).
    pub density_min: f64,

    /// Upper bound for density (prevents runaway compression).
    pub density_max: f64,

    /// Representative display colour [R, G, B].
    pub base_color: [u8; 3],
}

impl Material {
    /// Returns the parameter set associated with this material.
    pub fn props(self) -> MaterialProps {
        match self {
            // ── Planetary bodies ─────────────────────────────────────────── //
            Material::Rocky => MaterialProps {
                base_density: 5500.0,
                compressibility: 0.15,
                density_min: 3000.0,
                density_max: 13_000.0,
                base_color: [139, 90, 43],
            },

            Material::Icy => MaterialProps {
                base_density: 1500.0,
                compressibility: 0.08,
                density_min: 800.0,
                density_max: 3000.0,
                base_color: [180, 220, 240],
            },

            Material::Gas => MaterialProps {
                base_density: 700.0,
                compressibility: 0.35,
                density_min: 200.0,
                density_max: 5000.0,
                base_color: [210, 140, 60],
            },

            Material::IceGiant => MaterialProps {
                base_density: 1300.0,
                compressibility: 0.12,
                density_min: 900.0,
                density_max: 2500.0,
                base_color: [64, 164, 223],
            },

            // ── Small solar-system bodies ────────────────────────────────── //
            Material::Asteroid => MaterialProps {
                base_density: 2200.0,
                compressibility: 0.05,
                density_min: 1200.0,
                density_max: 4000.0,
                base_color: [80, 75, 68],
            },

            Material::Comet => MaterialProps {
                base_density: 500.0,
                compressibility: 0.03,
                density_min: 300.0,
                density_max: 900.0,
                base_color: [160, 190, 215],
            },

            Material::DustCloud => MaterialProps {
                base_density: 120.0,
                compressibility: 0.0,
                density_min: 20.0,
                density_max: 400.0,
                base_color: [194, 176, 148],
            },

            // ── Stellar objects ───────────────────────────────────────────── //
            Material::Star => MaterialProps {
                base_density: 1400.0,
                compressibility: 0.60,
                density_min: 500.0,
                density_max: 1.0e5,
                base_color: [255, 220, 100],
            },

            Material::BrownDwarf => MaterialProps {
                base_density: 50_000.0,
                compressibility: 0.50,
                density_min: 20_000.0,
                density_max: 1.0e5,
                base_color: [160, 60, 20],
            },

            Material::WhiteDwarf => MaterialProps {
                base_density: 1.0e7,
                compressibility: 0.80,
                density_min: 1.0e6,
                density_max: 1.0e9,
                base_color: [200, 220, 255],
            },
        }
    }

    /// Short human-readable name for UI display.
    pub fn display_name(self) -> &'static str {
        match self {
            Material::Rocky => "Rocky",
            Material::Icy => "Icy",
            Material::Gas => "Gas Giant",
            Material::IceGiant => "Ice Giant",
            Material::Asteroid => "Asteroid",
            Material::Comet => "Comet",
            Material::DustCloud => "Dust Cloud",
            Material::Star => "Star",
            Material::BrownDwarf => "Brown Dwarf",
            Material::WhiteDwarf => "White Dwarf",
        }
    }

    #[inline]
    pub fn is_diffuse(self) -> bool {
        matches!(self, Material::DustCloud)
    }
}

/// Computes the bulk density for a body of given material and mass.
///
/// Model:
/// ρ(m) = base_density · (1 + k · log10(m))
///
/// where k = compressibility.  The result is clamped to the material bounds.
///
/// Guarantees:
/// - ρ > 0
/// - monotonic (non-decreasing for m ≥ 1)
pub fn density(material: Material, mass: f64) -> f64 {
    let props = material.props();

    // avoid log(0) and negative masses
    let m = mass.abs().max(1e-12);

    let factor = 1.0 + props.compressibility * m.log10().max(0.0);
    let rho = props.base_density * factor;

    rho.clamp(props.density_min, props.density_max)
}

/// Radius from mass and density:
/// r = (3m / 4πρ)^(1/3)
///
/// This is the primary geometric relation used throughout the simulation.
#[inline]
pub fn radius_from_mass_density(mass: f64, density: f64) -> f64 {
    ((3.0 * mass) / (4.0 * PI * density.max(1e-30))).cbrt()
}
