//! Body construction presets — convenience layer, not core physics.
//!
//! A [`BodyPreset`] bundles the defaults a researcher reaches for when
//! spawning a body of a familiar physical class: density (or a
//! mass-dependent EOS), a representative colour, an optional
//! luminosity model, and the radiation-pressure receiver coefficient.
//! Presets are consumed once at construction time by
//! [`crate::domain::body::Body::from_preset`] and never referenced
//! again — [`Body`] owns its physical state directly, with no
//! ontological back-reference to the preset that produced it.
//!
//! This file is the **batteries-included** counterpart to the
//! low-level [`Body::new`](crate::domain::body::Body::new)
//! constructor; the two together form the project's two-level body
//! API:
//!
//! * Low-level: explicit physics, what the paper documents.
//!   ```ignore
//!   let earth = Body::new(3.003e-6, 5514.0)
//!       .at_3d(1.0, 0.0, 0.0)
//!       .with_velocity_3d(0.0, 6.283, 0.0);
//!   ```
//! * High-level: ergonomic preset for exploration, templates, tests.
//!   ```ignore
//!   let earth = Body::rocky(3.003e-6).at(1.0, 0.0);
//!   ```
//!
//! Removing the runtime [`Material`] enum (the older shape) means
//! presets cannot be inspected on a live body; if downstream code
//! needs to know which preset produced a body, the body's `name` is
//! the right channel — names are stable user data, presets are
//! construction-time defaults.
//!
//! # Density model
//!
//! Built-in presets that derive density from mass use the
//! Leinhardt & Stewart (2012) / Benz & Asphaug (1999)
//! parameterisation:
//!
//! ```text
//! ρ(m) = ρ₀ · (m / m₀)^α    clamped to [ρ_min, ρ_max]
//! ```
//!
//! [`DensitySource::Fixed`] skips the model and pins ρ to a literal
//! value — used for compact stellar remnants whose density is set by
//! degeneracy pressure rather than EOS extrapolation, and for
//! user-defined custom presets where the author knows the right ρ
//! directly.
//!
//! # Luminosity model
//!
//! Presets carry `luminosity: Option<LuminositySource>`:
//!
//! * `None` — non-luminous (planets, dust, cold remnants).
//! * `Some(Fixed(L))` — user-supplied bolometric luminosity (a
//!   pulsar with a measured X-ray luminosity, for example).
//! * `Some(Model(_))` — computed from mass and radius at
//!   construction time. **Static after construction**; if the user
//!   mutates the body's mass or density, luminosity does not
//!   recompute automatically. Most N-body simulations don't model
//!   stellar evolution; users who need that recompute manually.
//!
//! Luminosity models assume `mass` is in solar masses and `radius`
//! is in solar radii — match the canonical `solar_au` unit system.
//! For other unit systems, the caller is expected to override
//! `body.luminosity` after construction.
//!
//! # References
//!
//! * Leinhardt & Stewart (2012). *Collisions between gravity-dominated
//!   bodies I.* ApJ 745, 79.
//! * Benz & Asphaug (1999). *Catastrophic disruptions revisited.*
//!   Icarus 142, 5–20.
//! * Salaris & Cassisi (2005). *Evolution of Stars and Stellar
//!   Populations.* §5.3.
//! * Burrows et al. (1997). *A nongray theory of extrasolar giant
//!   planets and brown dwarfs.* ApJ 491, 856.
//! * Tout et al. (1996). *Zero-age main-sequence radii and
//!   luminosities as analytic functions of mass and metallicity.*
//!   MNRAS 281, 257.

use crate::domain::body::Body;

// ── Density ───────────────────────────────────────────────────────────────────

/// How a [`BodyPreset`] supplies bulk density to constructed bodies.
///
/// The choice is part of the preset definition; constructors never
/// fall back silently. A preset with [`DensitySource::Fixed`] always
/// produces a body with that exact density; a preset with
/// [`DensitySource::Model`] computes ρ from mass via the embedded
/// power law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DensitySource {
    /// Explicit density, independent of mass. Use when a single
    /// value is the right answer (compact remnants, dimensionless
    /// presets, custom presets where the author already knows ρ).
    Fixed(f64),
    /// Mass-dependent EOS via the Leinhardt-Stewart power law.
    Model(DensityModel),
}

impl DensitySource {
    /// Evaluate the source at the requested mass [simulation units].
    /// `Fixed` ignores `mass`; `Model` evaluates and clamps.
    #[inline]
    pub fn density_at(&self, mass: f64) -> f64 {
        match self {
            Self::Fixed(rho) => *rho,
            Self::Model(model) => model.density_at(mass),
        }
    }
}

/// Power-law equation of state used by [`DensitySource::Model`].
///
/// ```text
/// ρ(m) = ρ₀ · (m / m₀)^α    clamped to [ρ_min, ρ_max]
/// ```
///
/// All values in simulation units consistent with the project
/// convention (kg/m³ for ρ when masses are in solar units; treat as
/// dimensionless otherwise).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DensityModel {
    pub rho_0: f64,
    pub anchor_mass: f64,
    pub alpha: f64,
    pub rho_min: f64,
    pub rho_max: f64,
}

impl DensityModel {
    /// Evaluate ρ(m) with the published clamp. Negative or zero mass
    /// is treated as a tiny positive value so the power law stays
    /// finite (caller-side `mass.abs().max(1e-30)`).
    #[inline]
    pub fn density_at(&self, mass: f64) -> f64 {
        let m = mass.abs().max(1e-30);
        let rho = self.rho_0 * (m / self.anchor_mass).powf(self.alpha);
        rho.clamp(self.rho_min, self.rho_max)
    }
}

// ── Luminosity ────────────────────────────────────────────────────────────────

/// Bolometric luminosity supplier. `None` on the preset itself means
/// non-luminous; this enum is wrapped in `Option` precisely so the
/// caller can distinguish "no luminosity" from "fixed value zero".
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LuminositySource {
    /// Caller-supplied luminosity; preset just carries it through.
    Fixed(f64),
    /// Computed from mass (and sometimes radius) at construction.
    Model(LuminosityModel),
}

impl LuminositySource {
    /// Evaluate the source given mass [solar masses] and radius
    /// [solar radii]. `Fixed` ignores the inputs.
    #[inline]
    pub fn compute(&self, mass_solar: f64, radius_solar: f64) -> f64 {
        match self {
            Self::Fixed(l) => *l,
            Self::Model(model) => model.compute(mass_solar, radius_solar),
        }
    }
}

/// Selectable analytic luminosity laws. All return luminosity in
/// solar luminosities (`L_☉`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LuminosityModel {
    /// Main-sequence mass–luminosity relation with smooth knee at
    /// `M ≈ 0.43 M_☉`; α blends 2.3 → 3.5 across the
    /// fully-convective / radiative-core boundary.
    /// Reference: Salaris & Cassisi (2005) §5.3; Tout et al. (1996).
    MainSequence,
    /// Deuterium-burning luminosity for sub-stellar objects above
    /// 13 M_Jup. Smooth onset; zero below the threshold.
    /// Reference: Burrows et al. (1997).
    BrownDwarfBurrows,
    /// Stefan-Boltzmann cooling at `T_eff = 10 000 K` from the body's
    /// radius — appropriate for white-dwarf-class compact remnants.
    /// Reference: Koester & Chanmugam (1990).
    WhiteDwarfRadius,
}

impl LuminosityModel {
    /// Evaluate the model. Mass in solar masses, radius in solar
    /// radii. All return luminosity in solar luminosities.
    pub fn compute(&self, mass_solar: f64, radius_solar: f64) -> f64 {
        match self {
            Self::MainSequence => main_sequence_luminosity(mass_solar),
            Self::BrownDwarfBurrows => brown_dwarf_luminosity(mass_solar),
            Self::WhiteDwarfRadius => white_dwarf_luminosity(radius_solar),
        }
    }
}

// ── BodyClass ─────────────────────────────────────────────────────────────────

/// Visual / UX taxonomy for a [`Body`].
///
/// Distinct from [`BodyPreset`] — class is a coarse label the renderer
/// and inspector use to group bodies (filter trails by category, label
/// inspector sections), not a physics input. A single preset can be
/// instantiated under different classes: ICY anchors Europa
/// ([`Moon`](Self::Moon)) and Pluto ([`Asteroid`](Self::Asteroid))
/// equally well, and Earth's Moon is a ROCKY body classed as
/// [`Moon`](Self::Moon).
///
/// The class a [`BodyPreset`] suggests via its
/// [`default_class`](BodyPreset::default_class) field is just a
/// starting point; templates and user code can override it through
/// [`Body::with_class`](crate::domain::body::Body::with_class).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BodyClass {
    /// Stellar object: main-sequence star, brown dwarf, white dwarf,
    /// neutron star, black hole. Default for STAR / BROWN_DWARF /
    /// WHITE_DWARF presets.
    Star,
    /// Major planetary body (terrestrial, gas, ice giant). Default for
    /// ROCKY / GAS / ICE_GIANT presets.
    Planet,
    /// Natural satellite (icy or rocky). Templates set this on bodies
    /// orbiting a planet regardless of the underlying preset.
    Moon,
    /// Minor planetary body in a heliocentric (or stellar-centric)
    /// orbit. Default for ASTEROID; also covers KBOs and TNOs when
    /// instantiated from ICY.
    Asteroid,
    /// Volatile-rich body on a typically eccentric orbit. Default for
    /// COMET.
    Comet,
    /// Catch-all for hand-built bodies, test particles, and bodies
    /// loaded from snapshots that predate class persistence.
    Unknown,
}

impl BodyClass {
    /// Stable label suitable for UI controls and snapshot debugging.
    pub const fn display_name(self) -> &'static str {
        match self {
            Self::Star => "Star",
            Self::Planet => "Planet",
            Self::Moon => "Moon",
            Self::Asteroid => "Asteroid",
            Self::Comet => "Comet",
            Self::Unknown => "Unknown",
        }
    }

    /// One-byte codec for snapshot persistence. Stable across schema
    /// versions; new variants append.
    pub const fn to_u8(self) -> u8 {
        match self {
            Self::Star => 0,
            Self::Planet => 1,
            Self::Moon => 2,
            Self::Asteroid => 3,
            Self::Comet => 4,
            Self::Unknown => 5,
        }
    }

    /// Inverse of [`to_u8`](Self::to_u8). Unknown bytes round-trip to
    /// [`Self::Unknown`] so a forward-incompatible save never injects
    /// a wrong category.
    pub const fn from_u8(v: u8) -> Self {
        match v {
            0 => Self::Star,
            1 => Self::Planet,
            2 => Self::Moon,
            3 => Self::Asteroid,
            4 => Self::Comet,
            _ => Self::Unknown,
        }
    }
}

// ── BodyPreset ────────────────────────────────────────────────────────────────

/// Construction preset: the bundle of defaults
/// [`Body::from_preset`](crate::domain::body::Body::from_preset)
/// applies when stamping out a body of a familiar physical class.
///
/// Presets are consumed at construction and never persisted in the
/// resulting [`Body`]. Holding `&'static BodyPreset` references in
/// templates and UI dropdowns is cheap; storing by value works
/// equally well for user-defined presets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyPreset {
    /// Short human-readable label (UI, auto-naming, snapshot
    /// reconstruction). Not inspected by physics.
    pub display_name: &'static str,
    /// Default colour `[R, G, B]` — purely visual, applied iff the
    /// caller doesn't override via [`Body`]'s `with_color`-style
    /// builder. Never read by the physics path.
    pub default_color: [u8; 3],
    /// Radiation-pressure receiver coefficient `Q_pr`. Non-zero
    /// values mark the body as a radiation receiver consumed by
    /// [`crate::physics::radiation`]. Zero on stars and large
    /// planets; positive on dust grains, asteroids, comets.
    pub default_q_pr: f64,
    /// How density is supplied. Required choice — no silent default.
    pub density: DensitySource,
    /// Optional luminosity model; `None` for non-luminous classes.
    pub luminosity: Option<LuminositySource>,
    /// Default UX taxonomy applied when a [`Body`] is constructed from
    /// this preset. Templates can override on a per-body basis when
    /// the preset's natural class does not match the role (Earth's
    /// Moon is a [`ROCKY`] preset classed as
    /// [`Moon`](BodyClass::Moon)).
    pub default_class: BodyClass,

    /// Class-typical Bond albedo (`[0, 1]`, dimensionless), used by
    /// the photometry pipeline as a placeholder when the body's
    /// surface reflectivity isn't quoted directly. Real-body templates
    /// override via [`Body::with_albedo`] with the published value.
    /// Stars carry `0.0` — they emit, they don't reflect.
    pub default_albedo: f64,
}

impl BodyPreset {
    /// Default reference mass for this preset, used by spawn UIs that
    /// pre-fill a sensible mass when the user picks the preset.
    /// Pulled from the density model's anchor mass when the preset
    /// uses [`DensitySource::Model`]; otherwise returns `1.0`.
    #[inline]
    pub fn default_mass(&self) -> f64 {
        match self.density {
            DensitySource::Model(ref m) => m.anchor_mass,
            DensitySource::Fixed(_) => 1.0,
        }
    }

    /// Construct a [`Body`] from this preset at the requested mass.
    /// Convenience wrapper around [`Body::from_preset`] for callers
    /// who already hold the preset.
    #[inline]
    pub fn build(&self, mass: f64) -> Body {
        Body::from_preset(self, mass)
    }
}

// ── Built-in presets ──────────────────────────────────────────────────────────
//
// Reference densities and EOS exponents reproduce the values from
// the previous `Material` taxonomy bit-for-bit, so existing tests
// and templates keep their numeric outputs stable across the
// migration.

/// Dirty snowball: high porosity, volatile-dominated (comet nuclei).
/// Calibrated to 67P (Sierks et al. 2015: ~530 kg/m³).
pub const COMET: BodyPreset = BodyPreset {
    display_name: "Comet",
    default_color: [160, 190, 215],
    default_q_pr: 0.9,
    // 67P-class anchor: ≈10¹³ kg = 5e-18 M_☉. The α = 0.01 makes the
    // power law nearly flat across the cometary mass range, so the
    // exact anchor mostly sets the rho_min branch.
    density: DensitySource::Model(DensityModel {
        rho_0: 500.0,
        anchor_mass: 5e-18,
        alpha: 0.01,
        rho_min: 200.0,
        rho_max: 900.0,
    }),
    luminosity: None,
    default_class: BodyClass::Comet,
    // Halley-class nucleus Bond ≈ 0.04; comet nuclei are among the
    // darkest natural surfaces in the solar system.
    default_albedo: 0.04,
};

/// Rubble-pile or monolithic rock (C/S/M-type asteroids).
/// Calibrated to Carry (2012) bulk density compilation.
pub const ASTEROID: BodyPreset = BodyPreset {
    display_name: "Asteroid",
    default_color: [80, 75, 68],
    default_q_pr: 1.0,
    // Ceres-class anchor: 9.4 × 10²⁰ kg ≈ 4.7e-10 M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 2500.0,
        anchor_mass: 4.7e-10,
        alpha: 0.02,
        rho_min: 1200.0,
        rho_max: 5500.0,
    }),
    luminosity: None,
    default_class: BodyClass::Asteroid,
    // C-type / S-type asteroids span Bond 0.03–0.20; the default
    // sits at the population median (Tedesco & Veeder 1992).
    default_albedo: 0.10,
};

/// Silicate / iron body (terrestrial planets). ρ₀ = 5514 kg/m³
/// at 1 M_⊕ (Dziewonski & Anderson 1981).
/// Reproduces Moon, Mars, Earth, Venus densities to ~10 %.
pub const ROCKY: BodyPreset = BodyPreset {
    display_name: "Rocky",
    default_color: [139, 90, 43],
    default_q_pr: 0.0,
    // Earth anchor: 1 M_⊕ = 3.0034 × 10⁻⁶ M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 5514.0,
        anchor_mass: 3.0034e-6,
        alpha: 0.115,
        rho_min: 3000.0,
        rho_max: 13_000.0,
    }),
    luminosity: None,
    default_class: BodyClass::Planet,
    // Terrestrial-planet population median (Mercury 0.07, Mars
    // 0.25, Earth 0.31, Venus 0.77 — the spread is huge, so this
    // is a placeholder; named bodies override.
    default_albedo: 0.30,
};

/// Volatile-rich body (icy moons, ocean worlds, KBOs). Calibrated to
/// Ganymede (Schubert et al. 2004: 1936 kg/m³).
pub const ICY: BodyPreset = BodyPreset {
    display_name: "Icy",
    default_color: [180, 220, 240],
    default_q_pr: 0.7,
    // Ganymede anchor: 0.025 M_⊕ ≈ 7.45 × 10⁻⁸ M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 1936.0,
        anchor_mass: 7.45e-8,
        alpha: 0.05,
        rho_min: 800.0,
        rho_max: 3500.0,
    }),
    luminosity: None,
    default_class: BodyClass::Moon,
    // Icy moon Bond spans 0.10 (Callisto) to 0.99 (Enceladus);
    // 0.50 lands at the population mean for Galilean / Saturnian
    // satellites.
    default_albedo: 0.50,
};

/// Ice-dominated giant with small rocky core (Uranus/Neptune).
/// Calibrated to Neptune (Guillot 2005: ~1638 kg/m³).
pub const ICE_GIANT: BodyPreset = BodyPreset {
    display_name: "Ice Giant",
    default_color: [64, 164, 223],
    default_q_pr: 0.0,
    // Neptune anchor: 17.15 M_⊕ ≈ 5.151 × 10⁻⁵ M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 1638.0,
        anchor_mass: 5.151e-5,
        alpha: 0.12,
        rho_min: 900.0,
        rho_max: 3000.0,
    }),
    luminosity: None,
    default_class: BodyClass::Planet,
    // Uranus 0.30, Neptune 0.29 — population mean.
    default_albedo: 0.29,
};

/// Gas-dominated giant (Jupiter/Saturn/hot Jupiters). Calibrated to
/// Jupiter (1326 kg/m³); H/He envelope with steep compression.
pub const GAS: BodyPreset = BodyPreset {
    display_name: "Gas Giant",
    default_color: [210, 140, 60],
    default_q_pr: 0.0,
    // Jupiter anchor: 317.8 M_⊕ ≈ 9.5435 × 10⁻⁴ M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 1326.0,
        anchor_mass: 9.5435e-4,
        alpha: 0.18,
        rho_min: 200.0,
        rho_max: 8_000.0,
    }),
    luminosity: None,
    default_class: BodyClass::Planet,
    // Jupiter 0.50, Saturn 0.34 — population mean.
    default_albedo: 0.42,
};

/// Sub-stellar object below the hydrogen-burning limit (~13–80 M_Jup).
/// Density anchored at 40 M_Jup (Chabrier et al. 2009);
/// luminosity from deuterium burning.
pub const BROWN_DWARF: BodyPreset = BodyPreset {
    display_name: "Brown Dwarf",
    default_color: [160, 60, 20],
    default_q_pr: 0.0,
    // 40 M_Jup anchor (Chabrier et al. 2009): 40 × 9.5435 × 10⁻⁴
    // ≈ 0.0382 M_☉.
    density: DensitySource::Model(DensityModel {
        rho_0: 50_000.0,
        anchor_mass: 0.0382,
        alpha: 0.22,
        rho_min: 20_000.0,
        rho_max: 2.0e5,
    }),
    luminosity: Some(LuminositySource::Model(LuminosityModel::BrownDwarfBurrows)),
    default_class: BodyClass::Star,
    // Self-luminous body — emits its own flux. Reflected
    // contribution is negligible against the intrinsic output.
    default_albedo: 0.0,
};

/// Main-sequence star, F/G/K spectral type (Sun-like and warmer).
/// Density falls with mass on the upper main sequence (radiation
/// pressure wins). Luminosity from the smoothed Salaris & Cassisi
/// mass-luminosity relation.
///
/// For low-mass M-dwarfs (TRAPPIST-1, Proxima, Barnard's), use
/// [`RED_DWARF`] instead — different density-mass slope and
/// observably different colour.
pub const STAR: BodyPreset = BodyPreset {
    display_name: "Star",
    default_color: [255, 220, 100],
    default_q_pr: 0.0,
    // Solar anchor: 1 M_☉, ρ₀ = 1408 kg/m³ (Sun's bulk density).
    density: DensitySource::Model(DensityModel {
        rho_0: 1408.0,
        anchor_mass: 1.0,
        alpha: -0.35,
        rho_min: 100.0,
        rho_max: 1.0e5,
    }),
    luminosity: Some(LuminositySource::Model(LuminosityModel::MainSequence)),
    default_class: BodyClass::Star,
    default_albedo: 0.0,
};

/// Low-mass main-sequence star, M spectral type (TRAPPIST-1,
/// Proxima Centauri, Barnard's Star, ~75% of all stars by count).
///
/// Distinguished from [`STAR`]:
/// * Reddish colour (effective temperature 2300–3700 K vs 5000+ for
///   G/K stars).
/// * Density rises sharply at the bottom of the main sequence —
///   very-low-mass M dwarfs approach the brown-dwarf regime in
///   compactness (TRAPPIST-1 ≈ 50 000 kg/m³, Proxima ≈ 56 000).
/// * Anchor at 0.1 M_☉, where ρ ≈ 50 000 kg/m³ is well-measured.
///
/// Reference: Chabrier & Baraffe (2000) ARA&A 38; Mann et al. (2019).
pub const RED_DWARF: BodyPreset = BodyPreset {
    display_name: "Red Dwarf",
    default_color: [220, 100, 60],
    default_q_pr: 0.0,
    density: DensitySource::Model(DensityModel {
        rho_0: 50_000.0,
        anchor_mass: 0.1,
        alpha: -0.5,
        rho_min: 1_000.0,
        rho_max: 2.0e5,
    }),
    luminosity: Some(LuminositySource::Model(LuminosityModel::MainSequence)),
    default_class: BodyClass::Star,
    default_albedo: 0.0,
};

/// Degenerate stellar remnant supported by electron degeneracy.
/// Anchored to Sirius B (~3×10⁶ kg/m³); luminosity computed from
/// radius via Stefan-Boltzmann at 10 000 K.
pub const WHITE_DWARF: BodyPreset = BodyPreset {
    display_name: "White Dwarf",
    default_color: [200, 220, 255],
    default_q_pr: 0.0,
    // Sirius B anchor ≈ 1.018 M_☉; treat as 1.0 M_☉ for the EOS pivot.
    density: DensitySource::Model(DensityModel {
        rho_0: 3.0e6,
        anchor_mass: 1.0,
        // Observational fit; Chandrasekhar α = 2 diverges outside
        // [0.4, 1.2] M_☉ so the empirical 1.2 keeps the model in a
        // physically-bounded range across the merger / accretion
        // events common in N-body simulations.
        alpha: 1.20,
        rho_min: 1.0e5,
        rho_max: 1.0e9,
    }),
    luminosity: Some(LuminositySource::Model(LuminosityModel::WhiteDwarfRadius)),
    default_class: BodyClass::Star,
    default_albedo: 0.0,
};

/// Catalogue of built-in presets. Used by the spawn UI to populate
/// dropdowns and by snapshot loaders to reconstruct legacy material
/// references.
pub const ALL: &[&BodyPreset] = &[
    &COMET,
    &ASTEROID,
    &ROCKY,
    &ICY,
    &ICE_GIANT,
    &GAS,
    &BROWN_DWARF,
    &RED_DWARF,
    &STAR,
    &WHITE_DWARF,
];

// ── Luminosity model implementations ─────────────────────────────────────────

/// Logistic sigmoid: smooth step from 0 → 1 centred at `m0` with width `w`.
#[inline]
fn logistic(m: f64, m0: f64, w: f64) -> f64 {
    1.0 / (1.0 + ((m0 - m) / w).exp())
}

/// Main-sequence mass-luminosity with continuously differentiable α.
/// Returns L in solar luminosities. Reference: Salaris & Cassisi (2005);
/// Tout et al. (1996).
fn main_sequence_luminosity(mass_solar: f64) -> f64 {
    if mass_solar <= 0.0 {
        return 0.0;
    }
    let alpha = 2.3 + (3.5 - 2.3) * logistic(mass_solar, 0.43, 0.15);
    mass_solar.powf(alpha)
}

/// Deuterium-burning luminosity for sub-stellar objects.
/// Reference: Burrows et al. (1997).
fn brown_dwarf_luminosity(mass_solar: f64) -> f64 {
    if mass_solar <= 0.013 {
        return 0.0;
    }
    let onset = logistic(mass_solar, 0.013, 0.002);
    1e-3 * (mass_solar / 0.05).powi(2) * onset
}

/// Stefan-Boltzmann cooling luminosity for a 10 000 K white dwarf.
/// Reference: Koester & Chanmugam (1990).
fn white_dwarf_luminosity(radius_solar: f64) -> f64 {
    const T_EFF: f64 = 10_000.0;
    const T_SUN: f64 = 5_778.0;
    let t_ratio = T_EFF / T_SUN;
    radius_solar * radius_solar * t_ratio.powi(4)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, rel_tol: f64) -> bool {
        let denom = a.abs().max(b.abs()).max(1e-30);
        (a - b).abs() / denom < rel_tol
    }

    // ── Density model parity with the previous Material taxonomy ─────────────
    //
    // Mass arguments are in M_☉ — the canonical solar_au unit. Earlier
    // revisions of these tests passed M_⊕ because the anchor masses
    // were specified in M_⊕; the unit-mismatch fix moved every anchor
    // to M_☉ so the assertions had to follow.

    const M_EARTH_IN_SOLAR: f64 = 3.0034e-6;

    #[test]
    fn rocky_earth_density_recovers_observed_value() {
        let rho = ROCKY.density.density_at(M_EARTH_IN_SOLAR);
        assert!(approx_eq(rho, 5514.0, 0.10), "Earth ρ = {rho}");
    }

    #[test]
    fn rocky_moon_density_recovers_observed_value() {
        let rho = ROCKY.density.density_at(0.0123 * M_EARTH_IN_SOLAR);
        assert!(approx_eq(rho, 3346.0, 0.10), "Moon ρ = {rho}");
    }

    #[test]
    fn gas_jupiter_density_recovers_observed_value() {
        let rho = GAS.density.density_at(317.8 * M_EARTH_IN_SOLAR);
        assert!(approx_eq(rho, 1326.0, 0.10), "Jupiter ρ = {rho}");
    }

    #[test]
    fn star_sun_density_recovers_observed_value() {
        let rho = STAR.density.density_at(1.0);
        assert!(approx_eq(rho, 1408.0, 0.10), "Sun ρ = {rho}");
    }

    #[test]
    fn icy_pluto_density_recovers_observed_value() {
        let rho = ICY.density.density_at(0.0022 * M_EARTH_IN_SOLAR);
        assert!(approx_eq(rho, 1854.0, 0.10), "Pluto ρ = {rho}");
    }

    // ── DensitySource: Fixed bypasses the model ──────────────────────────────

    #[test]
    fn fixed_density_ignores_mass() {
        let src = DensitySource::Fixed(1234.5);
        assert_eq!(src.density_at(1e-12), 1234.5);
        assert_eq!(src.density_at(1e30), 1234.5);
    }

    // ── Density clamps respect bounds ────────────────────────────────────────

    #[test]
    fn density_model_respects_clamps() {
        for preset in ALL {
            if let DensitySource::Model(ref m) = preset.density {
                for &mass in &[1e-12_f64, 1e-4, 1.0, 1e4, 1e12] {
                    let rho = m.density_at(mass);
                    assert!(
                        rho >= m.rho_min && rho <= m.rho_max,
                        "{}: ρ({mass:.0e}) = {rho:.0} outside [{:.0}, {:.0}]",
                        preset.display_name,
                        m.rho_min,
                        m.rho_max,
                    );
                }
            }
        }
    }

    // ── Luminosity models match the previous helper functions ────────────────

    #[test]
    fn main_sequence_solar_luminosity_is_unity() {
        let l = LuminosityModel::MainSequence.compute(1.0, 1.0);
        assert!(approx_eq(l, 1.0, 0.01), "L(1 M_☉) = {l}");
    }

    #[test]
    fn main_sequence_increases_with_mass() {
        let l1 = LuminosityModel::MainSequence.compute(1.0, 1.0);
        let l2 = LuminosityModel::MainSequence.compute(2.0, 1.0);
        assert!(l2 > l1);
    }

    #[test]
    fn brown_dwarf_below_threshold_is_zero() {
        let l = LuminosityModel::BrownDwarfBurrows.compute(0.01, 1.0);
        assert_eq!(l, 0.0);
    }

    #[test]
    fn white_dwarf_radius_dependent_only() {
        let l_a = LuminosityModel::WhiteDwarfRadius.compute(0.6, 0.01);
        let l_b = LuminosityModel::WhiteDwarfRadius.compute(1.0, 0.01);
        // Ignores mass — same radius means same luminosity.
        assert_eq!(l_a, l_b);
    }

    // ── LuminositySource: Fixed bypasses the model ───────────────────────────

    #[test]
    fn fixed_luminosity_ignores_mass_and_radius() {
        let src = LuminositySource::Fixed(7.5);
        assert_eq!(src.compute(1.0, 1.0), 7.5);
        assert_eq!(src.compute(100.0, 0.01), 7.5);
    }

    // ── Preset metadata ──────────────────────────────────────────────────────

    #[test]
    fn presets_have_distinct_display_names() {
        let names: Vec<&str> = ALL.iter().map(|p| p.display_name).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), names.len(), "display names must be unique");
    }

    #[test]
    fn default_mass_pulls_from_density_anchor() {
        assert_eq!(ROCKY.default_mass(), M_EARTH_IN_SOLAR);
        assert_eq!(STAR.default_mass(), 1.0);
        assert_eq!(GAS.default_mass(), 9.5435e-4);
    }

    #[test]
    fn radiation_receivers_have_positive_q_pr() {
        for preset in [&COMET, &ASTEROID, &ICY] {
            assert!(preset.default_q_pr > 0.0, "{} should be a receiver", preset.display_name);
        }
    }

    #[test]
    fn massive_classes_have_zero_q_pr() {
        for preset in [&ROCKY, &ICE_GIANT, &GAS, &BROWN_DWARF, &STAR, &WHITE_DWARF] {
            assert_eq!(
                preset.default_q_pr, 0.0,
                "{} should not be a receiver",
                preset.display_name
            );
        }
    }

    // ── Built-in presets are usable in const context ─────────────────────────

    #[test]
    fn presets_are_const() {
        const _: BodyPreset = ROCKY;
        const _: BodyPreset = STAR;
        const _: &[&BodyPreset] = ALL;
    }
}
