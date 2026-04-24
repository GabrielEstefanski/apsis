//! [`TemplateKind`] — type-safe enumeration of the built-in scenario presets.
//!
//! Primary API for script consumers:
//!
//! ```ignore
//! use gravity_sim_core::core::system::System;
//! use gravity_sim_core::templates::TemplateKind;
//!
//! let mut sys = System::from_template(TemplateKind::SolarSystem);
//! sys.integrate_for(100.0);
//! ```
//!
//! Each variant maps 1:1 to a builder function in [`crate::templates::presets`].
//! Adding a preset means: add a variant here, wire its name + builder in the
//! match arms below, and write its `presets::*` fn. The runtime catalog and
//! the string-keyed lookup both derive from this enum — no duplicate tables
//! to keep in sync.

use crate::templates::Template;
use crate::templates::category::TemplateCategory;
use crate::templates::presets::{hierachical::simple_three_body, *};

/// Built-in scenario presets.
///
/// Type-safe primary identifier for [`System::from_template`](crate::core::system::System::from_template).
/// The string-keyed variant (for config files, CLI args, plugin registration)
/// is [`System::from_template_str`](crate::core::system::System::from_template_str) —
/// but prefer this enum in Rust code for exhaustive autocomplete and
/// compile-time typo rejection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TemplateKind {
    // ── Single bodies ────────────────────────────────────────────────────────
    Star,
    BrownDwarf,
    GasGiant,
    RockyPlanet,

    // ── Multi-body systems ───────────────────────────────────────────────────
    BinaryStars,
    StarWithCompanion,
    SolarSystem,
    Trappist1,
    Kepler36,
    AlphaCentauriAb,
    Hd80606,
    SunEarthLagrange,
    SunEarthUnstableLagrange,
    JupiterTrojans,
    Hierarchical,

    // ── Three-body problems ──────────────────────────────────────────────────
    ThreeBodyChaoticEjection,
    ThreeBodyFigureEight,
    ThreeBodyLagrangeTriangle,
}

/// Error returned by [`TemplateKind::from_name`] /
/// [`System::from_template_str`](crate::core::system::System::from_template_str)
/// when the given name does not match any built-in preset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownTemplate {
    pub name: String,
}

impl std::fmt::Display for UnknownTemplate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "unknown template \"{}\"; known templates: {}",
            self.name,
            TemplateKind::all().iter().map(|t| t.name()).collect::<Vec<_>>().join(", ")
        )
    }
}

impl std::error::Error for UnknownTemplate {}

impl TemplateKind {
    /// Every built-in preset, in canonical order. Useful for UI listings,
    /// test harnesses, and "iterate over all presets" scripts.
    pub fn all() -> &'static [TemplateKind] {
        &[
            TemplateKind::Star,
            TemplateKind::BrownDwarf,
            TemplateKind::GasGiant,
            TemplateKind::RockyPlanet,
            TemplateKind::BinaryStars,
            TemplateKind::StarWithCompanion,
            TemplateKind::SolarSystem,
            TemplateKind::Trappist1,
            TemplateKind::Kepler36,
            TemplateKind::AlphaCentauriAb,
            TemplateKind::Hd80606,
            TemplateKind::SunEarthLagrange,
            TemplateKind::SunEarthUnstableLagrange,
            TemplateKind::JupiterTrojans,
            TemplateKind::Hierarchical,
            TemplateKind::ThreeBodyChaoticEjection,
            TemplateKind::ThreeBodyFigureEight,
            TemplateKind::ThreeBodyLagrangeTriangle,
        ]
    }

    /// Stable display name for this preset. Same string used by
    /// [`from_name`](Self::from_name), the UI catalog, and config-file
    /// round-trips — renaming one breaks all three by design.
    pub fn name(self) -> &'static str {
        match self {
            TemplateKind::Star => "Star",
            TemplateKind::BrownDwarf => "Brown Dwarf",
            TemplateKind::GasGiant => "Gas Giant",
            TemplateKind::RockyPlanet => "Rocky Planet",
            TemplateKind::BinaryStars => "Binary Stars",
            TemplateKind::StarWithCompanion => "Star + Comp.",
            TemplateKind::SolarSystem => "Solar System",
            TemplateKind::Trappist1 => "TRAPPIST-1",
            TemplateKind::Kepler36 => "Kepler-36",
            TemplateKind::AlphaCentauriAb => "Alpha Centauri AB",
            TemplateKind::Hd80606 => "HD 80606 System",
            TemplateKind::SunEarthLagrange => "Sun–Earth L4/L5",
            TemplateKind::SunEarthUnstableLagrange => "Sun–Earth L1/L2/L3",
            TemplateKind::JupiterTrojans => "Jupiter Trojans",
            TemplateKind::Hierarchical => "Hierarchical",
            TemplateKind::ThreeBodyChaoticEjection => "3-Body Chaotic Ejection",
            TemplateKind::ThreeBodyFigureEight => "3-Body Figure Eight",
            TemplateKind::ThreeBodyLagrangeTriangle => "3-Body Lagrange Triangle",
        }
    }

    /// UI category this preset belongs to. Drives the template panel grouping
    /// in the interactive app; headless callers can ignore this.
    pub fn category(self) -> TemplateCategory {
        use TemplateCategory::*;
        match self {
            TemplateKind::Star
            | TemplateKind::BrownDwarf
            | TemplateKind::GasGiant
            | TemplateKind::RockyPlanet => Bodies,

            TemplateKind::BinaryStars
            | TemplateKind::StarWithCompanion
            | TemplateKind::SolarSystem
            | TemplateKind::Trappist1
            | TemplateKind::Kepler36
            | TemplateKind::AlphaCentauriAb
            | TemplateKind::Hd80606
            | TemplateKind::SunEarthLagrange
            | TemplateKind::SunEarthUnstableLagrange
            | TemplateKind::JupiterTrojans
            | TemplateKind::Hierarchical => Systems,

            TemplateKind::ThreeBodyChaoticEjection
            | TemplateKind::ThreeBodyFigureEight
            | TemplateKind::ThreeBodyLagrangeTriangle => ThreeBodyProblems,
        }
    }

    /// Build the [`Template`] (bodies + metadata) for this preset, using
    /// `seed` for any randomised placement (clusters, trojans, etc.).
    /// Deterministic presets ignore the seed.
    pub fn build(self, seed: u64) -> Template {
        match self {
            TemplateKind::Star => star(seed),
            TemplateKind::BrownDwarf => brown_dwarf(seed),
            TemplateKind::GasGiant => gas_giant(seed),
            TemplateKind::RockyPlanet => rocky_planet(seed),
            TemplateKind::BinaryStars => binary_star(seed),
            TemplateKind::StarWithCompanion => star_companion(seed),
            TemplateKind::SolarSystem => solar_system(seed),
            TemplateKind::Trappist1 => trappist_1(seed),
            TemplateKind::Kepler36 => kepler_36(seed),
            TemplateKind::AlphaCentauriAb => alpha_centauri_ab(seed),
            TemplateKind::Hd80606 => hd_80606(seed),
            TemplateKind::SunEarthLagrange => sun_earth_lagrange(seed),
            TemplateKind::SunEarthUnstableLagrange => sun_earth_unstable_lagrange(seed),
            TemplateKind::JupiterTrojans => jupiter_trojans(seed),
            TemplateKind::Hierarchical => simple_three_body(seed),
            TemplateKind::ThreeBodyChaoticEjection => three_body_chaotic_ejection(seed),
            TemplateKind::ThreeBodyFigureEight => three_body_figure_eight(seed),
            TemplateKind::ThreeBodyLagrangeTriangle => three_body_lagrange_triangle(seed),
        }
    }

    /// String-keyed lookup. Returns `Err` with the list of known names when
    /// the input does not match. Used by config-file loaders and the
    /// [`System::from_template_str`](crate::core::system::System::from_template_str)
    /// escape hatch.
    pub fn from_name(name: &str) -> Result<Self, UnknownTemplate> {
        Self::all()
            .iter()
            .copied()
            .find(|k| k.name() == name)
            .ok_or_else(|| UnknownTemplate { name: name.to_string() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_variants_have_unique_names() {
        let mut names: Vec<&str> = TemplateKind::all().iter().map(|t| t.name()).collect();
        names.sort();
        let len_before = names.len();
        names.dedup();
        assert_eq!(len_before, names.len(), "duplicate template names");
    }

    #[test]
    fn all_variants_have_a_category() {
        // Exhaustive match in `category()` means the compiler already guarantees
        // this, but calling it once verifies no accidental `unreachable!()`.
        for &t in TemplateKind::all() {
            let _ = t.category();
        }
    }

    #[test]
    fn from_name_roundtrips_every_variant() {
        for &t in TemplateKind::all() {
            assert_eq!(TemplateKind::from_name(t.name()).unwrap(), t);
        }
    }

    #[test]
    fn from_name_rejects_typos() {
        let err = TemplateKind::from_name("solarsystem").unwrap_err();
        assert_eq!(err.name, "solarsystem");
    }

    #[test]
    fn build_produces_nonempty_body_list_for_systems() {
        // Every "system" preset should have bodies; single-body presets have 1.
        for &t in TemplateKind::all() {
            let tpl = t.build(0);
            assert!(!tpl.bodies.is_empty(), "{:?} produced zero bodies", t);
        }
    }

    #[test]
    fn with_seed_rebuilds_randomised_preset() {
        use crate::core::system::System;

        // Jupiter Trojans uses the seed for cluster layout; different seeds
        // must produce different first-body positions.
        let sys0 = System::from_template(TemplateKind::JupiterTrojans);
        let sys1 = System::from_template(TemplateKind::JupiterTrojans).with_seed(42);

        assert_eq!(sys0.bodies().len(), sys1.bodies().len());
        let any_differ = sys0
            .bodies()
            .iter()
            .zip(sys1.bodies().iter())
            .any(|(a, b)| (a.x - b.x).abs() > 1e-12 || (a.y - b.y).abs() > 1e-12);
        assert!(any_differ, "with_seed must rebuild randomised preset");
    }

    #[test]
    fn manual_mutation_clears_template_source() {
        use crate::core::system::System;
        use crate::domain::body::Body;

        // After from_template + manual add, a subsequent with_seed must NOT
        // wipe the manually-added body. This is enforced by the auto-clear
        // of template_source inside add_body / add_named_bodies / etc.
        let mut sys = System::from_template(TemplateKind::JupiterTrojans);
        let n_before = sys.bodies().len();
        sys.add_body(Body::rocky(1.0).at(100.0, 100.0));
        assert_eq!(sys.bodies().len(), n_before + 1);

        let sys = sys.with_seed(99);
        // Bodies must still be the mutated count, not the fresh template count.
        assert_eq!(
            sys.bodies().len(),
            n_before + 1,
            "with_seed on a manually-mutated system must not rebuild"
        );
    }
}
