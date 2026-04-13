pub mod binary;
pub mod single_bodies;
pub mod systems;
pub mod threebodyproblems;

pub use alpha_centauri_ab::alpha_centauri_ab;
pub use binary::{binary_star, star_companion};
pub use hd_80606_b_system::hd_80606;
pub use jupiter_trojan::jupiter_trojans;
pub use kepler_36::kepler_36;
pub use single_bodies::{brown_dwarf, gas_giant, rocky_planet, star};
pub use solar_system::solar_system;
pub use sun_earth_lagrange::sun_earth_lagrange;
pub use sun_earth_unstable_lagrange::sun_earth_unstable_lagrange;
pub use three_body_chaotic_ejection::three_body_chaotic_ejection;
pub use three_body_figure_eight::three_body_figure_eight;
pub use three_body_lagrange_triangle::three_body_lagrange_triangle;
pub use trappist_one::trappist_1;

pub use systems::{
    alpha_centauri_ab, hd_80606_b_system, hierachical, jupiter_trojan, kepler_36, solar_system,
    sun_earth_lagrange, sun_earth_unstable_lagrange, trappist_one,
};

pub use threebodyproblems::{
    three_body_chaotic_ejection, three_body_figure_eight, three_body_lagrange_triangle,
};
