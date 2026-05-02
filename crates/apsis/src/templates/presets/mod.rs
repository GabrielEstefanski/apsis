pub mod binary;
pub mod single_bodies;
pub mod systems;
pub mod threebodyproblems;

pub use binary::{binary_star, star_companion};
pub use single_bodies::{brown_dwarf, gas_giant, rocky_planet, star};

pub use systems::{
    alpha_centauri_ab::alpha_centauri_ab, hd_80606_b_system::hd_80606, hot_jupiter::hot_jupiter,
    jupiter_trojan::jupiter_trojans, kepler_36::kepler_36, pluto_charon::pluto_charon,
    solar_system::solar_system, sun_earth_lagrange::sun_earth_lagrange,
    sun_earth_moon::sun_earth_moon, sun_earth_unstable_lagrange::sun_earth_unstable_lagrange,
    trappist_one::trappist_1,
};

pub use threebodyproblems::{
    three_body_chaotic_ejection::three_body_chaotic_ejection,
    three_body_euler_collinear::three_body_euler_collinear,
    three_body_figure_eight::three_body_figure_eight,
    three_body_lagrange_equilateral::three_body_lagrange_equilateral,
    three_body_pythagorean::three_body_pythagorean,
};
