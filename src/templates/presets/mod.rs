pub mod binary;
pub mod single_bodies;
pub mod systems;
pub mod threebodyproblems;

pub use binary::{binary_star, star_companion};
pub use jupiter_trojan::jupiter_trojans;
pub use single_bodies::{brown_dwarf, gas_giant, rocky_planet, star};
pub use solar_system::solar_system;
pub use sun_earth_lagrange::sun_earth_lagrange;
pub use sun_earth_unstable_lagrange::sun_earth_unstable_lagrange;

pub use three_body_chaotic_ejection::three_body_chaotic_ejection;
pub use three_body_figure_eight::three_body_figure_eight;
pub use three_body_lagrange_triangle::three_body_lagrange_triangle;

pub use systems::{
    hierachical, jupiter_trojan, solar_system, sun_earth_lagrange, sun_earth_unstable_lagrange,
};

pub use threebodyproblems::{
    three_body_chaotic_ejection, three_body_figure_eight, three_body_lagrange_triangle,
};
