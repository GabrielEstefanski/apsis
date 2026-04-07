pub mod binary;
pub mod single_bodies;
pub mod simple_system;
pub mod solar_system;

pub use binary::{binary_star, star_companion};
pub use single_bodies::{brown_dwarf, gas_giant, rocky_planet, star};
pub use simple_system::simple_system;
pub use solar_system::{inner_solar_system, star_gas_giant};
