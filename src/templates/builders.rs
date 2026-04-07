pub fn circular_orbit(center_mass: f64, radius: f64, phase: f64) -> ([f64; 2], [f64; 2]) {
    let x = radius * phase.cos();
    let y = radius * phase.sin();

    let v = (center_mass / radius).sqrt();

    let vx = -v * phase.sin();
    let vy = v * phase.cos();

    ([x, y], [vx, vy])
}
