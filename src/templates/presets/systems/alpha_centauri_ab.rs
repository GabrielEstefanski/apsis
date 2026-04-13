use crate::{
    core::materials::Material,
    templates::{Template, TemplateBody},
};

pub fn alpha_centauri_ab() -> Template {
    let mut bodies = Vec::with_capacity(3);

    let m1: f64 = 1.10;
    let m2: f64 = 0.91;
    let m3: f64 = 0.12;

    let total_mass_ab = m1 + m2;

    let a_ab: f64 = 23.4;
    let e_ab: f64 = 0.52;

    let r_peri_ab = a_ab * (1.0 - e_ab);
    let v_ab = (total_mass_ab * (1.0 + e_ab) / r_peri_ab).sqrt();

    let r1 = r_peri_ab * (m2 / total_mass_ab);
    let r2 = r_peri_ab * (m1 / total_mass_ab);

    let v1 = v_ab * (m2 / total_mass_ab);
    let v2 = v_ab * (m1 / total_mass_ab);

    bodies.push(TemplateBody {
        name: Some("Alpha Centauri A"),
        mass: m1,
        material: Material::Star,
        position: Some([-r1, 0.0]),
        velocity: [0.0, -v1],
        spin: 0.0,
    });

    bodies.push(TemplateBody {
        name: Some("Alpha Centauri B"),
        mass: m2,
        material: Material::Star,
        position: Some([r2, 0.0]),
        velocity: [0.0, v2],
        spin: 0.0,
    });

    let a_p: f64 = 13000.0;
    let e_p: f64 = 0.7;

    let x = a_p * (1.0 - e_p);
    let y = 2000.0;

    let r = (x * x + y * y).sqrt();

    let v = (total_mass_ab * (2.0 / r - 1.0 / a_p)).sqrt();

    let vx = -y / r * v;
    let vy = x / r * v;

    bodies.push(TemplateBody {
        name: Some("Proxima Centauri"),
        mass: m3,
        material: Material::Star,
        position: Some([x, y]),
        velocity: [vx, vy],
        spin: 0.0,
    });

    Template {
        name: "Alpha Centauri ABC",
        description: "Triple system",
        bodies,
        display_scale: 1.0,
        suggested_dt: Some(0.002),
    }
}
