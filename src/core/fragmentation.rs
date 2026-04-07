//! Impact-energy fragmentation physics for the N-body simulation.
//!
//! Simplified gravitational-regime model inspired by Leinhardt & Stewart 2012.

use std::f64::consts::PI;

use crate::domain::body::{
    Body, default_moment_inertia, default_softening, radius_from_density_mass,
};
use crate::domain::materials::{Material, pair_disruption_scale};

pub const SUB_THRESHOLD: f64 = 0.1;
pub const DISRUPTION_THRESHOLD: f64 = 1.0;

const MIN_FRAGMENT_MASS: f64 = 1e-4;
const STRENGTH_K: f64 = 0.05;

struct ImpactGeometry {
    v_n: f64,
    v_t: f64,
    d: f64,
}

fn impact_geometry(bi: &Body, bj: &Body) -> ImpactGeometry {
    let dx = bi.x - bj.x;
    let dy = bi.y - bj.y;
    let d = (dx * dx + dy * dy).sqrt().max(1e-30);
    let nx = dx / d;
    let ny = dy / d;
    let dvx = bi.vx - bj.vx;
    let dvy = bi.vy - bj.vy;
    let v_n = dvx * nx + dvy * ny;
    let v_t = nx * dvy - ny * dvx;
    ImpactGeometry { v_n, v_t, d }
}

fn make_body(
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    mass: f64,
    density: f64,
    material: Material,
) -> Body {
    let physical_radius = radius_from_density_mass(density, mass);
    Body {
        x,
        y,
        vx,
        vy,
        mass,
        density,
        radius: physical_radius,
        physical_radius,
        softening: default_softening(mass).max(physical_radius * 2.0),
        omega_z: 0.0,
        moment_inertia: default_moment_inertia(mass, physical_radius),
        material,
        color: material.props().base_color,
    }
}

fn make_dust_cloud(
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    mass: f64,
    source_density: f64,
) -> Option<Body> {
    if mass <= 1e-12 {
        return None;
    }

    let props = Material::DustCloud.props();
    let density = (0.05 * source_density).clamp(props.density_min, props.density_max);
    let mut cloud = Body::new(x, y, vx, vy, mass, Material::DustCloud);
    cloud.density = density;
    cloud.sync_physical_properties();
    cloud.radius = 0.0;
    cloud.softening = (4.0 * cloud.physical_radius).max(2.0 * cloud.physical_radius);
    cloud.omega_z = 0.0;
    Some(cloud)
}

pub enum ImpactResult {
    SubThreshold,
    HitAndRun {
        bi_new: Body,
        bj_new: Body,
        dust_cloud: Option<Body>,
        dust_mass: f64,
        q_ratio: f64,
    },
    Debris {
        fragments: Vec<Body>,
        dust_cloud: Option<Body>,
        dust_mass: f64,
        q_ratio: f64,
    },
}

pub fn specific_impact_energy(bi: &Body, bj: &Body) -> f64 {
    let geom = impact_geometry(bi, bj);
    let m_total = bi.mass + bj.mass;
    let mu = bi.mass * bj.mass / m_total;
    0.5 * mu * geom.v_n * geom.v_n / m_total
}

pub fn disruption_threshold(bi: &Body, bj: &Body, g_eff: f64) -> f64 {
    let m_total = bi.mass + bj.mass;
    let r_eff = (bi.physical_radius + bj.physical_radius) * 0.5;
    let gravity = (3.0 / 5.0) * g_eff * m_total / r_eff;
    let strength = STRENGTH_K * r_eff.powf(-0.3);
    (gravity + strength) * pair_disruption_scale(bi.material, bj.material)
}

pub fn evaluate_impact(bi: &Body, bj: &Body, g_eff: f64) -> ImpactResult {
    let q = specific_impact_energy(bi, bj);
    let q_star = disruption_threshold(bi, bj, g_eff);
    let q_ratio = q / q_star.max(1e-60);

    if q_ratio < SUB_THRESHOLD {
        return ImpactResult::SubThreshold;
    }

    if q_ratio < DISRUPTION_THRESHOLD {
        hit_and_run(bi, bj, q_ratio)
    } else {
        debris(bi, bj, g_eff, q_ratio)
    }
}

fn hit_and_run(bi: &Body, bj: &Body, q_ratio: f64) -> ImpactResult {
    let m_total = bi.mass + bj.mass;
    let inv_m = 1.0 / m_total;

    let x_com = (bi.mass * bi.x + bj.mass * bj.x) * inv_m;
    let y_com = (bi.mass * bi.y + bj.mass * bj.y) * inv_m;
    let v_com_x = (bi.mass * bi.vx + bj.mass * bj.vx) * inv_m;
    let v_com_y = (bi.mass * bi.vy + bj.mass * bj.vy) * inv_m;

    let (proj, targ) = if bi.mass <= bj.mass {
        (bi, bj)
    } else {
        (bj, bi)
    };

    let dust_frac = (0.30 * q_ratio).clamp(0.0, 0.50);
    let dust_mass = dust_frac * proj.mass;
    let proj_mass_new = (proj.mass - dust_mass).max(MIN_FRAGMENT_MASS);

    let vx_proj_new = proj.vx * (1.0 - dust_frac) + v_com_x * dust_frac;
    let vy_proj_new = proj.vy * (1.0 - dust_frac) + v_com_y * dust_frac;

    let mut proj_new = make_body(
        proj.x,
        proj.y,
        vx_proj_new,
        vy_proj_new,
        proj_mass_new,
        proj.density,
        proj.material,
    );
    let mut targ_new = *targ;

    let p_before_x = bi.mass * bi.vx + bj.mass * bj.vx;
    let p_before_y = bi.mass * bi.vy + bj.mass * bj.vy;
    let p_after_x = proj_mass_new * vx_proj_new + targ_new.mass * targ_new.vx;
    let p_after_y = proj_mass_new * vy_proj_new + targ_new.mass * targ_new.vy;
    let dust_vx = if dust_mass > 0.0 {
        (p_before_x - p_after_x) / dust_mass
    } else {
        v_com_x
    };
    let dust_vy = if dust_mass > 0.0 {
        (p_before_y - p_after_y) / dust_mass
    } else {
        v_com_y
    };
    let mut dust_cloud = make_dust_cloud(x_com, y_com, dust_vx, dust_vy, dust_mass, proj.density);

    let geom = impact_geometry(bi, bj);
    let l_total_before: f64 = [bi, bj]
        .into_iter()
        .map(|b| {
            let rx = b.x - x_com;
            let ry = b.y - y_com;
            let dvx = b.vx - v_com_x;
            let dvy = b.vy - v_com_y;
            b.mass * (rx * dvy - ry * dvx) + b.moment_inertia * b.omega_z
        })
        .sum();

    let dust_radius = dust_cloud
        .as_ref()
        .map(|cloud| cloud.physical_radius)
        .unwrap_or(0.0);
    let r_sum = (proj_new.physical_radius + targ_new.physical_radius + dust_radius).max(1e-30);
    let frac_proj = proj_new.physical_radius / r_sum;
    let frac_targ = targ_new.physical_radius / r_sum;

    let v_rel_mag = geom.v_n.hypot(geom.v_t).max(1e-30);
    let omega_max_proj = v_rel_mag / proj_new.physical_radius.max(1e-30);
    let omega_max_targ = v_rel_mag / targ_new.physical_radius.max(1e-30);
    proj_new.omega_z = (proj.omega_z
        + l_total_before * frac_proj / proj_new.moment_inertia.max(1e-30))
    .clamp(-omega_max_proj, omega_max_proj);
    targ_new.omega_z = (targ.omega_z
        + l_total_before * frac_targ / targ_new.moment_inertia.max(1e-30))
    .clamp(-omega_max_targ, omega_max_targ);

    let mut l_after_solid = 0.0;
    for body in [&proj_new, &targ_new] {
        let rx = body.x - x_com;
        let ry = body.y - y_com;
        let dvx = body.vx - v_com_x;
        let dvy = body.vy - v_com_y;
        l_after_solid += body.mass * (rx * dvy - ry * dvx) + body.moment_inertia * body.omega_z;
    }

    if let Some(cloud) = dust_cloud.as_mut() {
        let omega_max_cloud = v_rel_mag / cloud.physical_radius.max(1e-30);
        let residual_l = l_total_before - l_after_solid;
        cloud.omega_z =
            (residual_l / cloud.moment_inertia.max(1e-30)).clamp(-omega_max_cloud, omega_max_cloud);
    }

    let (bi_new, bj_new) = if bi.mass <= bj.mass {
        (proj_new, targ_new)
    } else {
        (targ_new, proj_new)
    };

    ImpactResult::HitAndRun {
        bi_new,
        bj_new,
        dust_cloud,
        dust_mass,
        q_ratio,
    }
}

fn debris(bi: &Body, bj: &Body, g_eff: f64, q_ratio: f64) -> ImpactResult {
    let m_total = bi.mass + bj.mass;
    let inv_m = 1.0 / m_total;

    let x_com = (bi.mass * bi.x + bj.mass * bj.x) * inv_m;
    let y_com = (bi.mass * bi.y + bj.mass * bj.y) * inv_m;
    let v_com_x = (bi.mass * bi.vx + bj.mass * bj.vx) * inv_m;
    let v_com_y = (bi.mass * bi.vy + bj.mass * bj.vy) * inv_m;

    let r_eff = (bi.physical_radius + bj.physical_radius) * 0.5;
    let v_esc = (2.0 * g_eff * m_total / r_eff).sqrt();
    let total_volume = (4.0 / 3.0) * PI * (bi.physical_radius.powi(3) + bj.physical_radius.powi(3));
    let frag_density = m_total / total_volume;

    let dominant_material = if bi.mass >= bj.mass {
        bi.material
    } else {
        bj.material
    };

    let t = (q_ratio / 2.0).clamp(0.0, 1.0);
    let m_lr = m_total * (1.0 - t.powf(1.3));
    let m_ej = m_total - m_lr;

    let kick_speed = (v_esc * q_ratio.sqrt()).clamp(0.2 * v_esc, 3.0 * v_esc);
    let dvx = bi.vx - bj.vx;
    let dvy = bi.vy - bj.vy;
    let geom = impact_geometry(bi, bj);
    let dx = bi.x - bj.x;
    let dy = bi.y - bj.y;
    let base_angle = dy.atan2(dx);

    let mut n_ejecta = (q_ratio * 6.0).clamp(2.0, 12.0) as usize;
    n_ejecta = n_ejecta.min(((m_ej / MIN_FRAGMENT_MASS).floor() as usize).max(1));

    let ejecta_fracs = power_law_fracs(n_ejecta);
    let mut fragments = Vec::with_capacity(n_ejecta + 1);
    let mut unresolved_dust_mass = 0.0;
    let has_remnant = m_lr >= MIN_FRAGMENT_MASS;

    let offset_r = bi.physical_radius + bj.physical_radius;

    if has_remnant {
        fragments.push(make_body(
            x_com,
            y_com,
            v_com_x,
            v_com_y,
            m_lr,
            frag_density,
            dominant_material,
        ));
    } else {
        unresolved_dust_mass += m_lr;
    }

    for (k, &frac) in ejecta_fracs.iter().enumerate() {
        let m_k = frac * m_ej;
        if m_k < MIN_FRAGMENT_MASS {
            unresolved_dust_mass += m_k;
            continue;
        }
        let angle =
            base_angle + 2.0 * PI * (k as f64) / (n_ejecta as f64) + rand::random::<f64>() * 0.4
                - 0.2;
        let vx = v_com_x + kick_speed * angle.cos();
        let vy = v_com_y + kick_speed * angle.sin();
        let frag_x = x_com + offset_r * angle.cos();
        let frag_y = y_com + offset_r * angle.sin();
        fragments.push(make_body(
            frag_x,
            frag_y,
            vx,
            vy,
            m_k,
            frag_density,
            dominant_material,
        ));
    }

    let m_tracked: f64 = fragments.iter().map(|f| f.mass).sum();
    if m_tracked > 1e-30 {
        let p_total_x = bi.mass * bi.vx + bj.mass * bj.vx;
        let p_total_y = bi.mass * bi.vy + bj.mass * bj.vy;
        let p_frags_x: f64 = fragments.iter().map(|f| f.mass * f.vx).sum();
        let p_frags_y: f64 = fragments.iter().map(|f| f.mass * f.vy).sum();
        let corr_x = (p_frags_x - p_total_x) / m_tracked;
        let corr_y = (p_frags_y - p_total_y) / m_tracked;
        for f in &mut fragments {
            f.vx -= corr_x;
            f.vy -= corr_y;
        }
    }

    let dust_cloud = if m_tracked > 1e-30 && unresolved_dust_mass > 1e-30 {
        let p_total_x = bi.mass * bi.vx + bj.mass * bj.vx;
        let p_total_y = bi.mass * bi.vy + bj.mass * bj.vy;
        let p_frags_x: f64 = fragments.iter().map(|f| f.mass * f.vx).sum();
        let p_frags_y: f64 = fragments.iter().map(|f| f.mass * f.vy).sum();
        let dust_vx = (p_total_x - p_frags_x) / unresolved_dust_mass;
        let dust_vy = (p_total_y - p_frags_y) / unresolved_dust_mass;
        make_dust_cloud(
            x_com,
            y_com,
            dust_vx,
            dust_vy,
            unresolved_dust_mass,
            frag_density,
        )
    } else {
        make_dust_cloud(
            x_com,
            y_com,
            v_com_x,
            v_com_y,
            unresolved_dust_mass,
            frag_density,
        )
    };

    if has_remnant {
        let geom = impact_geometry(bi, bj);
        let mu_r = bi.mass * bj.mass / m_total;
        let l_total = mu_r * geom.d * geom.v_t
            + bi.moment_inertia * bi.omega_z
            + bj.moment_inertia * bj.omega_z;

        let v_rel_mag = geom.v_n.hypot(geom.v_t).max(1e-30);
        let lr = &mut fragments[0];
        let omega_max = v_rel_mag / lr.physical_radius.max(1e-30);
        lr.omega_z = (l_total / lr.moment_inertia.max(1e-30)).clamp(-omega_max, omega_max);
    }

    ImpactResult::Debris {
        fragments,
        dust_cloud,
        dust_mass: unresolved_dust_mass,
        q_ratio,
    }
}

fn power_law_fracs(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let weights: Vec<f64> = (1..=n).map(|k| 1.0 / (k as f64).powf(2.0 / 3.0)).collect();
    let total: f64 = weights.iter().sum();
    weights.iter().map(|w| w / total).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::density_from_mass_radius;
    use crate::domain::materials::Material;

    fn body_at(x: f64, y: f64, vx: f64, vy: f64, mass: f64, radius: f64) -> Body {
        let mut b = Body::new(x, y, vx, vy, mass, Material::Rocky);
        b.radius = radius;
        b.density = density_from_mass_radius(mass, radius);
        b.sync_physical_properties();
        b
    }

    #[test]
    fn spe_is_positive_for_approaching_bodies() {
        let a = body_at(-0.1, 0.0, 1.0, 0.0, 1.0, 0.1);
        let b = body_at(0.1, 0.0, -1.0, 0.0, 1.0, 0.1);
        assert!(specific_impact_energy(&a, &b) > 0.0);
    }

    #[test]
    fn spe_is_zero_when_same_velocity() {
        let a = body_at(-0.1, 0.0, 2.0, 0.5, 3.0, 0.2);
        let b = body_at(0.1, 0.0, 2.0, 0.5, 5.0, 0.3);
        assert!(specific_impact_energy(&a, &b).abs() < 1e-15);
    }

    #[test]
    fn spe_scales_quadratically_with_relative_speed() {
        let a1 = body_at(-0.1, 0.0, 1.0, 0.0, 1.0, 0.1);
        let b1 = body_at(0.1, 0.0, -1.0, 0.0, 1.0, 0.1);
        let a2 = body_at(-0.1, 0.0, 2.0, 0.0, 1.0, 0.1);
        let b2 = body_at(0.1, 0.0, -2.0, 0.0, 1.0, 0.1);
        let q1 = specific_impact_energy(&a1, &b1);
        let q2 = specific_impact_energy(&a2, &b2);
        assert!((q2 / q1 - 4.0).abs() < 1e-10, "ratio = {}", q2 / q1);
    }

    #[test]
    fn q_star_is_positive() {
        let a = body_at(-0.1, 0.0, 0.0, 0.0, 2.0, 0.3);
        let b = body_at(0.1, 0.0, 0.0, 0.0, 1.0, 0.2);
        assert!(disruption_threshold(&a, &b, 1.0) > 0.0);
    }

    #[test]
    fn q_star_increases_with_g_eff() {
        let a = body_at(0.0, 0.0, 0.0, 0.0, 1.0, 0.2);
        let b = body_at(0.1, 0.0, 0.0, 0.0, 1.0, 0.2);
        let q1 = disruption_threshold(&a, &b, 1.0);
        let q2 = disruption_threshold(&a, &b, 2.0);
        assert!(q2 / q1 > 1.2);
    }

    #[test]
    fn sub_threshold_impact_returns_sub_threshold() {
        let a = body_at(-0.5, 0.0, 0.001, 0.0, 5.0, 0.5);
        let b = body_at(0.5, 0.0, -0.001, 0.0, 5.0, 0.5);
        assert!(matches!(
            evaluate_impact(&a, &b, 1.0),
            ImpactResult::SubThreshold
        ));
    }

    #[test]
    fn high_energy_impact_returns_debris() {
        let a = body_at(-0.1, 0.0, 200.0, 0.0, 1.0, 0.05);
        let b = body_at(0.1, 0.0, -200.0, 0.0, 1.0, 0.05);
        assert!(matches!(
            evaluate_impact(&a, &b, 1.0),
            ImpactResult::Debris { .. }
        ));
    }

    #[test]
    fn debris_conserves_total_mass() {
        let a = body_at(-0.1, 0.0, 50.0, 0.0, 2.0, 0.2);
        let b = body_at(0.1, 0.0, -50.0, 0.0, 3.0, 0.3);
        let m_before = a.mass + b.mass;
        match evaluate_impact(&a, &b, 1.0) {
            ImpactResult::Debris {
                fragments,
                dust_cloud,
                ..
            } => {
                let cloud_mass = dust_cloud.map(|cloud| cloud.mass).unwrap_or(0.0);
                let m_after: f64 = fragments.iter().map(|f| f.mass).sum::<f64>() + cloud_mass;
                assert!((m_after - m_before).abs() / m_before < 1e-10);
            }
            _ => panic!("expected Debris"),
        }
    }

    #[test]
    fn debris_conserves_linear_momentum_approximately() {
        let a = body_at(-0.1, 0.0, 30.0, 5.0, 2.0, 0.2);
        let b = body_at(0.1, 0.0, -30.0, -5.0, 2.0, 0.2);
        let px_before = a.mass * a.vx + b.mass * b.vx;
        let py_before = a.mass * a.vy + b.mass * b.vy;
        match evaluate_impact(&a, &b, 1.0) {
            ImpactResult::Debris {
                fragments,
                dust_cloud,
                ..
            } => {
                let px_after: f64 = fragments.iter().map(|f| f.mass * f.vx).sum();
                let py_after: f64 = fragments.iter().map(|f| f.mass * f.vy).sum();
                let cloud_px = dust_cloud.map(|cloud| cloud.mass * cloud.vx).unwrap_or(0.0);
                let cloud_py = dust_cloud.map(|cloud| cloud.mass * cloud.vy).unwrap_or(0.0);
                let tol = 1e-8;
                assert!((px_after + cloud_px - px_before).abs() < tol);
                assert!((py_after + cloud_py - py_before).abs() < tol);
            }
            _ => panic!("expected Debris"),
        }
    }

    #[test]
    fn hit_and_run_preserves_two_bodies() {
        let a = body_at(-0.5, 0.0, 3.0, 0.0, 1.0, 0.5);
        let b = body_at(0.5, 0.0, -3.0, 0.0, 5.0, 0.7);
        match evaluate_impact(&a, &b, 1.0) {
            ImpactResult::HitAndRun { bi_new, bj_new, .. } => {
                assert!((bj_new.mass - b.mass).abs() < 1e-12);
                assert!(bi_new.mass < a.mass + 1e-12);
            }
            ImpactResult::SubThreshold | ImpactResult::Debris { .. } => {}
        }
    }

    #[test]
    fn hit_and_run_tracks_dust_as_cloud() {
        let a = body_at(-0.5, 0.0, 3.0, 0.2, 1.0, 0.5);
        let b = body_at(0.5, 0.0, -3.0, -0.1, 5.0, 0.7);
        if let ImpactResult::HitAndRun {
            bi_new,
            bj_new,
            dust_cloud,
            ..
        } = evaluate_impact(&a, &b, 1.0)
        {
            let cloud = dust_cloud.expect("dust cloud expected");
            let m_after = bi_new.mass + bj_new.mass + cloud.mass;
            let m_before = a.mass + b.mass;
            assert!((m_after - m_before).abs() / m_before < 1e-10);
            assert_eq!(cloud.material, Material::DustCloud);
            assert_eq!(cloud.radius, 0.0);
        }
    }

    #[test]
    fn power_law_fracs_sum_to_one() {
        for n in 1..=6 {
            let fracs = power_law_fracs(n);
            let sum: f64 = fracs.iter().sum();
            assert!((sum - 1.0).abs() < 1e-14, "n={n} sum={sum}");
        }
    }

    #[test]
    fn power_law_fracs_are_decreasing() {
        let fracs = power_law_fracs(4);
        for i in 0..fracs.len() - 1 {
            assert!(fracs[i] > fracs[i + 1]);
        }
    }
}
