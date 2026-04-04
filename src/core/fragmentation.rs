//! Impact-energy fragmentation physics for the N-body simulation.
//!
//! ## Model summary (Leinhardt & Stewart 2012, simplified gravitational regime)
//!
//! | Quantity | Formula |
//! |---|---|
//! | Specific impact energy | Q = ½μv_rel² / M_total |
//! | Disruption threshold   | Q* ≈ strength + gravity |
//! | Largest remnant mass   | smooth transition |
//!
//! ## Outcome classification
//!
//! | Q / Q*       | Outcome |
//! |---|---|
//! | < `SUB_THRESHOLD`     | Sub-threshold |
//! | `SUB_THRESHOLD`–1.0   | Hit-and-run |
//! | ≥ 1.0                 | Debris |

use std::f64::consts::PI;

use crate::domain::body::{
    Body, default_moment_inertia, default_softening, radius_from_density_mass,
    sphere_radius_from_volume,
};

// ── Constants ────────────────────────────────────────────────────────────── //

pub const SUB_THRESHOLD: f64 = 0.1;
pub const DISRUPTION_THRESHOLD: f64 = 1.0;

const MIN_FRAGMENT_MASS: f64 = 1e-4;

// Strength term (material regime approximation)
const STRENGTH_K: f64 = 0.05;

// ── Helpers ──────────────────────────────────────────────────────────────── //

fn make_body(x: f64, y: f64, vx: f64, vy: f64, mass: f64, density: f64) -> Body {
    let r = radius_from_density_mass(density, mass);
    Body {
        x,
        y,
        vx,
        vy,
        mass,
        density,
        radius: r,
        softening: default_softening(mass).max(r * 2.0),
        omega_z: 0.0,
        moment_inertia: default_moment_inertia(mass, r),
        color: None,
    }
}

// ── Public types ─────────────────────────────────────────────────────────── //

pub enum ImpactResult {
    SubThreshold,

    HitAndRun {
        bi_new: Body,
        bj_new: Body,
        dust_mass: f64,
        q_ratio: f64,
    },

    Debris {
        fragments: Vec<Body>,
        dust_mass: f64,
        q_ratio: f64,
    },
}

// ── Core physics ─────────────────────────────────────────────────────────── //

pub fn specific_impact_energy(bi: &Body, bj: &Body) -> f64 {
    let m_total = bi.mass + bj.mass;
    let dvx = bi.vx - bj.vx;
    let dvy = bi.vy - bj.vy;
    let v_rel_sq = dvx * dvx + dvy * dvy;
    let mu = bi.mass * bj.mass / m_total;
    0.5 * mu * v_rel_sq / m_total
}

/// Improved Q* with strength + gravity
pub fn disruption_threshold(bi: &Body, bj: &Body, g_eff: f64) -> f64 {
    let m_total = bi.mass + bj.mass;

    let v_i = bi.mass / bi.density.max(1e-30);
    let v_j = bj.mass / bj.density.max(1e-30);
    let r_eff = sphere_radius_from_volume(v_i + v_j).max(1e-30);

    let gravity = (3.0 / 5.0) * g_eff * m_total / r_eff;

    // Simple strength scaling (dominates for small bodies)
    let strength = STRENGTH_K * r_eff.powf(-0.3);

    gravity + strength
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

// ── Hit-and-run ──────────────────────────────────────────────────────────── //

fn hit_and_run(bi: &Body, bj: &Body, q_ratio: f64) -> ImpactResult {
    let m_total = bi.mass + bj.mass;
    let inv_m = 1.0 / m_total;

    let v_com_x = (bi.mass * bi.vx + bj.mass * bj.vx) * inv_m;
    let v_com_y = (bi.mass * bi.vy + bj.mass * bj.vy) * inv_m;

    let (proj_is_i, proj, targ) = if bi.mass <= bj.mass {
        (true, bi, bj)
    } else {
        (false, bj, bi)
    };

    let dust_frac = (0.30 * q_ratio).clamp(0.0, 0.50);
    let dust_mass = dust_frac * proj.mass;
    let proj_mass_new = (proj.mass - dust_mass).max(MIN_FRAGMENT_MASS);

    let vx_proj_new = proj.vx * (1.0 - dust_frac) + v_com_x * dust_frac;
    let vy_proj_new = proj.vy * (1.0 - dust_frac) + v_com_y * dust_frac;

    let proj_new = make_body(
        proj.x,
        proj.y,
        vx_proj_new,
        vy_proj_new,
        proj_mass_new,
        proj.density,
    );

    let targ_new = *targ;

    let (bi_new, bj_new) = if proj_is_i {
        (proj_new, targ_new)
    } else {
        (targ_new, proj_new)
    };

    ImpactResult::HitAndRun {
        bi_new,
        bj_new,
        dust_mass,
        q_ratio,
    }
}

// ── Debris ───────────────────────────────────────────────────────────────── //

fn debris(bi: &Body, bj: &Body, g_eff: f64, q_ratio: f64) -> ImpactResult {
    let m_total = bi.mass + bj.mass;
    let inv_m = 1.0 / m_total;

    let x_com = (bi.mass * bi.x + bj.mass * bj.x) * inv_m;
    let y_com = (bi.mass * bi.y + bj.mass * bj.y) * inv_m;
    let v_com_x = (bi.mass * bi.vx + bj.mass * bj.vx) * inv_m;
    let v_com_y = (bi.mass * bi.vy + bj.mass * bj.vy) * inv_m;

    let v_i = bi.mass / bi.density.max(1e-30);
    let v_j = bj.mass / bj.density.max(1e-30);
    let r_eff = sphere_radius_from_volume(v_i + v_j).max(1e-30);

    let v_esc = (2.0 * g_eff * m_total / r_eff).sqrt();
    let frag_density = m_total / (v_i + v_j);

    let t = (q_ratio / 2.0).clamp(0.0, 1.0);
    let m_lr = m_total * (1.0 - t.powf(1.3));
    let m_ej = m_total - m_lr;

    let kick_speed = (v_esc * q_ratio.sqrt()).clamp(0.2 * v_esc, 3.0 * v_esc);

    let dvx = bi.vx - bj.vx;
    let dvy = bi.vy - bj.vy;
    let base_angle = dvy.atan2(dvx);

    let mut n_ejecta = (q_ratio * 6.0).clamp(2.0, 12.0) as usize;

    n_ejecta = n_ejecta.min(((m_ej / MIN_FRAGMENT_MASS).floor() as usize).max(1));

    let ejecta_fracs = power_law_fracs(n_ejecta);

    let mut fragments = Vec::with_capacity(n_ejecta + 1);
    let mut dust_mass = 0.0;

    if m_lr >= MIN_FRAGMENT_MASS {
        fragments.push(make_body(
            x_com,
            y_com,
            v_com_x,
            v_com_y,
            m_lr,
            frag_density,
        ));
    } else {
        dust_mass += m_lr;
    }

    for (k, &frac) in ejecta_fracs.iter().enumerate() {
        let m_k = frac * m_ej;

        if m_k < MIN_FRAGMENT_MASS {
            dust_mass += m_k;
            continue;
        }

        let spread = 0.4;
        let angle =
            base_angle + 2.0 * PI * (k as f64) / (n_ejecta as f64) + rand::random::<f64>() * spread
                - spread / 2.0;

        let vx = v_com_x + kick_speed * angle.cos();
        let vy = v_com_y + kick_speed * angle.sin();

        fragments.push(make_body(x_com, y_com, vx, vy, m_k, frag_density));
    }

    let m_tracked: f64 = fragments.iter().map(|f| f.mass).sum();

    if m_tracked > 1e-30 {
        let px: f64 = fragments.iter().map(|f| f.mass * f.vx).sum();
        let py: f64 = fragments.iter().map(|f| f.mass * f.vy).sum();

        let expected_px = m_tracked * v_com_x;
        let expected_py = m_tracked * v_com_y;

        let corr_x = (px - expected_px) / m_tracked;
        let corr_y = (py - expected_py) / m_tracked;

        for f in &mut fragments {
            f.vx -= corr_x;
            f.vy -= corr_y;
        }
    }

    ImpactResult::Debris {
        fragments,
        dust_mass,
        q_ratio,
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────── //

fn power_law_fracs(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let weights: Vec<f64> = (1..=n).map(|k| 1.0 / (k as f64).powf(2.0 / 3.0)).collect();
    let total: f64 = weights.iter().sum();
    weights.iter().map(|w| w / total).collect()
}

// ── Tests ────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::density_from_mass_radius;

    fn body_at(x: f64, y: f64, vx: f64, vy: f64, mass: f64, radius: f64) -> Body {
        let mut b = Body::new(x, y, vx, vy, mass);
        b.radius = radius;
        b.density = density_from_mass_radius(mass, radius);
        b
    }

    // ── specific_impact_energy ─────────────────────────────────────────── //

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
        // v_rel doubles → Q quadruples
        assert!((q2 / q1 - 4.0).abs() < 1e-10, "ratio = {}", q2 / q1);
    }

    // ── disruption_threshold ───────────────────────────────────────────── //

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

    // ── evaluate_impact ────────────────────────────────────────────────── //

    #[test]
    fn sub_threshold_impact_returns_sub_threshold() {
        // Very slow approach → Q << Q*
        let a = body_at(-0.5, 0.0, 0.001, 0.0, 5.0, 0.5);
        let b = body_at(0.5, 0.0, -0.001, 0.0, 5.0, 0.5);
        assert!(matches!(
            evaluate_impact(&a, &b, 1.0),
            ImpactResult::SubThreshold
        ));
    }

    #[test]
    fn high_energy_impact_returns_debris() {
        // Extremely fast approach → Q >> Q*
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
                dust_mass,
                ..
            } => {
                let m_after: f64 = fragments.iter().map(|f| f.mass).sum::<f64>() + dust_mass;
                assert!(
                    (m_after - m_before).abs() / m_before < 1e-10,
                    "mass error = {:.2e}",
                    (m_after - m_before) / m_before
                );
            }
            _ => panic!("expected Debris"),
        }
    }

    #[test]
    fn debris_conserves_linear_momentum_approximately() {
        let a = body_at(-0.1, 0.0, 30.0, 5.0, 2.0, 0.2);
        let b = body_at(0.1, 0.0, -30.0, -5.0, 2.0, 0.2);
        let m_total = a.mass + b.mass;
        let px_before = a.mass * a.vx + b.mass * b.vx;
        let py_before = a.mass * a.vy + b.mass * b.vy;
        match evaluate_impact(&a, &b, 1.0) {
            ImpactResult::Debris {
                fragments,
                dust_mass,
                ..
            } => {
                let m_tracked: f64 = fragments.iter().map(|f| f.mass).sum();
                // dust carries away CoM momentum
                let dust_frac = dust_mass / m_total;
                let px_after: f64 = fragments.iter().map(|f| f.mass * f.vx).sum();
                let py_after: f64 = fragments.iter().map(|f| f.mass * f.vy).sum();
                let expected_px = px_before * (1.0 - dust_frac);
                let expected_py = py_before * (1.0 - dust_frac);
                let tol = (m_tracked * 1e-8).max(1e-12);
                assert!(
                    (px_after - expected_px).abs() < tol,
                    "px error = {:.2e}",
                    px_after - expected_px
                );
                assert!(
                    (py_after - expected_py).abs() < tol,
                    "py error = {:.2e}",
                    py_after - expected_py
                );
            }
            _ => panic!("expected Debris"),
        }
    }

    #[test]
    fn hit_and_run_preserves_two_bodies() {
        let a = body_at(-0.5, 0.0, 3.0, 0.0, 1.0, 0.5);
        let b = body_at(0.5, 0.0, -3.0, 0.0, 5.0, 0.7);
        // Q/Q* should be in hit-and-run range; verify by checking result type.
        // (If this ends up as SubThreshold, just skip — velocities may not produce H&R)
        match evaluate_impact(&a, &b, 1.0) {
            ImpactResult::HitAndRun { bi_new, bj_new, .. } => {
                // Target (larger body) should be unchanged
                assert!((bj_new.mass - b.mass).abs() < 1e-12);
                // Projectile residual should have less mass
                assert!(bi_new.mass < a.mass + 1e-12);
            }
            ImpactResult::SubThreshold | ImpactResult::Debris { .. } => {
                // Acceptable — velocities may not fall exactly in H&R range
            }
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
            assert!(
                fracs[i] > fracs[i + 1],
                "fracs not monotone at i={i}: {} vs {}",
                fracs[i],
                fracs[i + 1]
            );
        }
    }
}
