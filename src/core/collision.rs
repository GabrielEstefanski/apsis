//! Collision detection and resolution for the N-body simulation.

use crate::core::fragmentation::{self, ImpactResult};
use crate::domain::body::{Body, default_moment_inertia, sphere_radius_from_volume, sphere_volume};
use crate::domain::materials::pair_restitution;
use std::collections::VecDeque;

/// CoR <= this value selects **astrophysics mode**:
/// - gravitationally *bound* contacts → merge
/// - gravitationally *unbound* contacts → pass-through (gravity already
///   handles the interaction; no impulse is injected)
///
/// CoR > this value selects **arcade mode**: impulse bounce for all contacts.
pub const MERGE_COR_THRESHOLD: f64 = 0.1;

/// Aggregate outcome of collision handling for a single step.
#[derive(Debug, Clone, Default)]
pub struct CollisionOutcome {
    pub merges: usize,
    pub bounces: usize,
    pub near_misses: usize,
    /// Number of individually tracked fragment bodies spawned this step.
    pub fragments_spawned: usize,
    /// Number of hit-and-run events this step.
    pub hit_and_runs: usize,
    /// Total ejecta mass below the fragment tracking threshold.
    pub total_dust_mass: f64,
    /// Visual feedback events; consumed and cleared by the renderer each frame.
    pub impact_events: Vec<ImpactEvent>,
}

/// Snapshot of a collision event for visual feedback.
#[derive(Debug, Clone, Copy)]
pub struct ImpactEvent {
    /// Centre-of-mass of the two bodies at the moment of contact (world coords).
    pub x: f64,
    pub y: f64,
    /// Outward normal from bj toward bi (unit vector).
    pub nx: f64,
    pub ny: f64,
    /// Relative speed |vi − vj| at contact.
    pub v_rel: f64,
}

/// Result of the normalized sub-step contact solver.
pub struct ContactEvent {
    pub t_frac: f64,
}

/// Earliest collision scheduled inside a time window.
pub struct ScheduledCollision {
    pub time: f64,
    pub i: usize,
    pub j: usize,
}

enum CollisionResponse {
    Bounce((f64, f64, f64, f64)),
    Merge,
    HitAndRun {
        bi_new: Body,
        bj_new: Body,
        dust_cloud: Option<Body>,
        dust_mass: f64,
    },
    Fragments {
        bodies: Vec<Body>,
        dust_cloud: Option<Body>,
        dust_mass: f64,
    },
    None,
}

// ── Public API ────────────────────────────────────────────────────────────── //

/// Find the first contact time between two bodies using a **linear** trajectory
/// approximation (end-points only — no acceleration).
///
/// Returns the normalized time `t_frac ∈ [0, 1]` at which the separation
/// first equals `contact_dist`, or `None` if no contact occurs in the window.
///
/// Kept as `pub` for unit tests; the simulation loop uses
/// [`find_earliest_contact`] which applies the quadratic (accelerated) solver.
pub fn find_contact_time(
    r_i0: (f64, f64),
    r_i1: (f64, f64),
    r_j0: (f64, f64),
    r_j1: (f64, f64),
    contact_dist: f64,
) -> Option<ContactEvent> {
    let dx0 = r_i0.0 - r_j0.0;
    let dy0 = r_i0.1 - r_j0.1;

    if dx0 * dx0 + dy0 * dy0 <= contact_dist * contact_dist {
        return Some(ContactEvent { t_frac: 0.0 });
    }

    let ddx = (r_i1.0 - r_i0.0) - (r_j1.0 - r_j0.0);
    let ddy = (r_i1.1 - r_i0.1) - (r_j1.1 - r_j0.1);

    let a = ddx * ddx + ddy * ddy;
    if a < 1e-30 {
        return None;
    }

    let b = 2.0 * (dx0 * ddx + dy0 * ddy);
    let c = dx0 * dx0 + dy0 * dy0 - contact_dist * contact_dist;
    let discriminant = b * b - 4.0 * a * c;

    if discriminant < 0.0 {
        return None;
    }

    let t_frac = (-b - discriminant.sqrt()) / (2.0 * a);

    if t_frac > 1.0 + 1e-10 || t_frac < -1e-10 {
        return None;
    }

    Some(ContactEvent {
        t_frac: t_frac.max(0.0),
    })
}

/// Find the earliest approaching contact inside `[0, max_dt]`.
///
/// Uses a **quadratic** trajectory approximation:
///
/// ```text
/// Δr(τ) = Δr₀ + Δv·τ + ½Δa·τ²
/// ```
///
/// where `Δa` comes from `acc` (the accelerations at the start of the window).
/// Pass `&[]` to fall back to a linear (velocity-only) trajectory.
///
/// The approach direction is verified at the contact time: pairs that are
/// already in contact but separating are ignored.
pub fn find_earliest_contact(
    bodies: &[Body],
    acc: &[(f64, f64)],
    max_dt: f64,
) -> Option<ScheduledCollision> {
    if max_dt <= 0.0 {
        return None;
    }

    let mut best: Option<ScheduledCollision> = None;

    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let bi = &bodies[i];
            let bj = &bodies[j];
            let contact = bi.radius + bj.radius;
            if contact <= 0.0 {
                continue;
            }

            let p = (bi.x - bj.x, bi.y - bj.y);
            let q = (bi.vx - bj.vx, bi.vy - bj.vy);

            let (ai, aj) = if i < acc.len() && j < acc.len() {
                (acc[i], acc[j])
            } else {
                ((0.0, 0.0), (0.0, 0.0))
            };
            // w = ½ · (a_i − a_j)
            let w = (0.5 * (ai.0 - aj.0), 0.5 * (ai.1 - aj.1));

            let Some(tau) = find_contact_time_quadratic(p, q, w, contact, max_dt) else {
                continue;
            };

            let replace = match &best {
                Some(current) => tau < current.time,
                None => true,
            };

            if replace {
                best = Some(ScheduledCollision { time: tau, i, j });
            }
        }
    }

    best
}

/// Resolve a collision for a pair that is already at contact.
///
/// `g_eff` is the effective gravitational constant (`G₀ · g_factor`) used when
/// deciding whether the pair is gravitationally bound (merge vs. pass-through).
pub fn resolve_contact(
    bodies: &mut Vec<Body>,
    trails: &mut Vec<VecDeque<(f64, f64)>>,
    i: usize,
    j: usize,
    cor: f64,
    g_eff: f64,
) -> CollisionOutcome {
    let mut outcome = CollisionOutcome::default();

    if i >= bodies.len() || j >= bodies.len() || i == j {
        return outcome;
    }

    // Snapshot impact geometry *before* bodies are modified or removed.
    let impact_event = {
        let bi = &bodies[i];
        let bj = &bodies[j];
        let m_total = bi.mass + bj.mass;
        let x = (bi.mass * bi.x + bj.mass * bj.x) / m_total;
        let y = (bi.mass * bi.y + bj.mass * bj.y) / m_total;
        let dx = bi.x - bj.x;
        let dy = bi.y - bj.y;
        let d = (dx * dx + dy * dy).sqrt().max(1e-30);
        let dvx = bi.vx - bj.vx;
        let dvy = bi.vy - bj.vy;
        ImpactEvent {
            x,
            y,
            nx: dx / d,
            ny: dy / d,
            v_rel: (dvx * dvx + dvy * dvy).sqrt(),
        }
    };

    match collision_response(&bodies[i], &bodies[j], cor, g_eff) {
        CollisionResponse::Bounce((vxi, vyi, vxj, vyj)) => {
            bodies[i].vx = vxi;
            bodies[i].vy = vyi;
            bodies[j].vx = vxj;
            bodies[j].vy = vyj;
            outcome.bounces = 1;
            outcome.impact_events.push(impact_event);
        }
        CollisionResponse::Merge => {
            let merged = merge_pair(bodies[i], bodies[j]);
            bodies.swap_remove(j);
            trails.swap_remove(j);
            bodies[i] = merged;
            outcome.merges = 1;
            outcome.impact_events.push(impact_event);
        }
        CollisionResponse::HitAndRun {
            bi_new,
            bj_new,
            dust_cloud,
            dust_mass,
        } => {
            bodies[i] = bi_new;
            bodies[j] = bj_new;
            if let Some(cloud) = dust_cloud {
                bodies.push(cloud);
                trails.push(VecDeque::new());
            }
            outcome.hit_and_runs = 1;
            outcome.total_dust_mass += dust_mass;
            outcome.impact_events.push(impact_event);
        }
        CollisionResponse::Fragments {
            bodies: frags,
            dust_cloud,
            dust_mass,
        } => {
            let n = frags.len() + usize::from(dust_cloud.is_some());

            let (lo, hi) = if i < j { (i, j) } else { (j, i) };
            bodies.swap_remove(hi);
            trails.swap_remove(hi);
            bodies.swap_remove(lo);
            trails.swap_remove(lo);
            for frag in frags {
                bodies.push(frag);
                trails.push(VecDeque::new());
            }
            if let Some(cloud) = dust_cloud {
                bodies.push(cloud);
                trails.push(VecDeque::new());
            }
            outcome.fragments_spawned = n;
            outcome.total_dust_mass += dust_mass;
            outcome.impact_events.push(impact_event);
        }
        CollisionResponse::None => {
            // Pass-through: unbound fly-by in astrophysics mode, or already
            // separating. Count as a near-miss so the dt-controller can tighten.
            outcome.near_misses = 1;
        }
    }

    outcome
}

/// Backwards-compatible event-driven helper used by tests.
///
/// Drifts bodies with constant velocity and resolves contacts inside `[0, dt]`.
/// Accelerations are not available here, so the solver falls back to a
/// linear (velocity-only) trajectory.
pub fn detect_and_resolve(
    bodies: &mut Vec<Body>,
    trails: &mut Vec<VecDeque<(f64, f64)>>,
    dt: f64,
    _pre_positions: &[(f64, f64)],
    cor: f64,
    g_eff: f64,
) -> CollisionOutcome {
    let mut outcome = CollisionOutcome::default();
    let mut remaining = dt.max(0.0);
    let mut iterations = 0;
    let max_iterations = 16;

    while remaining > 1e-8 && iterations < max_iterations {
        iterations += 1;

        let Some(event) = find_earliest_contact(bodies, &[], remaining) else {
            drift_linear(bodies, remaining);
            break;
        };

        if event.time > 0.0 {
            drift_linear(bodies, event.time);
            remaining -= event.time;
        }

        let step_outcome = resolve_contact(bodies, trails, event.i, event.j, cor, g_eff);
        outcome.merges += step_outcome.merges;
        outcome.bounces += step_outcome.bounces;
        outcome.near_misses += step_outcome.near_misses;
        outcome.fragments_spawned += step_outcome.fragments_spawned;
        outcome.hit_and_runs += step_outcome.hit_and_runs;
        outcome.total_dust_mass += step_outcome.total_dust_mass;

        if step_outcome.merges == 0
            && step_outcome.bounces == 0
            && step_outcome.fragments_spawned == 0
            && step_outcome.hit_and_runs == 0
        {
            let advance = remaining.min(1e-8);
            if advance <= 0.0 {
                break;
            }
            drift_linear(bodies, advance);
            remaining -= advance;
        }
    }

    outcome
}

/// Resolve a contact between two approaching bodies using the impulse method.
///
/// Returns the new `(vx_i, vy_i, vx_j, vy_j)`, or `None` when the contact
/// should be treated as a pass-through (separating pair, unbound fly-by in
/// astrophysics mode, or merge signal).
pub fn resolve_bounce(bi: &Body, bj: &Body, cor: f64, g_eff: f64) -> Option<(f64, f64, f64, f64)> {
    match collision_response(bi, bj, cor, g_eff) {
        CollisionResponse::Bounce(v) => Some(v),
        CollisionResponse::Merge
        | CollisionResponse::HitAndRun { .. }
        | CollisionResponse::Fragments { .. }
        | CollisionResponse::None => None,
    }
}

/// Merge two bodies into one, conserving mass, linear momentum, volume, and spin.
///
/// ## Density model
///
/// Each body has a bulk density ρ = m/V.  The merged body's volume is the sum
/// of the individual volumes, and its density is the mass-weighted average:
///
/// ```text
/// V_i = m_i / ρ_i
/// V'  = V_i + V_j
/// ρ'  = (m_i + m_j) / V'
/// r'  = (3 V' / 4π)^(1/3)
/// ```
///
/// Softening is derived from the merged volume as well (same ∛-scaling),
/// then clamped so `r' ≤ ε'/2` to preserve the Plummer-flatcore invariant.
pub fn merge_pair(bi: Body, bj: Body) -> Body {
    let m = bi.mass + bj.mass;
    debug_assert!(m > 0.0, "merge_pair: total mass must be positive");

    // ── Kinematics ────────────────────────────────────────────────────────── //
    let inv_m = 1.0 / m;
    let x = (bi.mass * bi.x + bj.mass * bj.x) * inv_m;
    let y = (bi.mass * bi.y + bj.mass * bj.y) * inv_m;
    let vx = (bi.mass * bi.vx + bj.mass * bj.vx) * inv_m;
    let vy = (bi.mass * bi.vy + bj.mass * bj.vy) * inv_m;

    // ── Volume & density ──────────────────────────────────────────────────── //
    // Volume from each body's own density (V = m / ρ)
    let v_i = bi.mass / bi.density.max(1e-30);
    let v_j = bj.mass / bj.density.max(1e-30);
    let total_volume = v_i + v_j;

    let density = m / total_volume;

    let physical_radius = sphere_radius_from_volume(total_volume);

    let v_soft_i = sphere_volume(bi.softening).max(0.0);
    let v_soft_j = sphere_volume(bj.softening).max(0.0);
    let softening = sphere_radius_from_volume(v_soft_i + v_soft_j) * 2.0;
    let softening = softening.max(physical_radius * 2.0);

    // Enforce Plummer-flatcore invariant: r ≤ ε/2
    let radius = physical_radius.min(softening * 0.5);
    let moment_inertia = default_moment_inertia(m, physical_radius);

    // ── Angular momentum ──────────────────────────────────────────────────── //
    //
    // L_orbital = μ · (r_j − r_i) × (v_j − v_i)
    //           = μ · d · v_t   (v_t = n × Δv, CCW positive)
    //
    // This equals Σ mᵢ (rᵢ − r_com) × vᵢ for the two-body system, so it IS
    // the orbital angular momentum relative to the COM.  No absolute positions
    // are used: (r_j − r_i) is a pure relative vector.
    let mu = (bi.mass * bj.mass) * inv_m;
    let r_ij_x = bj.x - bi.x;
    let r_ij_y = bj.y - bi.y;
    let v_rel_x = bj.vx - bi.vx;
    let v_rel_y = bj.vy - bi.vy;
    let l_int = mu * (r_ij_x * v_rel_y - r_ij_y * v_rel_x);

    let l_spin = bi.moment_inertia * bi.omega_z + bj.moment_inertia * bj.omega_z;
    let omega_z = (l_int + l_spin) / moment_inertia.max(1e-30);

    // Clamp: surface equatorial speed (ω · R) must not exceed the relative
    // impact speed.  Frontal impact (v_t = 0) → l_int = 0 → ω ≈ 0.
    // Oblique impact (v_t large) → l_int large → ω large, but still bounded.
    let v_rel_mag = (v_rel_x * v_rel_x + v_rel_y * v_rel_y).sqrt();
    let omega_max = v_rel_mag / physical_radius.max(1e-30);
    let omega_z = omega_z.clamp(-omega_max, omega_max);

    Body {
        x,
        y,
        vx,
        vy,
        mass: m,
        softening,
        radius,
        physical_radius,
        density,
        omega_z,
        moment_inertia,

        // Dominant body (heavier) determines material and base colour.
        material: if bi.mass >= bj.mass {
            bi.material
        } else {
            bj.material
        },
        color: if bi.mass >= bj.mass {
            bi.color
        } else {
            bj.color
        },
    }
}

// ── Private helpers ───────────────────────────────────────────────────────── //

/// Decide how to respond to a contact between `bi` and `bj`.
///
/// **Astrophysics mode** (`cor ≤ MERGE_COR_THRESHOLD`):
/// - Gravitationally *bound* (`e_orb < 0`): `Merge`
/// - Gravitationally *unbound* (`e_orb ≥ 0`): `None` (pass-through)
///   Gravity already accounts for the deflection; applying an impulse here
///   would double-count it and inject phantom energy.
///
/// **Arcade mode** (`cor > MERGE_COR_THRESHOLD`):
/// - Apply impulse `J = −(1 + CoR) μ v_n` for all approaching contacts.
fn collision_response(bi: &Body, bj: &Body, cor: f64, g_eff: f64) -> CollisionResponse {
    if bi.is_diffuse_cloud() || bj.is_diffuse_cloud() {
        return CollisionResponse::None;
    }

    let dx = bi.x - bj.x;
    let dy = bi.y - bj.y;
    let d = (dx * dx + dy * dy).sqrt();
    if d < 1e-30 {
        return CollisionResponse::None;
    }

    let nx = dx / d;
    let ny = dy / d;

    let dvx = bi.vx - bj.vx;
    let dvy = bi.vy - bj.vy;

    let v_n = dvx * nx + dvy * ny;
    if v_n >= 0.0 {
        return CollisionResponse::None;
    }

    let v_rel_sq = dvx * dvx + dvy * dvy;
    let e_orb = 0.5 * v_rel_sq - g_eff * (bi.mass + bj.mass) / d.max(1e-30);

    let is_bound = e_orb < 0.0;

    if cor <= MERGE_COR_THRESHOLD {
        let impact = fragmentation::evaluate_impact(bi, bj, g_eff);

        if is_bound {
            match impact {
                ImpactResult::Debris {
                    fragments,
                    dust_cloud,
                    dust_mass,
                    ..
                } => {
                    return CollisionResponse::Fragments {
                        bodies: fragments,
                        dust_cloud,
                        dust_mass,
                    };
                }

                _ => {
                    return CollisionResponse::Merge;
                }
            }
        } else {
            match impact {
                ImpactResult::SubThreshold => {
                    return CollisionResponse::None;
                }

                ImpactResult::HitAndRun {
                    bi_new,
                    bj_new,
                    dust_cloud,
                    dust_mass,
                    ..
                } => {
                    return CollisionResponse::HitAndRun {
                        bi_new,
                        bj_new,
                        dust_cloud,
                        dust_mass,
                    };
                }

                ImpactResult::Debris {
                    fragments,
                    dust_cloud,
                    dust_mass,
                    ..
                } => {
                    return CollisionResponse::Fragments {
                        bodies: fragments,
                        dust_cloud,
                        dust_mass,
                    };
                }
            }
        }
    }

    // Arcade mode: impulse bounce.
    // Use the larger of the global CoR slider and the pair's material-derived
    // restitution so the slider always acts as a lower-bound override while
    // material physics still kicks in when the slider is at minimum.
    let effective_cor = cor.max(pair_restitution(bi.material, bj.material));
    let mu = bi.mass * bj.mass / (bi.mass + bj.mass);
    let j_mag = -(1.0 + effective_cor) * mu * v_n;

    CollisionResponse::Bounce((
        bi.vx + j_mag * nx / bi.mass,
        bi.vy + j_mag * ny / bi.mass,
        bj.vx - j_mag * nx / bj.mass,
        bj.vy - j_mag * ny / bj.mass,
    ))
}

/// Find the earliest contact time for one pair using a quadratic trajectory.
///
/// Trajectory model: `Δr(τ) = p + q·τ + w·τ²`
/// - `p` = relative position at start of window
/// - `q` = relative velocity at start of window
/// - `w` = ½ · relative acceleration (`½Δa`)
///
/// Algorithm:
/// 1. Initial-overlap check — returns `Some(0)` if already in contact AND approaching.
/// 2. Linear seed — quadratic formula on `|p + q·τ|² = R²`.
/// 3. Newton-Raphson refinement — corrects for the `w·τ²` term.
/// 4. Approach check at the refined contact time.
fn find_contact_time_quadratic(
    p: (f64, f64),
    q: (f64, f64),
    w: (f64, f64),
    r_contact: f64,
    max_dt: f64,
) -> Option<f64> {
    let r2 = r_contact * r_contact;

    // ── Step 1: already in contact? ─────────────────────────────────────── //
    let p_sq = p.0 * p.0 + p.1 * p.1;
    if p_sq <= r2 {
        // Only accept if the pair is still approaching (p · q < 0)
        let v_n = p.0 * q.0 + p.1 * q.1;
        return if v_n < 0.0 { Some(0.0) } else { None };
    }

    // ── Step 2: linear seed ──────────────────────────────────────────────── //
    // Solve (q·q)τ² + 2(p·q)τ + (p·p − R²) = 0
    let aq = q.0 * q.0 + q.1 * q.1;
    if aq < 1e-30 {
        return None; // Zero relative velocity
    }
    let bq = 2.0 * (p.0 * q.0 + p.1 * q.1);
    let cq = p_sq - r2;
    let disc = bq * bq - 4.0 * aq * cq;
    if disc < 0.0 {
        return None; // Linear path misses contact sphere
    }

    let t_lin = (-bq - disc.sqrt()) / (2.0 * aq);
    if t_lin < -1e-10 || t_lin > max_dt + 1e-10 {
        return None; // Contact outside the integration window
    }

    // ── Step 3: Newton-Raphson refinement ───────────────────────────────── //
    // Skip NR when acceleration term is negligible
    let w_sq = w.0 * w.0 + w.1 * w.1;
    if w_sq < 1e-30 {
        let tau = t_lin.clamp(0.0, max_dt);
        let rx = p.0 + q.0 * tau;
        let ry = p.1 + q.1 * tau;
        // Approach check: Δr · Δv < 0
        return if rx * q.0 + ry * q.1 < 0.0 {
            Some(tau)
        } else {
            None
        };
    }

    // f(τ)  = |p + q·τ + w·τ²|² − R²
    // f′(τ) = 2 · (p + q·τ + w·τ²) · (q + 2w·τ)
    let mut tau = t_lin.clamp(0.0, max_dt);
    for _ in 0..8 {
        let tau2 = tau * tau;
        let rx = p.0 + q.0 * tau + w.0 * tau2;
        let ry = p.1 + q.1 * tau + w.1 * tau2;
        let f = rx * rx + ry * ry - r2;

        let vrx = q.0 + 2.0 * w.0 * tau;
        let vry = q.1 + 2.0 * w.1 * tau;
        let df = 2.0 * (rx * vrx + ry * vry);

        if df.abs() < 1e-30 {
            break;
        }
        let step = f / df;
        tau -= step;
        if step.abs() < 1e-13 {
            break;
        }
    }

    // ── Step 4: validate ────────────────────────────────────────────────── //
    if tau < -1e-10 || tau > max_dt + 1e-10 {
        return None;
    }
    let tau = tau.clamp(0.0, max_dt);

    let tau2 = tau * tau;
    let rx = p.0 + q.0 * tau + w.0 * tau2;
    let ry = p.1 + q.1 * tau + w.1 * tau2;
    let vrx = q.0 + 2.0 * w.0 * tau;
    let vry = q.1 + 2.0 * w.1 * tau;

    // Approach check at the refined contact time
    if rx * vrx + ry * vry >= 0.0 {
        return None;
    }

    Some(tau)
}

fn drift_linear(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────── //

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::density_from_mass_radius;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::new(
            x,
            y,
            vx,
            vy,
            mass,
            crate::domain::materials::Material::Rocky,
        )
    }

    /// Create a body at position `(x, y)` with the given mass and a
    /// consistent `density` derived from the supplied radius.
    fn body_with_radius(x: f64, y: f64, vx: f64, vy: f64, mass: f64, r: f64) -> Body {
        let mut b = Body::new(
            x,
            y,
            vx,
            vy,
            mass,
            crate::domain::materials::Material::Rocky,
        );
        b.radius = r;
        b.density = density_from_mass_radius(mass, r);
        b.sync_physical_properties();
        b
    }

    fn contact_pair(
        vx_a: f64,
        vy_a: f64,
        mass_a: f64,
        vx_b: f64,
        vy_b: f64,
        mass_b: f64,
        r: f64,
    ) -> (Body, Body) {
        let a = body_with_radius(-r, 0.0, vx_a, vy_a, mass_a, r);
        let b = body_with_radius(r, 0.0, vx_b, vy_b, mass_b, r);
        (a, b)
    }

    fn make_scene(a: Body, b: Body) -> (Vec<Body>, Vec<VecDeque<(f64, f64)>>) {
        (vec![a, b], vec![VecDeque::new(), VecDeque::new()])
    }

    // ── merge_pair ─────────────────────────────────────────────────────── //

    #[test]
    fn merge_conserves_mass() {
        let (a, b) = contact_pair(0.0, 0.0, 3.0, 0.0, 0.0, 5.0, 0.1);
        let m = merge_pair(a, b);
        assert!((m.mass - 8.0).abs() < 1e-12);
    }

    #[test]
    fn merge_conserves_linear_momentum_x() {
        let (a, b) = contact_pair(2.0, 0.0, 1.0, -4.0, 0.0, 2.0, 0.05);
        let px_before = a.mass * a.vx + b.mass * b.vx;
        let m = merge_pair(a, b);
        assert!((m.mass * m.vx - px_before).abs() < 1e-12);
    }

    #[test]
    fn merge_conserves_linear_momentum_y() {
        let (a, b) = contact_pair(0.0, 3.0, 2.0, 0.0, -1.0, 4.0, 0.05);
        let py_before = a.mass * a.vy + b.mass * b.vy;
        let m = merge_pair(a, b);
        assert!((m.mass * m.vy - py_before).abs() < 1e-12);
    }

    #[test]
    fn merge_dissipates_kinetic_energy() {
        let (a, b) = contact_pair(3.0, 0.0, 1.0, -3.0, 0.0, 1.0, 0.1);
        let ke_before = 0.5 * a.mass * (a.vx.powi(2) + a.vy.powi(2))
            + 0.5 * b.mass * (b.vx.powi(2) + b.vy.powi(2));
        let m = merge_pair(a, b);
        let ke_after = 0.5 * m.mass * (m.vx.powi(2) + m.vy.powi(2));
        assert!(ke_after <= ke_before + 1e-12);
    }

    #[test]
    fn merge_places_body_at_com() {
        let a = body_with_radius(0.0, 0.0, 0.0, 0.0, 2.0, 0.5);
        let b = body_with_radius(6.0, 0.0, 0.0, 0.0, 4.0, 0.5);
        let m = merge_pair(a, b);
        assert!((m.x - 4.0).abs() < 1e-12);
        assert!(m.y.abs() < 1e-12);
    }

    #[test]
    fn merge_radius_cannot_shrink_below_either_constituent() {
        let (a, b) = contact_pair(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.1);
        let m = merge_pair(a, b);
        assert!(m.radius >= a.radius.max(b.radius));
    }

    #[test]
    fn off_centre_approach_produces_nonzero_spin() {
        let a = body_with_radius(-0.05, 0.1, 1.0, 0.0, 1.0, 0.2);
        let b = body_with_radius(0.05, -0.1, -1.0, 0.0, 1.0, 0.2);
        let m = merge_pair(a, b);
        assert!(m.omega_z.abs() > 1e-12);
    }

    #[test]
    fn head_on_collision_produces_zero_spin() {
        let a = body_with_radius(-0.1, 0.0, 1.0, 0.0, 1.0, 0.2);
        let b = body_with_radius(0.1, 0.0, -1.0, 0.0, 1.0, 0.2);
        let m = merge_pair(a, b);
        assert!(m.omega_z.abs() < 1e-12);
    }

    // ── find_contact_time (linear, pub) ────────────────────────────────── //

    #[test]
    fn find_contact_detects_crossing_trajectories() {
        let event = find_contact_time((0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 0.0), 0.1);
        assert!(event.is_some());
        let t = event.unwrap().t_frac;
        assert!(t > 0.0 && t < 1.0, "t_frac = {t}");
    }

    #[test]
    fn find_contact_returns_none_for_parallel_miss() {
        let event = find_contact_time((0.0, 0.0), (1.0, 0.0), (0.0, 5.0), (1.0, 5.0), 0.1);
        assert!(event.is_none());
    }

    #[test]
    fn find_contact_returns_zero_for_initial_overlap() {
        let event = find_contact_time((0.0, 0.0), (0.5, 0.0), (0.0, 0.0), (0.5, 0.0), 0.1);
        assert!(event.is_some());
        assert_eq!(event.unwrap().t_frac, 0.0);
    }

    #[test]
    fn find_contact_returns_none_when_bodies_diverge_without_touching() {
        let event = find_contact_time((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (4.0, 0.0), 0.1);
        assert!(event.is_none());
    }

    // ── find_contact_time_quadratic ────────────────────────────────────── //

    #[test]
    fn quadratic_ccd_matches_linear_when_acceleration_is_zero() {
        // With w=0, quadratic NR should give the same result as the linear solver.
        let p = (-1.0, 0.0);
        let q = (2.0, 0.0); // approaching: Δv = +2 means i moves toward j faster
        // Actually: p = r_i - r_j = -1 means i is to the left of j
        // q = v_i - v_j = 2 means i moves right faster → approaching
        // Contact at p + q*t = 0, so t = 0.5 → but contact at |Δr| = 0.1
        // p + q*t = -1 + 2t = ±0.1 → t = 0.45 or t = 0.55; first is 0.45
        let tau = find_contact_time_quadratic(p, q, (0.0, 0.0), 0.1, 1.0);
        assert!(tau.is_some(), "should detect contact");
        let t = tau.unwrap();
        assert!((t - 0.45).abs() < 1e-10, "t = {t}");
    }

    #[test]
    fn quadratic_ccd_refines_timing_with_acceleration() {
        // Two bodies approaching on a curved path.
        // Pure linear would predict a slightly different time than quadratic.
        let p = (-1.0, 0.0);
        let q = (2.0, 0.0);
        let w = (-0.5, 0.0); // deceleration (½Δa)
        // Quadratic: Δx(τ) = -1 + 2τ - 0.5τ²
        // Contact when |Δx| = 0.1: -1 + 2τ - 0.5τ² = -0.1 → 0.5τ² - 2τ + 0.9 = 0
        // τ = (2 ± √(4 - 1.8)) / 1 = (2 ± √2.2) / 1
        // τ₁ = 2 - 1.4832... ≈ 0.5168
        let tau = find_contact_time_quadratic(p, q, w, 0.1, 1.0);
        assert!(tau.is_some(), "should detect contact");
        let t = tau.unwrap();
        // Expected: (2 - sqrt(2.2)) / 1.0 ≈ 0.5168
        let expected = 2.0 - (2.2f64).sqrt();
        assert!(
            (t - expected).abs() < 1e-8,
            "t = {t}, expected ≈ {expected}"
        );
    }

    #[test]
    fn quadratic_ccd_returns_none_for_already_separating_overlap() {
        // Bodies overlapping but moving apart
        let p = (-0.05, 0.0); // |p| = 0.05 < R = 0.1 → in contact
        let q = (-2.0, 0.0); // p · q = (-0.05)(-2) = 0.1 > 0 → separating
        let tau = find_contact_time_quadratic(p, q, (0.0, 0.0), 0.1, 1.0);
        assert!(tau.is_none(), "separating overlap should return None");
    }

    #[test]
    fn quadratic_ccd_returns_zero_for_overlapping_approaching() {
        let p = (-0.05, 0.0); // in contact
        let q = (2.0, 0.0); // approaching: p · q = (-0.05)(2) = -0.1 < 0
        let tau = find_contact_time_quadratic(p, q, (0.0, 0.0), 0.1, 1.0);
        assert_eq!(tau, Some(0.0));
    }

    // ── find_earliest_contact ──────────────────────────────────────────── //

    #[test]
    fn earliest_contact_ignores_separating_overlap() {
        let a = body_with_radius(-0.001, 0.0, -1.0, 0.0, 1.0, 0.01);
        let b = body_with_radius(0.001, 0.0, 1.0, 0.0, 1.0, 0.01);
        let bodies = vec![a, b];
        assert!(find_earliest_contact(&bodies, &[], 0.01).is_none());
    }

    // ── resolve_bounce / collision_response ───────────────────────────── //

    #[test]
    fn bounce_cor1_conserves_kinetic_energy() {
        let (a, b) = contact_pair(2.0, 0.0, 1.0, -2.0, 0.0, 1.0, 0.5);
        let ke_before = 0.5 * (a.vx.powi(2) + b.vx.powi(2));
        let (vxi, vyi, vxj, vyj) = resolve_bounce(&a, &b, 1.0, 1.0).unwrap();
        let ke_after = 0.5 * (vxi.powi(2) + vyi.powi(2)) + 0.5 * (vxj.powi(2) + vyj.powi(2));
        assert!((ke_after - ke_before).abs() < 1e-10);
    }

    #[test]
    fn bounce_conserves_linear_momentum() {
        let (a, b) = contact_pair(3.0, 1.0, 2.0, -1.0, -0.5, 3.0, 0.5);
        let px_before = a.mass * a.vx + b.mass * b.vx;
        let py_before = a.mass * a.vy + b.mass * b.vy;
        let (vxi, vyi, vxj, vyj) = resolve_bounce(&a, &b, 0.6, 1.0).unwrap();
        let px_after = a.mass * vxi + b.mass * vxj;
        let py_after = a.mass * vyi + b.mass * vyj;
        assert!((px_after - px_before).abs() < 1e-10);
        assert!((py_after - py_before).abs() < 1e-10);
    }

    #[test]
    fn bounce_leaves_bodies_separating() {
        let (a, b) = contact_pair(2.0, 0.0, 1.0, -2.0, 0.0, 1.0, 0.5);
        let (vxi, _, vxj, _) = resolve_bounce(&a, &b, 0.5, 1.0).unwrap();
        let v_n_after = (vxi - vxj) * (-1.0);
        assert!(v_n_after > 0.0);
    }

    #[test]
    fn bound_inelastic_collision_signals_merge() {
        let (a, b) = contact_pair(0.01, 0.0, 10.0, -0.01, 0.0, 10.0, 0.5);
        let result = resolve_bounce(&a, &b, 0.0, 1.0);
        assert!(result.is_none());
    }

    /// In astrophysics mode (cor=0), an unbound fly-by must **pass through**
    /// without any impulse.  Injecting a zero-restitution bounce here would
    /// double-count the gravitational interaction and create phantom energy.
    #[test]
    fn unbound_flyby_passthrough_in_astro_mode() {
        let (a, b) = contact_pair(100.0, 0.0, 1.0, -100.0, 0.0, 1.0, 0.5);
        let result = resolve_bounce(&a, &b, 0.0, 1.0);
        assert!(
            result.is_none(),
            "expected pass-through for unbound fly-by in astrophysics mode"
        );
    }

    /// In arcade mode (cor > threshold), an unbound fly-by receives a normal
    /// elastic/partial-restitution bounce impulse.
    #[test]
    fn unbound_flyby_bounces_in_arcade_mode() {
        let (a, b) = contact_pair(100.0, 0.0, 1.0, -100.0, 0.0, 1.0, 0.5);
        let result = resolve_bounce(&a, &b, 0.5, 1.0);
        assert!(
            result.is_some(),
            "expected bounce for unbound fly-by in arcade mode (cor=0.5)"
        );
    }

    #[test]
    fn bounce_conserves_angular_momentum_z() {
        let (a, b) = contact_pair(3.0, 1.0, 1.0, -3.0, 0.0, 1.0, 0.5);
        let lz_before = a.mass * (a.x * a.vy - a.y * a.vx) + b.mass * (b.x * b.vy - b.y * b.vx);
        let (vxi, vyi, vxj, vyj) = resolve_bounce(&a, &b, 0.7, 1.0).unwrap();
        let lz_after = a.mass * (a.x * vyi - a.y * vxi) + b.mass * (b.x * vyj - b.y * vxj);
        assert!((lz_after - lz_before).abs() < 1e-10);
    }

    // ── detect_and_resolve (integration) ──────────────────────────────── //

    #[test]
    fn fly_by_with_positive_orbital_energy_is_not_merged() {
        let a = body_with_radius(-0.001, 0.0, 100.0, 0.0, 1.0, 0.01);
        let b = body_with_radius(0.001, 0.0, -100.0, 0.0, 1.0, 0.01);
        let (mut bodies, mut trails) = make_scene(a, b);
        let out = detect_and_resolve(&mut bodies, &mut trails, 0.01, &[], 0.0, 1.0);
        // High-energy encounter: no merge, but may produce fragments
        assert_eq!(out.merges, 0);
    }

    #[test]
    fn bound_slow_approach_triggers_merge() {
        let a = body_with_radius(-0.001, 0.0, 0.001, 0.0, 10.0, 0.01);
        let b = body_with_radius(0.001, 0.0, -0.001, 0.0, 10.0, 0.01);
        let (mut bodies, mut trails) = make_scene(a, b);
        let out = detect_and_resolve(&mut bodies, &mut trails, 0.01, &[], 0.0, 1.0);
        assert_eq!(out.merges, 1);
        assert_eq!(bodies.len(), 1);
    }

    #[test]
    fn separating_overlapping_bodies_are_not_merged() {
        let a = body_with_radius(-0.001, 0.0, -1.0, 0.0, 1.0, 0.01);
        let b = body_with_radius(0.001, 0.0, 1.0, 0.0, 1.0, 0.01);
        let (mut bodies, mut trails) = make_scene(a, b);
        let out = detect_and_resolve(&mut bodies, &mut trails, 0.01, &[], 0.0, 1.0);
        assert_eq!(out.merges, 0);
    }

    #[test]
    fn merged_body_conserves_total_mass_of_system() {
        let a = body_with_radius(-0.001, 0.0, 0.001, 0.0, 3.0, 0.01);
        let b = body_with_radius(0.001, 0.0, -0.001, 0.0, 5.0, 0.01);
        let m_before = a.mass + b.mass;
        let (mut bodies, mut trails) = make_scene(a, b);
        detect_and_resolve(&mut bodies, &mut trails, 0.01, &[], 0.0, 1.0);
        let m_after: f64 = bodies.iter().map(|b| b.mass).sum();
        assert!((m_after - m_before).abs() < 1e-12);
    }

    #[test]
    fn merged_body_conserves_total_linear_momentum() {
        let a = body_with_radius(-0.001, 0.0, 0.5, 0.0, 2.0, 0.01);
        let b = body_with_radius(0.001, 0.0, -0.5, 0.0, 3.0, 0.01);
        let px_before = a.mass * a.vx + b.mass * b.vx;
        let (mut bodies, mut trails) = make_scene(a, b);
        detect_and_resolve(&mut bodies, &mut trails, 0.01, &[], 0.0, 1.0);
        let px_after: f64 = bodies.iter().map(|b| b.mass * b.vx).sum();
        assert!((px_after - px_before).abs() < 1e-12);
    }

    #[test]
    fn sub_step_detection_catches_tunnelling_bodies() {
        let a = body_with_radius(0.0, 0.0, 0.0, 0.0, 1.0, 0.1);
        let b = body_with_radius(0.3, 0.0, 0.0, 0.0, 1.0, 0.1);
        let (mut bodies, mut trails) = make_scene(a, b);
        let pre = vec![(-5.0, 0.0), (5.0, 0.0)];
        let out = detect_and_resolve(&mut bodies, &mut trails, 0.01, &pre, 0.0, 1.0);
        assert_eq!(out.merges, 0);
    }

    // ── density / volume conservation ─────────────────────────────────── //

    /// V' = V_i + V_j: merged volume equals the sum of constituent volumes.
    #[test]
    fn merge_conserves_volume() {
        let a = body_with_radius(-0.1, 0.0, 0.0, 0.0, 2.0, 0.3);
        let b = body_with_radius(0.1, 0.0, 0.0, 0.0, 3.0, 0.5);
        let v_before = a.mass / a.density + b.mass / b.density;
        let m = merge_pair(a, b);
        // Use density to recover volume: V' = m' / ρ'
        let v_after = m.mass / m.density;
        assert!(
            (v_after - v_before).abs() / v_before < 1e-10,
            "volume error = {:.2e}",
            (v_after - v_before) / v_before
        );
    }

    /// ρ' = (m_i + m_j) / (V_i + V_j): density is mass-to-volume ratio.
    #[test]
    fn merge_density_equals_total_mass_over_total_volume() {
        let a = body_with_radius(-0.1, 0.0, 0.0, 0.0, 1.0, 0.2);
        let b = body_with_radius(0.1, 0.0, 0.0, 0.0, 4.0, 0.4);
        let v_total = a.mass / a.density + b.mass / b.density;
        let expected_density = (a.mass + b.mass) / v_total;
        let m = merge_pair(a, b);
        assert!(
            (m.density - expected_density).abs() / expected_density < 1e-10,
            "density error = {:.2e}",
            (m.density - expected_density) / expected_density
        );
    }

    /// Denser constituent always produces a smaller merged radius for equal mass.
    #[test]
    fn denser_body_produces_smaller_merged_radius() {
        // Both with same mass but different densities
        let a_dense = body_with_radius(0.0, 0.0, 0.0, 0.0, 1.0, 0.1); // small radius = dense
        let b_dense = body_with_radius(1.0, 0.0, 0.0, 0.0, 1.0, 0.1);
        let a_loose = body_with_radius(0.0, 0.0, 0.0, 0.0, 1.0, 0.5); // large radius = light
        let b_loose = body_with_radius(1.0, 0.0, 0.0, 0.0, 1.0, 0.5);
        let m_dense = merge_pair(a_dense, b_dense);
        let m_loose = merge_pair(a_loose, b_loose);
        assert!(
            m_dense.radius < m_loose.radius,
            "denser pair should produce smaller merged body: {} vs {}",
            m_dense.radius,
            m_loose.radius
        );
    }
}
