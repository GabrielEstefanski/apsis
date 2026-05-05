//! Visual smoothing of osculating orbital elements for the predicted-orbit
//! overlay.
//!
//! # Why this exists
//!
//! In an N-body system, the *osculating* ellipse — the Keplerian orbit
//! consistent with a body's instantaneous (r, v) — is **not** invariant.
//! Indirect perturbations on the primary (e.g. Jupiter tugging the Sun)
//! make the heliocentric energy of every other body oscillate at the
//! perturber's period. With time-warp, those oscillations alias into
//! seconds-of-wallclock, and the predicted-orbit overlay visibly jitters
//! every frame. The physics is correct; the *display* lies about the
//! orbit being stable.
//!
//! This module is a pure render-side EMA filter on the **linear primitives**
//! of the orbit (ε, ex, ey, h — see [`apsis::physics::orbital::OrbitInvariants`]),
//! reconstructing (a, e, ω) only at the end. It never feeds back into the
//! integrator.
//!
//! # Why linear primitives, not (a, e, ω)
//!
//! * `a = −GM/(2ε)` is non-linear in ε; smoothing `a` directly biases
//!   high-eccentricity orbits and is unstable near ε ≈ 0 (parabolic).
//! * `ω = atan2(ey, ex)` has a wraparound at ±π; EMA on the angle would
//!   alias π → −π as a discontinuity. Smoothing the Cartesian
//!   eccentricity-vector components (ex, ey) avoids this without needing
//!   sin/cos splitting or unwrap heuristics.
//! * (ε, ex, ey, h) are all linear (or bilinear) in the state vector,
//!   making EMA a faithful low-pass on the underlying signal.
//!
//! # Time constant τ
//!
//! Per body, τ is set in **simulation time** (not wallclock) so smoothing
//! attenuates the same physical frequency regardless of time-warp. The
//! dominant perturber is identified by the magnitude of the differential
//! acceleration each sibling exerts on the body in the primary's frame:
//!
//! ```text
//! Δa_j = G·M_j · [(r_p − r_j)/|r_p − r_j|³  −  (r_i − r_j)/|r_i − r_j|³]
//! ```
//!
//! τ is the *weighted mean* of sibling orbital periods, weighted by
//! `|Δa_j|`. This subsumes hysteresis (no argmax flicker) and naturally
//! handles bi-perturber regimes (Saturn perturbed by both Jupiter and
//! Uranus) without picking a winner.
//!
//! Final τ = K · weighted_period, clamped to [0.01·P_self, 3·P_self] to
//! protect against pathological perturbers and to bound the filter's
//! memory.

use apsis::domain::body::Body;
use apsis::physics::orbital::{
    OrbitInvariants, OrbitalElements, compute_invariants, elements_anchored_to_body,
    elements_from_invariants,
};
use std::collections::HashMap;
use std::f64::consts::TAU;

/// Default smoothing strength: τ_perturber = K · weighted_period.
///
/// K = 1.5 attenuates the dominant perturber's frequency by ~−20 dB (≈ 90%
/// signal removal) for a first-order EMA, while leaving secular dynamics
/// (precession, long-period drift) essentially untouched. Values below 1.0
/// leak jitter visibly; values above 2.5 begin clipping real precession.
const DEFAULT_K: f64 = 1.5;

/// Lower clamp for τ_sim, expressed as a fraction of the body's own period.
/// Prevents oversmoothing if no perturber is found and forces a meaningful
/// floor in pathological cases.
const TAU_MIN_PERIOD_FRAC: f64 = 0.01;

/// Upper clamp for τ_sim, expressed as a multiple of the body's own period.
///
/// The cap exists only as a safety net against pathological perturbers
/// (e.g. a sibling whose period is many orders of magnitude longer than
/// the body's). It must NOT engage in normal solar-system regimes:
/// inner planets (P ≤ a few yr) need τ ≈ 18 yr to attenuate Jupiter's
/// 12-yr indirect forcing, so a tight cap (e.g. 3·P_self) under-filters
/// Earth/Venus/Mars and leaves visible jitter on the orbit overlay.
///
/// 50× P_self gives Earth τ_max = 50 yr — comfortably above the 18 yr
/// the weighted-mean wants — while still bounding pathological cases.
/// Real apsidal precession in the solar system is ~1°/century at most,
/// so even τ = 50 yr blurs only ~0.5° of secular signal: invisible at
/// interactive playback timescales.
const TAU_MAX_PERIOD_FRAC: f64 = 50.0;

/// Per-body smoothing state. Stored in the cache and updated each frame.
///
/// Identity is composed of `primary_key` and `gm`: any change (hierarchy
/// flip, mass change from a merge, GM change from constant retune)
/// invalidates the entry and resets the EMA from the new instantaneous
/// invariants.
///
/// `last_t_sim` records the simulation time of the last EMA update for
/// this body. Computing α from the per-body Δt rather than the per-frame
/// Δt keeps smoothing correct when overlay membership is intermittent
/// (e.g., a body falls out of top-N for several frames, then re-enters).
#[derive(Debug, Clone)]
struct SmoothState {
    /// Stable identity of the primary body. Cache invalidates on change.
    primary_key: String,
    /// Gravitational parameter at last update. Cache invalidates on change.
    gm: f64,
    /// Simulation time at the previous update for this body.
    last_t_sim: f64,
    /// Smoothed specific orbital energy.
    energy: f64,
    /// Smoothed specific angular momentum vector (3D).
    h_vec: apsis::math::Vec3,
    /// Smoothed eccentricity (Laplace–Runge–Lenz) vector (3D).
    e_vec: apsis::math::Vec3,
    /// Latest body-relative position. Not smoothed — captures the most
    /// recent observation so the [`elements_anchored_to_body`] fallback
    /// to [`elements_from_invariants`] still recovers anomalies on the
    /// rare unbound / near-circular paths.
    r_rel: apsis::math::Vec3,
    /// Latest body-relative velocity. Not smoothed; same rationale as
    /// `r_rel`.
    v_rel: apsis::math::Vec3,
}

impl SmoothState {
    fn from_invariants(inv: &OrbitInvariants, primary_key: String, gm: f64, t_sim: f64) -> Self {
        Self {
            primary_key,
            gm,
            last_t_sim: t_sim,
            energy: inv.energy,
            h_vec: inv.h_vec,
            e_vec: inv.e_vec,
            r_rel: inv.r_rel,
            v_rel: inv.v_rel,
        }
    }

    fn invariants(&self) -> OrbitInvariants {
        OrbitInvariants {
            energy: self.energy,
            h_vec: self.h_vec,
            e_vec: self.e_vec,
            r_rel: self.r_rel,
            v_rel: self.v_rel,
        }
    }
}

#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Stateful EMA cache for osculating orbital elements.
///
/// Identifies bodies by stable name (the simulation's authored name), so
/// `swap_remove` shuffles in the body list don't corrupt the cache.
/// Removed bodies are pruned on demand via [`Self::prune`].
pub struct OrbitSmoother {
    cache: HashMap<String, SmoothState>,
    /// Smoothing strength multiplier. K ∈ [0.5, 3.0] is the practical range;
    /// see [`DEFAULT_K`] for the default rationale.
    k: f64,
}

impl OrbitSmoother {
    pub fn new() -> Self {
        Self { cache: HashMap::new(), k: DEFAULT_K }
    }

    /// Drops cache entries whose body name no longer exists in the live
    /// simulation. O(N + M) where N = cache size, M = live names.
    pub fn prune(&mut self, live_names: &[String]) {
        if self.cache.is_empty() {
            return;
        }
        let alive: std::collections::HashSet<&str> =
            live_names.iter().map(String::as_str).collect();
        self.cache.retain(|k, _| alive.contains(k.as_str()));
    }

    /// Drops every cache entry. Use when topology changes catastrophically
    /// (system reset, restore from snapshot).
    pub fn clear(&mut self) {
        self.cache.clear();
    }

    /// Returns the smoothed osculating elements for body `idx`, advancing
    /// the cache by the per-body Δt (`t_sim` − last cached update).
    ///
    /// `siblings` should contain every body sharing `primary_idx` as its
    /// hierarchy primary; the smoother filters out `idx` itself. Pass an
    /// empty slice when no siblings exist (pure two-body system) — the
    /// smoother short-circuits in that case since there is no perturbation
    /// to filter.
    ///
    /// `t_sim` is **simulation time**, not wall-clock. Per-body Δt makes
    /// smoothing correct even when a body's overlay is intermittent (e.g.
    /// dropping in and out of top-N between frames): α is computed from
    /// the gap since *its own* last update, not since the previous frame.
    ///
    /// Returns `None` only when the underlying state is degenerate
    /// (coincident bodies, zero GM); callers should skip the overlay.
    pub fn smoothed(
        &mut self,
        bodies: &[Body],
        names: &[String],
        idx: usize,
        primary_idx: usize,
        siblings: &[usize],
        g_factor: f64,
        t_sim: f64,
    ) -> Option<OrbitalElements> {
        let inv = compute_invariants(bodies, idx, primary_idx, g_factor)?;
        let gm = g_factor * (bodies[idx].mass + bodies[primary_idx].mass);

        // Hyperbolic / parabolic: skip smoothing. These are transient
        // states; jitter is irrelevant and EMA across a regime change
        // (bound ↔ unbound) would lag visibly. Drop any cached entry so a
        // future re-capture in bound state starts cleanly.
        if inv.energy >= 0.0 {
            self.cache.remove(&names[idx]);
            return Some(elements_from_invariants(&inv, primary_idx, gm));
        }

        // Bound: a > 0, period finite.
        let a_self = -gm / (2.0 * inv.energy);
        let period_self = TAU * (a_self * a_self * a_self / gm).sqrt();

        let tau_sim =
            compute_tau(bodies, idx, primary_idx, siblings, g_factor, period_self, self.k);

        let key = &names[idx];
        let primary_key = &names[primary_idx];

        let invariants = match self.cache.get_mut(key) {
            Some(st) if st.primary_key == *primary_key && (st.gm - gm).abs() < 1e-15 => {
                let dt_sim = (t_sim - st.last_t_sim).max(0.0);
                let alpha = if tau_sim > 0.0 && dt_sim > 0.0 {
                    1.0 - (-dt_sim / tau_sim).exp()
                } else {
                    0.0
                };
                // Lerp each component of `h_vec` and `e_vec`
                // independently. These vectors carry magnitude as a
                // physical quantity (`|h_vec|` is the conserved
                // angular momentum, `|e_vec|` is the eccentricity), so
                // component-wise lerp is the correct EMA: the
                // recovered magnitudes and the recovered angles
                // (`atan2` of components) are well-defined for any
                // intermediate state. Renormalising after lerp would
                // invent magnitude information and break the
                // smoother's contract that the cached state is a
                // faithful weighted average of recent observations.
                //
                // Angular continuity (`π → −π` wrap of the recovered
                // `ω` / `Ω` between frames) is the responsibility of
                // any unwrap layer above this — see the
                // `physics::orbital` module-level docs.
                st.energy = lerp(st.energy, inv.energy, alpha);
                st.h_vec = apsis::math::Vec3::new(
                    lerp(st.h_vec.x, inv.h_vec.x, alpha),
                    lerp(st.h_vec.y, inv.h_vec.y, alpha),
                    lerp(st.h_vec.z, inv.h_vec.z, alpha),
                );
                st.e_vec = apsis::math::Vec3::new(
                    lerp(st.e_vec.x, inv.e_vec.x, alpha),
                    lerp(st.e_vec.y, inv.e_vec.y, alpha),
                    lerp(st.e_vec.z, inv.e_vec.z, alpha),
                );
                // r/v are snapshots, not EMA-smoothed: they change every
                // frame and the smoother only consumes them in fallback
                // paths that need a current observation.
                st.r_rel = inv.r_rel;
                st.v_rel = inv.v_rel;
                st.last_t_sim = t_sim;
                st.invariants()
            },
            _ => {
                // Cold start, primary changed, or GM changed: snap to
                // current invariants without lerp.
                let st = SmoothState::from_invariants(&inv, primary_key.clone(), gm, t_sim);
                let snapshot = st.invariants();
                self.cache.insert(key.clone(), st);
                snapshot
            },
        };

        // Anchor the displayed ellipse through the body's current
        // position: smoothed (a, e) for shape, ω recomputed from the
        // instantaneous state so the orbit polyline contains the body
        // exactly. Without this, a multi-body system shows the time-
        // averaged ellipse, which is offset from where the body
        // actually is by ~1% of `a` — visually wrong.
        Some(elements_anchored_to_body(
            &invariants,
            primary_idx,
            gm,
            &bodies[idx],
            &bodies[primary_idx],
        ))
    }
}

impl Default for OrbitSmoother {
    fn default() -> Self {
        Self::new()
    }
}

/// Computes τ_sim (simulation time) for the EMA filter on body `idx`.
///
/// Algorithm:
///   1. For each sibling j ≠ idx: score = |Δa_j|, the differential
///      acceleration j exerts in the primary's frame.
///   2. period_j = orbital period of j around its (shared) primary.
///   3. τ = K · Σ_j (score_j · period_j) / Σ_j score_j.
///   4. Clamp to [TAU_MIN_FRAC · P_self, TAU_MAX_FRAC · P_self].
///
/// When no sibling has finite period (e.g. all unbound) or the sibling
/// list is empty, τ falls back to TAU_MIN_FRAC · P_self — effectively no
/// smoothing, which is correct for a true two-body system.
fn compute_tau(
    bodies: &[Body],
    idx: usize,
    primary_idx: usize,
    siblings: &[usize],
    g_factor: f64,
    period_self: f64,
    k: f64,
) -> f64 {
    let tau_min = TAU_MIN_PERIOD_FRAC * period_self;
    let tau_max = TAU_MAX_PERIOD_FRAC * period_self;

    if siblings.is_empty() || !period_self.is_finite() {
        return tau_min;
    }

    let bp = &bodies[primary_idx];
    let bi = &bodies[idx];
    let r_p = (bp.x, bp.y);
    let r_i = (bi.x, bi.y);

    let mut total_score = 0.0_f64;
    let mut weighted_period_sum = 0.0_f64;

    for &j in siblings {
        if j == idx || j == primary_idx {
            continue;
        }
        let bj = &bodies[j];
        let r_j = (bj.x, bj.y);

        let dx_pj = r_p.0 - r_j.0;
        let dy_pj = r_p.1 - r_j.1;
        let r_pj = (dx_pj * dx_pj + dy_pj * dy_pj).sqrt();
        let dx_ij = r_i.0 - r_j.0;
        let dy_ij = r_i.1 - r_j.1;
        let r_ij = (dx_ij * dx_ij + dy_ij * dy_ij).sqrt();
        if r_pj < 1e-15 || r_ij < 1e-15 {
            continue;
        }
        let inv_pj3 = 1.0 / (r_pj * r_pj * r_pj);
        let inv_ij3 = 1.0 / (r_ij * r_ij * r_ij);
        let dax = dx_pj * inv_pj3 - dx_ij * inv_ij3;
        let day = dy_pj * inv_pj3 - dy_ij * inv_ij3;
        let score = g_factor * bj.mass * (dax * dax + day * day).sqrt();
        if !(score > 0.0) {
            continue;
        }

        // Period of j around the shared primary. Skip unbound siblings —
        // they don't contribute a meaningful period to the weighted mean.
        let inv_j = match compute_invariants(bodies, j, primary_idx, g_factor) {
            Some(inv) if inv.energy < 0.0 => inv,
            _ => continue,
        };
        let gm_j = g_factor * (bj.mass + bp.mass);
        let a_j = -gm_j / (2.0 * inv_j.energy);
        let period_j = TAU * (a_j * a_j * a_j / gm_j).sqrt();
        if !period_j.is_finite() {
            continue;
        }

        total_score += score;
        weighted_period_sum += score * period_j;
    }

    if total_score <= 0.0 {
        return tau_min;
    }

    let weighted_period = weighted_period_sum / total_score;
    (k * weighted_period).clamp(tau_min, tau_max)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        Body::rocky(mass).at(x, y).with_velocity(vx, vy)
    }

    fn circular_velocity(m: f64, r: f64) -> f64 {
        (m / r).sqrt()
    }

    /// Two-body circular system: smoother must converge exactly because
    /// (ε, ex, ey, h) are constant for an unperturbed Kepler orbit, so EMA
    /// on a constant signal is identity for any α.
    #[test]
    fn two_body_kepler_no_drift() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, v_c, 1e-10)];
        let names = vec!["sun".to_string(), "earth".to_string()];

        let mut s = OrbitSmoother::new();
        let el0 = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0).unwrap();

        // Repeated calls on identical state must be byte-identical.
        for k in 1..=100 {
            let t = k as f64 * 1e-3;
            let el = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, t).unwrap();
            assert_eq!(el.a, el0.a);
            assert_eq!(el.e, el0.e);
        }
    }

    /// EMA on a step input must converge to the new value at the rate set
    /// by α.
    #[test]
    fn ema_converges_to_new_state() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let mut bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, v_c, 1e-10)];
        let names = vec!["sun".to_string(), "earth".to_string()];

        let mut s = OrbitSmoother::new();
        // Cold-start at circular state.
        let el0 = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0).unwrap();
        assert!(el0.e < 1e-6);

        // Perturb to e ≈ 0.3 by boosting tangential velocity. With no
        // siblings τ = TAU_MIN · P_self, so each step has α very close to
        // 1; convergence is fast.
        bodies[1].vy *= 1.15;
        for k in 1..=200 {
            s.smoothed(&bodies, &names, 1, 0, &[], 1.0, k as f64);
        }
        let el_final = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 1000.0).unwrap();
        assert!(el_final.e > 0.05, "EMA should have converged to perturbed e, got {}", el_final.e);
    }

    /// Hyperbolic orbits skip the cache entirely and return raw elements.
    #[test]
    fn hyperbolic_skips_cache() {
        let m = 1e6;
        let r_peri = 10.0;
        let v_peri: f64 = 1.5 * (2.0_f64 * m / r_peri).sqrt();
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r_peri, 0.0, 0.0, v_peri, 1e-10)];
        let names = vec!["sun".to_string(), "comet".to_string()];

        let mut s = OrbitSmoother::new();
        let _ = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0);
        assert!(!s.cache.contains_key("comet"), "hyperbolic state must not enter cache");
    }

    /// Changing the primary key must reset the cache for that body —
    /// observable by the recorded `primary_key` and by the fact that GM
    /// changes when the primary's mass differs.
    #[test]
    fn primary_change_resets_cache() {
        let m_a = 1e6;
        let m_b = 4e6; // different mass → different GM → cache invalidates
        let r = 10.0;
        let v_c_a = circular_velocity(m_a, r);
        let bodies = vec![
            body(-r, 0.0, 0.0, 0.0, m_a),
            body(2.0 * r, 0.0, 0.0, 0.0, m_b),
            // Satellite at origin moving +y: bound to either anchor.
            body(0.0, 0.0, 0.0, v_c_a, 1e-10),
        ];
        let names = vec!["a".to_string(), "b".to_string(), "p".to_string()];

        let mut s = OrbitSmoother::new();
        let _ = s.smoothed(&bodies, &names, 2, 0, &[], 1.0, 0.0);
        let st_a = s.cache.get("p").expect("first call must cache");
        assert_eq!(st_a.primary_key, "a");
        let gm_a = st_a.gm;

        // Switch primary → cache slot under "p" must rebind, with the new
        // primary key and a different GM (because m_b ≠ m_a).
        let _ = s.smoothed(&bodies, &names, 2, 1, &[], 1.0, 0.1);
        let st_b = s.cache.get("p").expect("second call must re-cache under new primary");
        assert_eq!(st_b.primary_key, "b");
        assert_ne!(st_b.gm, gm_a);
    }

    /// `prune` removes entries whose name no longer exists.
    #[test]
    fn prune_drops_dead_bodies() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, v_c, 1e-10)];
        let names = vec!["sun".to_string(), "earth".to_string()];
        let mut s = OrbitSmoother::new();
        let _ = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0);
        assert!(s.cache.contains_key("earth"));

        s.prune(&[String::from("sun")]);
        assert!(!s.cache.contains_key("earth"));
    }

    /// After smoothing + anchoring, the body's actual position must lie
    /// on the displayed ellipse to within numerical tolerance. This is
    /// the contract that motivates the anchor refactor: under EMA-only,
    /// the body would sit ~1% of `a` off the displayed orbit.
    fn body_lies_on_orbit(el: &OrbitalElements, body: &Body, primary: &Body) -> f64 {
        // Geometric check: at the body's azimuth from focus, the orbit's
        // predicted distance must equal the body's actual distance.
        let rx = body.x - primary.x;
        let ry = body.y - primary.y;
        let r = (rx * rx + ry * ry).sqrt();
        let theta = ry.atan2(rx);
        let nu = theta - el.omega;
        let p = el.a * (1.0 - el.e * el.e);
        let r_orbit = p / (1.0 + el.e * nu.cos());
        (r - r_orbit).abs() / r.max(1e-15)
    }

    /// Prograde Kepler orbit: body must sit exactly on the anchored
    /// orbit even on the very first frame (cold start uses raw
    /// invariants, anchor still applies).
    #[test]
    fn anchored_orbit_contains_body_prograde() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        // Eccentric: tangential velocity 1.2× circular → e ≈ 0.44
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, 1.2 * v_c, 1e-10)];
        let names = vec!["sun".to_string(), "p".to_string()];
        let mut s = OrbitSmoother::new();
        let el = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0).unwrap();
        let err = body_lies_on_orbit(&el, &bodies[1], &bodies[0]);
        assert!(err < 1e-9, "prograde body off orbit by {}", err);
    }

    /// Retrograde Kepler orbit (h < 0): the sign(h·(r·v)) formulation
    /// must place the body on the same side of the orbit as it actually
    /// orbits. With the older sign(r·v) formula this would render the
    /// orbit flipped.
    #[test]
    fn anchored_orbit_contains_body_retrograde() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let bodies = vec![
            body(0.0, 0.0, 0.0, 0.0, m),
            // tangential velocity inverted → retrograde, eccentric
            body(r, 0.0, 0.0, -1.2 * v_c, 1e-10),
        ];
        let names = vec!["sun".to_string(), "p".to_string()];
        let mut s = OrbitSmoother::new();
        let el = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0).unwrap();
        let err = body_lies_on_orbit(&el, &bodies[1], &bodies[0]);
        assert!(err < 1e-9, "retrograde body off orbit by {}", err);
    }

    /// After many EMA updates with perturbed state, the displayed orbit
    /// must still contain the body. Demonstrates that anchoring is
    /// independent of how stale the smoothed shape is.
    #[test]
    fn anchored_orbit_contains_body_after_smoothing() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, 1.05 * v_c, 1e-10)];
        let names = vec!["sun".to_string(), "p".to_string()];
        let mut s = OrbitSmoother::new();
        // Run many updates with the same state — smoother converges,
        // then anchor must still place body on orbit exactly.
        for k in 0..50 {
            s.smoothed(&bodies, &names, 1, 0, &[], 1.0, k as f64 * 0.01);
        }
        let el = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.5).unwrap();
        let err = body_lies_on_orbit(&el, &bodies[1], &bodies[0]);
        assert!(err < 1e-9, "post-EMA body off orbit by {}", err);
    }

    /// Per-body Δt: a body whose overlay drops out for many frames and
    /// then re-enters must not be smoothed with a stale tiny α — the
    /// per-body last_t_sim correctly accumulates the elapsed sim time.
    #[test]
    fn intermittent_overlay_uses_per_body_dt() {
        let m = 1e6;
        let r = 10.0;
        let v_c = circular_velocity(m, r);
        let mut bodies = vec![body(0.0, 0.0, 0.0, 0.0, m), body(r, 0.0, 0.0, v_c, 1e-10)];
        let names = vec!["sun".to_string(), "earth".to_string()];

        let mut s = OrbitSmoother::new();
        // First update at t=0 cold-starts at circular state.
        let _ = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 0.0).unwrap();

        // Body changes drastically while overlay is *not* requested for
        // many sim seconds (no smoother calls). When overlay returns,
        // a single call with a long Δt should yield α ≈ 1 (essentially
        // snap to current state), not α ≈ small (frame-step assumption).
        bodies[1].vy *= 1.3;
        let el = s.smoothed(&bodies, &names, 1, 0, &[], 1.0, 1e6).unwrap();
        // After a Δt much larger than τ_min, EMA should have effectively
        // collapsed to the current invariants.
        let raw_inv = compute_invariants(&bodies, 1, 0, 1.0).unwrap();
        let gm = 1.0 * (bodies[1].mass + bodies[0].mass);
        let raw = elements_from_invariants(&raw_inv, 0, gm);
        let rel_err = (el.e - raw.e).abs() / raw.e.max(1e-9);
        assert!(rel_err < 1e-3, "long-Δt α should snap to raw, got rel_err {}", rel_err);
    }

    /// τ must lie within [TAU_MIN, TAU_MAX] · P_self regardless of perturber period.
    #[test]
    fn tau_clamp_respects_bounds() {
        let m_central = 1e6;
        let r_inner = 1.0;
        let r_outer = 100.0;
        let v_inner = circular_velocity(m_central, r_inner);
        let v_outer = circular_velocity(m_central, r_outer);
        let bodies = vec![
            body(0.0, 0.0, 0.0, 0.0, m_central),
            body(r_inner, 0.0, 0.0, v_inner, 1.0), // huge "inner perturber"
            body(r_outer, 0.0, 0.0, v_outer, 1e-10), // outer body
        ];
        let inv_outer = compute_invariants(&bodies, 2, 0, 1.0).unwrap();
        let gm = 1.0 * (bodies[2].mass + bodies[0].mass);
        let a_outer = -gm / (2.0 * inv_outer.energy);
        let period_outer = TAU * (a_outer * a_outer * a_outer / gm).sqrt();

        let tau = compute_tau(&bodies, 2, 0, &[1], 1.0, period_outer, DEFAULT_K);
        assert!(tau >= TAU_MIN_PERIOD_FRAC * period_outer - 1e-9);
        assert!(tau <= TAU_MAX_PERIOD_FRAC * period_outer + 1e-9);
    }
}
