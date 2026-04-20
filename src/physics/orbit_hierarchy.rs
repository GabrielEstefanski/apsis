//! Frame-coherent gravitational hierarchy.
//!
//! For each body, picks the single "primary" it orbits around and assigns
//! a level (0 = root primary, 1 = planetary, 2 = satellite, …). A body
//! with no bound primary is classified `Free` — ejected or hyperbolic to
//! every candidate in the system.
//!
//! # Scoring: SOI nesting + specific orbital energy
//!
//! Raw `G·m/r²` ranking flips chaotically near transitions because it
//! depends only on geometry; two candidates at similar distance give
//! near-identical scores. Raw specific orbital energy `ε = v_rel²/2 − μ/r`
//! does better (it includes velocity, so dynamical state pulls candidates
//! apart near saddles), but does not nest correctly: a moon orbiting a
//! planet is typically *more energetically bound* to the distant massive
//! star than to its nearby planet, because `μ_star` dominates `μ_planet`.
//!
//! The right abstraction is the patched-conics **sphere of influence**:
//! a body belongs to the primary whose SOI most tightly contains it. The
//! algorithm is:
//!
//! 1. Compute each body's Laplace SOI radius `r_SOI(j) ≈ r_{j,H} ·
//!    (m_j / m_H)^(2/5)` where `H` is the globally most massive body
//!    (treated as the cosmic frame). For `H` itself, SOI is infinite.
//! 2. For each body `i`, pick the primary `j` with the **smallest** SOI
//!    that still contains `i` (`r_ij < r_SOI(j)`) and to which `i` is
//!    energetically bound (`ε_ij < 0`). Ties broken by minimum `ε`
//!    (most-bound wins).
//! 3. Hysteresis: retain the previous primary unless the new winner is
//!    in a *meaningfully* smaller SOI layer (`SOI_new < SOI_prev × (1 −
//!    margin)`). Default margin 15 %. Near-equal candidates stick to
//!    whatever was chosen last frame; genuine layer changes still flip.
//!
//! # Cadence
//!
//! [`tick`](OrbitHierarchy::tick) is gated: the hierarchy refreshes at
//! most every `cadence` (default 500 ms). Topology changes (body count
//! delta) bypass the gate so newly-added bodies do not render stale data.
//!
//! # Known limitations
//!
//! * Equal-mass binaries: with two exactly-equal primaries the Laplace
//!   formula degenerates (each SOI covers the full separation). The
//!   index-based tiebreak yields a consistent choice but the physical
//!   ideal is "orbit around the barycenter" — a later pass.
//! * The global SOI reference (most massive body) is a single-frame
//!   approximation. For decentralised systems it is less principled but
//!   remains well-defined.

use crate::domain::body::Body;
use std::time::{Duration, Instant};

/// Cap on hierarchy depth. Deeper chains saturate here so a pathological
/// input cannot spin the level walker forever.
pub const MAX_LEVEL: u8 = 32;

/// Default cadence for [`OrbitHierarchy::tick`].
pub const DEFAULT_CADENCE: Duration = Duration::from_millis(500);

/// Default hysteresis margin. New primary must be in an SOI layer at
/// least this fraction smaller than the previous primary's SOI to win.
pub const DEFAULT_HYSTERESIS_MARGIN: f64 = 0.15;

/// Laplace SOI exponent: `r_SOI ≈ r × (m_body / m_ref)^EXPONENT`.
/// 2/5 is the classical patched-conics value.
const SOI_EXPONENT: f64 = 2.0 / 5.0;

/// Classification of a body in the current hierarchy snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitClass {
    /// Body is bound to a primary; `level` is its depth in the tree.
    /// 0 = root primary (no primary above it), 1 = planetary, 2 = satellite,
    /// 3+ = sub-satellite.
    Bound { level: u8 },
    /// No bound primary — `ε ≥ 0` for every candidate, or no candidate's
    /// SOI contains this body. Typical for ejected bodies and hyperbolic
    /// flybys.
    Free,
}

/// Snapshot of the gravitational hierarchy plus state for stable
/// frame-to-frame updates.
#[derive(Debug, Clone)]
pub struct OrbitHierarchy {
    /// `primary[i] = Some(j)` — body `i`'s dominant primary is `j`.
    /// `primary[i] = None`    — body `i` is classified [`OrbitClass::Free`]
    /// or is itself a root.
    primary: Vec<Option<usize>>,
    /// Cached depth in the hierarchy tree. Meaningless when class is Free.
    level: Vec<u8>,
    /// Cached class for O(1) lookup.
    class: Vec<OrbitClass>,
    /// SOI radius of each body under the last-used reference. Reused by
    /// the hysteresis check next tick.
    soi: Vec<f64>,
    /// Specific binding energy of `i` w.r.t. its current primary. Negative
    /// = bound, 0.0 when Free. Reused by tests and top-N scoring.
    bind_energy: Vec<f64>,

    last_computed: Option<Instant>,
    last_body_count: usize,

    /// Recompute cadence. `Duration::ZERO` = every tick.
    pub cadence: Duration,
    /// Fractional margin: new primary's SOI must be less than
    /// `prev_soi × (1 - margin)` to unseat the previous primary.
    pub hysteresis_margin: f64,
}

impl Default for OrbitHierarchy {
    fn default() -> Self {
        Self::new()
    }
}

impl OrbitHierarchy {
    pub fn new() -> Self {
        Self {
            primary: Vec::new(),
            level: Vec::new(),
            class: Vec::new(),
            soi: Vec::new(),
            bind_energy: Vec::new(),
            last_computed: None,
            last_body_count: 0,
            cadence: DEFAULT_CADENCE,
            hysteresis_margin: DEFAULT_HYSTERESIS_MARGIN,
        }
    }

    /// Current primary of body `idx`, or `None` if body is Free / index is
    /// out of range / hierarchy has never been computed.
    pub fn primary(&self, idx: usize) -> Option<usize> {
        self.primary.get(idx).copied().flatten()
    }

    /// Classification of body `idx`.
    pub fn class(&self, idx: usize) -> Option<OrbitClass> {
        self.class.get(idx).copied()
    }

    /// Convenience: level of body `idx`, or `None` if out of range or Free.
    pub fn level(&self, idx: usize) -> Option<u8> {
        match self.class(idx)? {
            OrbitClass::Bound { level } => Some(level),
            OrbitClass::Free => None,
        }
    }

    /// Specific binding energy at last recompute. 0.0 for Free bodies.
    pub fn bind_energy(&self, idx: usize) -> f64 {
        self.bind_energy.get(idx).copied().unwrap_or(0.0)
    }

    /// SOI radius used for body `idx` at last recompute. Infinite for the
    /// globally heaviest body (the SOI reference). 0.0 before first tick.
    pub fn soi(&self, idx: usize) -> f64 {
        self.soi.get(idx).copied().unwrap_or(0.0)
    }

    /// Cadence-gated recompute. Returns `true` if the hierarchy was
    /// refreshed this call. Topology changes (body count delta) force an
    /// immediate refresh.
    pub fn tick(&mut self, bodies: &[Body], g_factor: f64) -> bool {
        let topology_changed = self.last_body_count != bodies.len();
        let first_run = self.last_computed.is_none();
        let due = self
            .last_computed
            .map_or(true, |t| Instant::now().duration_since(t) >= self.cadence);
        if first_run || topology_changed || due {
            self.recompute(bodies, g_factor);
            true
        } else {
            false
        }
    }

    /// Unconditional recompute. Bypasses the cadence gate. Called by
    /// tests, by [`tick`](Self::tick) when gates permit, and by UI code
    /// that wants immediate feedback.
    pub fn recompute(&mut self, bodies: &[Body], g_factor: f64) {
        let n = bodies.len();

        // Snapshot previous state for hysteresis. Cloning is cheap and
        // keeps the decision logic trivially correct (compare against the
        // pre-recompute state, not against partial updates).
        let prev_primary: Vec<Option<usize>> = self.primary.clone();
        let prev_soi: Vec<f64> = self.soi.clone();

        self.primary.clear();
        self.primary.resize(n, None);
        self.bind_energy.clear();
        self.bind_energy.resize(n, 0.0);
        self.level.clear();
        self.level.resize(n, 0);
        self.class.clear();
        self.class.resize(n, OrbitClass::Free);
        self.soi.clear();
        self.soi.resize(n, 0.0);

        // Sparse / trivial system — nothing to classify.
        if n == 0 {
            self.last_computed = Some(Instant::now());
            self.last_body_count = 0;
            return;
        }
        if n == 1 {
            self.soi[0] = f64::INFINITY;
            self.class[0] = OrbitClass::Free;
            self.last_computed = Some(Instant::now());
            self.last_body_count = 1;
            return;
        }

        // ── SOI assignment ──────────────────────────────────────────────
        // Globally heaviest body is the cosmic frame (SOI = ∞). All other
        // SOIs are computed relative to it via Laplace's formula.
        let heaviest = heaviest_index(bodies);
        self.soi[heaviest] = f64::INFINITY;
        let h = &bodies[heaviest];
        let m_ref = h.mass.max(f64::MIN_POSITIVE);
        for j in 0..n {
            if j == heaviest {
                continue;
            }
            let bj = &bodies[j];
            let dx = bj.x - h.x;
            let dy = bj.y - h.y;
            let r = (dx * dx + dy * dy).sqrt();
            if !r.is_finite() || r == 0.0 {
                self.soi[j] = 0.0;
                continue;
            }
            let ratio = (bj.mass / m_ref).max(0.0);
            self.soi[j] = r * ratio.powf(SOI_EXPONENT);
        }

        // ── Primary selection per body ──────────────────────────────────
        for i in 0..n {
            let winner = pick_primary(i, bodies, &self.soi, g_factor);
            let chosen = apply_hysteresis(
                i,
                winner,
                prev_primary.get(i).copied().flatten(),
                &prev_soi,
                bodies,
                &self.soi,
                g_factor,
                self.hysteresis_margin,
            );
            if let Some((p, eps)) = chosen {
                self.primary[i] = Some(p);
                self.bind_energy[i] = eps;
            }
        }

        // ── Break any residual cycles (rare — possible at exact mass
        // ties where SOI formula and index-tiebreak interact) ───────────
        break_cycles(&mut self.primary, bodies);

        // ── Levels ──────────────────────────────────────────────────────
        for i in 0..n {
            if self.primary[i].is_some() {
                self.level[i] = walk_level(i, &self.primary);
            }
        }

        // ── Class finalisation ──────────────────────────────────────────
        let mut has_child = vec![false; n];
        for i in 0..n {
            if let Some(p) = self.primary[i] {
                if p < n {
                    has_child[p] = true;
                }
            }
        }
        for i in 0..n {
            self.class[i] = match self.primary[i] {
                Some(_) => OrbitClass::Bound { level: self.level[i] },
                None if has_child[i] => OrbitClass::Bound { level: 0 },
                None => OrbitClass::Free,
            };
        }

        self.last_computed = Some(Instant::now());
        self.last_body_count = n;
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

fn heaviest_index(bodies: &[Body]) -> usize {
    let mut best = 0;
    let mut best_mass = bodies[0].mass;
    for (i, b) in bodies.iter().enumerate().skip(1) {
        if b.mass > best_mass {
            best_mass = b.mass;
            best = i;
        }
    }
    best
}

/// Specific orbital energy of body `a` as a test particle in `b`'s field.
/// Returns `None` if the pair is coincident or any value is non-finite.
fn binding_energy(a: &Body, b: &Body, g_factor: f64) -> Option<f64> {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let r2 = dx * dx + dy * dy;
    if r2 == 0.0 || !r2.is_finite() {
        return None;
    }
    let r = r2.sqrt();
    let dvx = a.vx - b.vx;
    let dvy = a.vy - b.vy;
    let v2 = dvx * dvx + dvy * dvy;
    let mu = g_factor * b.mass;
    let eps = 0.5 * v2 - mu / r;
    if eps.is_finite() { Some(eps) } else { None }
}

/// Returns the pair-distance between two bodies. None if non-finite.
fn distance(a: &Body, b: &Body) -> Option<f64> {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let r2 = dx * dx + dy * dy;
    if !r2.is_finite() {
        return None;
    }
    Some(r2.sqrt())
}

/// Chooses body `i`'s primary without hysteresis.
///
/// Selection rule:
/// 1. Valid candidates must contain `i` inside their SOI, be bound to `i`
///    (`ε < 0`), and (by mass priority) satisfy `m_j > m_i`, with equal
///    masses broken by index.
/// 2. Among valid candidates, the one with the **smallest SOI** wins
///    (deepest layer). Ties on SOI are broken by smallest `ε`.
///
/// Returns `Some((primary_idx, eps))` or `None` if no candidate is valid.
fn pick_primary(
    i: usize,
    bodies: &[Body],
    soi: &[f64],
    g_factor: f64,
) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64, f64)> = None; // (j, soi_j, eps)
    for j in 0..bodies.len() {
        if j == i {
            continue;
        }
        // Mass priority with index tiebreak — enforces an acyclic ordering.
        let m_i = bodies[i].mass;
        let m_j = bodies[j].mass;
        let more_massive = m_j > m_i || (m_j == m_i && j < i);
        if !more_massive {
            continue;
        }
        let soi_j = soi[j];
        if !(soi_j > 0.0) {
            continue;
        }
        let Some(r) = distance(&bodies[i], &bodies[j]) else {
            continue;
        };
        if r >= soi_j {
            continue; // outside candidate's sphere of influence
        }
        let Some(eps) = binding_energy(&bodies[i], &bodies[j], g_factor) else {
            continue;
        };
        if eps >= 0.0 {
            continue; // unbound w.r.t. this candidate
        }
        match best {
            None => best = Some((j, soi_j, eps)),
            Some((_, cur_soi, cur_eps)) => {
                if soi_j < cur_soi || (soi_j == cur_soi && eps < cur_eps) {
                    best = Some((j, soi_j, eps));
                }
            }
        }
    }
    best.map(|(j, _soi, eps)| (j, eps))
}

/// Applies hysteresis to the frame's primary choice.
///
/// Keeps the previous primary when:
/// * the previous primary is still a valid candidate this frame (inside
///   SOI, bound), and
/// * the new winner's SOI is not strictly smaller than `prev_soi ×
///   (1 − margin)` — i.e. the new candidate has not moved into a
///   meaningfully deeper layer.
///
/// Switches unconditionally when the previous primary has become invalid
/// (for example, the body escaped its SOI).
fn apply_hysteresis(
    i: usize,
    winner: Option<(usize, f64)>,
    prev_primary: Option<usize>,
    prev_soi: &[f64],
    bodies: &[Body],
    soi: &[f64],
    g_factor: f64,
    margin: f64,
) -> Option<(usize, f64)> {
    let (w_idx, w_eps) = winner?;
    let Some(prev_j) = prev_primary else {
        return Some((w_idx, w_eps));
    };
    if prev_j >= bodies.len() || prev_j == i || prev_j == w_idx {
        return Some((w_idx, w_eps));
    }
    // Prev must still be a valid candidate.
    let prev_soi_val = soi.get(prev_j).copied().unwrap_or(0.0);
    if !(prev_soi_val > 0.0) {
        return Some((w_idx, w_eps));
    }
    let Some(r_prev) = distance(&bodies[i], &bodies[prev_j]) else {
        return Some((w_idx, w_eps));
    };
    if r_prev >= prev_soi_val {
        return Some((w_idx, w_eps));
    }
    let Some(eps_prev) = binding_energy(&bodies[i], &bodies[prev_j], g_factor) else {
        return Some((w_idx, w_eps));
    };
    if eps_prev >= 0.0 {
        return Some((w_idx, w_eps));
    }
    // Use the previous frame's SOI for the comparison when available; if
    // the body set changed, fall back to current-frame SOI.
    let reference_soi = prev_soi
        .get(prev_j)
        .copied()
        .filter(|s| *s > 0.0)
        .unwrap_or(prev_soi_val);
    let winner_soi = soi.get(w_idx).copied().unwrap_or(f64::INFINITY);
    // Retain previous unless the winner is in a strictly deeper SOI layer
    // (smaller SOI by at least the margin fraction).
    if winner_soi < reference_soi * (1.0 - margin) {
        Some((w_idx, w_eps))
    } else {
        Some((prev_j, eps_prev))
    }
}

/// Breaks cycles in the primary graph. With the mass-priority tiebreak
/// cycles are nearly impossible, but finite-precision edge cases can
/// still produce them — this pass guarantees the output is a forest.
fn break_cycles(primary: &mut [Option<usize>], bodies: &[Body]) {
    let n = primary.len();
    let mut state = vec![0u8; n]; // 0 unseen, 1 on-path, 2 finalised
    let mut path: Vec<usize> = Vec::with_capacity(16);

    for start in 0..n {
        if state[start] != 0 {
            continue;
        }
        path.clear();
        let mut cur = start;
        loop {
            match state[cur] {
                1 => {
                    let cycle_start = path.iter().position(|&x| x == cur).unwrap();
                    let cycle = &path[cycle_start..];
                    let heaviest = *cycle
                        .iter()
                        .max_by(|&&a, &&b| {
                            bodies[a]
                                .mass
                                .partial_cmp(&bodies[b].mass)
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap();
                    primary[heaviest] = None;
                    for &p in cycle {
                        state[p] = 2;
                    }
                    break;
                }
                2 => break,
                _ => {}
            }
            state[cur] = 1;
            path.push(cur);
            match primary[cur] {
                Some(next) if next < n => cur = next,
                _ => break,
            }
        }
        for &p in &path {
            state[p] = 2;
        }
    }
}

/// Walks the primary chain from `start` up to a root. Caller must have
/// ensured the graph is acyclic. Saturates at [`MAX_LEVEL`].
fn walk_level(start: usize, primary: &[Option<usize>]) -> u8 {
    let mut cur = start;
    let mut depth: u32 = 0;
    while let Some(next) = primary[cur] {
        cur = next;
        depth += 1;
        if depth >= MAX_LEVEL as u32 {
            return MAX_LEVEL;
        }
    }
    depth as u8
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::body::Body;
    use crate::domain::materials::Material;

    fn body(x: f64, y: f64, vx: f64, vy: f64, mass: f64) -> Body {
        let mut b = Body::new(x, y, vx, vy, mass, Material::Asteroid);
        b.physical_radius = 0.0;
        b
    }

    fn circular(r: f64, mu: f64) -> f64 {
        (mu / r).sqrt()
    }

    #[test]
    fn star_planet_moon_levels_0_1_2() {
        let g = 1.0;
        let m_star = 1.0e6;
        let m_planet = 1.0e3;
        let m_moon = 1.0;

        let r_p = 100.0;
        let v_p = circular(r_p, g * m_star);
        let r_m = 1.0;
        let v_m_rel = circular(r_m, g * m_planet);

        let star = body(0.0, 0.0, 0.0, 0.0, m_star);
        let planet = body(r_p, 0.0, 0.0, v_p, m_planet);
        let moon = body(r_p + r_m, 0.0, 0.0, v_p + v_m_rel, m_moon);
        let bodies = vec![star, planet, moon];

        let mut h = OrbitHierarchy::new();
        h.recompute(&bodies, g);

        assert_eq!(
            h.class(0),
            Some(OrbitClass::Bound { level: 0 }),
            "star is the root primary",
        );
        assert_eq!(
            h.class(1),
            Some(OrbitClass::Bound { level: 1 }),
            "planet orbits the star (level 1)",
        );
        assert_eq!(
            h.class(2),
            Some(OrbitClass::Bound { level: 2 }),
            "moon orbits the planet (level 2) — the whole point of SOI nesting",
        );
        assert_eq!(h.primary(1), Some(0));
        assert_eq!(h.primary(2), Some(1));
    }

    #[test]
    fn ejected_body_is_classified_free() {
        let g = 1.0;
        let m_star = 1.0e6;
        let star = body(0.0, 0.0, 0.0, 0.0, m_star);
        // Speed far above escape velocity from anything plausible.
        let rogue = body(10.0, 0.0, 1000.0, 0.0, 1.0);

        let mut h = OrbitHierarchy::new();
        h.recompute(&vec![star, rogue], g);

        assert_eq!(h.class(1), Some(OrbitClass::Free));
        assert_eq!(h.primary(1), None);
    }

    #[test]
    fn hysteresis_retains_primary_across_marginal_soi_flip() {
        // Moon between two planets with near-identical SOIs. Without
        // hysteresis, a tiny shift flips which SOI is smaller; with
        // hysteresis, the previous primary sticks.
        let g = 1.0;
        let m_star = 1.0e6;
        let m_planet = 1.0e3;
        // Laplace SOI ≈ r × (1e-3)^0.4 ≈ r × 0.0631
        // Place planets at nearly-equal distances so SOIs differ by <1%.

        // Frame 1 — A slightly closer, so SOI_A < SOI_B. Moon is in both.
        let star1 = body(0.0, 0.0, 0.0, 0.0, m_star);
        let pa1 = body(100.0, 0.0, 0.0, 0.0, m_planet);
        let pb1 = body(100.0, 10.0, 0.0, 0.0, m_planet);
        let moon1 = body(100.0, 5.0, 0.0, 0.0, 1.0);
        let bodies1 = vec![star1, pa1, pb1, moon1];

        let mut h = OrbitHierarchy::new();
        h.recompute(&bodies1, g);
        let initial = h.primary(3).expect("moon should have a primary");

        // Frame 2 — planets swap their SOI ordering by a fraction of a
        // percent. Without hysteresis, primary would flip; with 15 %
        // margin, the previous pick must stand.
        let star2 = body(0.0, 0.0, 0.0, 0.0, m_star);
        let pa2 = body(100.01, 0.0, 0.0, 0.0, m_planet);
        let pb2 = body(100.0, 10.0, 0.0, 0.0, m_planet);
        let moon2 = body(100.0, 5.0, 0.0, 0.0, 1.0);
        h.recompute(&vec![star2, pa2, pb2, moon2], g);

        assert_eq!(
            h.primary(3),
            Some(initial),
            "hysteresis must preserve the previous primary on a sub-margin flip",
        );
    }

    #[test]
    fn hysteresis_yields_when_new_candidate_is_dramatically_deeper() {
        // Same structure, but in frame 2 one planet's SOI is now much
        // smaller (moved close to the star) — a genuine layer change.
        let g = 1.0;
        let m_star = 1.0e6;
        let m_planet = 1.0e3;

        let star1 = body(0.0, 0.0, 0.0, 0.0, m_star);
        let pa1 = body(100.0, 0.0, 0.0, 0.0, m_planet);
        let pb1 = body(100.0, 10.0, 0.0, 0.0, m_planet);
        let moon1 = body(100.0, 5.0, 0.0, 0.0, 1.0);

        let mut h = OrbitHierarchy::new();
        h.recompute(&vec![star1, pa1, pb1, moon1], g);
        let initial = h.primary(3).unwrap();

        // Drop pb2 to a 10× smaller distance → SOI_B shrinks ~10× → new
        // winner. Keep the moon adjacent to B so it is still inside.
        let star2 = body(0.0, 0.0, 0.0, 0.0, m_star);
        let pa2 = body(100.0, 0.0, 0.0, 0.0, m_planet);
        let pb2 = body(10.0, 1.0, 0.0, 0.0, m_planet);
        let moon2 = body(10.0, 1.1, 0.0, 0.0, 1.0);
        h.recompute(&vec![star2, pa2, pb2, moon2], g);

        assert_ne!(
            h.primary(3),
            Some(initial),
            "a dramatic SOI-layer change must override hysteresis",
        );
    }

    #[test]
    fn topology_change_bypasses_cadence_gate() {
        let g = 1.0;
        let star = body(0.0, 0.0, 0.0, 0.0, 1.0e6);
        let planet = body(100.0, 0.0, 0.0, circular(100.0, 1.0e6), 1.0e3);

        let mut h = OrbitHierarchy::new();
        h.cadence = Duration::from_secs(3600); // prevent normal ticks

        assert!(h.tick(&vec![star, planet], g), "first tick always recomputes");
        assert!(
            !h.tick(&vec![star, planet], g),
            "cadence should gate unchanged topology",
        );

        let moon = body(101.0, 0.0, 0.0, circular(100.0, 1.0e6) + 1.0, 1.0);
        assert!(
            h.tick(&vec![star, planet, moon], g),
            "a body-count change must force an immediate recompute",
        );
        assert_eq!(h.class(2), Some(OrbitClass::Bound { level: 2 }));
    }

    #[test]
    fn binary_cycle_is_broken_heaviest_becomes_root() {
        let g = 1.0;
        let m_a = 1.0e6;
        let m_b = 1.0e3;
        let r = 10.0;
        let a = body(-r, 0.0, 0.0, 0.1, m_a);
        let b = body(r, 0.0, 0.0, -100.0, m_b);

        let mut h = OrbitHierarchy::new();
        h.recompute(&vec![a, b], g);

        assert_eq!(h.class(0), Some(OrbitClass::Bound { level: 0 }));
        assert_eq!(h.class(1), Some(OrbitClass::Bound { level: 1 }));
        assert_eq!(h.primary(1), Some(0));
        assert_eq!(h.primary(0), None);
    }

    #[test]
    fn empty_and_single_body_are_safe() {
        let mut h = OrbitHierarchy::new();
        h.recompute(&[], 1.0);
        assert_eq!(h.primary(0), None);

        let lone = body(0.0, 0.0, 0.0, 0.0, 1.0);
        h.recompute(&vec![lone], 1.0);
        assert_eq!(h.class(0), Some(OrbitClass::Free));
    }
}
