//! SoA layout perf harness — measures the production SoA force-eval path
//! against a faithful in-file shadow of the AoS path that lived in
//! `tree.rs` / `engine.rs` at commit `03cb1a6` (pre-SoA refactor).
//!
//! Lab notebook: `docs/experiments/2026-05-10-soa-layout.md`.
//!
//! Opt-in:
//!
//! ```text
//! cargo test --release -p apsis perf_soa_aos_vs_soa -- --ignored --nocapture
//! ```
//!
//! ## Why a shadow path inside the harness
//!
//! The Tier 2 hypothesis (`t_walk_AoS / t_walk_SoA ∈ [1.20, 1.50]` at
//! N = 10⁴) requires both wall-time numerator and denominator measured on
//! the same hardware, same compiler, same body distribution, same warm-cache
//! state. Production removed the AoS path entirely in commit 5 (`253abe9`),
//! so the harness reintroduces it locally as the [`aos_baseline`] sub-module
//! — a verbatim copy of the AoS implementation read out of git history at
//! commit `03cb1a6`. Same `Node` struct, same `Octree` shape, same
//! `bh_eval_body` arithmetic; only the data source differs (`&[Body]`
//! direct field reads instead of `&BodyArrays` indexed reads).
//!
//! The shadow is contained in this file — production `tree.rs` /
//! `engine.rs` stay SoA-only. When the SoA experiment closes in the
//! §Decision commit, this whole file goes (per the perf-2 / perf-4
//! closure pattern); the AoS shadow is removed with it, leaving the
//! lab notebook as the sole record.
//!
//! ## Frozen variables
//!
//! * Seeds: `0x6F637472`, `0x71756164`, `0x6D6F7274` (perf 2×2 / engine
//!   ceiling / MAC canonical set; cross-experiment comparability)
//! * Body distribution: sphere log-normal mass
//! * `θ`: `0.5` (production canonical)
//! * Multipole order: quadrupole always-on (perf 2×2 §Decision)
//! * MAC: classical `s/d < θ` (MAC §Decision)
//! * N grid: `{1_000, 5_000, 10_000}` (matches MAC notebook)
//! * Warmup runs (discarded): 3 per (cell, seed)
//! * Measured runs: 5 per (cell, seed); within-seed median is the cell value
//!
//! CSV output: `target/perf-soa/profile.csv`. Per-row schema in
//! [`write_header`]; one row per (n, seed, path).

#![allow(dead_code)]

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use crate::domain::body::Body;
use crate::domain::body_arrays::BodyArrays;
use crate::math::Vec3;
use crate::physics::gravity::{
    BarnesHutEngine, DEFAULT_LEAF, G, Kernel, NO_CHILD, Node, PlummerKernel, WalkCounters,
    pair_eps2,
};

const SEEDS: [u64; 3] = [0x6F637472, 0x71756164, 0x6D6F7274];
const N_VALUES: [usize; 3] = [1_000, 5_000, 10_000];
const THETA: f64 = 0.5;
const WARMUP_RUNS: usize = 3;
const MEASURED_RUNS: usize = 5;

// ── AoS baseline shadow — verbatim from commit 03cb1a6 ──────────────────────── //
//
// Faithful copy of the AoS Octree + walk + exact-eval that lived in
// `tree.rs` / `engine.rs` before the SoA refactor (PR-perf-5 commits 4-5).
// The shadow reuses production's `Node` struct so the per-node memory
// layout is bit-identical to what the AoS path would have observed in
// production; only the data-source code differs from production's SoA.
//
// Contained in this sub-module so production stays SoA-only and the
// shadow goes away when the harness file is removed at experiment closure.

mod aos_baseline {
    use super::{Body, DEFAULT_LEAF, G, Kernel, NO_CHILD, Node, Vec3, WalkCounters, pair_eps2};
    use rayon::prelude::*;

    const TREE_PAD: f64 = 1e-2;

    pub(super) struct Octree {
        nodes: Vec<Node<DEFAULT_LEAF>>,
        max_depth: usize,
    }

    impl Octree {
        pub(super) fn new(max_depth: usize) -> Self {
            Self { nodes: Vec::new(), max_depth }
        }

        pub(super) fn nodes(&self) -> &[Node<DEFAULT_LEAF>] {
            &self.nodes
        }

        pub(super) fn build(&mut self, bodies: &[Body]) {
            self.nodes.clear();

            if bodies.is_empty() {
                return;
            }

            let mut min_x = bodies[0].pos_x;
            let mut max_x = bodies[0].pos_x;
            let mut min_y = bodies[0].pos_y;
            let mut max_y = bodies[0].pos_y;
            let mut min_z = bodies[0].pos_z;
            let mut max_z = bodies[0].pos_z;
            for b in &bodies[1..] {
                min_x = min_x.min(b.pos_x);
                max_x = max_x.max(b.pos_x);
                min_y = min_y.min(b.pos_y);
                max_y = max_y.max(b.pos_y);
                min_z = min_z.min(b.pos_z);
                max_z = max_z.max(b.pos_z);
            }
            let cx = 0.5 * (min_x + max_x);
            let cy = 0.5 * (min_y + max_y);
            let cz = 0.5 * (min_z + max_z);
            let extent = (max_x - min_x).max(max_y - min_y).max(max_z - min_z);
            let mut half = 0.5 * extent;
            half = if half <= 0.0 { TREE_PAD } else { half * 1.0001 + TREE_PAD };

            self.nodes.push(Node::new(cx, cy, cz, half));

            for i in 0..bodies.len() {
                self.insert(0, i, bodies, 0);
            }

            self.aggregate_mass(0, bodies);
            self.aggregate_quadrupole(0, bodies);
        }

        fn insert(
            &mut self,
            mut node_idx: usize,
            body_idx: usize,
            bodies: &[Body],
            mut depth: usize,
        ) {
            loop {
                if depth > self.max_depth {
                    let node = &mut self.nodes[node_idx];
                    if (node.body_len as usize) < DEFAULT_LEAF {
                        node.bodies[node.body_len as usize] = body_idx as u32;
                        node.body_len += 1;
                    }
                    return;
                }

                if self.nodes[node_idx].is_leaf() {
                    let len = self.nodes[node_idx].body_len as usize;

                    if len < DEFAULT_LEAF || depth == self.max_depth {
                        if (self.nodes[node_idx].body_len as usize) < DEFAULT_LEAF {
                            self.nodes[node_idx].bodies[len] = body_idx as u32;
                            self.nodes[node_idx].body_len += 1;
                        }
                        return;
                    }

                    let existing_len = self.nodes[node_idx].body_len as usize;
                    let existing = self.nodes[node_idx].bodies;
                    self.nodes[node_idx].body_len = 0;

                    self.subdivide(node_idx);

                    for &bi in &existing[..existing_len] {
                        let child = self.child_octant(node_idx, bi as usize, bodies);
                        self.insert(child, bi as usize, bodies, depth + 1);
                    }
                }

                node_idx = self.child_octant(node_idx, body_idx, bodies);
                depth += 1;
            }
        }

        fn subdivide(&mut self, idx: usize) {
            let (cx, cy, cz, half) = {
                let n = &self.nodes[idx];
                (n.cx, n.cy, n.cz, n.half)
            };
            let h = half * 0.5;
            let mut children = [NO_CHILD; 8];
            for octant in 0..8 {
                let sx = if octant & 0b001 != 0 { h } else { -h };
                let sy = if octant & 0b010 != 0 { h } else { -h };
                let sz = if octant & 0b100 != 0 { h } else { -h };
                children[octant] = self.push_node(cx + sx, cy + sy, cz + sz, h) as u32;
            }
            self.nodes[idx].children = children;
        }

        fn push_node(&mut self, cx: f64, cy: f64, cz: f64, half: f64) -> usize {
            let idx = self.nodes.len();
            self.nodes.push(Node::new(cx, cy, cz, half));
            idx
        }

        fn child_octant(&self, node_idx: usize, body_idx: usize, bodies: &[Body]) -> usize {
            let n = &self.nodes[node_idx];
            let b = bodies[body_idx];
            let octant = ((b.pos_z >= n.cz) as usize) << 2
                | ((b.pos_y >= n.cy) as usize) << 1
                | (b.pos_x >= n.cx) as usize;
            self.nodes[node_idx].children[octant] as usize
        }

        fn aggregate_mass(&mut self, idx: usize, bodies: &[Body]) -> (f64, f64, f64, f64) {
            if self.nodes[idx].is_leaf() {
                let len = self.nodes[idx].body_len as usize;
                let mut m = 0.0_f64;
                let mut wx = 0.0_f64;
                let mut wy = 0.0_f64;
                let mut wz = 0.0_f64;

                for k in 0..len {
                    let b = bodies[self.nodes[idx].bodies[k] as usize];
                    m += b.mass;
                    wx += b.mass * b.pos_x;
                    wy += b.mass * b.pos_y;
                    wz += b.mass * b.pos_z;
                }

                self.nodes[idx].body_count = len as u32;
                self.nodes[idx].mass = m;

                if m > 0.0 {
                    self.nodes[idx].com_x = wx / m;
                    self.nodes[idx].com_y = wy / m;
                    self.nodes[idx].com_z = wz / m;
                    return (
                        m,
                        self.nodes[idx].com_x,
                        self.nodes[idx].com_y,
                        self.nodes[idx].com_z,
                    );
                }
                return (0.0, 0.0, 0.0, 0.0);
            }

            let children = self.nodes[idx].children;
            let mut m = 0.0_f64;
            let mut wx = 0.0_f64;
            let mut wy = 0.0_f64;
            let mut wz = 0.0_f64;
            let mut cnt = 0u32;

            for &c in &children {
                if c == NO_CHILD {
                    continue;
                }
                let (cm, cx, cy, cz) = self.aggregate_mass(c as usize, bodies);
                m += cm;
                wx += cm * cx;
                wy += cm * cy;
                wz += cm * cz;
                cnt += self.nodes[c as usize].body_count;
            }

            self.nodes[idx].body_count = cnt;
            self.nodes[idx].mass = m;
            if m > 0.0 {
                self.nodes[idx].com_x = wx / m;
                self.nodes[idx].com_y = wy / m;
                self.nodes[idx].com_z = wz / m;
            }

            (
                self.nodes[idx].mass,
                self.nodes[idx].com_x,
                self.nodes[idx].com_y,
                self.nodes[idx].com_z,
            )
        }

        fn aggregate_quadrupole(&mut self, idx: usize, bodies: &[Body]) {
            if self.nodes[idx].is_leaf() {
                let cmx = self.nodes[idx].com_x;
                let cmy = self.nodes[idx].com_y;
                let cmz = self.nodes[idx].com_z;
                let len = self.nodes[idx].body_len as usize;

                let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);

                for k in 0..len {
                    let b = bodies[self.nodes[idx].bodies[k] as usize];
                    let dx = b.pos_x - cmx;
                    let dy = b.pos_y - cmy;
                    let dz = b.pos_z - cmz;
                    let d2 = dx * dx + dy * dy + dz * dz;
                    q_xx += b.mass * (3.0 * dx * dx - d2);
                    q_xy += b.mass * 3.0 * dx * dy;
                    q_xz += b.mass * 3.0 * dx * dz;
                    q_yy += b.mass * (3.0 * dy * dy - d2);
                    q_yz += b.mass * 3.0 * dy * dz;
                }

                let n = &mut self.nodes[idx];
                n.q_xx = q_xx;
                n.q_xy = q_xy;
                n.q_xz = q_xz;
                n.q_yy = q_yy;
                n.q_yz = q_yz;
                return;
            }

            let children = self.nodes[idx].children;
            for &c in &children {
                if c != NO_CHILD {
                    self.aggregate_quadrupole(c as usize, bodies);
                }
            }

            let pcom_x = self.nodes[idx].com_x;
            let pcom_y = self.nodes[idx].com_y;
            let pcom_z = self.nodes[idx].com_z;

            let (mut q_xx, mut q_xy, mut q_xz, mut q_yy, mut q_yz) = (0.0, 0.0, 0.0, 0.0, 0.0);

            for &c in &children {
                if c == NO_CHILD {
                    continue;
                }
                let child = &self.nodes[c as usize];
                let dx = child.com_x - pcom_x;
                let dy = child.com_y - pcom_y;
                let dz = child.com_z - pcom_z;
                let d2 = dx * dx + dy * dy + dz * dz;
                let m = child.mass;

                q_xx += child.q_xx + m * (3.0 * dx * dx - d2);
                q_xy += child.q_xy + m * 3.0 * dx * dy;
                q_xz += child.q_xz + m * 3.0 * dx * dz;
                q_yy += child.q_yy + m * (3.0 * dy * dy - d2);
                q_yz += child.q_yz + m * 3.0 * dy * dz;
            }

            let n = &mut self.nodes[idx];
            n.q_xx = q_xx;
            n.q_xy = q_xy;
            n.q_xz = q_xz;
            n.q_yy = q_yy;
            n.q_yz = q_yz;
        }
    }

    pub(super) fn evaluate_profile(
        octree: &Octree,
        bodies: &[Body],
        theta: f64,
        exact_threshold: usize,
        kernel: &dyn Kernel,
        acc: &mut [Vec3],
    ) -> (f64, WalkCounters) {
        let n = bodies.len();
        acc.fill(Vec3::ZERO);

        if n == 0 {
            return (0.0, WalkCounters::default());
        }

        if n <= exact_threshold {
            return (exact_eval(bodies, kernel, acc), WalkCounters::default());
        }

        let nodes = octree.nodes();

        // Parallel walk via rayon — matches the pre-PR-perf-5 production
        // walk (engine.rs evaluate_profile used the same `into_par_iter`
        // pattern at commit 03cb1a6). Sequential AoS would conflate the
        // layout effect with a parallelism effect that production never
        // had on this code path.
        let results: Vec<(Vec3, f64, WalkCounters)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                bh_eval_body(nodes, i, &bodies[i], bodies, theta, kernel, &mut stack)
            })
            .collect();

        let mut potential = 0.0_f64;
        let mut counters = WalkCounters::default();
        for (i, (a, phi, c)) in results.into_iter().enumerate() {
            acc[i] = a;
            potential += bodies[i].mass * phi;
            counters.merge(&c);
        }

        (0.5 * potential, counters)
    }

    fn exact_eval(bodies: &[Body], kernel: &dyn Kernel, acc: &mut [Vec3]) -> f64 {
        let n = bodies.len();
        let mut potential = 0.0_f64;

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = bodies[j].pos_x - bodies[i].pos_x;
                let dy = bodies[j].pos_y - bodies[i].pos_y;
                let dz = bodies[j].pos_z - bodies[i].pos_z;
                let eps2 = pair_eps2(bodies[i].softening, bodies[j].softening);
                let r_sq = dx * dx + dy * dy + dz * dz;

                let fac = G * kernel.acceleration_factor(r_sq, eps2);

                acc[i].x += bodies[j].mass * dx * fac;
                acc[i].y += bodies[j].mass * dy * fac;
                acc[i].z += bodies[j].mass * dz * fac;
                acc[j].x -= bodies[i].mass * dx * fac;
                acc[j].y -= bodies[i].mass * dy * fac;
                acc[j].z -= bodies[i].mass * dz * fac;

                let phi_ij = -G * bodies[j].mass * kernel.potential(r_sq, eps2);
                potential += bodies[i].mass * phi_ij;
            }
        }

        potential
    }

    #[inline(always)]
    fn bh_eval_body(
        nodes: &[Node<DEFAULT_LEAF>],
        body_idx: usize,
        body: &Body,
        bodies: &[Body],
        theta: f64,
        kernel: &dyn Kernel,
        stack: &mut Vec<u32>,
    ) -> (Vec3, f64, WalkCounters) {
        let mut a = Vec3::ZERO;
        let mut phi = 0.0_f64;
        let mut counters = WalkCounters::default();

        stack.clear();
        if !nodes.is_empty() {
            stack.push(0);
        }

        while let Some(raw) = stack.pop() {
            counters.n_node_visits += 1;
            let node = &nodes[raw as usize];
            if node.mass <= 0.0 {
                continue;
            }

            if node.is_leaf() {
                for k in 0..node.body_len as usize {
                    let bi = node.bodies[k] as usize;
                    if bi == body_idx {
                        continue;
                    }
                    let other = bodies[bi];
                    let dx = other.pos_x - body.pos_x;
                    let dy = other.pos_y - body.pos_y;
                    let dz = other.pos_z - body.pos_z;
                    let eps2 = pair_eps2(body.softening, other.softening);
                    let r_sq = dx * dx + dy * dy + dz * dz;

                    let fac = G * other.mass * kernel.acceleration_factor(r_sq, eps2);
                    a.x += dx * fac;
                    a.y += dy * fac;
                    a.z += dz * fac;
                    phi += -G * other.mass * kernel.potential(r_sq, eps2);
                    counters.n_leaf_interactions += 1;
                }
                continue;
            }

            let dx = node.com_x - body.pos_x;
            let dy = node.com_y - body.pos_y;
            let dz = node.com_z - body.pos_z;
            let eps2 = body.softening * body.softening;
            let d = (dx * dx + dy * dy + dz * dz + eps2).sqrt();

            if node.size() / d < theta {
                let r_sq = dx * dx + dy * dy + dz * dz;
                let fac = G * node.mass * kernel.acceleration_factor(r_sq, eps2);
                a.x += dx * fac;
                a.y += dy * fac;
                a.z += dz * fac;
                phi += -G * node.mass * kernel.potential(r_sq, eps2);

                let q_zz = -(node.q_xx + node.q_yy);
                let qr_x = node.q_xx * dx + node.q_xy * dy + node.q_xz * dz;
                let qr_y = node.q_xy * dx + node.q_yy * dy + node.q_yz * dz;
                let qr_z = node.q_xz * dx + node.q_yz * dy + q_zz * dz;
                let rqr = dx * qr_x + dy * qr_y + dz * qr_z;

                let inv_r2 = 1.0 / (r_sq + eps2);
                let inv_r5 = fac / node.mass * inv_r2;
                let inv_r7 = inv_r5 * inv_r2;

                let coef_qr = -G * inv_r5;
                let coef_r = 2.5 * G * rqr * inv_r7;

                a.x += coef_qr * qr_x + coef_r * dx;
                a.y += coef_qr * qr_y + coef_r * dy;
                a.z += coef_qr * qr_z + coef_r * dz;
                phi += -0.5 * G * rqr * inv_r5;
                counters.n_bh_accepted += 1;
            } else {
                for &c in &node.children {
                    if c != NO_CHILD {
                        stack.push(c);
                    }
                }
            }
        }

        (a, phi, counters)
    }
}

// ── Harness — AoS vs SoA at canonical seeds ─────────────────────────────────── //

#[test]
#[ignore = "perf experiment; opt-in via cargo test --release perf_soa_aos_vs_soa -- --ignored --nocapture"]
fn perf_soa_aos_vs_soa() {
    let out_dir = perf_soa_output_dir();
    fs::create_dir_all(&out_dir).expect("create perf output dir");
    let csv_path = out_dir.join("profile.csv");
    let mut writer = fs::File::create(&csv_path).expect("create profile.csv");
    write_header(&mut writer);

    eprintln!("[perf-soa] AoS shadow vs SoA production, N={:?}, seeds={:?}", N_VALUES, SEEDS);
    eprintln!(
        "[perf-soa]   protocol: per (N, seed): {} warmup + {} measured runs, within-seed median",
        WARMUP_RUNS, MEASURED_RUNS
    );

    let t_total = Instant::now();
    for &n in &N_VALUES {
        for &seed in &SEEDS {
            let bodies = sphere_distribution_lognormal(n, seed);
            let aos = measure_aos(&bodies);
            let soa = measure_soa(&bodies);

            write_row(&mut writer, n, seed, "AoS_baseline", &aos);
            write_row(&mut writer, n, seed, "SoA_production", &soa);
            print_pair(n, seed, &aos, &soa);
        }
    }
    eprintln!("[perf-soa] runtime: {:.1}s", t_total.elapsed().as_secs_f64());
    eprintln!("[perf-soa] wrote {}", csv_path.display());
}

#[derive(Debug, Clone, Copy)]
struct Row {
    t_total_median_ms: f64,
    t_build_median_ms: f64,
    t_walk_median_ms: f64,
    t_pack_median_ms: f64,
    n_node_visits: u64,
    n_bh_accepted: u64,
    n_leaf_interactions: u64,
}

fn measure_aos(bodies: &[Body]) -> Row {
    let kernel = PlummerKernel::new();
    let mut octree = aos_baseline::Octree::new(16);
    let mut acc = vec![Vec3::ZERO; bodies.len()];

    // Warmup
    for _ in 0..WARMUP_RUNS {
        octree.build(bodies);
        let _ = aos_baseline::evaluate_profile(&octree, bodies, THETA, 1, &kernel, &mut acc);
    }

    let mut total_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut build_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut walk_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut last_counters = (0u64, 0u64, 0u64);
    for _ in 0..MEASURED_RUNS {
        let t_total = Instant::now();
        let t_build = Instant::now();
        octree.build(bodies);
        let build = t_build.elapsed().as_secs_f64() * 1000.0;

        let t_walk = Instant::now();
        let (_, c) = aos_baseline::evaluate_profile(&octree, bodies, THETA, 1, &kernel, &mut acc);
        let walk = t_walk.elapsed().as_secs_f64() * 1000.0;

        total_ms.push(t_total.elapsed().as_secs_f64() * 1000.0);
        build_ms.push(build);
        walk_ms.push(walk);
        last_counters = (c.n_node_visits, c.n_bh_accepted, c.n_leaf_interactions);
    }
    sort_in_place(&mut total_ms);
    sort_in_place(&mut build_ms);
    sort_in_place(&mut walk_ms);

    Row {
        t_total_median_ms: total_ms[MEASURED_RUNS / 2],
        t_build_median_ms: build_ms[MEASURED_RUNS / 2],
        t_walk_median_ms: walk_ms[MEASURED_RUNS / 2],
        t_pack_median_ms: 0.0,
        n_node_visits: last_counters.0,
        n_bh_accepted: last_counters.1,
        n_leaf_interactions: last_counters.2,
    }
}

fn measure_soa(bodies: &[Body]) -> Row {
    let mut engine = BarnesHutEngine::new(16);
    engine.set_exact_threshold(1);
    let mut arrays = BodyArrays::with_capacity(bodies.len());
    let mut acc = vec![Vec3::ZERO; bodies.len()];

    // Warmup
    for _ in 0..WARMUP_RUNS {
        arrays.pack_from(bodies);
        engine.build(&arrays);
        let _ = engine.evaluate_profile(&arrays, THETA, &mut acc);
    }

    let mut total_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut build_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut walk_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut pack_ms = Vec::with_capacity(MEASURED_RUNS);
    let mut last_counters = (0u64, 0u64, 0u64);
    for _ in 0..MEASURED_RUNS {
        let t_total = Instant::now();
        let t_pack = Instant::now();
        arrays.pack_from(bodies);
        let pack = t_pack.elapsed().as_secs_f64() * 1000.0;

        let t_build = Instant::now();
        engine.build(&arrays);
        let build = t_build.elapsed().as_secs_f64() * 1000.0;

        let t_walk = Instant::now();
        let (_, c) = engine.evaluate_profile(&arrays, THETA, &mut acc);
        let walk = t_walk.elapsed().as_secs_f64() * 1000.0;

        total_ms.push(t_total.elapsed().as_secs_f64() * 1000.0);
        build_ms.push(build);
        walk_ms.push(walk);
        pack_ms.push(pack);
        last_counters = (c.n_node_visits, c.n_bh_accepted, c.n_leaf_interactions);
    }
    sort_in_place(&mut total_ms);
    sort_in_place(&mut build_ms);
    sort_in_place(&mut walk_ms);
    sort_in_place(&mut pack_ms);

    Row {
        t_total_median_ms: total_ms[MEASURED_RUNS / 2],
        t_build_median_ms: build_ms[MEASURED_RUNS / 2],
        t_walk_median_ms: walk_ms[MEASURED_RUNS / 2],
        t_pack_median_ms: pack_ms[MEASURED_RUNS / 2],
        n_node_visits: last_counters.0,
        n_bh_accepted: last_counters.1,
        n_leaf_interactions: last_counters.2,
    }
}

fn sort_in_place(v: &mut [f64]) {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

// ── CSV / printing ─────────────────────────────────────────────────────────── //

fn write_header(writer: &mut fs::File) {
    writeln!(
        writer,
        "n,seed,path,t_total_ms,t_build_ms,t_walk_ms,t_pack_ms,\
         n_node_visits,n_bh_accepted,n_leaf_interactions,n_runs,warmup_runs"
    )
    .unwrap();
}

fn write_row(writer: &mut fs::File, n: usize, seed: u64, path: &str, r: &Row) {
    writeln!(
        writer,
        "{},0x{:X},{},{:.6},{:.6},{:.6},{:.6},{},{},{},{},{}",
        n,
        seed,
        path,
        r.t_total_median_ms,
        r.t_build_median_ms,
        r.t_walk_median_ms,
        r.t_pack_median_ms,
        r.n_node_visits,
        r.n_bh_accepted,
        r.n_leaf_interactions,
        MEASURED_RUNS,
        WARMUP_RUNS
    )
    .unwrap();
}

fn print_pair(n: usize, seed: u64, aos: &Row, soa: &Row) {
    let walk_ratio = aos.t_walk_median_ms / soa.t_walk_median_ms;
    let total_ratio = aos.t_total_median_ms / soa.t_total_median_ms;
    let build_ratio = aos.t_build_median_ms / soa.t_build_median_ms;
    let pack_pct = soa.t_pack_median_ms / soa.t_total_median_ms.max(1e-12) * 100.0;

    eprintln!(
        "[perf-soa] N={:>5} seed=0x{:X}  AoS: total={:>7.3}ms (build={:.3} walk={:.3})",
        n, seed, aos.t_total_median_ms, aos.t_build_median_ms, aos.t_walk_median_ms
    );
    eprintln!(
        "[perf-soa] N={:>5} seed=0x{:X}  SoA: total={:>7.3}ms (pack={:.3} build={:.3} walk={:.3})",
        n,
        seed,
        soa.t_total_median_ms,
        soa.t_pack_median_ms,
        soa.t_build_median_ms,
        soa.t_walk_median_ms
    );
    eprintln!(
        "[perf-soa] N={:>5} seed=0x{:X}  ratios AoS/SoA: walk={:.3}  build={:.3}  total={:.3}  | pack={:.2}% of SoA total",
        n, seed, walk_ratio, build_ratio, total_ratio, pack_pct
    );
}

fn perf_soa_output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/perf-soa")
}

// ── Body distribution (matches engine_ceiling.rs / perf 2×2 / MAC) ─────────── //

fn sphere_distribution_lognormal(n: usize, seed: u64) -> Vec<Body> {
    let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
    let mut next_u64 = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        state
    };
    let mut next_unit = || (next_u64() >> 11) as f64 / (1u64 << 53) as f64;

    let mut bodies = Vec::with_capacity(n);
    while bodies.len() < n {
        let x = 2.0 * next_unit() - 1.0;
        let y = 2.0 * next_unit() - 1.0;
        let z = 2.0 * next_unit() - 1.0;
        if x * x + y * y + z * z > 1.0 {
            continue;
        }
        let u1 = next_unit().max(1e-12);
        let u2 = next_unit();
        let normal = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        let mass = normal.exp();
        let mut b = Body::rocky(mass).at(x, y).with_velocity(0.0, 0.0);
        b.pos_z = z;
        bodies.push(b);
    }
    bodies
}
