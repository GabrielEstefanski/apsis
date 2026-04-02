use crate::domain::body::Body;
use rayon::prelude::*;

pub const G: f64 = 1.0;
pub const EPS: f64 = 1e-4;

const LEAF_CAPACITY: usize = 8;
const EXACT_THRESHOLD: usize = 64;
const NO_CHILD: u32 = u32::MAX;

#[derive(Clone, Copy)]
struct Node {
    cx: f64,
    cy: f64,
    half: f64,
    mass: f64,
    com_x: f64,
    com_y: f64,
    body_count: u32,
    body_len: u8,
    bodies: [u32; LEAF_CAPACITY],
    children: [u32; 4],
}

impl Node {
    fn new(cx: f64, cy: f64, half: f64) -> Self {
        Self {
            cx,
            cy,
            half,
            mass: 0.0,
            com_x: 0.0,
            com_y: 0.0,
            body_count: 0,
            body_len: 0,
            bodies: [0u32; LEAF_CAPACITY],
            children: [NO_CHILD; 4],
        }
    }

    #[inline]
    fn is_leaf(&self) -> bool {
        self.children[0] == NO_CHILD
    }

    #[inline]
    fn size(&self) -> f64 {
        self.half * 2.0
    }
}

pub struct BarnesHutEngine {
    max_depth: usize,
    nodes: Vec<Node>,
}

impl BarnesHutEngine {
    pub fn new(max_depth: usize) -> Self {
        Self {
            max_depth,
            nodes: Vec::new(),
        }
    }

    pub fn build(&mut self, bodies: &[Body]) {
        self.nodes.clear();
        if bodies.is_empty() {
            return;
        }

        let mut min_x = bodies[0].x;
        let mut max_x = bodies[0].x;
        let mut min_y = bodies[0].y;
        let mut max_y = bodies[0].y;

        for b in &bodies[1..] {
            min_x = min_x.min(b.x);
            max_x = max_x.max(b.x);
            min_y = min_y.min(b.y);
            max_y = max_y.max(b.y);
        }

        let cx = 0.5 * (min_x + max_x);
        let cy = 0.5 * (min_y + max_y);
        let mut half = 0.5 * (max_x - min_x).max(max_y - min_y);
        half = if half <= 0.0 {
            EPS.sqrt()
        } else {
            half * 1.0001 + EPS.sqrt()
        };

        self.nodes.push(Node::new(cx, cy, half));

        for i in 0..bodies.len() {
            self.insert_body(0, i, bodies);
        }

        self.aggregate_mass(0, bodies);
    }

    pub fn evaluate(&self, bodies: &[Body], theta: f64, acc: &mut [(f64, f64)]) -> f64 {
        let n = bodies.len();
        acc.fill((0.0, 0.0));

        if n == 0 {
            return 0.0;
        }

        if n <= EXACT_THRESHOLD {
            return exact_eval(bodies, acc);
        }

        let nodes = &self.nodes[..];

        let results: Vec<(f64, f64, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut stack = Vec::with_capacity(128);
                bh_eval_body(nodes, i, &bodies[i], bodies, theta, &mut stack)
            })
            .collect();

        let mut potential_sum = 0.0_f64;
        for (i, (ax, ay, phi)) in results.into_iter().enumerate() {
            acc[i] = (ax, ay);
            potential_sum += bodies[i].mass * phi;
        }
        0.5 * potential_sum
    }

    pub fn theta_error_proxy(&self, body_idx: usize, bodies: &[Body], theta: f64) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }

        let body = &bodies[body_idx];
        let eps2 = EPS * EPS;
        let mut violation_sum = 0.0_f64;
        let mut weight_sum = 0.0_f64;

        let mut stack: Vec<u32> = Vec::with_capacity(64);
        stack.push(0);

        while let Some(raw) = stack.pop() {
            let node = &self.nodes[raw as usize];
            if node.mass <= 0.0 || node.is_leaf() {
                continue;
            }

            let dx = node.com_x - body.x;
            let dy = node.com_y - body.y;
            let d = (dx * dx + dy * dy + eps2).sqrt();
            let ratio = node.size() / d;

            if ratio < theta {
                violation_sum += node.mass * ratio * ratio;
                weight_sum += node.mass;
            } else {
                for &c in &node.children {
                    if c != NO_CHILD {
                        stack.push(c);
                    }
                }
            }
        }

        if weight_sum > 0.0 {
            (violation_sum / weight_sum).sqrt()
        } else {
            0.0
        }
    }

    fn insert_body(&mut self, root: usize, body_idx: usize, bodies: &[Body]) {
        self.insert_body_at(root, body_idx, bodies, 0);
    }

    fn insert_body_at(
        &mut self,
        mut node_idx: usize,
        body_idx: usize,
        bodies: &[Body],
        mut depth: usize,
    ) {
        loop {
            if depth > self.max_depth {
                let node = &mut self.nodes[node_idx];
                if (node.body_len as usize) < LEAF_CAPACITY {
                    node.bodies[node.body_len as usize] = body_idx as u32;
                    node.body_len += 1;
                }
                return;
            }

            if self.nodes[node_idx].is_leaf() {
                let len = self.nodes[node_idx].body_len as usize;

                if len < LEAF_CAPACITY || depth == self.max_depth {
                    if (self.nodes[node_idx].body_len as usize) < LEAF_CAPACITY {
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
                    let child = self.child_for_body(node_idx, bi as usize, bodies);
                    self.insert_body_at(child, bi as usize, bodies, depth + 1);
                }
            }

            node_idx = self.child_for_body(node_idx, body_idx, bodies);
            depth += 1;
        }
    }

    fn subdivide(&mut self, idx: usize) {
        let (cx, cy, half) = {
            let n = &self.nodes[idx];
            (n.cx, n.cy, n.half)
        };
        let h = half * 0.5;
        let c = [
            self.push_node(cx - h, cy - h, h),
            self.push_node(cx + h, cy - h, h),
            self.push_node(cx - h, cy + h, h),
            self.push_node(cx + h, cy + h, h),
        ];
        self.nodes[idx].children = [c[0] as u32, c[1] as u32, c[2] as u32, c[3] as u32];
    }

    fn push_node(&mut self, cx: f64, cy: f64, half: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node::new(cx, cy, half));
        idx
    }

    fn child_for_body(&self, node_idx: usize, body_idx: usize, bodies: &[Body]) -> usize {
        let n = &self.nodes[node_idx];
        let b = bodies[body_idx];
        let q = match (b.x >= n.cx, b.y >= n.cy) {
            (false, false) => 0,
            (true, false) => 1,
            (false, true) => 2,
            (true, true) => 3,
        };
        self.nodes[node_idx].children[q] as usize
    }

    fn aggregate_mass(&mut self, idx: usize, bodies: &[Body]) -> (f64, f64, f64) {
        if self.nodes[idx].is_leaf() {
            let len = self.nodes[idx].body_len as usize;
            let mut m = 0.0_f64;
            let mut wx = 0.0_f64;
            let mut wy = 0.0_f64;

            for k in 0..len {
                let bi = self.nodes[idx].bodies[k] as usize;
                let b = bodies[bi];
                m += b.mass;
                wx += b.mass * b.x;
                wy += b.mass * b.y;
            }

            self.nodes[idx].body_count = len as u32;
            self.nodes[idx].mass = m;

            if m > 0.0 {
                self.nodes[idx].com_x = wx / m;
                self.nodes[idx].com_y = wy / m;
                return (m, self.nodes[idx].com_x, self.nodes[idx].com_y);
            }
            self.nodes[idx].com_x = 0.0;
            self.nodes[idx].com_y = 0.0;
            return (0.0, 0.0, 0.0);
        }

        let children = self.nodes[idx].children;
        let mut m = 0.0_f64;
        let mut wx = 0.0_f64;
        let mut wy = 0.0_f64;
        let mut cnt = 0u32;

        for &c in &children {
            if c == NO_CHILD {
                continue;
            }
            let (cm, cx, cy) = self.aggregate_mass(c as usize, bodies);
            m += cm;
            wx += cm * cx;
            wy += cm * cy;
            cnt += self.nodes[c as usize].body_count;
        }

        self.nodes[idx].body_count = cnt;
        self.nodes[idx].mass = m;

        if m > 0.0 {
            self.nodes[idx].com_x = wx / m;
            self.nodes[idx].com_y = wy / m;
        }

        (
            self.nodes[idx].mass,
            self.nodes[idx].com_x,
            self.nodes[idx].com_y,
        )
    }
}

fn exact_eval(bodies: &[Body], acc: &mut [(f64, f64)]) -> f64 {
    let n = bodies.len();
    let eps2 = EPS * EPS;
    let mut potential = 0.0_f64;

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = bodies[j].x - bodies[i].x;
            let dy = bodies[j].y - bodies[i].y;

            let d2 = dx * dx + dy * dy + eps2;
            let inv_r = d2.sqrt().recip();
            let inv_r3 = inv_r * inv_r * inv_r;
            let fac = G * inv_r3;

            let fx = dx * fac;
            let fy = dy * fac;

            acc[i].0 += bodies[j].mass * fx;
            acc[i].1 += bodies[j].mass * fy;
            acc[j].0 -= bodies[i].mass * fx;
            acc[j].1 -= bodies[i].mass * fy;

            potential += -G * bodies[i].mass * bodies[j].mass * inv_r;
        }
    }

    potential
}

fn bh_eval_body(
    nodes: &[Node],
    body_idx: usize,
    body: &Body,
    bodies: &[Body],
    theta: f64,
    stack: &mut Vec<u32>,
) -> (f64, f64, f64) {
    let eps2 = EPS * EPS;
    let mut ax = 0.0_f64;
    let mut ay = 0.0_f64;
    let mut phi = 0.0_f64;

    stack.clear();
    if !nodes.is_empty() {
        stack.push(0);
    }

    while let Some(raw) = stack.pop() {
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
                let dx = other.x - body.x;
                let dy = other.y - body.y;
                let d2 = dx * dx + dy * dy + eps2;
                let inv_r = d2.sqrt().recip();
                let gm = G * other.mass;
                let fac = gm * inv_r * inv_r * inv_r;
                ax += dx * fac;
                ay += dy * fac;
                phi += -gm * inv_r;
            }
            continue;
        }

        let dx = node.com_x - body.x;
        let dy = node.com_y - body.y;
        let d2 = dx * dx + dy * dy + eps2;
        let d = d2.sqrt();

        if node.size() / d < theta {
            let inv_r = d.recip();
            let gm = G * node.mass;
            let fac = gm * inv_r * inv_r * inv_r;
            ax += dx * fac;
            ay += dy * fac;
            phi += -gm * inv_r;
        } else {
            for &c in &node.children {
                if c != NO_CHILD {
                    stack.push(c);
                }
            }
        }
    }

    (ax, ay, phi)
}
