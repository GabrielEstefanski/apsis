use apsis::domain::body::Body;

const BELT_MIN_MEMBERS: usize = 8;
const BELT_BAND_FRACTION: f64 = 0.15;
const TRAIL_MASS_RATIO: f64 = 1e-6;

pub struct BodyRenderHints {
    pub show_trail: bool,
    pub belt_ring: Option<BeltRing>,
}

pub struct BeltRing {
    pub belt_id: usize,
    pub radius_world: f32,
    pub width_world: f32,
}

pub fn compute_render_hints(bodies: &[Body]) -> Vec<BodyRenderHints> {
    let dominant = bodies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.mass.partial_cmp(&b.1.mass).unwrap_or(std::cmp::Ordering::Equal));

    let Some((dom_idx, dom)) = dominant else {
        return bodies
            .iter()
            .map(|_| BodyRenderHints { show_trail: true, belt_ring: None })
            .collect();
    };

    let mut hints: Vec<BodyRenderHints> = bodies
        .iter()
        .map(|b| {
            let mass_ratio = b.mass / dom.mass;
            BodyRenderHints { show_trail: mass_ratio > TRAIL_MASS_RATIO, belt_ring: None }
        })
        .collect();

    detect_belts(bodies, dom_idx, &mut hints);
    hints
}

pub fn detect_belts(bodies: &[Body], dom_idx: usize, hints: &mut [BodyRenderHints]) {
    let dom = &bodies[dom_idx];

    let candidates: Vec<(usize, f64)> = bodies
        .iter()
        .enumerate()
        .filter(|(i, _)| !hints[*i].show_trail)
        .map(|(i, b)| {
            let dx = b.x - dom.x;
            let dy = b.y - dom.y;
            (i, (dx * dx + dy * dy).sqrt())
        })
        .collect();

    if candidates.len() < BELT_MIN_MEMBERS {
        return;
    }

    let mut sorted = candidates;
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let mut belt_id = 0usize;
    let mut i = 0;
    while i < sorted.len() {
        let anchor = sorted[i].1;
        let band_max = anchor * (1.0 + BELT_BAND_FRACTION);
        let group: Vec<_> = sorted[i..].iter().take_while(|(_, r)| *r <= band_max).collect();

        if group.len() >= BELT_MIN_MEMBERS {
            let mean_r = group.iter().map(|(_, r)| r).sum::<f64>() / group.len() as f64;
            let width = (group.last().unwrap().1 - group.first().unwrap().1).max(mean_r * 2.0);

            for (body_idx, _) in &group {
                hints[*body_idx].belt_ring = Some(BeltRing {
                    belt_id,
                    radius_world: mean_r as f32,
                    width_world: width as f32,
                });
            }
            belt_id += 1;
            i += group.len();
        } else {
            i += 1;
        }
    }
}
