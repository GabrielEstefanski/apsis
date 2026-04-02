use crate::domain::body::Body;
use crate::physics::gravity::BarnesHutEngine;

pub struct StepResult<'a> {
    pub acc1: &'a [(f64, f64)],
    pub potential: f64,
}

pub fn step<'a>(
    bodies: &mut [Body],
    dt: f64,
    theta: f64,
    engine: &mut BarnesHutEngine,
    scratch_acc: &'a mut Vec<(f64, f64)>,
) -> StepResult<'a> {
    let n = bodies.len();

    if scratch_acc.len() != n {
        scratch_acc.resize(n, (0.0, 0.0));
    }

    engine.build(bodies);
    engine.evaluate(bodies, theta, scratch_acc);

    for i in 0..n {
        let b = &mut bodies[i];
        b.vx += 0.5 * scratch_acc[i].0 * dt;
        b.vy += 0.5 * scratch_acc[i].1 * dt;

        b.x += b.vx * dt;
        b.y += b.vy * dt;
    }

    engine.build(bodies);
    let potential = engine.evaluate(bodies, theta, scratch_acc);

    for i in 0..n {
        let b = &mut bodies[i];
        b.vx += 0.5 * scratch_acc[i].0 * dt;
        b.vy += 0.5 * scratch_acc[i].1 * dt;
    }

    StepResult {
        acc1: scratch_acc.as_slice(),
        potential,
    }
}
