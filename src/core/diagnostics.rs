use crate::core::body::Body;

#[derive(Debug, Clone, Copy, Default)]
pub struct SimulationDiagnostics {
    pub max_acc: f64,
    pub jerk: f64,
    pub max_vel: f64,
}

#[derive(Debug, Clone)]
pub struct DiagnosticsComputer {
    last_acc: Vec<(f64, f64)>,
}

impl DiagnosticsComputer {
    pub fn new() -> Self {
        Self { last_acc: Vec::new() }
    }

    pub fn compute(
        &mut self,
        acc: &[(f64, f64)],
        bodies: &[Body],
        current_dt: f64,
    ) -> SimulationDiagnostics {
        let n = acc.len();
        let has_last = self.last_acc.len() == n && current_dt > 0.0;

        let mut max_acc = 0.0_f64;
        let mut max_jerk = 0.0_f64;
        let mut max_vel = 0.0_f64;

        for i in 0..n {
            let (ax, ay) = acc[i];
            let a2 = ax * ax + ay * ay;
            max_acc = max_acc.max(a2.sqrt());

            if has_last {
                let (lax, lay) = self.last_acc[i];
                let dax = ax - lax;
                let day = ay - lay;
                let j = (dax * dax + day * day).sqrt() / current_dt;
                max_jerk = max_jerk.max(j);
            }
        }

        for b in bodies {
            let v2 = b.vx * b.vx + b.vy * b.vy;
            max_vel = max_vel.max(v2.sqrt());
        }

        if self.last_acc.len() != n {
            self.last_acc.resize(n, (0.0, 0.0));
        }
        self.last_acc.copy_from_slice(acc);

        SimulationDiagnostics { max_acc, jerk: max_jerk, max_vel }
    }
}
