use crate::domain::body::Body;
use crate::math::Vec3;

#[derive(Debug, Clone, Copy, Default)]
pub struct SimulationDiagnostics {
    pub max_acc: f64,
    pub jerk: f64,
    pub max_vel: f64,
}

#[derive(Debug, Clone)]
pub struct DiagnosticsComputer {
    last_acc: Vec<Vec3>,
}

impl Default for DiagnosticsComputer {
    fn default() -> Self {
        Self::new()
    }
}

impl DiagnosticsComputer {
    pub fn new() -> Self {
        Self { last_acc: Vec::new() }
    }

    pub fn compute(
        &mut self,
        acc: &[Vec3],
        bodies: &[Body],
        current_dt: f64,
    ) -> SimulationDiagnostics {
        let n = acc.len();
        let has_last = self.last_acc.len() == n && current_dt > 0.0;

        let mut max_acc = 0.0_f64;
        let mut max_jerk = 0.0_f64;
        let mut max_vel = 0.0_f64;

        for i in 0..n {
            let a = acc[i];
            max_acc = max_acc.max(a.length());

            if has_last {
                let da = a - self.last_acc[i];
                let j = da.length() / current_dt;
                max_jerk = max_jerk.max(j);
            }
        }

        for b in bodies {
            let v2 = b.vel_x * b.vel_x + b.vel_y * b.vel_y + b.vel_z * b.vel_z;
            max_vel = max_vel.max(v2.sqrt());
        }

        if self.last_acc.len() != n {
            self.last_acc.resize(n, Vec3::ZERO);
        }
        self.last_acc.copy_from_slice(acc);

        SimulationDiagnostics { max_acc, jerk: max_jerk, max_vel }
    }
}
