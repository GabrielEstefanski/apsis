//! Structured comparison of two records. See ADR-012 §"Semantic diff API".

use crate::records::frame::Snapshot;
use crate::records::header::Header;
use crate::records::reader::{Record, RecordError};

/// Result of comparing two records: a categorised list of header
/// differences and a summary of the frame stream.
#[derive(Debug, Clone, PartialEq)]
pub struct RecordDiff {
    pub header: Vec<HeaderChange>,
    pub frames: FrameStreamDiff,
}

impl RecordDiff {
    /// `true` when both header and frame stream are byte-identical
    /// (no header changes, trailer hashes match, final-snapshot rms = 0).
    pub fn is_empty(&self) -> bool {
        self.header.is_empty()
            && self.frames.trailer_blake3_match
            && self.frames.trajectory_rms_at_final == Some(0.0)
    }
}

/// One categorised difference between the two headers. Field order is
/// `(before, after)` throughout.
#[derive(Debug, Clone, PartialEq)]
pub enum HeaderChange {
    OperatorAdded { name: String, version: String, crate_hash: String },
    OperatorRemoved { name: String, version: String },
    OperatorVersionChanged { name: String, before: String, after: String },
    OperatorCrateHashChanged { name: String, before: String, after: String },
    IntegratorKindChanged { before: String, after: String },
    IntegratorDtModeChanged { before: String, after: String },
    IntegratorInitialDtChanged { before: f64, after: f64 },
    IntegratorParamsChanged,
    KernelVariantChanged { before: String, after: String },
    KernelSofteningChanged { before: Option<f64>, after: Option<f64> },
    SeedChanged { before: u64, after: u64 },
    UnitSystemChanged { before: String, after: String },
    ApsisVersionChanged { before: String, after: String },
    RustcVersionChanged { before: String, after: String },
    CargoLockChanged { before: String, after: String },
    BodyCountChanged { before: usize, after: usize },
}

/// Frame-stream summary. `trajectory_rms_at_final` is `None` when the
/// records carry different body counts (no point-to-point pairing).
#[derive(Debug, Clone, PartialEq)]
pub struct FrameStreamDiff {
    pub event_count: (usize, usize),
    pub diagnostic_count: (usize, usize),
    pub snapshot_count: (usize, usize),
    pub trajectory_rms_at_final: Option<f64>,
    pub trailer_blake3_match: bool,
    pub trailer_step_count: (u64, u64),
}

impl Record {
    /// Categorised diff against another record.
    pub fn diff(&self, other: &Self) -> Result<RecordDiff, RecordError> {
        let header = diff_headers(self.header(), other.header());
        let frames = diff_frames(self, other)?;
        Ok(RecordDiff { header, frames })
    }
}

fn diff_headers(a: &Header, b: &Header) -> Vec<HeaderChange> {
    let mut out = Vec::new();

    if a.apsis.version != b.apsis.version {
        out.push(HeaderChange::ApsisVersionChanged {
            before: a.apsis.version.clone(),
            after: b.apsis.version.clone(),
        });
    }
    if a.apsis.rustc_version != b.apsis.rustc_version {
        out.push(HeaderChange::RustcVersionChanged {
            before: a.apsis.rustc_version.clone(),
            after: b.apsis.rustc_version.clone(),
        });
    }
    if a.reproducibility.cargo_lock_blake3 != b.reproducibility.cargo_lock_blake3 {
        out.push(HeaderChange::CargoLockChanged {
            before: a.reproducibility.cargo_lock_blake3.clone(),
            after: b.reproducibility.cargo_lock_blake3.clone(),
        });
    }
    if a.reproducibility.seed != b.reproducibility.seed {
        out.push(HeaderChange::SeedChanged {
            before: a.reproducibility.seed,
            after: b.reproducibility.seed,
        });
    }
    if a.unit_system != b.unit_system {
        out.push(HeaderChange::UnitSystemChanged {
            before: format!("{:?}", a.unit_system),
            after: format!("{:?}", b.unit_system),
        });
    }
    if a.integrator.kind != b.integrator.kind {
        out.push(HeaderChange::IntegratorKindChanged {
            before: a.integrator.kind.clone(),
            after: b.integrator.kind.clone(),
        });
    }
    if a.integrator.dt_mode != b.integrator.dt_mode {
        out.push(HeaderChange::IntegratorDtModeChanged {
            before: a.integrator.dt_mode.clone(),
            after: b.integrator.dt_mode.clone(),
        });
    }
    if a.integrator.initial_dt != b.integrator.initial_dt {
        out.push(HeaderChange::IntegratorInitialDtChanged {
            before: a.integrator.initial_dt,
            after: b.integrator.initial_dt,
        });
    }
    if a.integrator.params != b.integrator.params {
        out.push(HeaderChange::IntegratorParamsChanged);
    }
    if a.kernel.variant != b.kernel.variant {
        out.push(HeaderChange::KernelVariantChanged {
            before: a.kernel.variant.clone(),
            after: b.kernel.variant.clone(),
        });
    }
    if a.kernel.softening != b.kernel.softening {
        out.push(HeaderChange::KernelSofteningChanged {
            before: a.kernel.softening,
            after: b.kernel.softening,
        });
    }
    if a.bodies.count != b.bodies.count {
        out.push(HeaderChange::BodyCountChanged { before: a.bodies.count, after: b.bodies.count });
    }

    out.extend(diff_operators(a, b));
    out
}

fn diff_operators(a: &Header, b: &Header) -> Vec<HeaderChange> {
    let mut out = Vec::new();
    let a_names: std::collections::BTreeMap<&str, &crate::records::header::OperatorMeta> =
        a.operators.iter().map(|op| (op.name.as_str(), op)).collect();
    let b_names: std::collections::BTreeMap<&str, &crate::records::header::OperatorMeta> =
        b.operators.iter().map(|op| (op.name.as_str(), op)).collect();

    for (name, op_b) in &b_names {
        if !a_names.contains_key(name) {
            out.push(HeaderChange::OperatorAdded {
                name: op_b.name.clone(),
                version: op_b.version.clone(),
                crate_hash: op_b.crate_hash.clone(),
            });
        }
    }
    for (name, op_a) in &a_names {
        match b_names.get(name) {
            None => {
                out.push(HeaderChange::OperatorRemoved {
                    name: op_a.name.clone(),
                    version: op_a.version.clone(),
                });
            },
            Some(op_b) => {
                if op_a.version != op_b.version {
                    out.push(HeaderChange::OperatorVersionChanged {
                        name: op_a.name.clone(),
                        before: op_a.version.clone(),
                        after: op_b.version.clone(),
                    });
                }
                if op_a.crate_hash != op_b.crate_hash {
                    out.push(HeaderChange::OperatorCrateHashChanged {
                        name: op_a.name.clone(),
                        before: op_a.crate_hash.clone(),
                        after: op_b.crate_hash.clone(),
                    });
                }
            },
        }
    }
    out
}

fn diff_frames(a: &Record, b: &Record) -> Result<FrameStreamDiff, RecordError> {
    let events_a = a.events()?.count();
    let events_b = b.events()?.count();
    let diag_a = a.diagnostics()?.count();
    let diag_b = b.diagnostics()?.count();
    let snap_a = a.dense()?.count();
    let snap_b = b.dense()?.count();

    let (_, final_a) = a.bookends()?;
    let (_, final_b) = b.bookends()?;
    let rms = final_snapshot_rms(&final_a, &final_b);

    let trailer_blake3_match = a.trailer().blake3 == b.trailer().blake3;
    let trailer_step_count = (a.trailer().step_count, b.trailer().step_count);

    Ok(FrameStreamDiff {
        event_count: (events_a, events_b),
        diagnostic_count: (diag_a, diag_b),
        snapshot_count: (snap_a, snap_b),
        trajectory_rms_at_final: rms,
        trailer_blake3_match,
        trailer_step_count,
    })
}

fn final_snapshot_rms(a: &Snapshot, b: &Snapshot) -> Option<f64> {
    if a.bodies.len() != b.bodies.len() || a.bodies.is_empty() {
        return None;
    }
    let n = a.bodies.len() as f64;
    let sum: f64 = a
        .bodies
        .iter()
        .zip(&b.bodies)
        .map(|(ba, bb)| {
            let dx = ba.pos[0] - bb.pos[0];
            let dy = ba.pos[1] - bb.pos[1];
            let dz = ba.pos[2] - bb.pos[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    Some((sum / n).sqrt())
}
