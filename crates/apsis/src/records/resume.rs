//! Mid-run resume from a recorded snapshot. See ADR-012
//! §"Mid-run snapshot resume".
//!
//! v0.2 ships a `restore_into` API: the caller already holds a
//! `System` constructed with the same bodies, integrator type, units,
//! and operators as the recorded run; `restore_into` reloads the
//! dynamic state (body positions/velocities + integrator scratch) from
//! the n-th captured snapshot/`ResumeState` pair.
//!
//! Full reconstruction (`Record::resume_from(idx) -> System`) requires
//! mapping `BodyMeta` back to constructor presets and walking the
//! operator registry; deferred to a follow-up alongside the operator
//! re-instantiation protocol.

use crate::core::system::System;
use crate::physics::integrator::traits::ResumeError;
use crate::records::frame::Frame;
use crate::records::reader::{Record, RecordError};

/// Errors specific to `restore_into`. Wraps both record-side I/O and
/// integrator-side payload-validation failures.
#[derive(Debug)]
pub enum RestoreError {
    Record(RecordError),
    /// `snapshot_idx` is out of range for this record.
    SnapshotIndexOutOfRange {
        idx: usize,
        available: usize,
    },
    /// Snapshot and `ResumeState` frames are not paired one-to-one. The
    /// record must be written with
    /// [`RecordHook::with_resume_capture(true)`](crate::records::RecordHook::with_resume_capture).
    MissingResumeState {
        snapshot_idx: usize,
    },
    /// Body count in the snapshot does not match the System's
    /// `bodies.len()`.
    BodyCountMismatch {
        expected: usize,
        found: usize,
    },
    /// Integrator kind recorded in the header does not match the active
    /// integrator on the System.
    IntegratorMismatch {
        expected: String,
        found: String,
    },
    Integrator(ResumeError),
}

impl std::fmt::Display for RestoreError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Record(e) => write!(f, "record read: {e}"),
            Self::SnapshotIndexOutOfRange { idx, available } => {
                write!(f, "snapshot index {idx} out of range (record has {available})")
            },
            Self::MissingResumeState { snapshot_idx } => {
                write!(
                    f,
                    "no ResumeState frame paired with snapshot {snapshot_idx} \
                     (record written without with_resume_capture)"
                )
            },
            Self::BodyCountMismatch { expected, found } => {
                write!(f, "body count mismatch: System has {expected}, snapshot has {found}")
            },
            Self::IntegratorMismatch { expected, found } => {
                write!(f, "integrator mismatch: record header says {expected}, System is {found}")
            },
            Self::Integrator(e) => write!(f, "integrator restore: {e}"),
        }
    }
}

impl std::error::Error for RestoreError {}

impl From<RecordError> for RestoreError {
    fn from(e: RecordError) -> Self {
        Self::Record(e)
    }
}

impl From<ResumeError> for RestoreError {
    fn from(e: ResumeError) -> Self {
        Self::Integrator(e)
    }
}

/// Restore `sys` to the dynamic state captured at the n-th Snapshot
/// frame of `record`. The System's body count, integrator kind, and
/// units must already match the record; body positions/velocities,
/// integrator scratch, `t`, and `steps` are mutated.
///
/// **Diagnostic baseline:** energy/Lz baselines are taken from the
/// post-restore state, so subsequent `rel_energy_error` /
/// `rel_angular_momentum_error` readings measure drift from
/// `t = snapshot.t`, not from `t = 0` of the original record. Compare
/// against the `Diagnostic` frames in the source record if you need
/// continuity with the original timeline.
pub fn restore_into(
    sys: &mut System,
    record: &Record,
    snapshot_idx: usize,
) -> Result<(), RestoreError> {
    let recorded_kind = record.header().integrator.kind.as_str();
    let active_kind = sys.integrator_kind().label();
    if recorded_kind != active_kind {
        return Err(RestoreError::IntegratorMismatch {
            expected: recorded_kind.to_string(),
            found: active_kind.to_string(),
        });
    }

    let mut snapshots = Vec::new();
    let mut resume_states = Vec::new();
    for f in frames_until_trailer(record)? {
        match f? {
            Frame::Snapshot(s) => snapshots.push(s),
            Frame::ResumeState(r) => resume_states.push(r),
            _ => {},
        }
    }
    if snapshot_idx >= snapshots.len() {
        return Err(RestoreError::SnapshotIndexOutOfRange {
            idx: snapshot_idx,
            available: snapshots.len(),
        });
    }
    if snapshots.len() != resume_states.len() {
        return Err(RestoreError::MissingResumeState { snapshot_idx });
    }

    let snap = &snapshots[snapshot_idx];
    let resume = &resume_states[snapshot_idx];

    if snap.bodies.len() != sys.bodies().len() {
        return Err(RestoreError::BodyCountMismatch {
            expected: sys.bodies().len(),
            found: snap.bodies.len(),
        });
    }

    for (b, s) in sys.bodies.iter_mut().zip(snap.bodies.iter()) {
        b.pos_x = s.pos[0];
        b.pos_y = s.pos[1];
        b.pos_z = s.pos[2];
        b.vel_x = s.vel[0];
        b.vel_y = s.vel[1];
        b.vel_z = s.vel[2];
    }
    sys.t = snap.t;
    sys.steps = resume.step_count;
    sys.integrator.restore_resume_state(&resume.bytes)?;
    sys.refresh_energy_diagnostics();
    Ok(())
}

fn frames_until_trailer(record: &Record) -> Result<Vec<Result<Frame, RecordError>>, RecordError> {
    let snaps: Vec<_> = record.dense()?.map(|r| r.map(Frame::Snapshot)).collect();
    let res: Vec<_> = record.resume_states()?.map(|r| r.map(Frame::ResumeState)).collect();
    let mut out = Vec::with_capacity(snaps.len() + res.len());
    out.extend(snaps);
    out.extend(res);
    Ok(out)
}
